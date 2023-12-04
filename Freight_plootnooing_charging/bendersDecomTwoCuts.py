#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
bendersDecomp21FebWithoutZ
notes:
    1. Both feasibility and benders cut
    2. Dual will be printed just before adding feasibility cut
    3. from model n and z is ommited
"""

# start counting time
import time
tic = time.perf_counter()

# import librarires for formulation and solving optimization model
from pyomo.environ import *
from pyomo.opt import SolverFactory, SolverStatus, TerminationCondition
import random

# Argument parsing strategy through terminal 
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("l", type=float,
                    help="input as route length")
parser.add_argument("cs", type=int,
                    help="input as number of charging station")
parser.add_argument("v", type=int,
                    help="input as number of vehicles")
parser.add_argument("p", type=int,
                    help="input as number of platoons")
args = parser.parse_args()

#======================= Data preparation ======================
# define environment: route length, charging station, vehicle and platoon number
routeLength = args.l# Mile 
numberOfCS = args.cs
csID = range(1,numberOfCS+1)
csIDwithO = range(0, numberOfCS+1)
csIDwithOD = range(0, numberOfCS+2) 
vehicleID = range(0, args.v)
platoonID = range(0, args.p)

# find distance between station of index i and j
distanceBetStation = dict()
for i in csIDwithOD:
	for j in csIDwithOD:
		if i==j:
			distanceBetStation[i,j] = 0
		else:
			distanceBetStation[i,j] = abs((routeLength/numberOfCS)*(i-j))

# collect vehcile configuration from manufacturer sites/ report		
# eCascadia: 475kwh, 1.9kwh/mile  https://www.autoweek.com/news/green-cars/a36506185/electric-big-rig-semi-trucks/
# eM2: 315kwh, 1.4kwh/mile 
# MT50e: 220kwh, 1.76kwh/mile
efficiencyVehicle = 1.5 #{k:random.uniform(1.4,1.9) for k in vehicleID}
batteryCapVehicle = 300 #{k:random.randrange(220,400) for k in vehicleID}


#different departure and tentative arrival time for different vehicles
diffDepTime=dict()
random.seed(a=1)
for i in vehicleID:
    diffDepTime[i]=random.randint(0, 4)
    
diffArrTime=dict()
for i in vehicleID:
    diffArrTime[i]=diffDepTime[i]+20
    

# define solver
solvername = 'cplex'
solver=SolverFactory(solvername)

import csv

# create class from where call different models/functions into main class 
class Models:
    # create fun to furnish master problem and solve    
    def masterProblem(self):
        # Sets
        mp.S = Set(initialize=csID,doc="set of cs")
        mp.J = Set(initialize=vehicleID,doc="set of vehicles")
        mp.I = Set(initialize=platoonID,doc="set of platoons")
        mp.swithO = Set(initialize=csIDwithO,doc="s union O") 
        mp.swithOD = Set(initialize=csIDwithOD,doc="s union OD") 
    
        # Variables
        mp.y = Var(mp.I, mp.J, mp.swithO, within=Binary, doc = "vehicle $j$ joins platoon $i$ at station $s$")
        #mp.n = Var(mp.I, mp.swithO, within = NonNegativeIntegers, doc = "number of vehicle in platoon i at station s")
        mp.w = Var(mp.I,mp.swithO, within =Binary, doc = "indicator variables")
        
        # Additional MP variable
        mp.zLower = Var(within=Reals)
        
        # Parameters
        mp.d = Param(mp.swithOD, mp.swithOD, initialize=distanceBetStation, doc = "distance between s1 and s2")
        mp.m = Param(mp.J, initialize =efficiencyVehicle, doc = "energy efficiency of vehicle kWh/ mile")
        mp.delta = Param(mp.I, mp.swithO, initialize = .13, doc = "energy saving % of platoon i depart from station s")
        mp.c_energy = Param(mp.J, initialize = 0.1165 , doc = "energy cost of vehicle j https://www.eia.gov/state/print.php?sid=FL")
        mp.N_max = Param(initialize = len(vehicleID), doc = "maximum vehicle in a platoon")
        
        energyCost  = sum(mp.c_energy[j]*mp.d[s,s+1]*mp.m[j]*(1-sum(mp.delta[i,s]*mp.y[i,j,s] for i in mp.I)) for j in mp.J for s in mp.swithO)
        
        # Define constraints                                  
        def cons1Rule(mp):
            return mp.zLower>=energyCost
        mp.cons1=Constraint(rule=cons1Rule,doc="")
        
        
        def joinPlatoonRule1(mp,j,s):
        	return -sum(mp.y[i,j,s] for i in mp.I) >= -1
        mp.joinPlatoon1= Constraint(mp.J, mp.swithO, rule = joinPlatoonRule1, doc = "vehicle j join at most 1 platoon at cs s")
        
        
        def joinPlatoonRule2(mp,i,s):
            return sum (mp.y[i,j,s] for j in mp.J) >= 0
        mp.joinPlatoon2 = Constraint(mp.I, mp.swithO, rule = joinPlatoonRule2, doc = "platoon i takes at most 1 vehicle at cs s")
        
# =============================================================================
#         def vehicleNoinPlatoonRule(mp,i,s):
#             return mp.n[i,s] == sum(mp.y[i,j,s] for j in mp.J) 
#         mp.vehicleNoinPlatoon = Constraint(mp.I,mp.swithO,rule=vehicleNoinPlatoonRule)
# =============================================================================
        
        def platoonNoLimitRule1(mp,i,s):
        	#return -mp.n[i,s] >= -mp.w[i,s]*mp.N_max
            return -sum(mp.y[i,j,s] for j in mp.J) >= -mp.w[i,s]*mp.N_max
        mp.platoonNoLimit1 = Constraint(mp.I, mp.swithO, rule = platoonNoLimitRule1)
        
        def platoonNoLimitRule2(mp,i,s):
        	return sum(mp.y[i,j,s] for j in mp.J) >= mp.w[i,s]*2
        mp.platoonNoLimit2 = Constraint(mp.I, mp.swithO, rule = platoonNoLimitRule2)
        
        mp.cuts = ConstraintList()                                       
        
        #define objective function
        def objectiveRule(mp):
            return mp.zLower
                
        mp.objective=Objective(rule=objectiveRule,sense=minimize,doc="objective function")
        
        global resultsMp
        resultsMp = solver.solve(mp)
        
        return mp
        
    # create fun to furnish subproblem and solve    
    def subProblem(self):    
        
        # Sets
        sp.S = Set(initialize=csID,doc="set of cs")
        sp.J = Set(initialize=vehicleID,doc="set of vehicles")
        sp.I = Set(initialize=platoonID,doc="set of platoons")
        sp.swithO = Set(initialize=csIDwithO,doc="s union O") 
        sp.swithOD = Set(initialize=csIDwithOD,doc="s union OD") 
    
    
        # Variables
        sp.t_jdep = Var(sp.J, sp.swithOD, within=NonNegativeReals, doc = "time of vehicle j depart from station s")
        sp.t_jarr = Var(sp.J, sp.swithOD, within = NonNegativeReals, doc = "time of vehicle j arrive at station s")
        sp.t_idep = Var(sp.I, sp.swithOD, within = NonNegativeReals, doc = "time of platoon i depart from station s") #EMPTY
        sp.t_iarr = Var(sp.I, sp.swithOD, within = NonNegativeReals, doc = "time of platoon i arrive at station s") #EMPTY
        sp.soc_arr = Var(sp.J,sp.swithOD, within = NonNegativeReals, doc = "soc of vehicle j at station s when arrive")
        sp.soc_dep = Var(sp.J,sp.swithOD, within = NonNegativeReals, doc = "soc of vehicle j at station s when departs")
        sp.t_cha = Var(sp.J,sp.S, within = NonNegativeReals, doc = "charging time of vehicle j at station s")
        sp.t_wai = Var(sp.J,sp.S, within = NonNegativeReals, doc = "waiting time (after charge) of vehicle j at station s")
        sp.dummy = Var(sp.J,within = NonNegativeReals, doc = "dummy var for delay cost calc")
        
        # Parameters
        sp.d = Param(sp.swithOD, sp.swithOD, initialize=distanceBetStation, doc = "distance between s1 and s2")
        sp.m = Param(sp.J, initialize =efficiencyVehicle, doc = "energy efficiency of vehicle kWh/ mile")
        sp.v_j = Param(sp.J, initialize = 50, doc = "average speed of vehicle j")
        sp.v_i = Param(sp.I, initialize = 50, doc = "average speed of platoon i")
        sp.r_cha = Param(sp.S, initialize = 100, doc = "charging speed kw of station s")
        sp.t_depmin = Param(sp.J, initialize = 0, doc = "earliest depart time of vehicle j")
        sp.t_arrmax = Param(sp.J, initialize = diffArrTime, doc = "latest arrive time hr of vehicle j")
        sp.soc_depmin = Param(sp.J, initialize = 1, doc = "minimum departure soc of vehicle j")
        sp.soc_arrmin = Param(sp.J, initialize = 0.2, doc = "minimum arrival soc of vehicle j")
        sp.soc_min = Param(sp.J, initialize = 0.2 , doc = "minimim soc allowed for vehicle j")
        sp.soc_max = Param(sp.J, initialize = 1,doc = "maximum soc allowed for vehicle j")
        sp.M = Param(initialize = 25, doc = "big number")
        sp.delta = Param(sp.I, sp.swithO, initialize = .13, doc = "energy saving % of platoon i depart from station s")
        sp.c_cha = Param(sp.S, initialize =.35 , doc = "charging cost per kwh at station s per kWh https://freewiretech.com/difference-between-ev-charging-levels/")
        sp.c_delay = Param(sp.J, initialize = 50 , doc = "delay cost of vehicle j ")
        sp.c_energy = Param(sp.J, initialize = 0.1165 , doc = "energy cost of vehicle j https://www.eia.gov/state/print.php?sid=FL")
        sp.b = Param(sp.J,initialize = batteryCapVehicle, doc = "battery capacity of vehicle j")
        sp.N_max = Param(initialize = len(vehicleID), doc = "maximum vehicle in a platoon")
        sp.t_dep = Param(sp.J, initialize = diffDepTime, doc = "different departure time for vehicles")
    
                                     
        # Constraints
        # SOC 
        def initialDepartSOCRule(sp,j):
        	return sp.soc_dep[j,0] == sp.soc_depmin[j]
        sp.initialSOC = Constraint(sp.J, rule = initialDepartSOCRule, doc = "init depart SOC = depat soc from origin")
    
        def lastArrivalSOCRule(sp,j):
        	return sp.soc_arr[j,list(sp.swithOD.data())[-1]]>=sp.soc_arrmin[j]
        sp.lastArrivalSOC = Constraint(sp.J, rule = lastArrivalSOCRule, doc = "last arrival SOC = arrival minimum soc")
    
        def departSOCfromCSRule(sp,j,s):
        	return sp.soc_dep[j,s] -sp.soc_arr[j,s] == (sp.r_cha[s]/sp.b[j])*sp.t_cha[j,s] 
        sp.departSOCfromCS = Constraint(sp.J, sp.S, rule = departSOCfromCSRule)
    
        def arriveSOCatCSRule_ineq1(sp,j,s):
            return sp.soc_arr[j,s+1]-sp.soc_dep[j,s] == - ((sp.d[s,s+1]*sp.m[j])/sp.b[j])*(1 - sum(sp.delta[i,s]*mp.y[i,j,s].value for i in sp.I))
    
        sp.arriveSOCatCS_ineq1 = Constraint(sp.J,sp.swithO, rule = arriveSOCatCSRule_ineq1)
            
        def arriveSOCMinRule(sp,j,s):
        	return sp.soc_arr[j,s] >= sp.soc_min[j]
        sp.arriveSOCMin = Constraint(sp.J, sp.S, rule = arriveSOCMinRule)
    
        def departSOCMaxRule(sp,j,s):
        	return -sp.soc_dep[j,s] >= -sp.soc_max[j] 
        sp.departSOCMax = Constraint(sp.J, sp.S, rule = departSOCMaxRule)
        
        # Time
        def initialDepartTimeRule(sp,j):
            return sp.t_jdep[j,0]>=sp.t_depmin[j] 
        sp.initialDepartTime = Constraint(sp.J, rule = initialDepartTimeRule)
    
        def departTimeRule_ineq1(sp,j,s):
            return sp.t_jdep[j,s]-sp.t_jarr[j,s]-sp.t_cha[j,s]- sp.t_wai[j,s] == 0
        sp.departTime_ineq1 = Constraint(sp.J, sp.S, rule = departTimeRule_ineq1)
            
        def arriveTimeRule(sp,j,s):
        	return sp.t_jarr[j,s+1] - sp.t_jdep[j,s] == sp.d[s,s+1]/sp.v_j[j]
        sp.arriveTime = Constraint(sp.J, sp.swithO, rule = arriveTimeRule)
    
        def platoonDepartRule(sp,i,s):
            return sp.t_idep[i,s] == sp.t_iarr[i,s]
        sp.platoonDepart = Constraint(sp.I,sp.S,rule=platoonDepartRule,doc="platton arrive time and depart time same")
    
        def platoonArriveRule(sp,i,s):
            return sp.t_iarr[i,s+1] - sp.t_idep[i,s] == sp.d[s,s+1]/sp.v_i[i]
        sp.platoonArriveTime = Constraint(sp.I, sp.swithO, rule = platoonArriveRule)
    
        def platoonVehicleDepRule1(sp,i,j,s):
        	return sp.t_jdep[j,s] - sp.t_idep[i,s] >= -sp.M*(1-mp.y[i,j,s].value)  
        sp.platoonVehicleDep1 = Constraint(sp.I, sp.J, sp.swithO, rule = platoonVehicleDepRule1)
    
        def platoonVehicleDepRule2(sp,i,j,s):
        	return  -sp.t_jdep[j,s] + sp.t_idep[i,s] >= -sp.M*(1-mp.y[i,j,s].value) 
        sp.platoonVehicleDep2 = Constraint(sp.I, sp.J, sp.swithO, rule = platoonVehicleDepRule2)
    
        def nonNegativeDelayRule1(sp,j):
            return sp.dummy[j] - sp.t_jarr[j,list(sp.swithOD.data())[-1]]>=-sp.t_arrmax[j] 
        sp.nonNegativeDelay1 = Constraint(sp.J,rule = nonNegativeDelayRule1)
                                                       
        #define objective function
        def objectiveRuleSp1(sp):
            chargingCost = sum(sp.c_cha[s]*sp.r_cha[s]*sp.t_cha[j,s] for j in sp.J for s in sp.S)
            delayCost = sum(sp.c_delay[j]*sp.dummy[j] for j in sp.J)
            
            return (chargingCost+delayCost)
                
        sp.objective=Objective(rule=objectiveRuleSp1,sense=minimize,doc="objective function")
        sp.dual = Suffix(direction = Suffix.IMPORT)
        
        global resultsSp
        resultsSp = solver.solve(sp)
        
        return sp
    
    # create fun to furnish new subproblem (when initial subproblem is infeasible) and solve    
    def newSubProblemwSlack(self):
        
        # Sets
        newSp1.S = Set(initialize=csID,doc="set of cs")
        newSp1.J = Set(initialize=vehicleID,doc="set of vehicles")
        newSp1.I = Set(initialize=platoonID,doc="set of platoons")
        newSp1.swithO = Set(initialize=csIDwithO,doc="s union O") 
        newSp1.swithOD = Set(initialize=csIDwithOD,doc="s union OD") 
        
        # Define variables
        # slack variable
        newSp1.sPlatVehDep1 = Var(newSp1.I, newSp1.J, newSp1.swithO, within = NonNegativeReals)
        newSp1.sPlatVehDep2 = Var(newSp1.I, newSp1.J, newSp1.swithO, within = NonNegativeReals)
        newSp1.sArriveSOCatCS_ineq1 = Var(newSp1.J, newSp1.swithO, within = NonNegativeReals, doc = "slack for equality with ycap")
        newSp1.sArriveSOCatCS_ineq2 = Var(newSp1.J, newSp1.swithO, within = NonNegativeReals, doc = "slack for equality with ycap")
        
        # variables similar to initial subproblem
        newSp1.t_jdep = Var(newSp1.J, newSp1.swithOD, within=NonNegativeReals, doc = "time of vehicle j depart from station s")
        newSp1.t_jarr = Var(newSp1.J, newSp1.swithOD, within = NonNegativeReals, doc = "time of vehicle j arrive at station s")
        newSp1.t_idep = Var(newSp1.I, newSp1.swithOD, within = NonNegativeReals, doc = "time of platoon i depart from station s") #Emodels.masterProblem()TY
        newSp1.t_iarr = Var(newSp1.I, newSp1.swithOD, within = NonNegativeReals, doc = "time of platoon i arrive at station s") #EMPTY
        newSp1.soc_arr = Var(newSp1.J,newSp1.swithOD, within = NonNegativeReals, doc = "soc of vehicle j at station s when arrive")
        newSp1.soc_dep = Var(newSp1.J,newSp1.swithOD, within = NonNegativeReals, doc = "soc of vehicle j at station s when departs")
        newSp1.t_cha = Var(newSp1.J,newSp1.S, within = NonNegativeReals, doc = "charging time of vehicle j at station s")
        newSp1.t_wai = Var(newSp1.J,newSp1.S, within = NonNegativeReals, doc = "waiting time (after charge) of vehicle j at station s")
        newSp1.dummy = Var(newSp1.J,within = NonNegativeReals, doc = "dummy var for delay cost calc")
        
        # Define parameters
        newSp1.d = Param(newSp1.swithOD, newSp1.swithOD, initialize=distanceBetStation, doc = "distance between s1 and s2")
        newSp1.m = Param(newSp1.J, initialize =efficiencyVehicle, doc = "energy efficiency of vehicle kWh/ mile")
        newSp1.v_j = Param(newSp1.J, initialize = 50, doc = "average newSp1eed of vehicle j")
        newSp1.v_i = Param(newSp1.I, initialize = 50, doc = "average newSp1eed of platoon i")
        newSp1.r_cha = Param(newSp1.S, initialize = 100, doc = "charging newSp1eed kw of station s")
        newSp1.t_que = Param(newSp1.J, newSp1.S, initialize = .25, doc = "waiting time hr before charge at s")
        newSp1.t_depmin = Param(newSp1.J, initialize = 0, doc = "earliest depart time of vehicle j")
        newSp1.t_arrmax = Param(newSp1.J, initialize = diffArrTime, doc = "latest arrive time hr of vehicle j")
        newSp1.soc_depmin = Param(newSp1.J, initialize = 1, doc = "minimum departure soc of vehicle j")
        newSp1.soc_arrmin = Param(newSp1.J, initialize = 0.2, doc = "minimum arrival soc of vehicle j")
        newSp1.soc_min = Param(newSp1.J, initialize = 0.2 , doc = "minimim soc allowed for vehicle j")
        newSp1.soc_max = Param(newSp1.J, initialize = 1,doc = "maximum soc allowed for vehicle j")
        newSp1.M = Param(initialize = 25, doc = "big number")
        newSp1.delta = Param(newSp1.I, newSp1.swithO, initialize = .13, doc = "energy saving % of platoon i depart from station s")
        newSp1.c_cha = Param(newSp1.S, initialize =.35 , doc = "charging cost per kwh at station s per kWh https://freewiretech.com/difference-between-ev-charging-levels/")
        newSp1.c_delay = Param(newSp1.J, initialize = 50 , doc = "delay cost of vehicle j ")
        newSp1.c_energy = Param(newSp1.J, initialize = 0.1165 , doc = "energy cost of vehicle j https://www.eia.gov/state/print.php?sid=FL")
        newSp1.b = Param(newSp1.J,initialize = batteryCapVehicle, doc = "battery capacity of vehicle j")
        newSp1.N_max = Param(initialize = len(vehicleID), doc = "maximum vehicle in a platoon")
        
       

        # Eequality constraints
        def initialDepartSOCRule(newSp1,j):
        	return newSp1.soc_dep[j,0] ==  newSp1.soc_depmin[j]
        newSp1.initialSOC = Constraint(newSp1.J, rule = initialDepartSOCRule, doc = "init depart SOC = depat soc from origin")
        
            
        def departSOCfromCSRule(newSp1,j,s):
        	return newSp1.soc_dep[j,s] -newSp1.soc_arr[j,s] == (newSp1.r_cha[s]/newSp1.b[j])*newSp1.t_cha[j,s] 
        newSp1.departSOCfromCS = Constraint(newSp1.J, newSp1.S, rule = departSOCfromCSRule)
        
        # slack added
        def arriveSOCatCSRule_ineq1(newSp1,j,s):
            return newSp1.soc_arr[j,s+1]-newSp1.soc_dep[j,s]+newSp1.sArriveSOCatCS_ineq1[j,s] == -((newSp1.d[s,s+1]*newSp1.m[j])/newSp1.b[j])*(1 - sum(newSp1.delta[i,s]*mp.y[i,j,s].value for i in newSp1.I))
        newSp1.arriveSOCatCS_ineq1 = Constraint(newSp1.J,newSp1.swithO, rule = arriveSOCatCSRule_ineq1)
        
        def departTimeRule_ineq1(newSp1,j,s):
            return newSp1.t_jdep[j,s]-newSp1.t_jarr[j,s]-newSp1.t_cha[j,s]-newSp1.t_wai[j,s] == 0
        newSp1.departTime_ineq1 = Constraint(newSp1.J, newSp1.S, rule = departTimeRule_ineq1)
        
        def arriveTimeRule(newSp1,j,s):
        	return newSp1.t_jarr[j,s+1] -newSp1.t_jdep[j,s] == newSp1.d[s,s+1]/newSp1.v_j[j]
        newSp1.arriveTime = Constraint(newSp1.J, newSp1.swithO, rule = arriveTimeRule)
        
        def platoonDepartRule(newSp1,i,s):
            return newSp1.t_idep[i,s] == newSp1.t_iarr[i,s]
        newSp1.platoonDepart = Constraint(newSp1.I,newSp1.S,rule=platoonDepartRule,doc="platton arrive time and depart time same")
    
        def platoonArriveRule(newSp1,i,s):
            return newSp1.t_iarr[i,s+1]-newSp1.t_idep[i,s] == newSp1.d[s,s+1]/newSp1.v_i[i]
        newSp1.platoonArriveTime = Constraint(newSp1.I, newSp1.swithO, rule = platoonArriveRule)
        
        # Constraints
        def lastArrivalSOCRule(newSp1,j):
        	return newSp1.soc_arr[j,list(newSp1.swithOD.data())[-1]]>=newSp1.soc_arrmin[j]
        newSp1.lastArrivalSOC = Constraint(newSp1.J, rule = lastArrivalSOCRule, doc = "last arrival SOC = arrival minimum soc")
         
        def arriveSOCMinRule(newSp1,j,s):
        	return newSp1.soc_arr[j,s]>=newSp1.soc_min[j]
        newSp1.arriveSOCMin = Constraint(newSp1.J, newSp1.S, rule = arriveSOCMinRule)
         
        def departSOCMaxRule(newSp1,j,s):
        	return -newSp1.soc_dep[j,s] >=-newSp1.soc_max[j]
        newSp1.departSOCMax = Constraint(newSp1.J, newSp1.S, rule = departSOCMaxRule)
         
        # Time
        def initialDepartTimeRule(newSp1,j):
            return newSp1.t_jdep[j,0] >= newSp1.t_depmin[j]   
        newSp1.initialDepartTime = Constraint(newSp1.J, rule = initialDepartTimeRule)
       
        # slack added
        def platoonVehicleDepRule1(newSp1,i,j,s):
            return newSp1.t_jdep[j,s]-newSp1.t_idep[i,s]+newSp1.sPlatVehDep1[i,j,s] >= -newSp1.M*(1-mp.y[i,j,s].value)
        newSp1.platoonVehicleDep1 = Constraint(newSp1.I, newSp1.J, newSp1.swithO, rule = platoonVehicleDepRule1)
        
        # slack added
        def platoonVehicleDepRule2(newSp1,i,j,s):
            return  newSp1.t_idep[i,s]-newSp1.t_jdep[j,s]  + newSp1.sPlatVehDep2[i,j,s]>= -(newSp1.M*(1-mp.y[i,j,s].value))# FIND DUAL
        newSp1.platoonVehicleDep2 = Constraint(newSp1.I, newSp1.J, newSp1.swithO, rule = platoonVehicleDepRule2)
    
        def nonNegativeDelayRule1(newSp1,j):
            return newSp1.dummy[j] - newSp1.t_jarr[j,list(newSp1.swithOD.data())[-1]]>= -newSp1.t_arrmax[j] 
        newSp1.nonNegativeDelay1 = Constraint(newSp1.J,rule = nonNegativeDelayRule1)
        
        
        # Objective
        def objectiveRuleNewSp1(newSp1):
            slackSum1 = sum(newSp1.sPlatVehDep1[i,j,s] for i in newSp1.I for j in newSp1.J for s in newSp1.S)
            slackSum2 = sum(newSp1.sPlatVehDep2[i,j,s] for i in newSp1.I for j in newSp1.J for s in newSp1.S)
            slackSum3 = sum(newSp1.sArriveSOCatCS_ineq1[j,s] for j in newSp1.J for s in newSp1.swithO)
            
            return slackSum1+slackSum2+slackSum3       
        newSp1.objective = Objective(rule = objectiveRuleNewSp1, sense = minimize)
        
        newSp1.dual = Suffix(direction = Suffix.IMPORT)
        
        global resultsnewSp1
        resultsnewSp1 = solver.solve(newSp1)
        
        return newSp1 
    
# Bender's decomposition algortihm
# create object of class Models
models = Models()
# create master problem concrete model
mp = ConcreteModel()
# solve master problem
models.masterProblem()
# estimate lower bound
LB = value(mp.objective)
print(f'Initial LB = {LB}')
# assume upper bound initially
import math
UB = math.inf
print(f'Initial UB = {UB}')

# import library to write file
import json

# loop continue until the gap between upper and lower bound is less than 0.05 default value: .00001,\
    # but here, i am easygoing.

# For small testing, use itearate
iterationNo = 0
while (UB-LB >= 0.05): #and iterationNo<10):  
    # initialize empty dictionary for dual
    dual_platoonVehicleDep1 = {}
    dual_platoonVehicleDep2 = {}
    dual_arriveSOCatCS = {}
    dual_iniSOC = {}
    dual_lastArrivalSOC = {}
    dual_initialDepartTime = {}
    dual_arriveTime = {}
    dual_arriveSOCMin = {}
    dual_departSOCMax = {}
    dual_platoonArriveTime = {}
    dual_nonNegativeDelay1 = {}
    # crete subproblem
    sp = ConcreteModel()
    # solev subproblem
    models.subProblem()
   
    # Ops! follow tasks when primal subproblem is infeasible i.e., dual of suproblem is unbounded
    if (resultsSp.solver.status == SolverStatus.ok) and (resultsSp.solver.termination_condition == TerminationCondition.infeasible):
        
        # create new subproblem
        newSp1 = ConcreteModel()
        # solve new subproblem
        models.newSubProblemwSlack()
        
        # estimate dual values and fill in empty dictionary
        for j in newSp1.J:
            for s in newSp1.S:
                dual_arriveSOCMin[j,s] = newSp1.dual[newSp1.arriveSOCMin[j,s]]
                dual_departSOCMax[j,s] = newSp1.dual[newSp1.departSOCMax[j,s]]
            
        for i in newSp1.I:
            for j in newSp1.J:
                for s in newSp1.swithO:
                    dual_platoonVehicleDep1[i,j,s] = newSp1.dual[newSp1.platoonVehicleDep1[i,j,s]]
                    dual_platoonVehicleDep2[i,j,s] = newSp1.dual[newSp1.platoonVehicleDep2[i,j,s]]
                    
        for j in newSp1.J:
            for s in newSp1.swithO:
                dual_arriveSOCatCS[j,s] = newSp1.dual[newSp1.arriveSOCatCS_ineq1[j,s]]
                dual_arriveTime[j,s] = newSp1.dual[newSp1.arriveTime[j,s]]

        for j in newSp1.J:
            dual_iniSOC[j] = newSp1.dual[newSp1.initialSOC[j]]
            dual_lastArrivalSOC[j] = newSp1.dual[newSp1.lastArrivalSOC[j]]
            dual_nonNegativeDelay1[j] = newSp1.dual[newSp1.nonNegativeDelay1[j]]
            dual_initialDepartTime[j] = newSp1.dual[newSp1.initialDepartTime[j]]

        for i in newSp1.I:
            for s in newSp1.swithO:
            	dual_platoonArriveTime[i,s] = newSp1.dual[newSp1.platoonArriveTime[i,s]]
       
        
        # test if dual is changing at each iteartion
# =============================================================================
#         nDual_platoonVehicleDep1 = {str(key): value for key, value in dual_platoonVehicleDep1.items()} 
#         nDual_platoonVehicleDep2 = {str(key): value for key, value in dual_platoonVehicleDep2.items()}
#         nDual_arriveSOCatCS = {str(key): value for key, value in dual_arriveSOCatCS.items()}
#         nDual_iniSOC = {str(key): value for key, value in dual_iniSOC.items()}
#         nDual_lastArrivalSOC = {str(key): value for key, value in dual_lastArrivalSOC.items()}
#         nDual_nonNegativeDelay1 = {str(key): value for key, value in dual_nonNegativeDelay1.items()}
#         nDual_initialDepartTime = {str(key): value for key, value in dual_initialDepartTime.items()}
#         nDual_arriveTime = {str(key): value for key, value in  dual_arriveTime.items()}
#         nDual_arriveSOCMin = {str(key): value for key, value in dual_arriveSOCMin.items()}
#         nDual_departSOCMax = {str(key): value for key, value in dual_departSOCMax.items()}
#         nDual_platoonArriveTime = {str(key): value for key, value in dual_platoonArriveTime.items()}
# =============================================================================
        
        
# =============================================================================
#         with open('nDual_platoonVehicleDep1_iterationNo'+str(iterationNo)+'.json', 'w') as fd:
#             fd.write(json.dumps(nDual_platoonVehicleDep1))    
#         with open('nDual_platoonVehicleDep2_iterationNo'+str(iterationNo)+'.json', 'w') as fd:
#             fd.write(json.dumps(nDual_platoonVehicleDep2)) 
#         with open('nDual_arriveSOCatCS_iterationNo'+str(iterationNo)+'.json', 'w') as fd:
#             fd.write(json.dumps(nDual_arriveSOCatCS)) 
#         with open('nDual_iniSOC_iterationNo'+str(iterationNo)+'.json', 'w') as fd:
#             fd.write(json.dumps(nDual_iniSOC)) 
#         with open('nDual_lastArrivalSOC_iterationNo'+str(iterationNo)+'.json','w') as fd:
#             fd.write(json.dumps(nDual_lastArrivalSOC)) 
#         with open('nDual_nonNegativeDelay1_iterationNo'+str(iterationNo)+'.json', 'w') as fd:
#             fd.write(json.dumps(nDual_nonNegativeDelay1))
#         with open('nDual_initialDepartTime_iterationNo'+str(iterationNo)+'.json', 'w') as fd:
#             fd.write(json.dumps(nDual_initialDepartTime))
#         with open('nDual_arriveTime_iterationNo'+str(iterationNo)+'.json', 'w') as fd:
#             fd.write(json.dumps(nDual_arriveTime))
#         with open('nDual_arriveSOCMin_iterationNo'+str(iterationNo)+'.json', 'w') as fd:
#             fd.write(json.dumps(nDual_arriveSOCMin))
#         with open('nDual_departSOCMax_iterationNo'+str(iterationNo)+'.json', 'w') as fd:
#             fd.write(json.dumps(nDual_departSOCMax))
#         with open('nDual_platoonArriveTime_iterationNo'+str(iterationNo)+'.json', 'w') as fd:
#             fd.write(json.dumps(nDual_platoonArriveTime))
# =============================================================================
            
        # generate feasibility cut            
        def feasibilityCutRule(mp):
             return sum(-newSp1.M*(1-mp.y[i,j,s])*dual_platoonVehicleDep1[i,j,s] for i in newSp1.I for j in newSp1.J for s in newSp1.swithO)\
                +sum(-(newSp1.M*(1-mp.y[i,j,s]))*dual_platoonVehicleDep2[i,j,s] for i in newSp1.I for j in newSp1.J for s in newSp1.swithO)\
                +sum(-((newSp1.d[s,s+1]*newSp1.m[j])/newSp1.b[j])*(1 - sum(newSp1.delta[i,s]*mp.y[i,j,s] for i in newSp1.I))*dual_arriveSOCatCS[j,s] for j in newSp1.J for s in newSp1.swithO)\
                +sum(newSp1.soc_arrmin[j]*dual_lastArrivalSOC[j] for j in newSp1.J)\
                +sum(newSp1.soc_depmin[j]*dual_iniSOC[j] for j in newSp1.J)\
                +sum((-newSp1.t_arrmax[j])*dual_nonNegativeDelay1[j] for j in newSp1.J)\
                +sum(newSp1.t_depmin[j]*dual_initialDepartTime[j] for j in newSp1.J)\
                +sum((newSp1.d[s,s+1]/newSp1.v_j[j])*dual_arriveTime[j,s] for j in newSp1.J for s in newSp1.swithO)\
                +sum(newSp1.soc_min[j]*dual_arriveSOCMin[j,s] for j in newSp1.J for s in newSp1.S)\
                +sum((-newSp1.soc_max[j])*dual_departSOCMax[j,s] for j in newSp1.J for s in newSp1.S)\
                +sum((newSp1.d[s,s+1]/newSp1.v_i[i])*dual_platoonArriveTime[i,s] for i in newSp1.I for s in newSp1.swithO)<= 0
        
        # add feasibility cut in master problem
        mp.cuts.add(feasibilityCutRule(mp))
        
    # if initial subproblem is feasible, do as follows.    
    if (resultsSp.solver.status == SolverStatus.ok) and (resultsSp.solver.termination_condition == TerminationCondition.optimal):
        # estimate the original problem's objective  
        energyCostWithYcap  = sum(mp.c_energy[j]*mp.d[s,s+1]*mp.m[j]*(1-sum(mp.delta[i,s]*mp.y[i,j,s].value for i in mp.I)) for j in mp.J for s in mp.swithO)
        # update the upper bound
        UB = min(UB, (value(sp.objective)+energyCostWithYcap))
        # print UB when it is updated
        print(f'UB = {UB}')
        # fill in empty dictionary with dual of constraints from subproblems
        for j in sp.J:
            for s in sp.S:
                dual_arriveSOCMin[j,s] = sp.dual[sp.arriveSOCMin[j,s]]
                dual_departSOCMax[j,s] = sp.dual[sp.departSOCMax[j,s]]
            
        
        for i in sp.I:
            for j in sp.J:
                for s in sp.swithO:
                    dual_platoonVehicleDep1[i,j,s] = sp.dual[sp.platoonVehicleDep1[i,j,s]]
                    dual_platoonVehicleDep2[i,j,s] = sp.dual[sp.platoonVehicleDep2[i,j,s]]
                    
        for j in sp.J:
            for s in sp.swithO:
                dual_arriveSOCatCS[j,s] = sp.dual[sp.arriveSOCatCS_ineq1[j,s]]
                #dual4_2[j,s] = sp.dual[sp.arriveSOCatCS_ineq2[j,s]]
                dual_arriveTime[j,s] = sp.dual[sp.arriveTime[j,s]]

        for j in sp.J:
            dual_iniSOC[j] = sp.dual[sp.initialSOC[j]]
            dual_lastArrivalSOC[j] = sp.dual[sp.lastArrivalSOC[j]]
            dual_nonNegativeDelay1[j] = sp.dual[sp.nonNegativeDelay1[j]]
            dual_initialDepartTime[j] = sp.dual[sp.initialDepartTime[j]]

        for i in sp.I:
            for s in sp.swithO:
            	dual_platoonArriveTime[i,s] = sp.dual[sp.platoonArriveTime[i,s]]
        
        # generate optimality cuts
        def optimalityCutRule(mp):
            return mp.zLower>=sum(mp.c_energy[j]*mp.d[s,s+1]*mp.m[j]*(1-sum(mp.delta[i,s]*mp.y[i,j,s] for i in mp.I)) for j in mp.J for s in mp.swithO)\
                +sum(-sp.M*(1-mp.y[i,j,s])*dual_platoonVehicleDep1[i,j,s] for i in sp.I for j in sp.J for s in sp.swithO)\
                +sum(-(sp.M*(1-mp.y[i,j,s]))*dual_platoonVehicleDep2[i,j,s] for i in sp.I for j in sp.J for s in sp.swithO)\
                +sum(-((sp.d[s,s+1]*sp.m[j])/sp.b[j])*(1 - sum(sp.delta[i,s]*mp.y[i,j,s] for i in sp.I))*dual_arriveSOCatCS[j,s] for j in sp.J for s in sp.swithO)\
                +sum(sp.soc_arrmin[j]*dual_lastArrivalSOC[j] for j in sp.J)\
                +sum(sp.soc_depmin[j]*dual_iniSOC[j] for j in sp.J)\
                +sum((-sp.t_arrmax[j])*dual_nonNegativeDelay1[j] for j in sp.J)\
                +sum(sp.t_depmin[j]*dual_initialDepartTime[j] for j in sp.J)\
                +sum((sp.d[s,s+1]/sp.v_j[j])*dual_arriveTime[j,s] for j in sp.J for s in sp.swithO)\
                +sum(sp.soc_min[j]*dual_arriveSOCMin[j,s] for j in sp.J for s in sp.S)\
                +sum((-sp.soc_max[j])*dual_departSOCMax[j,s] for j in sp.J for s in sp.S)\
                +sum((sp.d[s,s+1]/sp.v_i[i])*dual_platoonArriveTime[i,s] for i in sp.I for s in sp.swithO)
        
        # add optimality cuts in master problem
        mp.cuts.add(optimalityCutRule(mp))
        
        
    # solve the master problem
    resultsMp=solver.solve(mp)
    
# =============================================================================
#     with open('masterProblemVar'+str(iterationNo)+'.text', 'w') as fd:
#         for v in mp.component_objects(Var):
#             fd.write(str(v))
#             fd.write(str(v.get_values()))
#     
# =============================================================================
    LB = value(mp.objective)
    print(f'LB = {LB}\n')
    #iterationNo = iterationNo+1
    
toc = time.perf_counter()
print(f'Processing time = {toc-tic}\n')
    


# =============================================================================
#     for c in mp.component_objects(Constraint, active=True):
#         print ("   Constraint",c)
#         for index in c:
#             print ("      ", index, c[index])
# =============================================================================
 
    #iteration += 1
    
    
    
