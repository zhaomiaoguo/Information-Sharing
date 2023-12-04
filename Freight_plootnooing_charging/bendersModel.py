#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 23 14:10:01 2022

@author: rakibulalam
"""

# Library
from pyomo.environ import *
from pyomo.opt import SolverFactory, SolverStatus, TerminationCondition
import random

routeLength = 500 # Mile
#stationGap = 150 # Mile
numberOfCS = 15
#csID = range(1,int(routeLength/stationGap)-1)
csID = range(1,numberOfCS+1)
#csIDwithO = range(0, int(routeLength/stationGap)-1)
csIDwithO = range(0, numberOfCS+1)
#csIDwithOD = range(0,int(routeLength/stationGap)) 
csIDwithOD = range(0, numberOfCS+2) 
vehicleID = range(0,3)
platoonID = range(0,3)

distanceBetStation = dict()
for i in csIDwithOD:
	for j in csIDwithOD:
		if i==j:
			distanceBetStation[i,j] = 0
		else:
			distanceBetStation[i,j] = abs((routeLength/numberOfCS)*(i-j))


		
# eCascadia: 475kwh, 1.9kwh/mile  https://www.autoweek.com/news/green-cars/a36506185/electric-big-rig-semi-trucks/
# eM2: 315kwh, 1.4kwh/mile 
# MT50e: 220kwh, 1.76kwh/mile
efficiencyVehicle = 1.5#{k:random.uniform(1.4,1.9) for k in vehicleID}
batteryCapVehicle = 300#{k:random.randrange(220,400) for k in vehicleID}


#different departure time for differnt vehicles
diffDepTime=dict()
for i in vehicleID:
    diffDepTime[i]=random.randint(0, 4)
    
diffArrTime=dict()
for i in vehicleID:
    diffArrTime[i]=diffDepTime[i]+11.16 
# =============================================================================
# # Model
# model = ConcreteModel()
# 
# # Sets
# model.S = Set(initialize=csID,doc="set of cs")
# model.J = Set(initialize=vehicleID,doc="set of vehicles")
# model.I = Set(initialize=platoonID,doc="set of platoons")
# model.swithO = Set(initialize=csIDwithO,doc="s union O") 
# model.swithOD = Set(initialize=csIDwithOD,doc="s union OD") 
# 
# 
# # Variables
# model.y = Var(model.I, model.J, model.swithO, within=Binary, doc = "vehicle $j$ joins platoon $i$ at station $s$")
# model.t_jdep = Var(model.J, model.swithOD, within=NonNegativeReals, doc = "time of vehicle j depart from station s")
# model.t_jarr = Var(model.J, model.swithOD, within = NonNegativeReals, doc = "time of vehicle j arrive at station s")
# model.t_idep = Var(model.I, model.swithOD, within = NonNegativeReals, doc = "time of platoon i depart from station s") #EMPTY
# model.t_iarr = Var(model.I, model.swithOD, within = NonNegativeReals, doc = "time of platoon i arrive at station s") #EMPTY
# model.tau = Var(model.I,model.swithOD, within = NonNegativeReals, doc="time of platoon i arrive/depart at/from station s")
# model.soc_arr = Var(model.J,model.swithOD, within = NonNegativeReals, doc = "soc of vehicle j at station s when arrive")
# model.soc_dep = Var(model.J,model.swithOD, within = NonNegativeReals, doc = "soc of vehicle j at station s when departs")
# model.t_cha = Var(model.J,model.S, within = NonNegativeReals, doc = "charging time of vehicle j at station s")
# model.t_wai = Var(model.J,model.S, within = NonNegativeReals, doc = "waiting time (after charge) of vehicle j at station s")
# model.n = Var(model.I, model.swithO, within = NonNegativeIntegers, doc = "number of vehicle in platoon i at station s")
# #model.delta = Var(model.I, model.swithO, doc = "energy saving % of platoon i depart from station s")
# model.z = Var(model.J, model.S, within = Binary, doc = "indicator variable")
# model.dummy = Var(model.J,within = Reals, doc = "dummy var for delay cost calc")
# model.w = Var(model.I,model.swithO, within =Binary, doc = "indicator variables")
# # to address, KeyError: "Index '(0, 0)' is not valid for indexed component 'delta'", model.delta index updated
# # Parameters
# model.d = Param(model.swithOD, model.swithOD, initialize=distanceBetStation, doc = "distance between s1 and s2")
# model.m = Param(model.J, initialize =efficiencyVehicle, doc = "energy efficiency of vehicle kWh/ mile")
# model.v_j = Param(model.J, initialize = 50, doc = "average speed of vehicle j")
# model.v_i = Param(model.I, initialize = 50, doc = "average speed of platoon i")
# model.r_cha = Param(model.S, initialize = 100, doc = "charging speed kw of station s")
# model.t_que = Param(model.J, model.S, initialize = .25, doc = "waiting time hr before charge at s")
# model.t_depmin = Param(model.J, initialize = 0, doc = "earliest depart time of vehicle j")
# model.t_arrmax = Param(model.J, initialize = diffArrTime, doc = "latest arrive time hr of vehicle j")
# model.soc_depmin = Param(model.J, initialize = 1, doc = "minimum departure soc of vehicle j")
# model.soc_arrmin = Param(model.J, initialize = 0.2, doc = "minimum arrival soc of vehicle j")
# model.soc_min = Param(model.J, initialize = 0.2 , doc = "minimim soc allowed for vehicle j")
# model.soc_max = Param(model.J, initialize = 1,doc = "maximum soc allowed for vehicle j")
# model.M = Param(initialize = 10, doc = "big number")
# model.delta = Param(model.I, model.swithO, initialize = .13, doc = "energy saving % of platoon i depart from station s")
# model.c_cha = Param(model.S, initialize =.35 , doc = "charging cost per kwh at station s per kWh https://freewiretech.com/difference-between-ev-charging-levels/")
# model.c_delay = Param(model.J, initialize = 50 , doc = "delay cost of vehicle j ")
# model.c_energy = Param(model.J, initialize = 0.1165 , doc = "energy cost of vehicle j https://www.eia.gov/state/print.php?sid=FL")
# model.b = Param(model.J,initialize = batteryCapVehicle, doc = "battery capacity of vehicle j")
# model.N_max = Param(initialize = len(vehicleID), doc = "maximum vehicle in a platoon")
# model.t_dep = Param(model.J, initialize = diffDepTime, doc = "different departure time for vehicles")
# 
# 
# 
# # Objetive
# def objectiveRule(model):
#     chargingCost = sum(model.c_cha[s]*model.r_cha[s]*model.t_cha[j,s] for j in model.J for s in model.S)
#     #delayCost = sum(model.c_delay[j]*(model.t_jarr[j,list(model.swithOD.value)[-1]]-model.t_arrmax[j]) for j in model.J)
#     delayCost = sum(model.c_delay[j]*model.dummy[j] for j in model.J)
#     energyCost  = sum(model.c_energy[j]*model.d[s,s+1]*model.m[j]*(1-sum(model.delta[i,s]*model.y[i,j,s] for i in model.I)) for j in model.J for s in model.swithO)
#     return chargingCost + delayCost + energyCost
# model.objective = Objective(rule = objectiveRule, sense = minimize)
# 
# 
# 
# 
# solver = SolverFactory("cplex")
# #from datetime import datetime
# #t1 = datetime.now()
# results=solver.solve(model)
# #t2 = datetime.now()
# #print(f'solver time: {(t2.second)-(t1.second)}')
# 
# results.write()
# =============================================================================



# Master problem

# Create model
mp = ConcreteModel()

# Sets
mp.S = Set(initialize=csID,doc="set of cs")
mp.J = Set(initialize=vehicleID,doc="set of vehicles")
mp.I = Set(initialize=platoonID,doc="set of platoons")
mp.swithO = Set(initialize=csIDwithO,doc="s union O") 
mp.swithOD = Set(initialize=csIDwithOD,doc="s union OD") 


# Variables
mp.y = Var(mp.I, mp.J, mp.swithO, within=Binary, doc = "vehicle $j$ joins platoon $i$ at station $s$")
mp.n = Var(mp.I, mp.swithO, within = NonNegativeIntegers, doc = "number of vehicle in platoon i at station s")
#mp.z = Var(mp.J, mp.S, within = Binary, doc = "indicator variable")
mp.w = Var(mp.I,mp.swithO, within =Binary, doc = "indicator variables")
mp.z = Var(mp.J, mp.S, initialize =0, within = Binary, doc = "indicator variable") 

# Parameters
mp.d = Param(mp.swithOD, mp.swithOD, initialize=distanceBetStation, doc = "distance between s1 and s2")
mp.m = Param(mp.J, initialize =efficiencyVehicle, doc = "energy efficiency of vehicle kWh/ mile")
#mp.v_j = Param(mp.J, initialize = 50, doc = "average speed of vehicle j")
#mp.v_i = Param(mp.I, initialize = 50, doc = "average speed of platoon i")
#mp.r_cha = Param(mp.S, initialize = 100, doc = "charging speed kw of station s")
#mp.t_que = Param(mp.J, mp.S, initialize = .25, doc = "waiting time hr before charge at s")
#mp.t_depmin = Param(mp.J, initialize = 0, doc = "earliest depart time of vehicle j")
#mp.t_arrmax = Param(mp.J, initialize = diffArrTime, doc = "latest arrive time hr of vehicle j")
#mp.soc_depmin = Param(mp.J, initialize = 1, doc = "minimum departure soc of vehicle j")
#mp.soc_arrmin = Param(mp.J, initialize = 0.2, doc = "minimum arrival soc of vehicle j")
#mp.soc_min = Param(mp.J, initialize = 0.2 , doc = "minimim soc allowed for vehicle j")
#mp.soc_max = Param(mp.J, initialize = 1,doc = "maximum soc allowed for vehicle j")
#mp.M = Param(initialize = 10, doc = "big number")
mp.delta = Param(mp.I, mp.swithO, initialize = .13, doc = "energy saving % of platoon i depart from station s")
#mp.c_cha = Param(mp.S, initialize =.35 , doc = "charging cost per kwh at station s per kWh https://freewiretech.com/difference-between-ev-charging-levels/")
#mp.c_delay = Param(mp.J, initialize = 50 , doc = "delay cost of vehicle j ")
mp.c_energy = Param(mp.J, initialize = 0.1165 , doc = "energy cost of vehicle j https://www.eia.gov/state/print.php?sid=FL")
#mp.b = Param(mp.J,initialize = batteryCapVehicle, doc = "battery capacity of vehicle j")
mp.N_max = Param(initialize = len(vehicleID), doc = "maximum vehicle in a platoon")
#mp.t_dep = Param(mp.J, initialize = diffDepTime, doc = "different departure time for vehicles")


# Additional MP variable
mp.zLower = Var(within=Reals)


energyCost  = sum(mp.c_energy[j]*mp.d[s,s+1]*mp.m[j]*(1-sum(mp.delta[i,s]*mp.y[i,j,s] for i in mp.I)) for j in mp.J for s in mp.swithO)

# Define constraints                                  
def cons1Rule(mp):
    return mp.zLower>=energyCost
mp.cons1=Constraint(rule=cons1Rule,doc="")


# Energy
def joinPlatoonRule1(mp,j,s):
	return sum(mp.y[i,j,s] for i in mp.I) <= 1
mp.joinPlatoon1= Constraint(mp.J, mp.swithO, rule = joinPlatoonRule1, doc = "vehicle j join at most 1 platoon at cs s")


def joinPlatoonRule2(mp,i,s):
    return sum (mp.y[i,j,s] for j in mp.J) >= 0
mp.joinPlatoon2 = Constraint(mp.I, mp.swithO, rule = joinPlatoonRule2, doc = "platoon i takes at most 1 vehicle at cs s")

def vehicleNoinPlatoonRule(mp,i,s):
    return mp.n[i,s] == sum(mp.y[i,j,s] for j in mp.J) 
mp.vehicleNoinPlatoon = Constraint(mp.I,mp.swithO,rule=vehicleNoinPlatoonRule)

def platoonNoLimitRule1(mp,i,s):
	return mp.n[i,s] <= mp.w[i,s]*mp.N_max
mp.platoonNoLimit1 = Constraint(mp.I, mp.swithO, rule = platoonNoLimitRule1)

def platoonNoLimitRule2(mp,i,s):
	return mp.n[i,s] >= mp.w[i,s]*2
mp.platoonNoLimit2 = Constraint(mp.I, mp.swithO, rule = platoonNoLimitRule2)

                                       
#define objective function
def objectiveRule(mp):
    return mp.zLower
        
mp.objective=Objective(rule=objectiveRule,sense=minimize,doc="objective function")

solvername = 'cplex'
solver=SolverFactory(solvername)
resultsMP = solver.solve(mp)
resultsMP.write()
    
import math
LB = value(mp.objective)
UB = math.inf

iteration = 1
#while (UB-LB >= 0.05):
while iteration <= 1:   
    
    print(f'LB = {LB}, UB = {UB}')
    
    # Primal
    sp=ConcreteModel()
    
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
    sp.tau = Var(sp.I,sp.swithOD, within = NonNegativeReals, doc="time of platoon i arrive/depart at/from station s")
    sp.soc_arr = Var(sp.J,sp.swithOD, within = NonNegativeReals, doc = "soc of vehicle j at station s when arrive")
    sp.soc_dep = Var(sp.J,sp.swithOD, within = NonNegativeReals, doc = "soc of vehicle j at station s when departs")
    sp.t_cha = Var(sp.J,sp.S, within = NonNegativeReals, doc = "charging time of vehicle j at station s")
    sp.t_wai = Var(sp.J,sp.S, within = NonNegativeReals, doc = "waiting time (after charge) of vehicle j at station s")

    #sp.delta = Var(sp.I, sp.swithO, doc = "energy saving % of platoon i depart from station s")
    sp.dummy = Var(sp.J,within = Reals, doc = "dummy var for delay cost calc")
    
    # Parameters
    sp.d = Param(sp.swithOD, sp.swithOD, initialize=distanceBetStation, doc = "distance between s1 and s2")
    sp.m = Param(sp.J, initialize =efficiencyVehicle, doc = "energy efficiency of vehicle kWh/ mile")
    sp.v_j = Param(sp.J, initialize = 50, doc = "average speed of vehicle j")
    sp.v_i = Param(sp.I, initialize = 50, doc = "average speed of platoon i")
    sp.r_cha = Param(sp.S, initialize = 100, doc = "charging speed kw of station s")
    sp.t_que = Param(sp.J, sp.S, initialize = .25, doc = "waiting time hr before charge at s")
    sp.t_depmin = Param(sp.J, initialize = 0, doc = "earliest depart time of vehicle j")
    sp.t_arrmax = Param(sp.J, initialize = diffArrTime, doc = "latest arrive time hr of vehicle j")
    sp.soc_depmin = Param(sp.J, initialize = 1, doc = "minimum departure soc of vehicle j")
    sp.soc_arrmin = Param(sp.J, initialize = 0.2, doc = "minimum arrival soc of vehicle j")
    sp.soc_min = Param(sp.J, initialize = 0.2 , doc = "minimim soc allowed for vehicle j")
    sp.soc_max = Param(sp.J, initialize = 1,doc = "maximum soc allowed for vehicle j")
    sp.M = Param(initialize = 10, doc = "big number")
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
    	return sp.soc_dep[j,0] ==  sp.soc_depmin[j]
    sp.initialSOC = Constraint(sp.J, rule = initialDepartSOCRule, doc = "init depart SOC = depat soc from origin")

    def lastArrivalSOCRule(sp,j):
    	return sp.soc_arr[j,list(sp.swithOD.value)[-1]] >= sp.soc_arrmin[j] 
    sp.lastArrivalSOC = Constraint(sp.J, rule = lastArrivalSOCRule, doc = "last arrival SOC = arrival minimum soc")

    def departSOCfromCSRule(sp,j,s):
    	return sp.soc_dep[j,s] == sp.soc_arr[j,s] + (sp.r_cha[s]/sp.b[j])*sp.t_cha[j,s] 
    sp.departSOCfromCS = Constraint(sp.J, sp.S, rule = departSOCfromCSRule)

    def arriveSOCatCSRule(sp,j,s):
        return sp.soc_arr[j,s+1] == sp.soc_dep[j,s]- ((sp.d[s,s+1]*sp.m[j])/sp.b[j])*(1 - sum(sp.delta[i,s]*mp.y[i,j,s].value for i in sp.I))

    sp.arriveSOCatCS = Constraint(sp.J,sp.swithO, rule = arriveSOCatCSRule)

    def arriveSOCMinRule(sp,j,s):
    	return sp.soc_arr[j,s] >= sp.soc_min[j] 
    sp.arriveSOCMin = Constraint(sp.J, sp.S, rule = arriveSOCMinRule)

    def departSOCMaxRule(sp,j,s):
    	return sp.soc_dep[j,s] <= sp.soc_max[j] #
    sp.departSOCMax = Constraint(sp.J, sp.S, rule = departSOCMaxRule)

    # Time
    # ra 11/4
    def initialDepartTimeRule(sp,j):
        return sp.t_jdep[j,0] >= sp.t_depmin[j]
    	#return sp.t_jdep[j,0] == sp.t_dep[j]
    sp.initialDepartTime = Constraint(sp.J, rule = initialDepartTimeRule)
    #sp.initialDepartTime.pprint()

    def departTimeRule(sp,j,s):
        return sp.t_jdep[j,s] == sp.t_jarr[j,s] + sp.t_que[j,s]*mp.z[j,s].value + sp.t_cha[j,s] + sp.t_wai[j,s]
    sp.departTime = Constraint(sp.J, sp.S, rule = departTimeRule)

    def relChaAndQueRule1(sp,j,s):
        return -sp.t_cha[j,s] >=-sp.M*mp.z[j,s].value # FIND DUAL
    sp.relChaAndQueRel1 = Constraint(sp.J,sp.S, rule = relChaAndQueRule1)

    def relChaAndQueRule2(sp,j,s):
        return sp.t_cha[j,s] >= 0 
    sp.relChaAndQueRel2 = Constraint(sp.J,sp.S, rule = relChaAndQueRule2)

    def arriveTimeRule(sp,j,s):
    	return sp.t_jarr[j,s+1] == sp.t_jdep[j,s] + sp.d[s,s+1]/sp.v_j[j]
    sp.arriveTime = Constraint(sp.J, sp.swithO, rule = arriveTimeRule)

    def platoonDepartRule(sp,i,s):
        return sp.t_idep[i,s] == sp.t_iarr[i,s]
    sp.platoonDepart = Constraint(sp.I,sp.S,rule=platoonDepartRule,doc="platton arrive time and depart time same")

    def platoonArriveRule(sp,i,s):
        return sp.t_iarr[i,s+1] == sp.t_idep[i,s] + sp.d[s,s+1]/sp.v_i[i]
    sp.platoonArriveTime = Constraint(sp.I, sp.swithO, rule = platoonArriveRule)

    def platoonVehicleDepRule1(sp,i,j,s):
    	return -sp.M*(1-mp.y[i,j,s].value) <= sp.t_jdep[j,s] - sp.t_idep[i,s]  # FIND DUAL
    sp.platoonVehicleDep1 = Constraint(sp.I, sp.J, sp.swithO, rule = platoonVehicleDepRule1)

    def platoonVehicleDepRule2(sp,i,j,s):
    	return  sp.t_jdep[j,s] - sp.t_idep[i,s] <= sp.M*(1-mp.y[i,j,s].value) # FIND DUAL
    sp.platoonVehicleDep2 = Constraint(sp.I, sp.J, sp.swithO, rule = platoonVehicleDepRule2)

    def nonNegativeDelayRule1(sp,j):
        return sp.dummy[j] >= sp.t_jarr[j,list(sp.swithOD.value)[-1]]-sp.t_arrmax[j] 
    sp.nonNegativeDelay1 = Constraint(sp.J,rule = nonNegativeDelayRule1)

    def nonNegativeDelayRule2(sp,j):
        return sp.dummy[j] >= 0 
    sp.nonNegativeDelay2 = Constraint(sp.J,rule = nonNegativeDelayRule2)

                                               
    #define objective function
    def objectiveRuleSp1(sp):
        chargingCost = sum(sp.c_cha[s]*sp.r_cha[s]*sp.t_cha[j,s] for j in sp.J for s in sp.S)
        delayCost = sum(sp.c_delay[j]*sp.dummy[j] for j in sp.J)
        
        return (chargingCost+delayCost)
            
    sp.objective=Objective(rule=objectiveRuleSp1,sense=minimize,doc="objective function")
    sp.dual = Suffix(direction = Suffix.IMPORT)
    resultsSp = solver.solve(sp)
    
    dual1 = {}
    dual2 = {}
    dual3 = {}

    if (resultsSp.solver.status == SolverStatus.ok) and (resultsSp.solver.termination_condition == TerminationCondition.infeasible):
        newSp1 = ConcreteModel()
        
        # Sets
        newSp1.S = Set(initialize=csID,doc="set of cs")
        newSp1.J = Set(initialize=vehicleID,doc="set of vehicles")
        newSp1.I = Set(initialize=platoonID,doc="set of platoons")
        newSp1.swithO = Set(initialize=csIDwithO,doc="s union O") 
        newSp1.swithOD = Set(initialize=csIDwithOD,doc="s union OD") 
        
        # Define variables
        newSp1.sChaQue= Var(newSp1.J, newSp1.S, within=NonNegativeReals) 
        newSp1.sPlatVehDep = Var(newSp1.I, newSp1.J, newSp1.swithO, within = NonNegativeReals)
        newSp1.sPlatVehDep2 = Var(newSp1.I, newSp1.J, newSp1.swithO, within = NonNegativeReals)

        newSp1.t_jdep = Var(newSp1.J, newSp1.swithOD, within=NonNegativeReals, doc = "time of vehicle j depart from station s")
        newSp1.t_jarr = Var(newSp1.J, newSp1.swithOD, within = NonNegativeReals, doc = "time of vehicle j arrive at station s")
        newSp1.t_idep = Var(newSp1.I, newSp1.swithOD, within = NonNegativeReals, doc = "time of platoon i depart from station s") #EMPTY
        newSp1.t_iarr = Var(newSp1.I, newSp1.swithOD, within = NonNegativeReals, doc = "time of platoon i arrive at station s") #EMPTY
        newSp1.tau = Var(newSp1.I,newSp1.swithOD, within = NonNegativeReals, doc="time of platoon i arrive/depart at/from station s")
        newSp1.soc_arr = Var(newSp1.J,newSp1.swithOD, within = NonNegativeReals, doc = "soc of vehicle j at station s when arrive")
        newSp1.soc_dep = Var(newSp1.J,newSp1.swithOD, within = NonNegativeReals, doc = "soc of vehicle j at station s when departs")
        newSp1.t_cha = Var(newSp1.J,newSp1.S, within = NonNegativeReals, doc = "charging time of vehicle j at station s")
        newSp1.t_wai = Var(newSp1.J,newSp1.S, within = NonNegativeReals, doc = "waiting time (after charge) of vehicle j at station s")
        newSp1.z = Var(mp.J, mp.S, within = Binary, doc = "indicator variable") # THIS, WHEN IN MP, GIVES NONE. NONE CANNOT MULTIPLY WITH FLOAT IN CONSTRAINT ERROE SHOWS UP
        newSp1.dummy = Var(sp.J,within = Reals, doc = "dummy var for delay cost calc")
        
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
        newSp1.M = Param(initialize = 10, doc = "big number")
        newSp1.delta = Param(newSp1.I, newSp1.swithO, initialize = .13, doc = "energy saving % of platoon i depart from station s")
        newSp1.c_cha = Param(newSp1.S, initialize =.35 , doc = "charging cost per kwh at station s per kWh https://freewiretech.com/difference-between-ev-charging-levels/")
        newSp1.c_delay = Param(newSp1.J, initialize = 50 , doc = "delay cost of vehicle j ")
        newSp1.c_energy = Param(newSp1.J, initialize = 0.1165 , doc = "energy cost of vehicle j https://www.eia.gov/state/print.php?sid=FL")
        newSp1.b = Param(newSp1.J,initialize = batteryCapVehicle, doc = "battery capacity of vehicle j")
        newSp1.N_max = Param(initialize = len(vehicleID), doc = "maximum vehicle in a platoon")
        newSp1.t_dep = Param(newSp1.J, initialize = diffDepTime, doc = "different departure time for vehicles")

        # Constraints
        def lastArrivalSOCRule(newSp1,j):
        	return newSp1.soc_arr[j,list(newSp1.swithOD.value)[-1]] >= newSp1.soc_arrmin[j] 
        newSp1.lastArrivalSOC = Constraint(newSp1.J, rule = lastArrivalSOCRule, doc = "last arrival SOC = arrival minimum soc")

        
        
        def arriveSOCMinRule(newSp1,j,s):
        	return newSp1.soc_arr[j,s] >= newSp1.soc_min[j]
        newSp1.arriveSOCMin = Constraint(newSp1.J, newSp1.S, rule = arriveSOCMinRule)

        def departSOCMaxRule(newSp1,j,s):
        	return newSp1.soc_dep[j,s] <= newSp1.soc_max[j]
        newSp1.departSOCMax = Constraint(newSp1.J, newSp1.S, rule = departSOCMaxRule)

        # Time
        # ra 11/4
        def initialDepartTimeRule(newSp1,j):
            return newSp1.t_jdep[j,0]  >= newSp1.t_depmin[j]  
        	#return newSp1.t_jdep[j,0] == newSp1.t_dep[j]
        newSp1.initialDepartTime = Constraint(newSp1.J, rule = initialDepartTimeRule)
        #newSp1.initialDepartTime.pprint()

       
        def relChaAndQueRule1(newSp1,j,s):
            return -newSp1.t_cha[j,s] + newSp1.sChaQue[j,s] >=-(newSp1.M*newSp1.z[j,s]) 
        newSp1.relChaAndQueRel1 = Constraint(newSp1.J,newSp1.S, rule = relChaAndQueRule1)

        def relChaAndQueRule2(newSp1,j,s):
            return newSp1.t_cha[j,s]  >= 0 
        newSp1.relChaAndQueRel2 = Constraint(newSp1.J,newSp1.S, rule = relChaAndQueRule2)

        
        def platoonVehicleDepRule1(newSp1,i,j,s):
        	#return -newSp1.M*(1-mp.y[i,j,s].value) + newSp1.sPlatVehDep[i,j,s] <= newSp1.t_jdep[j,s] - newSp1.t_idep[i,s] 
            return newSp1.t_jdep[j,s]-newSp1.t_idep[i,s]+newSp1.sPlatVehDep[i,j,s] >= -newSp1.M*(1-mp.y[i,j,s].value)
        newSp1.platoonVehicleDep1 = Constraint(newSp1.I, newSp1.J, newSp1.swithO, rule = platoonVehicleDepRule1)

        def platoonVehicleDepRule2(newSp1,i,j,s):
        	#return  newSp1.t_jdep[j,s] - newSp1.t_idep[i,s] + newSp1.sPlatVehDep2[i,j,s]<= newSp1.M*(1-mp.y[i,j,s].value) 
            return  newSp1.t_idep[i,s]-newSp1.t_jdep[j,s]  + newSp1.sPlatVehDep2[i,j,s]>= -(newSp1.M*(1-mp.y[i,j,s].value))# FIND DUAL
        newSp1.platoonVehicleDep2 = Constraint(newSp1.I, newSp1.J, newSp1.swithO, rule = platoonVehicleDepRule2)

        def nonNegativeDelayRule1(newSp1,j):
            return newSp1.dummy[j]  >= newSp1.t_jarr[j,list(newSp1.swithOD.value)[-1]]-newSp1.t_arrmax[j] 
        newSp1.nonNegativeDelay1 = Constraint(newSp1.J,rule = nonNegativeDelayRule1)

        def nonNegativeDelayRule2(newSp1,j):
            return newSp1.dummy[j] >= 0
        newSp1.nonNegativeDelay2 = Constraint(newSp1.J,rule = nonNegativeDelayRule2)

        # Objective
        def objectiveRuleNewSp1(newSp1):
            slackSum1 = sum(newSp1.sChaQue[j,s] for j in newSp1.J for s in newSp1.S)
            slackSum2 = sum(newSp1.sPlatVehDep[i,j,s] for i in newSp1.I for j in newSp1.J for s in newSp1.S)
            slackSum3 = sum(newSp1.sPlatVehDep[i,j,s] for i in newSp1.I for j in newSp1.J for s in newSp1.S)
            return slackSum1+slackSum2+slackSum3
        newSp1.objective = Objective(rule = objectiveRuleNewSp1, sense = minimize)
        
        newSp1.dual = Suffix(direction = Suffix.IMPORT)
        resultsNewSp1 = solver.solve(newSp1)
        
        for j in newSp1.J:
            for s in newSp1.S:
                dual1 [j,s] = newSp1.dual[newSp1.relChaAndQueRel1[j,s]] 
        for i in newSp1.I:
            for j in newSp1.J:
                for s in newSp1.swithO:
                    dual2 [i,j,s] = newSp1.dual[newSp1.platoonVehicleDep1[i,j,s]]
                    dual3 [i,j,s] = newSp1.dual[newSp1.platoonVehicleDep2[i,j,s]]
         
        def infeasibilityCutRule(mp):
             return -(newSp1.M*newSp1.z[j,s])*dual1[j,s]-\
                 newSp1.M*(1-mp.y[i,j,s].value)*dual2 [i,j,s]-(newSp1.M*(1-mp.y[i,j,s].value))*dual3[i,j,s] <= 0
        mp.infesibility=Constraint(rule=cinfeasibilityCutRule,doc="")
        
        
    if (resultsSp1.solver.status == SolverStatus.ok) and (resultsSp1.solver.termination_condition == TerminationCondition.optimal):
        energyCostWithYcap  = sum(mp.c_energy[j]*mp.d[s,s+1]*mp.m[j]*(1-sum(mp.delta[i,s]*mp.y[i,j,s].value for i in mp.I)) for j in mp.J for s in mp.swithO)
       
        UB = min(UB, (value(sp.objective)+energyCostWithYcap))
        
        #print(UB)
        for j in sp.j:
            for s in sp.s:
                dual1[j,s] = sp.dual[sp.relChaAndQue1[j,s]]
        for i in sp1.I:
            for j in sp1.J:
                for s in sp1.swithO:
                    dual2 [i,j,s] = sp.dual[sp.platoonVehicleDep1[i,j,s]]
                    dual3 [i,j,s] = sp.dual[sp.platoonVehicleDep2[i,j,s]]
        
        def feasibilityCutRule(mp):
            return mp.zLower>=sum(mp.c_energy[j]*mp.d[s,s+1]*mp.m[j]*(1-sum(mp.delta[i,s]*mp.y[i,j,s].value for i in mp.I)) for j in mp.J for s in mp.swithO) +\
                sum(sp.M*mp.z[j,s]*dual1[j,s] for j in mp.J for s in mp.S)+sum ((-sp.M)*(1-mp.y[i,j,s]) for i in mp.I for j in mp.J for s in mp.swithO)+\
                sum(sp.M*(1-mp.y[i,j,s]) for i in mp.I for j in mp.J for s in mp.swithO)
        mp.feasibility=Constraint(rule=feasibilityCutRule,doc="")
    
    
    resultsMp=solver.solve(mp)
    LB = value(mp.objective)
    print(f'y = {mp.y.value}')
    
    iteration += 1
    
    
    
