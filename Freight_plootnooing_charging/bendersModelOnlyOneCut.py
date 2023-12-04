
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 23 22:36:47 2022

@author: rakibulalam
"""

# z is ommitted

# Library
import time
tic = time.perf_counter()

from pyomo.environ import *
from pyomo.opt import SolverFactory, SolverStatus, TerminationCondition
import random

#======================= Data preparation ======================
routeLength = 500# Mile 
numberOfCS = 20
csID = range(1,numberOfCS+1)
csIDwithO = range(0, numberOfCS+1)
csIDwithOD = range(0, numberOfCS+2) 
vehicleID = range(0,50)
platoonID = range(0,50)

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
batteryCapVehicle = 300 #{k:random.randrange(220,400) for k in vehicleID}


#different departure time for differnt vehicles
diffDepTime=dict()
random.seed(a=1)
for i in vehicleID:
    diffDepTime[i]=random.randint(0, 4)
    
diffArrTime=dict()
for i in vehicleID:
    diffArrTime[i]=diffDepTime[i]+11.16 
    

#======================= Genearal information solver and model =================
solvername = 'cplex'
solver=SolverFactory(solvername)


def delConsSubProb(): # we may need later in benders decomposition
    sp.del_component(sp.arriveSOCatCS_ineq1)
    sp.del_component(sp.arriveSOCatCS_ineq1_index)
    sp.del_component(sp.platoonVehicleDep1)
    sp.del_component(sp.platoonVehicleDep2)
    sp.del_component(sp.platoonVehicleDep1_index)
    sp.del_component(sp.platoonVehicleDep2_index)
    sp.del_component(sp.dual)

class Models:    
    def masterProblem(self):
    
        # Sets
        mp.S = Set(initialize=csID,doc="set of cs")
        mp.J = Set(initialize=vehicleID,doc="set of vehicles")
        mp.I = Set(initialize=platoonID,doc="set of platoons")
        mp.swithO = Set(initialize=csIDwithO,doc="s union O") 
        mp.swithOD = Set(initialize=csIDwithOD,doc="s union OD") 
    
        # Variables
        mp.y = Var(mp.I, mp.J, mp.swithO, within=Binary, doc = "vehicle $j$ joins platoon $i$ at station $s$")
        mp.n = Var(mp.I, mp.swithO, within = NonNegativeIntegers, doc = "number of vehicle in platoon i at station s")
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
        
        
        # Energy
        def joinPlatoonRule1(mp,j,s):
        	return -sum(mp.y[i,j,s] for i in mp.I) >= -1
        mp.joinPlatoon1= Constraint(mp.J, mp.swithO, rule = joinPlatoonRule1, doc = "vehicle j join at most 1 platoon at cs s")
        
        def joinPlatoonRule2(mp,i,s):
            return sum (mp.y[i,j,s] for j in mp.J) >= 0
        mp.joinPlatoon2 = Constraint(mp.I, mp.swithO, rule = joinPlatoonRule2, doc = "platoon i takes at most 1 vehicle at cs s")
        
        def vehicleNoinPlatoonRule(mp,i,s):
            return mp.n[i,s] == sum(mp.y[i,j,s] for j in mp.J) 
        mp.vehicleNoinPlatoon = Constraint(mp.I,mp.swithO,rule=vehicleNoinPlatoonRule)
        
        def platoonNoLimitRule1(mp,i,s):
        	return -mp.n[i,s] >= -mp.w[i,s]*mp.N_max
        mp.platoonNoLimit1 = Constraint(mp.I, mp.swithO, rule = platoonNoLimitRule1)
        
        def platoonNoLimitRule2(mp,i,s):
        	return mp.n[i,s] >= mp.w[i,s]*2
        mp.platoonNoLimit2 = Constraint(mp.I, mp.swithO, rule = platoonNoLimitRule2)
        
        mp.cuts = ConstraintList()                                       
        
        #define objective function
        def objectiveRule(mp):
            return mp.zLower
                
        mp.objective=Objective(rule=objectiveRule,sense=minimize,doc="objective function")
        
        global resultsMp
        resultsMp = solver.solve(mp)
        
        return mp
        
    
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
        sp.tau = Var(sp.I,sp.swithOD, within = NonNegativeReals, doc="time of platoon i arrive/depart at/from station s")
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
            return sp.soc_arr[j,s+1]-sp.soc_dep[j,s] == -((sp.d[s,s+1]*sp.m[j])/sp.b[j])*(1 - sum(sp.delta[i,s]*mp.y[i,j,s].value for i in sp.I))
    
        sp.arriveSOCatCS_ineq1 = Constraint(sp.J,sp.swithO, rule = arriveSOCatCSRule_ineq1)
        
# =============================================================================
#         def arriveSOCatCSRule_ineq2(sp,j,s):
#             return -sp.soc_arr[j,s+1]-sp.soc_dep[j,s] >= ((sp.d[s,s+1]*sp.m[j])/sp.b[j])*(1 - sum(sp.delta[i,s]*mp.y[i,j,s].data() for i in sp.I))
#     
#         sp.arriveSOCatCS_ineq2 = Constraint(sp.J,sp.swithO, rule = arriveSOCatCSRule_ineq2)
# =============================================================================
    
        def arriveSOCMinRule(sp,j,s):
        	return sp.soc_arr[j,s] >= sp.soc_min[j]
        sp.arriveSOCMin = Constraint(sp.J, sp.S, rule = arriveSOCMinRule)
    
        def departSOCMaxRule(sp,j,s):
        	return -sp.soc_dep[j,s] >= -sp.soc_max[j] 
        sp.departSOCMax = Constraint(sp.J, sp.S, rule = departSOCMaxRule)
        
        # Time
        def initialDepartTimeRule(sp,j):
            return sp.t_jdep[j,0]>=sp.t_depmin[j] 
        	#return sp.t_jdep[j,0] == sp.t_dep[j]
        sp.initialDepartTime = Constraint(sp.J, rule = initialDepartTimeRule)
    
        def departTimeRule_ineq1(sp,j,s):
            #return sp.t_jdep[j,s]-sp.t_jarr[j,s]-sp.t_cha[j,s]- sp.t_wai[j,s] == sp.t_que[j,s]*mp.z[j,s].data()
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
        	return sp.t_jdep[j,s] - sp.t_idep[i,s] >= -sp.M*(1-mp.y[i,j,s].value)  # FIND DUAL
        sp.platoonVehicleDep1 = Constraint(sp.I, sp.J, sp.swithO, rule = platoonVehicleDepRule1)
    
        def platoonVehicleDepRule2(sp,i,j,s):
        	return  -sp.t_jdep[j,s] + sp.t_idep[i,s] >= -sp.M*(1-mp.y[i,j,s].value) # FIND DUAL
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
    
    
    
# Bender's decomposition algortihm
    
import math
models = Models()
mp = ConcreteModel()
models.masterProblem()

LB = value(mp.objective)
print(f'Initial LB = {LB}')
UB = math.inf
print(f'Initial UB = {UB}')

while (UB-LB >= .001):   
    dual1 = {}
    dual2 = {}
    dual3 = {}
    dual4_1 = {}
    #dual5_1 = {}
    
    # newly added
    dual_iniSOC = {}
    dual_lastArrivalSOC = {}
    dual_nonNegativeDelay1 = {}
    dual_initialDepartTime = {}
    dual_arriveTime = {}
    dual_arriveSOCMin = {}
    dual_departSOCMax = {}
    dual_platoonArriveTime = {}
    
    
    sp = ConcreteModel()
    models.subProblem()
   
    
    if (resultsSp.solver.status == SolverStatus.ok) and (resultsSp.solver.termination_condition == TerminationCondition.infeasible):
        
        # Add slack variables
        sp.sArriveSOCatCS_ineq1 = Var(sp.J, sp.swithO, within = NonNegativeReals, doc = "slack for equality with ycap")
        sp.sPlatVehDep1 = Var(sp.I, sp.J, sp.swithO, within = NonNegativeReals)
        sp.sPlatVehDep2 = Var(sp.I, sp.J, sp.swithO, within = NonNegativeReals)
     
        delConsSubProb()
    
        def arriveSOCatCSRule_ineq1(sp,j,s):#slack
            return sp.soc_arr[j,s+1]-sp.soc_dep[j,s]+sp.sArriveSOCatCS_ineq1[j,s] == -((sp.d[s,s+1]*sp.m[j])/sp.b[j])*(1 - sum(sp.delta[i,s]*mp.y[i,j,s].value for i in sp.I))
        sp.arriveSOCatCS_ineq1 = Constraint(sp.J,sp.swithO, rule = arriveSOCatCSRule_ineq1)
    
        def platoonVehicleDepRule1(sp,i,j,s):
            return sp.t_jdep[j,s]-sp.t_idep[i,s]+sp.sPlatVehDep1[i,j,s]>= -sp.M*(1-mp.y[i,j,s].value)
        sp.platoonVehicleDep1 = Constraint(sp.I, sp.J, sp.swithO, rule = platoonVehicleDepRule1)
    
        def platoonVehicleDepRule2(sp,i,j,s):
            return sp.t_idep[i,s]-sp.t_jdep[j,s]+sp.sPlatVehDep2[i,j,s]>= -sp.M*(1-mp.y[i,j,s].value)                                                  
        sp.platoonVehicleDep2 = Constraint(sp.I, sp.J, sp.swithO, rule = platoonVehicleDepRule2)
        
        sp.del_component(sp.objective)
        
        chargingCost = sum(sp.c_cha[s]*sp.r_cha[s]*sp.t_cha[j,s] for j in sp.J for s in sp.S)
        delayCost = sum(sp.c_delay[j]*sp.dummy[j] for j in sp.J)
        
        def objectiveRuleSp1(sp):
           slackSum1 = sum(sp.sArriveSOCatCS_ineq1[j,s] for j in sp.J for s in sp.swithO)
           slackSum6 = sum(sp.sPlatVehDep1[i,j,s] for i in sp.I for j in sp.J for s in sp.swithO)
           slackSum7 = sum(sp.sPlatVehDep2[i,j,s] for i in sp.I for j in sp.J for s in sp.swithO)
           
           return chargingCost+delayCost+1000000*(slackSum1+slackSum6+slackSum7)
               
        sp.objective=Objective(rule=objectiveRuleSp1,sense=minimize,doc="objective function")
        sp.dual = Suffix(direction = Suffix.IMPORT)
        
        resultsSp=solver.solve(sp)
        
    if (resultsSp.solver.status == SolverStatus.ok) and (resultsSp.solver.termination_condition == TerminationCondition.optimal):
        energyCostWithYcap  = sum(sp.c_energy[j]*sp.d[s,s+1]*sp.m[j]*(1-sum(sp.delta[i,s]*mp.y[i,j,s].value for i in sp.I)) for j in sp.J for s in sp.swithO)
        UB = min(UB, (value(sp.objective)+energyCostWithYcap))
        #UB = min(UB, (data()(sp.objective)+sum(mp.y[i,j,s].data() for i in mp.I for j in mp.J for s in mp.swithO)))
        print(f'UB = {UB}')
        for j in sp.J:
            for s in sp.S:
                dual_arriveSOCMin[j,s] = sp.dual[sp.arriveSOCMin[j,s]]#zero in final iter
                dual_departSOCMax[j,s] = sp.dual[sp.departSOCMax[j,s]]#zero in final iter
                #dual5_1[j,s] = sp.dual[sp.departTime_ineq1[j,s]]  
            
        for i in sp.I:
            for j in sp.J:
                for s in sp.swithO:
                    dual2[i,j,s] = sp.dual[sp.platoonVehicleDep1[i,j,s]]#zero in final iter
                    dual3[i,j,s] = sp.dual[sp.platoonVehicleDep2[i,j,s]]#zero in final iter
                    
        for j in sp.J:
            for s in sp.swithO:
                dual4_1[j,s] = sp.dual[sp.arriveSOCatCS_ineq1[j,s]]#-105 in final iter
                dual_arriveTime[j,s] = sp.dual[sp.arriveTime[j,s]]#zero in final iter

        for j in sp.J:
            dual_iniSOC[j] = sp.dual[sp.initialSOC[j]]#-105 in final iter
            dual_lastArrivalSOC[j] = sp.dual[sp.lastArrivalSOC[j]]#105 in final iter
            dual_nonNegativeDelay1[j] = sp.dual[sp.nonNegativeDelay1[j]]#zero in final iter
            dual_initialDepartTime[j] = sp.dual[sp.initialDepartTime[j]]#zero in final iter

        for i in sp.I:
            for s in sp.swithO:
            	dual_platoonArriveTime[i,s] = sp.dual[sp.platoonArriveTime[i,s]]#zero in final iter
        
        def optimalityCutRule(mp):
            return mp.zLower>=sum(mp.c_energy[j]*mp.d[s,s+1]*mp.m[j]*(1-sum(mp.delta[i,s]*mp.y[i,j,s] for i in mp.I)) for j in mp.J for s in mp.swithO)\
                +sum(-sp.M*(1-mp.y[i,j,s])*dual2[i,j,s] for i in sp.I for j in sp.J for s in sp.swithO)\
                +sum(-(sp.M*(1-mp.y[i,j,s]))*dual3[i,j,s] for i in sp.I for j in sp.J for s in sp.swithO)\
                +sum(-((sp.d[s,s+1]*sp.m[j])/sp.b[j])*(1 - sum(sp.delta[i,s]*mp.y[i,j,s] for i in sp.I))*dual4_1[j,s] for j in sp.J for s in sp.swithO)\
                +sum(sp.soc_arrmin[j]*dual_lastArrivalSOC[j] for j in sp.J)\
                +sum(sp.soc_depmin[j]*dual_iniSOC[j] for j in sp.J)\
                +sum((-sp.t_arrmax[j])*dual_nonNegativeDelay1[j] for j in sp.J)\
                +sum(sp.t_depmin[j]*dual_initialDepartTime[j] for j in sp.J)\
                +sum((sp.d[s,s+1]/sp.v_j[j])*dual_arriveTime[j,s] for j in sp.J for s in sp.swithO)\
                +sum(sp.soc_min[j]*dual_arriveSOCMin[j,s] for j in sp.J for s in sp.S)\
                +sum((-sp.soc_max[j])*dual_departSOCMax[j,s] for j in sp.J for s in sp.S)\
                +sum((sp.d[s,s+1]/sp.v_i[i])*dual_platoonArriveTime[i,s] for i in sp.I for s in sp.swithO)
        #+sum(0*dual5_1[j,s] for j in sp.J for s in sp.S)\
        mp.cuts.add(optimalityCutRule(mp))
        
    
    resultsMp=solver.solve(mp)
    LB = value(mp.objective)
    print(f'LB = {LB}\n')
    
    
toc = time.perf_counter()
print(f'Processing time = {toc-tic}\n')
    






































# =============================================================================
# 
# 
# # original problem
# # Model
# originalModel = ConcreteModel()
# 
# # Sets
# originalModel.S = Set(initialize=csID,doc="set of cs")
# originalModel.J = Set(initialize=vehicleID,doc="set of vehicles")
# originalModel.I = Set(initialize=platoonID,doc="set of platoons")
# originalModel.swithO = Set(initialize=csIDwithO,doc="s union O") 
# originalModel.swithOD = Set(initialize=csIDwithOD,doc="s union OD") 
# 
# 
# # Variables
# originalModel.t_jdep = Var(originalModel.J, originalModel.swithOD, within=NonNegativeReals, doc = "time of vehicle j depart from station s")
# originalModel.t_jarr = Var(originalModel.J, originalModel.swithOD, within = NonNegativeReals, doc = "time of vehicle j arrive at station s")
# originalModel.t_idep = Var(originalModel.I, originalModel.swithOD, within = NonNegativeReals, doc = "time of platoon i depart from station s") #EMPTY
# originalModel.t_iarr = Var(originalModel.I, originalModel.swithOD, within = NonNegativeReals, doc = "time of platoon i arrive at station s") #EMPTY
# originalModel.tau = Var(originalModel.I,originalModel.swithOD, within = NonNegativeReals, doc="time of platoon i arrive/depart at/from station s")
# originalModel.soc_arr = Var(originalModel.J,originalModel.swithOD, within = NonNegativeReals, doc = "soc of vehicle j at station s when arrive")
# originalModel.soc_dep = Var(originalModel.J,originalModel.swithOD, within = NonNegativeReals, doc = "soc of vehicle j at station s when departs")
# originalModel.t_cha = Var(originalModel.J,originalModel.S, within = NonNegativeReals, doc = "charging time of vehicle j at station s")
# originalModel.t_wai = Var(originalModel.J,originalModel.S, within = NonNegativeReals, doc = "waiting time (after charge) of vehicle j at station s")
# originalModel.dummy = Var(originalModel.J,within = Reals, doc = "dummy var for delay cost calc")
# 
# # Parameters
# originalModel.d = Param(originalModel.swithOD, originalModel.swithOD, initialize=distanceBetStation, doc = "distance between s1 and s2")
# originalModel.m = Param(originalModel.J, initialize =efficiencyVehicle, doc = "energy efficiency of vehicle kWh/ mile")
# originalModel.v_j = Param(originalModel.J, initialize = 50, doc = "average speed of vehicle j")
# originalModel.v_i = Param(originalModel.I, initialize = 50, doc = "average speed of platoon i")
# originalModel.r_cha = Param(originalModel.S, initialize = 100, doc = "charging speed kw of station s")
# originalModel.t_depmin = Param(originalModel.J, initialize = 0, doc = "earliest depart time of vehicle j")
# originalModel.t_arrmax = Param(originalModel.J, initialize = diffArrTime, doc = "latest arrive time hr of vehicle j")
# originalModel.soc_depmin = Param(originalModel.J, initialize = 1, doc = "minimum departure soc of vehicle j")
# originalModel.soc_arrmin = Param(originalModel.J, initialize = 0.2, doc = "minimum arrival soc of vehicle j")
# originalModel.soc_min = Param(originalModel.J, initialize = 0.2 , doc = "minimim soc allowed for vehicle j")
# originalModel.soc_max = Param(originalModel.J, initialize = 1,doc = "maximum soc allowed for vehicle j")
# originalModel.M = Param(initialize = 25, doc = "big number")
# originalModel.delta = Param(originalModel.I, originalModel.swithO, initialize = .13, doc = "energy saving % of platoon i depart from station s")
# originalModel.c_cha = Param(originalModel.S, initialize =.35 , doc = "charging cost per kwh at station s per kWh https://freewiretech.com/difference-between-ev-charging-levels/")
# originalModel.c_delay = Param(originalModel.J, initialize = 50 , doc = "delay cost of vehicle j ")
# originalModel.c_energy = Param(originalModel.J, initialize = 0.1165 , doc = "energy cost of vehicle j https://www.eia.gov/state/print.php?sid=FL")
# originalModel.b = Param(originalModel.J,initialize = batteryCapVehicle, doc = "battery capacity of vehicle j")
# originalModel.N_max = Param(initialize = len(vehicleID), doc = "maximum vehicle in a platoon")
# originalModel.t_dep = Param(originalModel.J, initialize = diffDepTime, doc = "different departure time for vehicles")
# 
# 
# # Constraints
# # SOC 
# def initialDepartSOCRule(originalModel,j):
# 	return originalModel.soc_dep[j,0] ==  originalModel.soc_depmin[j]
# originalModel.initialSOC = Constraint(originalModel.J, rule = initialDepartSOCRule, doc = "init depart SOC = depat soc from origin")
# 
# def lastArrivalSOCRule(originalModel,j):
# 	return originalModel.soc_arr[j,list(originalModel.swithOD.data())[-1]] >= originalModel.soc_arrmin[j] 
# originalModel.lastArrivalSOC = Constraint(originalModel.J, rule = lastArrivalSOCRule, doc = "last arrival SOC = arrival minimum soc")
# 
# def departSOCfromCSRule(originalModel,j,s):
# 	return originalModel.soc_dep[j,s] == originalModel.soc_arr[j,s] + (originalModel.r_cha[s]/originalModel.b[j])*originalModel.t_cha[j,s]
# originalModel.departSOCfromCS = Constraint(originalModel.J, originalModel.S, rule = departSOCfromCSRule)
# 
# def arriveSOCatCSRule(originalModel,j,s):
#     return originalModel.soc_arr[j,s+1] == originalModel.soc_dep[j,s]- ((originalModel.d[s,s+1]*originalModel.m[j])/originalModel.b[j])*(1 - sum(originalModel.delta[i,s]*mp.y[i,j,s].data() for i in originalModel.I))
# 
# originalModel.arriveSOCatCS = Constraint(originalModel.J,originalModel.swithO, rule = arriveSOCatCSRule)
# 
# def arriveSOCMinRule(originalModel,j,s):
# 	return originalModel.soc_arr[j,s] >= originalModel.soc_min[j]
# originalModel.arriveSOCMin = Constraint(originalModel.J, originalModel.S, rule = arriveSOCMinRule)
# 
# def departSOCMaxRule(originalModel,j,s):
# 	return originalModel.soc_dep[j,s] <= originalModel.soc_max[j]
# originalModel.departSOCMax = Constraint(originalModel.J, originalModel.S, rule = departSOCMaxRule)
# 
# # Time
# # ra 11/4
# def initialDepartTimeRule(originalModel,j):
#     return originalModel.t_jdep[j,0] >= originalModel.t_depmin[j]  
# originalModel.initialDepartTime = Constraint(originalModel.J, rule = initialDepartTimeRule)
# 
# def departTimeRule(originalModel,j,s):
#     return originalModel.t_jdep[j,s] == originalModel.t_jarr[j,s] + originalModel.t_cha[j,s] + originalModel.t_wai[j,s]
# originalModel.departTime = Constraint(originalModel.J, originalModel.S, rule = departTimeRule)
# 
# # =============================================================================
# # def relChaAndQueRule1(model,j,s):
# #     return model.t_cha[j,s] <=model.M*model.z[j,s]
# # model.relChaAndQueRel1 = Constraint(model.J,model.S, rule = relChaAndQueRule1)
# # 
# # =============================================================================
# # =============================================================================
# # def relChaAndQueRule2(model,j,s):
# #     return model.t_cha[j,s] >= 0
# # model.relChaAndQueRel2 = Constraint(model.J,model.S, rule = relChaAndQueRule2)
# # 
# # =============================================================================
# def arriveTimeRule(originalModel,j,s):
# 	return originalModel.t_jarr[j,s+1] == originalModel.t_jdep[j,s] + originalModel.d[s,s+1]/originalModel.v_j[j]
# originalModel.arriveTime = Constraint(originalModel.J, originalModel.swithO, rule = arriveTimeRule)
# 
# def platoonDepartRule(originalModel,i,s):
#     return originalModel.t_idep[i,s] == originalModel.t_iarr[i,s]
# originalModel.platoonDepart = Constraint(originalModel.I,originalModel.S,rule=platoonDepartRule,doc="platton arrive time and depart time same")
# 
# def platoonArriveRule(originalModel,i,s):
#     return originalModel.t_iarr[i,s+1] == originalModel.t_idep[i,s] + originalModel.d[s,s+1]/originalModel.v_i[i]
# originalModel.platoonArriveTime = Constraint(originalModel.I, originalModel.swithO, rule = platoonArriveRule)
# 
# def platoonVehicleDepRule1(originalModel,i,j,s):
# 	return -originalModel.M*(1-mp.y[i,j,s].data()) <= originalModel.t_jdep[j,s] - originalModel.t_idep[i,s] 
# originalModel.platoonVehicleDep1 = Constraint(originalModel.I, originalModel.J, originalModel.swithO, rule = platoonVehicleDepRule1)
# 
# def platoonVehicleDepRule2(originalModel,i,j,s):
# 	return  originalModel.t_jdep[j,s] - originalModel.t_idep[i,s] <= originalModel.M*(1-mp.y[i,j,s].data())
# originalModel.platoonVehicleDep2 = Constraint(originalModel.I, originalModel.J, originalModel.swithO, rule = platoonVehicleDepRule2)
# 
# # =============================================================================
# # # Energy
# # def joinPlatoonRule1(originalModel,j,s):
# # 	return sum(mp.y[i,j,s].data() for i in originalModel.I) <= 1
# # originalModel.joinPlatoon1= Constraint(originalModel.J, originalModel.swithO, rule = joinPlatoonRule1, doc = "vehicle j join at most 1 platoon at cs s")
# # 
# # =============================================================================
# 
# # =============================================================================
# # def joinPlatoonRule2(originalModel,i,s):
# #     return sum (mp.y[i,j,s].data() for j in originalModel.J) >= 0
# # originalModel.joinPlatoon2 = Constraint(originalModel.I, originalModel.swithO, rule = joinPlatoonRule2, doc = "platoon i takes at most 1 vehicle at cs s")
# # 
# # def vehicleNoinPlatoonRule(originalModel,i,s):
# #     return mp.n[i,s].data() == sum(mp.y[i,j,s].data() for j in originalModel.J) #IMPORTANT 
# # originalModel.vehicleNoinPlatoon = Constraint(originalModel.I,originalModel.swithO,rule=vehicleNoinPlatoonRule)
# # 
# # 
# # def platoonNoLimitRule1(originalModel,i,s):
# # 	return mp.n[i,s].data() <= mp.w[i,s].data()*originalModel.N_max
# # originalModel.platoonNoLimit1 = Constraint(originalModel.I, originalModel.swithO, rule = platoonNoLimitRule1)
# # 
# # def platoonNoLimitRule2(originalModel,i,s):
# # 	return mp.n[i,s].data() >= mp.w[i,s].data()*2
# # originalModel.platoonNoLimit2 = Constraint(originalModel.I, originalModel.swithO, rule = platoonNoLimitRule2)
# # 
# # =============================================================================
# 
# 
# # def energySaveRule(model,i,s):
# # 	return model.delta[i,s] == .13#(.1+.7*(model.n[i,s]-2)+.13)/model.n[i,s]       nonlinear term error shows up 
# # model.energySave = Constraint(model.I, model.swithO, rule = energySaveRule)
# 
# 
# 
# def nonNegativeDelayRule1(originalModel,j):
#     return originalModel.dummy[j] >= originalModel.t_jarr[j,list(originalModel.swithOD.data())[-1]]-originalModel.t_arrmax[j]
# originalModel.nonNegativeDelay1 = Constraint(originalModel.J,rule = nonNegativeDelayRule1)
# 
# def nonNegativeDelayRule2(originalModel,j):
#     return originalModel.dummy[j] >= 0
# originalModel.nonNegativeDelay2 = Constraint(originalModel.J,rule = nonNegativeDelayRule2)
# 
# # Objetive
# def objectiveRule(originalModel):
#     chargingCost = sum(originalModel.c_cha[s]*originalModel.r_cha[s]*originalModel.t_cha[j,s] for j in originalModel.J for s in originalModel.S)
#     delayCost = sum(originalModel.c_delay[j]*originalModel.dummy[j] for j in originalModel.J)
#     energyCost  = sum(originalModel.c_energy[j]*originalModel.d[s,s+1]*originalModel.m[j]*(1-sum(originalModel.delta[i,s]*mp.y[i,j,s].data() for i in originalModel.I)) for j in originalModel.J for s in originalModel.swithO)
#     return chargingCost + delayCost + energyCost
# originalModel.objective = Objective(rule = objectiveRule, sense = minimize)
# 
# solvername = 'cplex'
# 
# solver=SolverFactory(solvername)
# 
# results=solver.solve(originalModel)
# 
# 
# results.write()
# 
# toc = time.perf_counter()
# print(f'Processing time = {toc-tic}\n')
# 
# # =============================================================================
# #     for c in mp.component_objects(Constraint, active=True):
# #         print ("   Constraint",c)
# #         for index in c:
# #             print ("      ", index, c[index])
# # =============================================================================
#  
#     #iteration += 1
# =============================================================================
    
    
    
