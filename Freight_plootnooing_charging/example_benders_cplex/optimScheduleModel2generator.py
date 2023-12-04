#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 12 21:45:40 2022

@author: rakibulalam
"""

# define timer
import time
tic = time.perf_counter()

#Import pyomo libraries
from pyomo.environ import *
from pyomo.opt import SolverFactory, SolverStatus, TerminationCondition
import random

# Argument parsing strategy through terminal
#import argparse
#arser = argparse.ArgumentParser()
#arser.add_argument("l", type=float,
#                   help="input as route length")
#arser.add_argument("cs", type=int,
#                   help="input as number of charging station")
#arser.add_argument("v", type=int,
#                   help="input as number of vehicles")
#arser.add_argument("p", type=int,
#                   help="input as number of platoons")
#args = parser.parse_args()


# uncomment argument values when above argument parser is uncommented    
routeLength = 200#args.l # Mile
numberOfCS = 20#args.cs
csID = range(1,numberOfCS+1)
csIDwithO = range(0, numberOfCS+1)
csIDwithOD = range(0, numberOfCS+2) 
vehicleID = range(0,5)#args.v)
platoonID = range(0,5)#args.p)

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
random.seed(a=1)
for i in vehicleID:
    diffDepTime[i]=random.randint(0, 4)
    
diffArrTime=dict()
for i in vehicleID:
    diffArrTime[i]=18#diffDepTime[i]+11.16 
# Model
model = ConcreteModel()

# Sets
model.S = Set(initialize=csID,doc="set of cs")
model.J = Set(initialize=vehicleID,doc="set of vehicles")
model.I = Set(initialize=platoonID,doc="set of platoons")
model.swithO = Set(initialize=csIDwithO,doc="s union O") 
model.swithOD = Set(initialize=csIDwithOD,doc="s union OD") 


# Variables
model.y = Var(model.I, model.J, model.swithO, within=Binary, doc = "vehicle $j$ joins platoon $i$ at station $s$")
model.t_jdep = Var(model.J, model.swithOD, within=NonNegativeReals, doc = "time of vehicle j depart from station s")
model.t_jarr = Var(model.J, model.swithOD, within = NonNegativeReals, doc = "time of vehicle j arrive at station s")
model.t_idep = Var(model.I, model.swithOD, within = NonNegativeReals, doc = "time of platoon i depart from station s") #EMPTY
model.t_iarr = Var(model.I, model.swithOD, within = NonNegativeReals, doc = "time of platoon i arrive at station s") #EMPTY
model.tau = Var(model.I,model.swithOD, within = NonNegativeReals, doc="time of platoon i arrive/depart at/from station s")
model.soc_arr = Var(model.J,model.swithOD, within = NonNegativeReals, doc = "soc of vehicle j at station s when arrive")
model.soc_dep = Var(model.J,model.swithOD, within = NonNegativeReals, doc = "soc of vehicle j at station s when departs")
model.t_cha = Var(model.J,model.S, within = NonNegativeReals, doc = "charging time of vehicle j at station s")
model.t_wai = Var(model.J,model.S, within = NonNegativeReals, doc = "waiting time (after charge) of vehicle j at station s")
model.dummy = Var(model.J,within = Reals, doc = "dummy var for delay cost calc")
model.w = Var(model.I,model.swithO, within =Binary, doc = "indicator variables")

# Parameters
model.d = Param(model.swithOD, model.swithOD, initialize=distanceBetStation, doc = "distance between s1 and s2")
model.m = Param(model.J, initialize =efficiencyVehicle, doc = "energy efficiency of vehicle kWh/ mile")
model.v_j = Param(model.J, initialize = 50, doc = "average speed of vehicle j")
model.v_i = Param(model.I, initialize = 50, doc = "average speed of platoon i")
model.r_cha = Param(model.S, initialize = 100, doc = "charging speed kw of station s")
#model.t_que = Param(model.J, model.S, initialize = .25, doc = "waiting time hr before charge at s")
model.t_depmin = Param(model.J, initialize = 0, doc = "earliest depart time of vehicle j")
model.t_arrmax = Param(model.J, initialize = diffArrTime, doc = "latest arrive time hr of vehicle j")
model.soc_depmin = Param(model.J, initialize = 1, doc = "minimum departure soc of vehicle j")
model.soc_arrmin = Param(model.J, initialize = 0.2, doc = "minimum arrival soc of vehicle j")
model.soc_min = Param(model.J, initialize = 0.2 , doc = "minimim soc allowed for vehicle j")
model.soc_max = Param(model.J, initialize = 1,doc = "maximum soc allowed for vehicle j")
model.M = Param(initialize = 25, doc = "big number")
model.delta = Param(model.I, model.swithO, initialize = .13, doc = "energy saving % of platoon i depart from station s")
model.c_cha = Param(model.S, initialize =.35 , doc = "charging cost per kwh at station s per kWh https://freewiretech.com/difference-between-ev-charging-levels/")
model.c_delay = Param(model.J, initialize = 50 , doc = "delay cost of vehicle j ")
model.c_energy = Param(model.J, initialize = 0.1165 , doc = "energy cost of vehicle j https://www.eia.gov/state/print.php?sid=FL")
model.b = Param(model.J,initialize = batteryCapVehicle, doc = "battery capacity of vehicle j")
model.N_max = Param(initialize = len(vehicleID), doc = "maximum vehicle in a platoon")
model.t_dep = Param(model.J, initialize = diffDepTime, doc = "different departure time for vehicles")


# Constraints
# SOC 
def initialDepartSOCRule(model,j):
	return model.soc_dep[j,0] ==  model.soc_depmin[j]
model.initialSOC = Constraint(model.J, rule = initialDepartSOCRule, doc = "init depart SOC = depat soc from origin")

def lastArrivalSOCRule(model,j):
	return model.soc_arr[j,list(model.swithOD.data())[-1]] >= model.soc_arrmin[j] 
model.lastArrivalSOC = Constraint(model.J, rule = lastArrivalSOCRule, doc = "last arrival SOC = arrival minimum soc")

def departSOCfromCSRule(model,j,s):
	return model.soc_dep[j,s] == model.soc_arr[j,s] + (model.r_cha[s]/model.b[j])*model.t_cha[j,s]
model.departSOCfromCS = Constraint(model.J, model.S, rule = departSOCfromCSRule)

def arriveSOCatCSRule(model,j,s):
    return model.soc_arr[j,s+1] == model.soc_dep[j,s]- ((model.d[s,s+1]*model.m[j])/model.b[j])*(1 - sum(model.delta[i,s]*model.y[i,j,s] for i in model.I))

model.arriveSOCatCS = Constraint(model.J,model.swithO, rule = arriveSOCatCSRule)

def arriveSOCMinRule(model,j,s):
	return model.soc_arr[j,s] >= model.soc_min[j]
model.arriveSOCMin = Constraint(model.J, model.S, rule = arriveSOCMinRule)

def departSOCMaxRule(model,j,s):
	return model.soc_dep[j,s] <= model.soc_max[j]
model.departSOCMax = Constraint(model.J, model.S, rule = departSOCMaxRule)

# Time
# ra 11/4
def initialDepartTimeRule(model,j):
    return model.t_jdep[j,0] >= model.t_depmin[j]  
	#return model.t_jdep[j,0] == model.t_dep[j]
model.initialDepartTime = Constraint(model.J, rule = initialDepartTimeRule)
#model.initialDepartTime.pprint()

def departTimeRule(model,j,s):
    return model.t_jdep[j,s] == model.t_jarr[j,s] + model.t_cha[j,s] + model.t_wai[j,s]
model.departTime = Constraint(model.J, model.S, rule = departTimeRule)

# =============================================================================
# def relChaAndQueRule1(model,j,s):
#     return model.t_cha[j,s] <=model.M*model.z[j,s]
# model.relChaAndQueRel1 = Constraint(model.J,model.S, rule = relChaAndQueRule1)
# 
# =============================================================================
# =============================================================================
# def relChaAndQueRule2(model,j,s):
#     return model.t_cha[j,s] >= 0
# model.relChaAndQueRel2 = Constraint(model.J,model.S, rule = relChaAndQueRule2)
# 
# =============================================================================
def arriveTimeRule(model,j,s):
	return model.t_jarr[j,s+1] == model.t_jdep[j,s] + model.d[s,s+1]/model.v_j[j]
model.arriveTime = Constraint(model.J, model.swithO, rule = arriveTimeRule)

def platoonDepartRule(model,i,s):
    return model.t_idep[i,s] == model.t_iarr[i,s]
model.platoonDepart = Constraint(model.I,model.S,rule=platoonDepartRule,doc="platton arrive time and depart time same")

def platoonArriveRule(model,i,s):
    return model.t_iarr[i,s+1] == model.t_idep[i,s] + model.d[s,s+1]/model.v_i[i]
model.platoonArriveTime = Constraint(model.I, model.swithO, rule = platoonArriveRule)

def platoonVehicleDepRule1(model,i,j,s):
	return -model.M*(1-model.y[i,j,s]) <= model.t_jdep[j,s] - model.t_idep[i,s] 
model.platoonVehicleDep1 = Constraint(model.I, model.J, model.swithO, rule = platoonVehicleDepRule1)

def platoonVehicleDepRule2(model,i,j,s):
	return  model.t_jdep[j,s] - model.t_idep[i,s] <= model.M*(1-model.y[i,j,s])
model.platoonVehicleDep2 = Constraint(model.I, model.J, model.swithO, rule = platoonVehicleDepRule2)

# Energy
def joinPlatoonRule1(model,j,s):
	return sum(model.y[i,j,s] for i in model.I) <= 1
model.joinPlatoon1= Constraint(model.J, model.swithO, rule = joinPlatoonRule1, doc = "vehicle j join at most 1 platoon at cs s")


def joinPlatoonRule2(model,i,s):
    return sum (model.y[i,j,s] for j in model.J) >= 0
model.joinPlatoon2 = Constraint(model.I, model.swithO, rule = joinPlatoonRule2, doc = "platoon i takes at most 1 vehicle at cs s")

def platoonNoLimitRule1(model,i,s):
	return sum(model.y[i,j,s] for j in model.J) <= model.w[i,s]*model.N_max
model.platoonNoLimit1 = Constraint(model.I, model.swithO, rule = platoonNoLimitRule1)

def platoonNoLimitRule2(model,i,s):
	return sum(model.y[i,j,s] for j in model.J) >= model.w[i,s]*2
model.platoonNoLimit2 = Constraint(model.I, model.swithO, rule = platoonNoLimitRule2)


def nonNegativeDelayRule1(model,j):
    return model.dummy[j] >= model.t_jarr[j,list(model.swithOD.data())[-1]]-model.t_arrmax[j]
model.nonNegativeDelay1 = Constraint(model.J,rule = nonNegativeDelayRule1)

def nonNegativeDelayRule2(model,j):
    return model.dummy[j] >= 0
model.nonNegativeDelay2 = Constraint(model.J,rule = nonNegativeDelayRule2)

# Objetive
def objectiveRule(model):
    chargingCost = sum(model.c_cha[s]*model.r_cha[s]*model.t_cha[j,s] for j in model.J for s in model.S)
    delayCost = sum(model.c_delay[j]*model.dummy[j] for j in model.J)
    energyCost  = sum(model.c_energy[j]*model.d[s,s+1]*model.m[j]*(1-sum(model.delta[i,s]*model.y[i,j,s] for i in model.I)) for j in model.J for s in model.swithO)
    return chargingCost + delayCost + energyCost
model.objective = Objective(rule = objectiveRule, sense = minimize)

# ignore following lines
# =============================================================================
# solvername = 'cplex'
# #model.dual = Suffix(direction = Suffix.IMPORT)
# solver=SolverFactory(solvername)
# ##solver.options['mip tolerances mipgap'] = .15
# results=solver.solve(model, tee = True)
# #
# results.write()
# =============================================================================
#
## write in docplex to call benders in CPPLEX new version
#import docplex.mp
#from docplex.mp import model_reader
#m = model_reader.ModelReader.read('optimCScheduleModel2.lp', ignore_names=True)
#m.parameters.benders.strategy = 3
#m.print_information()
#msol = m.solve(clean_before_solve=True)
#assert msol is not None, "model can't solve"
#m.report()
#toc = time.perf_counter()
#print(f'Processing time = {toc-tic}\n')

# =============================================================================
# for c in model.component_objects(Constraint, active=True):
#     print ("   Constraint",c)
#     for index in c:
#         print ("      ", index, model.dual[c[index]])
# 
# =============================================================================

#model.pprint()
#model.objective.pprint()
#model.y.pprint()
#model.t_jdep.pprint()
#model.t_jarr.pprint()
#model.t_idep.pprint()
# model.t_iarr = Var(model.I, model.swithOD, within = NonNegativeReals, doc = "time of platoon i arrive at station s")
# model.tau = Var(model.I,model.swithOD, within = NonNegativeReals, doc="time of platoon i arrive/depart at/from station s")
# model.soc_arr = Var(model.J,model.swithOD, within = NonNegativeReals, doc = "soc of vehicle j at station s when arrive")
# model.soc_dep = Var(model.J,model.swithOD, within = NonNegativeReals, doc = "soc of vehicle j at station s when departs")
# model.t_cha = Var(model.J,model.S, within = NonNegativeReals, doc = "charging time of vehicle j at station s")
# model.t_wai = Var(model.J,model.S, within = NonNegativeReals, doc = "waiting time (after charge) of vehicle j at station s")
# model.n

# plot y value of vehicle id 2 in platoon 1  and 2 for each cs
# =============================================================================
# =============================================================================
# y_0_2_s=[]
# y_1_2_s=[]
# #y_2_2_s=[]
# for s in model.S:
#     y_0_2_s.append(model.y[0,2,s].value)
#     y_1_2_s.append(model.y[1,2,s].value)
#     #y_2_2_s.append(model.y[2,2,s].value)
# import matplotlib.pyplot as plt
# fig, ax = plt.subplots()
# ax.set_xlim(1,18)
# ax.plot(list(model.S.value),y_0_2_s, label ="i=0")
# ax.plot(list(model.S.value),y_1_2_s, label = "i=1")
# #ax.plot(list(model.S.value),y_2_2_s, label = "i=2")
# ax.set(xlabel = "cs s", ylabel = "y value", title = "y value of j=2 for both i=0 (blue) and 1(orange)")
# plt.show()
# =============================================================================
#fig.savefig("plot1.png")
# 
# # 
# 
# for s in model.swithO:
#    print(model.y[0,0,s].value)
# 
# for s in model.swithO:
#    print(model.y[1,0,s].value)
# 
# for s in model.swithO:
#    print(model.y[2,0,s].value)
# 
# for s in model.swithO:
#    print(model.y[3,0,s].value)
# 
# for s in model.swithO:
#    print(model.y[4,0,s].value)
# 
# 
# 
# 
# for s in model.swithOD:
# 	print(model.t_jarr[0,s].value)
# 
# for s in model.S:
# 	print(model.t_cha[0,s].value)
# 
# 
# for s in model.swithOD:
# 	print(model.t_jdep[0,s].value)
#     
# 
# 
# for s in model.S:
# 	print(model.t_wai[0,s].value)
# 
# 
# 
# for s in model.swithOD:
# 	print(model.t_iarr[0,s].value)
# 
# for s in model.swithOD:
# 	print(model.t_idep[0,s].value)
# 
# for s in model.swithOD:
# 	print(model.t_iarr[1,s].value)
# 
# for s in model.swithOD:
# 	print(model.t_idep[1,s].value)
# 
# for s in model.swithOD:
# 	print(model.t_iarr[2,s].value)
# 
# for s in model.swithOD:
# 	print(model.t_idep[2,s].value)
# 
# for s in model.swithOD:
# 	print(model.t_iarr[4,s].value)
# 
# for s in model.swithOD:
# 	print(model.t_idep[4,s].value)
# =============================================================================


# find cost partwise
#print(f'charging cost: {sum(model.c_cha[s]*model.r_cha[s]*model.t_cha[j,s].value for j in model.J for s in model.S)}\n'
#       f'delay cost: {sum(model.c_delay[j]*model.dummy[j].value for j in model.J)}\n'
#       f'energy cost: {sum(model.c_energy[j]*model.d[s,s+1]*model.m[j]*(1-sum(model.delta[i,s]*model.y[i,j,s].value for i in model.I)) for j in model.J for s in model.swithO)}\n'
#       f'total cost: {model.objective.value()}\n')

      

# =============================================================================
# for i in model.I:
#     for s in model.S:
#        	print(f'i: {i} s: {s} n: {model.n[i,s].value}')
# 
# for j in model.J:
# 	print(model.t_jdep[j,0].value)
# 
# for i in model.I:
#  	print(model.t_idep[i,0].value)
#      
#      
# col_j=[]
# col_s=[]
# col_arrTime=[]
# col_chaTime=[]
# 
# for j in model.J:
#     for s in model.S:
#         col_j.append(j)
#         col_s.append(s)
#         col_arrTime.append(model.t_jarr[j,s].value)
#         col_chaTime.append(model.t_cha[j,s].value)
#         
# with open("chargingTimeSummaryVeh5CS20.csv","w")as file:
#     for elem1,elem2,elem3,elem4 in zip(col_j,col_s,col_arrTime,col_chaTime):
#         file.write(str(elem1)+"\t"+str(elem2)+"\t"+str(elem3)+"\t"+str(elem4)+"\n")
# =============================================================================
