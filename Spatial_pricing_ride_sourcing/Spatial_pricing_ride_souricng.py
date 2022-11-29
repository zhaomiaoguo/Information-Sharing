#  ___________________________________________________________________________
#
#  Pyomo: Python Optimization Modeling Objects
#  Copyright 2017 National Technology and Engineering Solutions of Sandia, LLC
#  Under the terms of Contract DE-NA0003525 with National Technology and
#  Engineering Solutions of Sandia, LLC, the U.S. Government retains certain
#  rights in this software.
#  This software is distributed under the 3-clause BSD License.
#  ___________________________________________________________________________

# Example 5.1.1 from
#
# Practical Bilevel Optimization: Algorithms and Applications
#   Jonathan Bard
from __future__ import division
from pyomo.environ import *
import random
from pyomo.opt import SolverFactory,SolverStatus,TerminationCondition
from pyutilib.misc.timing import tic,toc
import numpy as np
import matplotlib.pyplot as plt
import time
import os
#from time import clock, strftime,localtime
from time import strftime,localtime
from pyomo.environ import log as pyolog
from pyomo.environ import exp as pyoexp
from pyomo.bilevel import *

quadsol = 'cplex'
nlpsol = 'ipopt'
minlpsol = 'scipampl'
mipsol = 'cplex'
bisol = 'bilevel_blp_global'
infinity = float('inf')
class Lib:
    def Network_topo(N, A, R, S):
        # incidence matrix nodes to link
        mA = {}
        for n in N:
            for (r,s) in A:
                if n == r:
                    mA[n,r,s]=1.0
                elif n == s:
                    mA[n,r,s]=-1.0
                else:
                    mA[n,r,s]=0.0
        # incidance matrix between nodes and origin/destination
        mE = {}
        for n in N:
            for r in R:
                for s in S:
                    mE[n,r,s]=0.0
                    if n == r:
                        mE[n,r,s]=1.0
                    if n == s:
                        mE[n,r,s]=-1.0
                    if r == s:
                        mE[n,r,s]=0.0
        def NodesIn_init(n):
            retval = []
            for (i,j) in A:
                if j == n:
                    retval.append(i)
            return retval

        NodesIn = {}
        for n in N:
            NodesIn[n] = NodesIn_init(n)

        def NodesOut_init(n):
            retval = []
            for (i,j) in A:
                if i == n:
                    retval.append(j)
            return retval
        NodesOut = {}
        for n in N:
            NodesOut[n] = NodesOut_init(n)

        return mA, mE, NodesOut, NodesIn

    def Network_init(N, A, R, S, S2, K, alpha, cap, tff, b, supply, demand_cap, demand_slope, b0, b1, b2):
        mA,mE,NodesOut,NodesIn = Lib.Network_topo(N,A,R,S)
        Ntw = Network(nodes=N, links=A, origin=R, destination=S, origin2 = S2, destination2 = K, mA=mA, mE=mE, alpha=alpha, cap=cap, tff=tff, b=b, supply=supply, demand_cap=demand_cap, demand_slope = demand_slope, b0=b0, b1=b1, b2= b2, NodesOut=NodesOut, NodesIn = NodesIn)

        return Ntw

class IO():
    v = {}
    q = {}
    rho = {}
    def write_file(f1, dict1):
        for i in dict1.keys():
            f1.write("%i, %f," % (i, dict1[i]))
        f1.write("\n")
  
    def three_node():
        S = [2,3]
        R = [1]
        S2 = [2,3]
        K = [1]
        N = set(range(1,4))
        A = [(1,2), (1,3), (2,1), (2,3), (3,1), (3,2)]
        # utility parameters for charging facility choice
        b0 = {}
        for s in S:
            b0[s] = 0.0 # locational attractiveness
        b1 = 1 # travel time
        b2 = 3 # price

        # total origin
        supply = {} # existing drivers initial location
        demand_cap = {} # maximum number of riders when price equal to zero
        demand_slope = {} # slope of demand function: demand = demand - demand_slope*rho
        for n in N:
            supply[n] = 0
            if n in R:
                supply[n] = 50
        total_supply = sum(supply[r] for r in R)
        for n in N:
            demand_cap[n] = 0
            demand_slope[n] = 0
            if n in S:
                demand_cap[n] = 300
                demand_slope[n] =5
        #demand_cap[3] = 100
        #supply = {2:600, 11:100, 13: 600, 19:400, 21: 300}
        #demand = {1:300, 7:400, 14:600, 20:100, 24:600}

        # travel cost function: use BPR function
        alpha = {}
        cap = {}
        tff = {}
        b = {}
        for (r,s) in A:
                alpha[r,s] = 2.0
                b[r,s] = 0.15
                cap[r,s] = 10
                tff[r,s] = 10
        cap[1,3] = 5
        cap[3,1] = 5
        return Lib.Network_init(N, A, R, S, S2, K,alpha, cap, tff, b, supply, demand_cap, demand_slope, b0, b1, b2)
    def four_node():
        S = [2,4]
        R = [1,3]
        N = set(range(1,5))
        A = [(1,2), (1,3), (2,3), (3,4), (2,4), (3,2)]
        # utility parameters for charging facility choice
        b0 = {}
        for s in S:
            b0[s] = 0.0 # locational attractiveness
        b1 = 1 # travel time
        b2 = 0.06 # price

        # total origin
        supply = {} # existing drivers initial location
        demand_cap = {} # maximum number of riders when price equal to zero
        demand_slope = {} # slope of demand function: demand = demand - demand_slope*rho
        for n in N:
            supply[n] = 0
            if n in R:
                supply[n] =25
        total_supply = sum(supply[r] for r in R)
        for n in N:
            demand_cap[n] = 0
            demand_slope[n] = 0
            if n in S:
                demand_cap[n] = 30
                demand_slope[n] = 1

        #supply = {2:600, 11:100, 13: 600, 19:400, 21: 300}
        #demand = {1:300, 7:400, 14:600, 20:100, 24:600}

        # travel cost function: use BPR function
        alpha = {}
        cap = {}
        tff = {}
        b = {}
        for (r,s) in A:
                alpha[r,s] = 4.0
                b[r,s] = 0.14
                cap[r,s] = 10
                tff[r,s] = 10
        cap[1,3] = 5
        return Lib.Network_init(N, A, R, S, alpha, cap, tff, b, supply, demand_cap, demand_slope, b0, b1, b2)
    def Sioux_Fall():
        #S = [5, 15,20]
        #S2 = [5, 15,20]
        #R = [3, 13, 23]
        #K = [3, 13, 23]
        S = [1,3,5,7,9,11,13,15,17,19, 21,23]
        S2 = [1,3,5,7,9,11,13,15,17,19, 21,23]
        R = [2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24]
        K = [2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24]
        #R = [2, 4]
        #K = [2, 4]
        #S = [11,23]
        #S2 = [11,23]
        N = set(range(1,25))
        A = [(1,2),(1,3),
                (2,1),(2,6),
                (3,1),(3,4),(3,12),
                (4,3),(4,5),(4,11),
                (5,4),(5,6),(5,9),
                (6,2),(6,5),(6,8),
                (7,8),(7,18),
                (8,6),(8,7),(8,9),(8,16),
                (9,5),(9,8),(9,10),
                (10,9),(10,11),(10,15),(10,16),(10,17),
                (11,4),(11,10),(11,12),(11,14),
                (12,3),(12,11),(12,13),
                (13,12),(13,24),
                (14,11),(14,15),(14,23),
                (15,10),(15,14),(15,19),(15,22),
                (16,8),(16,10),(16,17),(16,18),
                (17,10),(17,16),(17,19),
                (18,7),(18,16),(18,20),
                (19,15),(19,17),(19,20),
                (20,18),(20,19),(20,21),(20,22),
                (21,20),(21,22),(21,24),
                (22,15),(22,20),(22,21),(22,23),
                (23,14),(23,22),(23,24),
                (24,13),(24,21),(24,23)]
        # utility parameters for charging facility choice
        b0 = {}
        for s in S:
            b0[s] = 0.0 # locational attractiveness
        b1 = 1.0 # travel time
        b2 = 3.0 # price

        # total origin
        supply = {} # existing drivers initial location
        demand_cap = {} # maximum number of riders when price equal to zero
        demand_slope = {} # slope of demand function: demand = demand - demand_slope*rho
        for n in N:
            supply[n] = 0
            if n in R:
                #supply[n] = 20 + 10* random.uniform(-1, 1)
                supply[n] = 100
        total_supply = sum(supply[r] for r in R)
        for n in N:
            demand_cap[n] = 0
            demand_slope[n] = 0
            if n in S:
                #demand_cap[n] = 10 + 10* random.uniform(-1, 1)
                demand_cap[n] = 250
                demand_slope[n] = 1
        #supply = {2:600, 11:100, 13: 600, 19:400, 21: 300}
        #demand = {1:300, 7:400, 14:600, 20:100, 24:600}

        # travel cost function: use BPR function
        alpha = {}
        cap = {}
        tff = {}
        b = {}
        for (r,s) in A:
            alpha[r,s] = 4.0
            b[r,s] = 0.15
        def init_cap():
            cap[1,2]=259.0020064
            cap[1,3]=234.0347319
            cap[2,1]=259.0020064
            cap[2,6]=49.58180928
            cap[3,1]=234.0347319
            cap[3,4]=171.1052372
            cap[3,12]=234.0347319
            cap[4,3]=171.1052372
            cap[4,5]=177.827941
            cap[4,11]=49.0882673
            cap[5,4]=177.827941
            cap[5,6]=49.47995469
            cap[5,9]=100
            cap[6,2]=49.58180928
            cap[6,5]=49.47995469
            cap[6,8]=48.98587646
            cap[7,8]=78.4181131
            cap[7,18]=234.0347319
            cap[8,6]=48.98587646
            cap[8,7]=78.4181131
            cap[8,9]=50.50193156
            cap[8,16]=50.45822583
            cap[9,5]=100
            cap[9,8]=50.50193156
            cap[9,10]=139.1578842
            cap[10,9]=139.1578842
            cap[10,11]=100
            cap[10,15]=135.1200155
            cap[10,16]=48.54917717
            cap[10,17]=49.93510694
            cap[11,4]=49.0882673
            cap[11,10]=100
            cap[11,12]=49.0882673
            cap[11,14]=48.76508287
            cap[12,3]=234.0347319
            cap[12,11]=49.0882673
            cap[12,13]=259.0020064
            cap[13,12]=259.0020064
            cap[13,24]=50.91256152
            cap[14,11]=48.76508287
            cap[14,15]=51.27526119
            cap[14,23]=49.24790605
            cap[15,10]=135.1200155
            cap[15,14]=51.27526119
            cap[15,19]=145.6475315
            cap[15,22]=95.99180565
            cap[16,8]=50.45822583
            cap[16,10]=48.54917717
            cap[16,17]=52.29910063
            cap[16,18]=196.7989671
            cap[17,10]=49.93510694
            cap[17,16]=52.29910063
            cap[17,19]=48.23950831
            cap[18,7]=234.0347319
            cap[18,16]=196.7989671
            cap[18,20]=234.0347319
            cap[19,15]=145.6475315
            cap[19,17]=48.23950831
            cap[19,20]=50.02607563
            cap[20,18]=234.0347319
            cap[20,19]=50.02607563
            cap[20,21]=50.5991234
            cap[20,22]=50.75697193
            cap[21,20]=50.5991234
            cap[21,22]=52.29910063
            cap[21,24]=48.85357564
            cap[22,15]=95.99180565
            cap[22,20]=50.75697193
            cap[22,21]=52.29910063
            cap[22,23]=50
            cap[23,14]=49.24790605
            cap[23,22]=50
            cap[23,24]=50.78508436
            cap[24,13]=50.91256152
            cap[24,21]=48.85357564
            cap[24,23]=50.78508436
        def init_cap2():
            cap[1,2]=776.682805381695
            cap[1,3]=701.812139046214
            cap[2,1]=776.682805381695
            cap[2,6]=148.683553701964
            cap[3,1]=701.812139046214
            cap[3,4]=513.102185929618
            cap[3,12]=701.812139046214
            cap[4,3]=513.102185929618
            cap[4,5]=533.261907932196
            cap[4,11]=147.203543662937
            cap[5,4]=533.261907932196
            cap[5,6]=148.378117038357
            cap[5,9]=299.875207986689
            cap[6,2]=148.683553701964
            cap[6,5]=148.378117038357
            cap[6,8]=146.896498918528
            cap[7,8]=235.156479757862
            cap[7,18]=701.812139046214
            cap[8,6]=146.896498918528
            cap[8,7]=235.156479757862
            cap[8,9]=151.442772302845
            cap[8,16]=151.311709654105
            cap[9,5]=299.875207986689
            cap[9,8]=151.442772302845
            cap[9,10]=417.299994674626
            cap[10,9]=417.299994674626
            cap[10,11]=299.875207986689
            cap[10,15]=405.191427512271
            cap[10,16]=145.586946014363
            cap[10,17]=149.743005794701
            cap[11,4]=147.203543662937
            cap[11,10]=299.875207986689
            cap[11,12]=147.203543662937
            cap[11,14]=146.234393681294
            cap[12,3]=701.812139046214
            cap[12,11]=147.203543662937
            cap[12,13]=776.682805381695
            cap[13,12]=776.682805381695
            cap[13,24]=152.674149749451
            cap[14,11]=146.234393681294
            cap[14,15]=153.761796139231
            cap[14,23]=147.682260696526
            cap[15,10]=405.191427512271
            cap[15,14]=153.761796139231
            cap[15,19]=436.760838013103
            cap[15,22]=287.855626843116
            cap[16,8]=151.311709654105
            cap[16,10]=145.586946014363
            cap[16,17]=156.83203678938
            cap[16,18]=590.15131190678
            cap[17,10]=149.743005794701
            cap[17,16]=156.83203678938
            cap[17,19]=144.658325876369
            cap[18,7]=701.812139046214
            cap[18,16]=590.15131190678
            cap[18,20]=701.812139046214
            cap[19,15]=436.760838013103
            cap[19,17]=144.658325876369
            cap[19,20]=150.015798343041
            cap[20,18]=701.812139046214
            cap[20,19]=150.015798343041
            cap[20,21]=151.734226535192
            cap[20,22]=152.207575142833
            cap[21,20]=151.734226535192
            cap[21,22]=156.83203678938
            cap[21,24]=146.499761559384
            cap[22,15]=287.855626843116
            cap[22,20]=152.207575142833
            cap[22,21]=156.83203678938
            cap[22,23]=149.937603993344
            cap[23,14]=147.682260696526
            cap[23,22]=149.937603993344
            cap[23,24]=152.291877350765
            cap[24,13]=152.674149749451
            cap[24,21]=146.499761559384
            cap[24,23]=152.291877350765
        init_cap()

        def init_tff():
            tff[1,2]=6
            tff[1,3]=4
            tff[2,1]=6
            tff[2,6]=5
            tff[3,1]=4
            tff[3,4]=4
            tff[3,12]=4
            tff[4,3]=4
            tff[4,5]=2
            tff[4,11]=6
            tff[5,4]=2
            tff[5,6]=4
            tff[5,9]=5
            tff[6,2]=5
            tff[6,5]=4
            tff[6,8]=2
            tff[7,8]=3
            tff[7,18]=2
            tff[8,6]=2
            tff[8,7]=3
            tff[8,9]=10
            tff[8,16]=5
            tff[9,5]=5
            tff[9,8]=10
            tff[9,10]=3
            tff[10,9]=3
            tff[10,11]=5
            tff[10,15]=6
            tff[10,16]=4
            tff[10,17]=8
            tff[11,4]=6
            tff[11,10]=5
            tff[11,12]=6
            tff[11,14]=4
            tff[12,3]=4
            tff[12,11]=6
            tff[12,13]=3
            tff[13,12]=3
            tff[13,24]=4
            tff[14,11]=4
            tff[14,15]=5
            tff[14,23]=4
            tff[15,10]=6
            tff[15,14]=5
            tff[15,19]=3
            tff[15,22]=3
            tff[16,8]=5
            tff[16,10]=4
            tff[16,17]=2
            tff[16,18]=3
            tff[17,10]=8
            tff[17,16]=2
            tff[17,19]=2
            tff[18,7]=2
            tff[18,16]=3
            tff[18,20]=4
            tff[19,15]=3
            tff[19,17]=2
            tff[19,20]=4
            tff[20,18]=4
            tff[20,19]=4
            tff[20,21]=6
            tff[20,22]=5
            tff[21,20]=6
            tff[21,22]=2
            tff[21,24]=3
            tff[22,15]=3
            tff[22,20]=5
            tff[22,21]=2
            tff[22,23]=4
            tff[23,14]=4
            tff[23,22]=4
            tff[23,24]=2
            tff[24,13]=4
            tff[24,21]=3
            tff[24,23]=2
        def init_tff2():
            tff[1,2]=12
            tff[1,3]=8
            tff[2,1]=12
            tff[2,6]=10
            tff[3,1]=8
            tff[3,4]=8
            tff[3,12]=8
            tff[4,3]=8
            tff[4,5]=4
            tff[4,11]=12
            tff[5,4]=4
            tff[5,6]=8
            tff[5,9]=10
            tff[6,2]=10
            tff[6,5]=8
            tff[6,8]=4
            tff[7,8]=6
            tff[7,18]=4
            tff[8,6]=4
            tff[8,7]=6
            tff[8,9]=20
            tff[8,16]=10
            tff[9,5]=10
            tff[9,8]=20
            tff[9,10]=6
            tff[10,9]=6
            tff[10,11]=10
            tff[10,15]=12
            tff[10,16]=8
            tff[10,17]=16
            tff[11,4]=12
            tff[11,10]=10
            tff[11,12]=12
            tff[11,14]=8
            tff[12,3]=8
            tff[12,11]=12
            tff[12,13]=6
            tff[13,12]=6
            tff[13,24]=8
            tff[14,11]=8
            tff[14,15]=10
            tff[14,23]=8
            tff[15,10]=12
            tff[15,14]=10
            tff[15,19]=6
            tff[15,22]=6
            tff[16,8]=10
            tff[16,10]=8
            tff[16,17]=4
            tff[16,18]=6
            tff[17,10]=16
            tff[17,16]=4
            tff[17,19]=4
            tff[18,7]=4
            tff[18,16]=6
            tff[18,20]=8
            tff[19,15]=6
            tff[19,17]=4
            tff[19,20]=8
            tff[20,18]=8
            tff[20,19]=8
            tff[20,21]=12
            tff[20,22]=10
            tff[21,20]=12
            tff[21,22]=4
            tff[21,24]=6
            tff[22,15]=6
            tff[22,20]=10
            tff[22,21]=4
            tff[22,23]=8
            tff[23,14]=8
            tff[23,22]=8
            tff[23,24]=4
            tff[24,13]=8
            tff[24,21]=6
            tff[24,23]=4
        init_tff()
        return Lib.Network_init(N, A, R, S,S2,K, alpha, cap, tff, b, supply, demand_cap, demand_slope, b0, b1, b2)

class Network:
    def __init__(self, nodes, links, origin, destination,origin2, destination2,mA, mE, alpha, cap, tff, b, supply, demand_cap, demand_slope, b0, b1, b2, NodesOut, NodesIn):
        self.N = nodes
        self.A = links
        self.R = origin
        self.S = destination
        self.S2 = origin2
        self.K =  destination2
        self.mA = mA
        self.mE = mE
        self.alpha = alpha
        self.cap = cap
        self.tff = tff
        self.b = b
        self.supply = supply
        self.demand_cap = demand_cap
        self.demand_slope = demand_slope
        self.b0 = b0
        self.b1 = b1
        self.b2 = b2
        self.NodesOut = NodesOut
        self.NodesIn = NodesIn
        self.epsilon = 0.0001
        self.congestion =True 
        self.capacity_scaler =1
        self.opt_gap = 0.001
        self.background_trip = True
    def combine(self, transaction_ini):
        if not self.congestion:
            for (i,j) in self.A:
                self.b[i,j] = 0
        else:
            for (i,j) in self.A:
                self.cap[i,j] = self.cap[i,j] * self.capacity_scaler

        model = AbstractModel('CDA')

        def transaction_init_rule(model, s):
            return transaction_ini[s]
        model.transaction =  Var(self.S, within=NonNegativeReals, initialize =transaction_init_rule)

        def q_bounds_rule(model,r,s):
            return (0.0000001, sum (self.supply[r] for r in self.R))
        model.q = Var(self.R, self.S, bounds = q_bounds_rule, initialize=10)#=1e-1)
        model.d = Var(self.S2, self.K, within=NonNegativeReals, initialize=50)#=1e-1)
        model.v_rs = Var(self.A,self.R, self.S, within=NonNegativeReals, initialize=10.0) # traffic flow
        model.v_s2k = Var(self.A,self.S2, self.K, within=NonNegativeReals, initialize=10.0) # traffic flow
        model.v = Var(self.A, within=NonNegativeReals, initialize=0.0) # traffic flow

        print("Objective...")
        def obj_rule(model):
            exp = self.b1/self.b2*sum(self.tff[r,s]*(model.v[r,s]+(self.b[r,s]/(self.alpha[r,s]+1.0))*(model.v[r,s]**(self.alpha[r,s]+1))/(self.cap[r,s]**(self.alpha[r,s]))) for (r,s) in self.A)
            exp1 = 1.0/self.b2*( sum( model.q[r,s]*( pyolog(model.q[r,s]) - 1.0 -self.b0[s]) for s in self.S for r in self.R))
            exp2 = sum(0.5/self.demand_slope[s] * model.transaction[s]**2 - 1/ self.demand_slope[s]*self.demand_cap[s]*model.transaction[s] for s in self.S)
            return exp+exp1+exp2
        model.obj = Objective(rule=obj_rule, sense=minimize)

        print("Flow conservation...")
        def con_v_rs_rule(model, n, r, s):
            if n == r and n ==s:
                return (0,sum(model.v_rs[n,j,r,s] for j in self.NodesOut[n]) - sum(model.v_rs[i,n,r,s] for i in self.NodesIn[n]), 0)
            elif n ==r:
                return (0, sum(model.v_rs[n,j,r,s] for j in self.NodesOut[n]) - sum(model.v_rs[i,n,r,s] for i in self.NodesIn[n]) - model.q[r,s], 0)
            elif n == s:
                return (0, sum(model.v_rs[n,j,r,s] for j in self.NodesOut[n]) - sum(model.v_rs[i,n,r,s] for i in self.NodesIn[n]) + model.q[r,s], 0)
            else:
                return (0,sum(model.v_rs[n,j,r,s] for j in self.NodesOut[n]) - sum(model.v_rs[i,n,r,s] for i in self.NodesIn[n]), 0)
        model.con_v_rs = Constraint(self.N, self.R, self.S, rule=con_v_rs_rule)
        def con_v_s2k_rule(model, n, s2, k):
            if n == s2 and n == k:
                return (0,sum(model.v_s2k[n,j,s2,k] for j in self.NodesOut[n]) - sum(model.v_s2k[i,n,s2,k] for i in self.NodesIn[n]), 0)
            elif n ==s2:
                return (0, sum(model.v_s2k[n,j,s2,k] for j in self.NodesOut[n]) - sum(model.v_s2k[i,n,s2,k] for i in self.NodesIn[n]) - model.d[s2,k], 0)
            elif n == k:
                return (0, sum(model.v_s2k[n,j,s2,k] for j in self.NodesOut[n]) - sum(model.v_s2k[i,n,s2,k] for i in self.NodesIn[n]) + model.d[s2,k], 0)
            else:
                return (0,sum(model.v_s2k[n,j,s2,k] for j in self.NodesOut[n]) - sum(model.v_s2k[i,n,s2,k] for i in self.NodesIn[n]), 0)
        model.con_v_s2k = Constraint(self.N, self.S2, self.K, rule=con_v_s2k_rule)

        def con_v_rule(model, i, j):
            return (0,sum(model.v_rs[i,j,r,s] for r in self.R for s in self.S)+sum(model.v_s2k[i,j,s2,k] for s2 in self.S2 for k in self.K) - model.v[i,j], 0)
        model.con_v = Constraint(self.A, rule=con_v_rule)

        print("Drivers origin...")
        def con_driver_rule(model, r):
            return (0, sum( model.q[r,s] for s in self.S) - self.supply[r], 0)
        model.con_driver = Constraint(self.R, rule=con_driver_rule)
        print("Riders origin-destination...")
        def con_riders_rule(model, s2, k):
            if self.background_trip:
                return (0, model.d[s2,k] -self.demand_cap[s2]/(len(self.K)), 0) # assume equaly distributed to each k
            else:
                return (0, model.d[s2,k],0)
        model.con_rider = Constraint(self.S2, self.K, rule=con_riders_rule)

        print("Balance...")
        def con_bal_rule(model, s):
            return (0, model.transaction[s] - sum(model.q[r,s] for r in self.R), 0)
        model.con_bal= Constraint(self.S, rule=con_bal_rule)

        print("Select solver...")
        optsolver = SolverFactory(nlpsol)
        print("Create instance...")
        inst = model.create_instance(report_timing=False)
        inst.pprint()
        #print("Dual direction...")
        inst.dual = Suffix(direction=Suffix.IMPORT)

        print("Start solving...")
        results = optsolver.solve(inst, tee=False, load_solutions=False)

        inst.solutions.load_from(results)

        v={}
        print("CDA: print v...")
        for (i,j) in self.A:
            v[i,j]=value(inst.v[i,j])
            print ("%i, %i : %f" % (i, j, v[i,j]))
        q={}
        print("CDA: print q...")
        for r in self.R:
            for s in self.S:
                q[r,s]=value(inst.q[r,s])
                print ("%i, %i : %f" % (r, s, q[r,s]))
        v_rs={}
        print("CDA: print v_rs...")
        for (i,j) in self.A:
            for r in self.R:
                for s in self.S:
                    v_rs[i,j,r,s] = value(inst.v_rs[i,j,r,s])
                    print ("%i, %i, %i, %i: %f" % (i,j,r,s, v_rs[i,j, r, s]))

        lamda_rs={}
        print("CDA: print lamda_rs...")
        inst.dual.display()
        for n in self.N:
            for r in self.R:
                for s in self.S:
                    lamda_rs[n,r,s] = -inst.dual[inst.con_v_rs[n,r,s]]
                    print ("%i, %i, %i: %f" % (n,r,s, lamda_rs[n, r, s]))
        gamma={}
        print("CDA: print gamma...")
        for r in self.R:
            gamma[r] = -inst.dual[inst.con_driver[r]]
            print ("%i: %f" % (r, gamma[r]))
        rho={}
        print("CDA: print rho...")
        for (s) in self.S:
            rho[s]= -inst.dual[inst.con_bal[s]]
            print ("%i: %f" % (s, rho[s]))

        transaction={}
        print("CDA: print transaction...")
        for (s) in self.S:
            transaction[s]=value(inst.transaction[s])
            print ("%i: %f" % (s, transaction[s]))
        return  q, v, lamda_rs, gamma,transaction, rho

    def problem(self,demand_cap, supply):
        for s in self.S:
            self.demand_cap[s] = demand_cap
        for r in self.R:
            self.supply[r] = supply
        rho_ini = {}
        transaction_ini = {}
        for s in self.S:
            rho_ini[s] = self.demand_cap[s]/self.demand_slope[s]/2.0
            transaction_ini[s] = self.demand_cap[s] - self.demand_slope[s] * rho_ini[s]
        # solve for the initialization
        q = {}
        v = {}
        lamda_rs = {}
        gamma = {}
        transaction = {}
        rho = {}
        Maxit = 1
        for iter in range(Maxit):
            time_bq[iter] = time.time()
            print ('Start iteration %i\t' % iter)
            q, v, lamda_rs, gamma, transaction, rho = self.combine(transaction_ini)
            print(rho)
            tt = 0
            ttf = 0
            for (r,s) in self.A:
                ttf = ttf + self.tff[r,s]*v[r,s]

                tt = tt+ v[r,s]*self.tff[r,s]*(1+self.b[r,s]*(v[r,s]**(self.alpha[r,s]))/(self.cap[r,s]**(self.alpha[r,s])))
            tt_ttf = tt/ttf
            print("congestion: %f" %tt_ttf)
            print("free flow time: %f" %ttf)
            print("total travel time: %f" %tt)
        return  q, v, rho

if __name__ == "__main__":
    #Ntw = IO.Sioux_Fall_6()
    Ntw = IO.Sioux_Fall()
    #Ntw = IO.three_node()
    #Ntw = IO.four_node()
    time_bq = {}
    start = time.time()
    tol = 0.1
    print ('Stopping critieria %f' % tol)
    #20
    Maxit = 1
    aux_str=strftime("%Y%m%d_%H%M", localtime())
    res_path=''+aux_str
    os.system('mkdir '+res_path)
    file_name = res_path+'/rho_bal.txt'
    file1 = open(file_name,"w") 
    for iter in range(Maxit):
        time_bq[iter] = time.time()
        print ('Start iteration %i\t' % iter)
        #demand_cap = 20+20*iter
        demand_cap = 300
        #supply = 20+20*iter
        supply = 100
        IO.q, IO.v, IO.rho = Ntw.problem(demand_cap, supply)
        IO.write_file(file1, IO.rho) 
    file1.close() #to change file access modes 

    end = time.time()
    el_time = end - start
    print ('Elapsed time %f' % el_time)
