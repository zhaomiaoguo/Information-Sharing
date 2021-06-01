from __future__ import division
from pyomo.environ import *
from pyomo.opt import SolverFactory,SolverStatus,TerminationCondition
from pyutilib.misc.timing import tic,toc
import numpy as np
import time
import os
from time import perf_counter, strftime,localtime
from pyomo.environ import log as pyolog

import random
quadsol = 'cplex'
nlsol = 'ipopt'

class Network:
    def __init__(self, nodes, links, origin, destination, facility, scenarios):
        self.N = nodes
        self.A = links
        self.R = origin
        self.S = destination
        self.K = facility
        self.Scn = scenarios
        self.I = {}
        self.C = {}

    def step_1(self, ADMM):
        print ('Starting Step 1:')
        print ('Solving investor profit maximization')
        c = {}
        g = {}
        for u in self.Scn.U:
            self.I.model[u].del_component('obj')
            def obj_rule_utmax(model):
                exp0 = -self.Scn.pr[u] * sum(ADMM.rho[u,k]*model.g[k] - (self.I.ca*model.c[k]**2 + self.I.cb*model.c[k]) - (self.I.ga*model.g[k]**2+self.I.gb*model.g[k]) for k in self.K)
                exp1 = ADMM.r_1*self.Scn.pr[u]/2.0 * sum( ( sum(ADMM.q[u,r,s,k] for r in self.R for s in self.S) - model.g[k])** 2 for k in self.K )
                exp2 = sum( ADMM.lbda[u,k]*(model.c[k] - ADMM.znu[k]) for k in self.K )
                exp3 = ADMM.r_2/2.0 * sum((model.c[k] - ADMM.znu[k])**2 for k in self.K)
                return exp0 + exp1 + exp2 + exp3
            self.I.model[u].add_component('obj', Objective(rule=obj_rule_utmax, sense=minimize))
            c_u, g_u = self.I.profit_maximization(ADMM,u)
            for k in self.K:
                c[u,k] = c_u[k]
                g[u,k] = g_u[k]
        ADMM.c = c
        ADMM.g = g
        print ('Solving combined distribution and assignment problem')
        q = {}
        x = {}
        v = {}
        for u in self.Scn.U:
            self.C.model[u].del_component('obj')
            def obj_rule_tap(model):
                exp0 = sum(self.C.tff[r,s]*(model.v[r,s]+(self.C.b[r,s]/(self.C.alpha[r,s]+1.0))*(model.v[r,s]**(self.C.alpha[r,s]+1))/(self.C.cap[r,s]**(self.C.alpha[r,s]))) for (r,s) in self.C.A)
                exp1 = 1.0/self.C.b1*( sum( model.q[r,s,k]*( pyolog(model.q[r,s,k]) - 1.0 + self.C.b3*ADMM.rho[u,k] -self.C.b0[k] - self.C.b2*ADMM.c[u,k] ) for k in self.K for s in self.S for r in self.R))
                exp2 =  sum((sum(self.C.e[r,s]*model.q[r, s, k] for r in self.R for s in self.S ) - ADMM.g[u, k])** 2 for k in self.K )
                return self.C.b1/self.C.b3*(exp0 + exp1) + ADMM.r_1/2.0 * exp2
            self.C.model[u].add_component('obj', Objective(rule=obj_rule_tap, sense=minimize))
            q_u, x_u, v_u = self.C.traffic_problem(ADMM,u)
            for (i,j) in self.A:
                v[u,i,j] = v_u[i,j]
            for r in self.R:
                for s in self.S:
                    for k in self.K:
                        q[u,r,s,k] = q_u[r,s,k]
                        for (i,j) in self.A:
                            x[u, i,j, r,s,k] = x_u[i,j,r,s,k]
        return q, x, v, c, g

    def step_2(self, ADMM):
        print ('Starting Step 2:')
        print ('Updating z, rho, lambda')
        z_new = {}
        rho_new = {}
        la_new = {}
        for k in self.K:
            z_new[k] = sum(ADMM.c[u, k] for u in self.Scn.U)/len(self.Scn.U)
        for u in self.Scn.U:
            for k in self.K:
                rho_new[u,k] = ADMM.rho[u,k] + ADMM.r_1 * self.Scn.pr[u]*(sum(sum(ADMM.q[u,r,s,k] for r in self.R) for s in self.S) - ADMM.g[u,k])
                la_new[u,k] = ADMM.lbda[u,k] + ADMM.r_2 * (ADMM.c[u,k]-z_new[k] )
        return z_new, rho_new, la_new

    def centralized_problem(self):
        model = AbstractModel('Convex Reformulation')
        model.c = Var(self.K, self.Scn.U, within=NonNegativeReals,initialize=1000.0)
        model.z = Var(self.K, within=NonNegativeReals,initialize=1000.0)
        model.g = Var(self.K, self.Scn.U,  within=NonNegativeReals,initialize=500.0)
        model.q = Var(self.R, self.S, self.K, self.Scn.U, within=NonNegativeReals, initialize=500.0)#=1e-1)
        model.v = Var(self.A, self.Scn.U, within=NonNegativeReals,initialize=0.0) # traffic flow
        model.x = Var(self.A, self.R, self.S, self.K, self.Scn.U, within=NonNegativeReals,initialize=0.0)   # traffic flow on link a associated with O-D pairs (r,s) and passing k
        # break up x into two pars: x1 and x2
        model.x1 = Var(self.A, self.R, self.S, self.K, self.Scn.U, within=NonNegativeReals,initialize=0.0)
        model.x2 = Var(self.A, self.R, self.S, self.K, self.Scn.U, within=NonNegativeReals,initialize=0.0)

        def obj_rule(model):
            exp0 = sum(self.Scn.pr[u]*sum(self.I.ca*model.c[k,u]**2+self.I.cb*model.c[k,u] for k in self.K) for u in self.Scn.U)
            exp1 = sum(self.Scn.pr[u]*sum(self.I.ga*model.g[k, u]**2+self.I.gb*model.g[k, u] for k in self.K) for u in self.Scn.U)
            exp20 = sum(self.Scn.pr[u]*sum(self.C.tff[r,s]*(model.v[r,s,u]+(self.C.b[r,s]/(self.C.alpha[r,s]+1.0))*(model.v[r,s,u]**(self.C.alpha[r,s]+1))/(self.C.cap[r,s]**(self.C.alpha[r,s]))) for (r,s) in self.A) for u in self.Scn.U)
            exp21 = sum(self.Scn.pr[u]*sum( model.q[r,s,k,u]*( pyolog(model.q[r,s,k,u]) - 1.0 - self.C.b0[k]) for k in self.K for s in self.S for r in self.R) for u in self.Scn.U)
            return exp0 + exp1 + self.C.b1/self.C.b3 * ( exp20+ 1.0/self.C.b1*exp21)
        model.obj = Objective(rule=obj_rule, sense=minimize)

        def market_clear_rule(model, k, u):
            return sum( self.C.e[r,s]*model.q[r, s, k,u] for r in self.R for s in self.S) - model.g[k,u] == 0
        model.market_clear = Constraint(self.K, self.Scn.U, rule = market_clear_rule)

        def capacity_rule(model, k, u):
            return model.g[k,u]-model.c[k,u] <= 0
        model.capacity = Constraint(self.K, self.Scn.U, rule=capacity_rule)

        def non_anti_rule(model, k, u):
            return model.c[k,u] == model.z[k]
        model.non_anti = Constraint(self.K,self.Scn.U, rule=non_anti_rule)

        def conb_rule(model, i, j, u):
            return sum( model.x[i,j,r,s,k,u] for k in self.K for s in self.S for r in self.R) == model.v[i,j,u]
        model.conb = Constraint(self.A, self.Scn.U, rule=conb_rule)

        def conb2_rule(model, i, j,r,s,k,u):
            return model.x1[i,j,r,s,k,u] + model.x2[i,j,r,s,k,u] == model.x[i,j,r,s,k,u]
        model.conb2 = Constraint(self.A, self.R, self.S, self.K, self.Scn.U, rule=conb2_rule)

        def conc2_rule(model, r, s, k,n, u):
            return sum(self.C.mA[n,i,j]*model.x1[i,j,r,s,k,u] for (i,j) in self.A) == model.q[r,s,k,u]*self.C.mE[n,r,k]
        model.conc2 = Constraint(self.R, self.S, self.K, self.N,self.Scn.U, rule=conc2_rule)

        def conc3_rule(model, r, s, k,n,u):
            return sum(self.C.mA[n,i,j]*model.x2[i,j,r,s,k,u] for (i,j) in self.A) == model.q[r,s,k,u]*self.C.mE[n,k,s]
        model.conc3 = Constraint(self.R, self.S, self.K, self.N, self.Scn.U, rule=conc3_rule)

        def cond_rule(model, r, s, u):
            return sum( model.q[r,s,k,u] for k in self.K) == self.C.d[r,s]*self.Scn.growth[u]
        model.cond = Constraint(self.R, self.S, self.Scn.U, rule=cond_rule)
        optsolver = SolverFactory(nlsol)
        inst = model.create_instance()
        inst.dual = Suffix(direction=Suffix.IMPORT)
        # initialization

#       inst.c[3].value = 585.8735090176859
#       inst.c[6].value = 228.81446723201296
#       inst.c[12].value = 472.7288981299494
#       inst.c[17].value = 440.00420118220916
#       inst.c[22].value = 772.5789256320681

        results = optsolver.solve(inst, load_solutions=False)
        inst.solutions.load_from(results)
        rho_ans = {}
        x = {}
        v = {}
        q = {}
        g = {}
        c = {}
        z = {}
        for k in self.K:
            z[k] = value(inst.z[k])

        for u in self.Scn.U:
            for k in self.K:
                rho_ans[u,k] = inst.dual[inst.market_clear[k,u]]
                print('Locational price at facility %i at scenario %i: %f' % (k,u, rho_ans[u,k]))
            for (i,j) in self.A:
                v[i,j,u]=value(inst.v[i,j,u])
                for r in self.R:
                    for s in self.S:
                        for k in self.K:
                            x[i,j,r,s,k,u]=value(inst.x[i,j,r,s,k,u])
            for r in self.R:
                for s in self.S:
                    for k in self.K:
                        q[r,s,k,u]=value(inst.q[r,s,k,u])
            for k in self.K:
                g[k,u]=value(inst.g[k,u])
                c[k,u]=value(inst.c[k,u])
        return g, c, q, x, v, rho_ans, z

    def init_ADMM(self):
        rho_0={}
        y = 500
        for u in self.Scn.U:
            for k in self.K:
                rho_0[u,k]=y-400.0*np.random.uniform(low = -1, high = 1)
        q_0  = {}
        for u in self.Scn.U:
            for r in self.R:
                for s in self.S:
                    for k in self.K:
                        q_0[u,r,s,k] = self.C.d[r,s]*self.Scn.growth[u]/len(self.K)
        z_0 = {}
        for k in Ntw.K:
            z_0[k] = sum(self.Scn.growth[u]*sum(self.C.d[r,s] for r in self.R for s in self.S) for u in self.Scn.U)/len(self.Scn.U)/len(self.K)
        la_0 = {}
        for u in self.Scn.U:
            for k in self.K:
                la_0[u,k] = 0
        c_0 = {}
        g_0 = {}
        for u in self.Scn.U:
            for k in self.K:
                c_0[u,k] = 0.0
                g_0[u,k] = 0.0
        Algo = ADMM(rho=rho_0,
                    q=q_0,
                    r1=100,#ADMM
                    r2=1,#PH
#                   pr=,
                    la=la_0,
                    z=z_0,
                    c=c_0,
                    g=g_0)#,
#                   growth)
        return Algo

    def exsu(self, ADMM):
        ES={}
        for u in self.Scn.U:
            for k in self.K:
                ES[u,k] = ADMM.g[k,u] - sum( self.C.e[r,s]*ADMM.q[r,s,k,u] for r in self.R for s in self.S)
        return ES

    def scen_diff(self, ADMM):
        SD={}
        for u in self.Scn.U:
            for k in Ntw.K:
                SD[u, k] = ADMM.c[k,u] - ADMM.znu[k]
        return SD

    # all the output can be saved in this function.
    def write_evs(self,EvES,EvSD,EvR,res_path):
        f_exsu = open((res_path+'/Resulting_exsu.csv'),'w')
        f_scdi = open((res_path+'/Resulting_scendiff.csv'),'w')
        f_prices = open((res_path+'/Resulting_prices.csv'),'w')
        for u in self.Scn.U:
            for k in self.K:
                aux_p = 'Scenario%i,Node%i,'%(u,k)
                aux_s = 'Scenario%i,Node%i,'%(u,k)
                aux_sd = 'Scenario%i,Node%i,'%(u,k)
                for nu in EvES:
                    aux_p = aux_p+('%f'%(EvR[nu][u,k]))+','
                    aux_s = aux_s+('%f'%(EvES[nu][u,k]))+','
                    aux_sd = aux_sd+('%f'%(EvSD[nu][u,k]))+','
                f_prices.write(aux_p+',\n')
                f_exsu.write(aux_s+',\n')
                f_scdi.write(aux_sd+',\n')
        f_prices.close()
        f_exsu.close()
        f_prices.close()

class Investors(Network):
    def __init__(self, nodes, links, origin, destination, facility, costca, costcb, costga, costgb, scenarios): #Locations,
#       self.K = Locations
        self.ca = costca
        self.cb = costcb
        self.ga = costga
        self.gb = costgb
        self.model = {}
        super().__init__(nodes, links, origin, destination, facility, scenarios)

    def profit_maximization_model(self, ADMM, u):
        model = AbstractModel('profit_max')
        model.c = Var(self.K, within=NonNegativeReals,initialize=1e5)
        model.g = Var(self.K, within=NonNegativeReals,initialize=1e5)
        def obj_rule(model):
            exp0 = -self.Scn.pr[u] * sum(ADMM.rho[u,k]*model.g[k] - (self.ca*model.c[k]**2 + self.cb*model.c[k]) - (self.ga*model.g[k]**2+self.gb*model.g[k]) for k in self.K)
            exp1 = ADMM.r_1*self.Scn.pr[u]/2.0 * sum( ( sum(ADMM.q[u,r,s,k] for r in self.R for s in self.S) - model.g[k])** 2 for k in self.K )
            exp2 = sum( ADMM.lbda[u,k]*(model.c[k] - ADMM.znu[k]) for k in self.K )
            exp3 = ADMM.r_2/2.0 * sum((model.c[k] - ADMM.znu[k])**2 for k in self.K)
            return exp0 + exp1 + exp2 + exp3
        model.obj = Objective(rule=obj_rule, sense=minimize)

        def capacity_rule(model, k):
            return model.g[k]-model.c[k] <= 0
        model.capacity = Constraint(self.K, rule=capacity_rule)
        inst = model.create_instance()
        return inst

    def profit_maximization(self,ADMM,u):
        optsolver = SolverFactory(quadsol)
        results = optsolver.solve(self.model[u], load_solutions=False)
        self.model[u].solutions.load_from(results)
        g={}
        c={}
        for k in self.K:
            g[k]=value(self.model[u].g[k])
            c[k]=value(self.model[u].c[k])
        return c, g

class Consumers(Network):
    def __init__(self,nodes, links, origin, destination, facility, b0, b1, b2, b3, inc, tff, b, alpha, cap, d, mA, mA_k, mE, mE_k, e, scenarios):
        self.b0 = b0
        self.b1 = b1
        self.b2 = b2
        self.b3 = b3
        self.inc = inc
        self.tff = tff
        self.b = b
        self.alpha = alpha
        self.cap = cap
        self.d = d
        self.mA = mA
        self.mA_k = mA_k
        self.mE = mE
        self.mE_k = mE_k
        self.e = e
        self.model = {}
        super().__init__(nodes, links, origin, destination, facility, scenarios)

    def traffic_problem_model(self, ADMM, u):
        model = AbstractModel('CDA')
        model.q = Var(self.R, self.S, self.K, within=NonNegativeReals, initialize=500.0)#=1e-1)
        model.v = Var(self.A, within=NonNegativeReals,initialize=0.0) # traffic flow
        model.x = Var(self.A, self.R, self.S, self.K, within=NonNegativeReals,initialize=0.0)   # traffic flow on link a associated with O-D pairs (r,s) and passing k
        # break up x into two pars: x1 and x2
        model.x1 = Var(self.A, self.R, self.S, self.K, within=NonNegativeReals,initialize=0.0)
        model.x2 = Var(self.A, self.R, self.S, self.K, within=NonNegativeReals,initialize=0.0)

        def objrule(model):
            exp0 = sum(self.tff[r,s]*(model.v[r,s]+(self.b[r,s]/(self.alpha[r,s]+1.0))*(model.v[r,s]**(self.alpha[r,s]+1))/(self.cap[r,s]**(self.alpha[r,s]))) for (r,s) in self.A)
            exp1 = 1.0/self.b1*( sum( model.q[r,s,k]*( pyolog(model.q[r,s,k]) - 1.0 + self.b3*ADMM.rho[u,k] -self.b0[k] - self.b2*ADMM.c[u,k] ) for k in self.K for s in self.S for r in self.R))
            exp2 =  sum((sum(self.e[r,s]*model.q[r, s, k] for r in self.R for s in self.S ) - ADMM.g[u, k])** 2 for k in self.K )
            return self.b1/self.b3*(exp0 + exp1) + ADMM.r_1/2.0 * exp2
        model.obj = Objective(rule=objrule)

        def conb_rule(model, i, j):
            return sum( model.x[i,j,r,s,k] for k in self.K for s in self.S for r in self.R) == model.v[i,j]
        model.conb = Constraint(self.A, rule=conb_rule)

        def conb2_rule(model, i, j,r,s,k):
            return model.x1[i,j,r,s,k] + model.x2[i,j,r,s,k] == model.x[i,j,r,s,k]
        model.conb2 = Constraint(self.A, self.R, self.S, self.K, rule=conb2_rule)

        def conc2_rule(model, r, s, k,n):
            return sum(self.mA[n,i,j]*model.x1[i,j,r,s,k] for (i,j) in self.A) == model.q[r,s,k]*self.mE[n,r,k]
        model.conc2 = Constraint(self.R, self.S, self.K, self.N, rule=conc2_rule)

        def conc3_rule(model, r, s, k,n):
            return sum(self.mA[n,i,j]*model.x2[i,j,r,s,k] for (i,j) in self.A) == model.q[r,s,k]*self.mE[n,k,s]
        model.conc3 = Constraint(self.R, self.S, self.K, self.N, rule=conc3_rule)

        def cond_rule(model, r, s):
            return sum( model.q[r,s,k] for k in self.K) == self.d[r,s] * self.Scn.growth[u]
        model.cond = Constraint(self.R, self.S, rule=cond_rule)
        inst = model.create_instance()
        return inst

    def traffic_problem(self, ADMM, u):
        x={}
        v={}
        q={}
        optsolver = SolverFactory(nlsol)
        results = optsolver.solve(self.model[u], load_solutions=False)
        self.model[u].solutions.load_from(results)
        for (i,j) in self.A:
            v[i,j]=value(self.model[u].v[i,j])
            for r in self.R:
                for s in self.S:
                    for k in self.K:
                        x[i,j,r,s,k]=value(self.model[u].x[i,j,r,s,k])
        for r in self.R:
            for s in self.S:
                for k in self.K:
                    q[r,s,k]=value(self.model[u].q[r,s,k])
        return q, x, v

class ADMM:
    def __init__(self, rho, q, r1, r2, la, z, c, g):#, growth):#pr,
        self.rho = rho
        self.q = q
        self.c = c
        self.g = g
        self.x = {}
        self.v = {}
        self.r_1 = r1
        self.r_2 = r2
        self.lbda = la
        self.znu = z
#       self.growth = growth

class Scenarios:
    def __init__(self, growth, U, pr ):
        self.growth = growth
        self.U = U
        self.pr = pr

def Example():
    S = [1,7,14,20,24]
    R = [2,11,13,19,21]
    K = [3,6,12,17,22] # [3,4, 5, 6, 8]
    U = set(range(1,21))
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
    I = [1]

    # utility parameters for charging facility choice
    b0 = {}
    for k in K:
        b0[k] = 0.0 # locational attractiveness
    b1 = 1 # travel time
    b2 = 0.0 # capacity
    b3 = 0.06 # price

    # income parameters for orgin and destination
    inc = {}
    for r in R:
        for s in S:
            inc[r,s]=1.0

    # travel demand between od
    d = {}
    for r in R:
        for s in S:
            d[r,s] = 100

    # EV adoption rate
    random.seed(1)
    growth = {}
    for u in U:
        growth[u] = random.uniform(1,3)
    
    # Probablity
    pr = {}
    for u in U:
        pr[u] = random.uniform(0,1) 
    pr_sum = sum(pr[u] for u in U)
    for u in U:
        pr[u] = pr[u]/pr_sum

    Scen = {}
    Scen = Scenarios(growth, U, pr)

    # travel cost function: use BPR function

    alpha = {}
    cap = {}
    tff = {}
    b = {}
    for (r,s) in A:
        alpha[r,s] = 4.0
        b[r,s] = 0.15

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

    # incidence matrix nodes
    mA_k={}
    for n in N:
        for (r,s) in A:
            if n == r:
                mA_k[n,r,s] = 1.0
            else:
                mA_k[n,r,s] = 0.0

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
    # incidence matrix
        for r in R:
            for k in K:
                mE[n,r,k]=0.0
                if n == r:
                    mE[n,r,k]=1.0
                if n == k:
                    mE[n,r,k]=-1.0
                if r == k:
                    mE[n,r,k]=0.0
        for k in K:
            for s in S:
                mE[n,k,s]=0.0
                if n == k:
                    mE[n,k,s]=1.0
                if n == s:
                    mE[n,k,s]=-1.0
                if s == k:
                    mE[n,k,s]=0.0

    # incidence matrix indicate nodes and charging
    mE_k={}
    for n in N:
        for r in R:
            for s in S:
                for k in K:
                    mE_k[n,r,s,k]=0.0
                    if n == k:
                        mE_k[n,r,s,k]=1.0

    # average energy demand between orgin and destination
    e = {}
    for r in R:
        for s in S:
            e[r,s]=1.0

    Ntw = Network(nodes=N,
                  links=A,
                  origin=R,
                  destination=S,
                  facility=K,
                  scenarios= Scen)
#(N, A, R, S, K, I, b0, b1, b2, b3, inc, tff, b, alpha, cap, d, mA, mA_k, mE, mE_k, e)
    costca = 0.100
    costcb = 170
    costga = 0.100
    costgb = 130

    Ntw.I = Investors(nodes=N,
                  links=A,
                  origin=R,
                  destination=S,
                  facility=K,
                  scenarios= Scen,
                  costca=costca,
                  costcb=costcb,
                  costga=costga,
                  costgb=costgb)

    Ntw.C = Consumers(nodes=N,
                  links=A,
                  origin=R,
                  destination=S,
                  facility=K,
                  scenarios=Scen,
                  b0=b0,
                  b1=b1,
                  b2=b2,
                  b3=b3,
                  inc=inc,
                  tff=tff,
                  b=b,
                  alpha=alpha,
                  cap=cap,
                  d=d,
                  mA=mA,
                  mA_k=mA_k,
                  mE=mE,
                  mE_k=mE_k,
                  e=e)

    return Ntw

if __name__ == "__main__":
    Ntw = Example()
    Algo = Ntw.init_ADMM()
    time_bq = {}
    start = time.time()
    SD_tol = 0.1
    ES_tol = 0.1
    print ('Stopping critieria %f' % ES_tol)
    Maxit = 1
    EE = {} # excess supply
    SS = {} # scenario difference
    RR = {} # prices
    RR[0] = Algo.rho
    EM = {} # max excess supply for all the locations
    SM = {} # max scenario difference for all the locations
    h_max = 4
    for u in Ntw.Scn.U:
        Ntw.I.model[u] = Ntw.I.profit_maximization_model(Algo,u)
        Ntw.C.model[u] = Ntw.C.traffic_problem_model(Algo,u)
    for iter in range(Maxit):
        time_bq[iter] = time.time()
        print ('Start iteration %i\t' % iter)
        # this are iteration for ADMM algorithm
        #Algo.q, Algo.x, Algo.v, Algo.c, Algo.g = Ntw.step_1(Algo)
        #Algo.znu, Algo.rho, Algo.lbda = Ntw.step_2(Algo)
       
        Algo.g, Algo.c, Algo.q, Algo.x, Algo.v, Algo.rho, Algo.znu = Ntw.centralized_problem()
        ES = Ntw.exsu(Algo)
        SD = Ntw.scen_diff(Algo)
        maxee = max( abs(ES[u, k]) for k in Ntw.K for u in Ntw.Scn.U)
        maxsd = max( abs(SD[u, k]) for k in Ntw.K for u in Ntw.Scn.U)
        RR[iter] = Algo.rho
        EE[iter] = ES
        SS[iter] = SD
        EM[iter] = maxee
        SM[iter] = maxsd
        print ('Iteration %i ES_{max} %f, SD_{max} %f'% (iter,maxee, maxsd))
        if iter > 5:
            if maxee <= ES_tol and maxsd <= SD_tol:
                if EM[iter-1] <= 1*ES_tol and EM[iter-2] <= 2*ES_tol and EM[iter-3] <= 3*ES_tol and SM[iter-1] <= 1*SD_tol and SM[iter-2] <= 2*SD_tol and SM[iter-3] <= 3*SD_tol:
                    print ('Equilibrium found!')
                    print ('Iteration %i ES_{max} %f, SD_{max} %f'% (iter,maxee, maxsd))
                    break
        if time.time()-start > h_max*3600:
            print ('More than %i hours!'%h_max)
            break
    end = time.time()
    el_time = end - start
    print ('Elapsed time %f' % el_time)
    aux_str=strftime("%Y%m%d_%H%M", localtime())
    pname='Results_'+aux_str
    os.system('mkdir '+pname)
    Ntw.write_evs(EvES=EE,
                    EvSD=SS,
                    EvR=RR,
                    res_path=pname)
    os.system('cp -r '+pname+'/ Results') # this line creates a copy of results so that we don't need to change code in plotting. TODA, a more efficient way is to write code in python to visulize the results directly.
