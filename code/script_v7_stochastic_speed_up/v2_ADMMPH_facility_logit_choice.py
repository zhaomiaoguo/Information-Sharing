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
import math
# from Transportaion_test_systems import import_matrix, transportation_network_topo
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

    def step_0(self, ADMM):
        print ('Starting Step 1:')
        print ('Solving investor profit maximization')
        c = {}
        g = {}
        for u in self.Scn.U:
            self.I.model[u].del_component('obj')
            def obj_rule_utmax(model):
                exp0 = - sum(ADMM.rho[u,k]*model.g[k] - (self.I.ca*model.c[k]**2 + self.I.cb*model.c[k]) - (self.I.ga*model.g[k]**2+self.I.gb*model.g[k]) for k in self.K)
                exp1 = ADMM.r_1/2.0 * sum( ( sum(ADMM.q[u,r,s,k] for r in self.R for s in self.S) - model.g[k])** 2 for k in self.K )
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

        print ('Calculate facility location choice:')
        q = {}
        for u in self.Scn.U:
            q_u = self.C.logit_facility_choice(ADMM, u)
            for r in self.R:
                for s in self.S:
                    for k in self.K:
                        q[u,r,s,k] = q_u[r,s,k]
        return q, c, g 

    def step_1(self, ADMM, ipoptObjs):
        print ('Starting Step 1:')
        print ('Solving investor profit maximization')
        c = {}
        g = {}
        for u in self.Scn.U:
            self.I.model[u].del_component('obj')
            
            def obj_rule_utmax(model):
                exp0 = - sum(ADMM.rho[u,k]*model.g[k] - (self.I.ca*model.c[k]**2 + self.I.cb*model.c[k]) - (self.I.ga*model.g[k]**2+self.I.gb*model.g[k]) for k in self.K)
                exp1 = ADMM.r_1/2.0 * sum( ( sum(ADMM.q[u,r,s,k] for r in self.R for s in self.S) - model.g[k])** 2 for k in self.K )
                exp2 = sum( ADMM.lbda[u,k]*(model.c[k] - ADMM.znu[k]) for k in self.K )
                exp3 = ADMM.r_2/2.0 * sum((model.c[k] - ADMM.znu[k])**2 for k in self.K)
                return exp0 + exp1 + exp2 + exp3
            
            def obj_rule_tap(model):
                exp0 = sum(self.C.tff[r,s]*(model.v[r,s]+(self.C.b[r,s]/(self.C.alpha[r,s]+1.0))*(model.v[r,s]**(self.C.alpha[r,s]+1))/(self.C.cap[r,s]**(self.C.alpha[r,s]))) for (r,s) in self.C.A)
                exp1 = 1.0/self.C.b1*( sum( model.q[r,s,k]*(pyolog(model.q[r,s,k]) - 1.0 -self.C.b0[k]) for k in self.K for s in self.S for r in self.R))
#                 exp2 = -sum(ADMM.rho[u,k] * sum(self.C.e[r,s]*model.q[r,s,k] for r in self.R for s in self.S) for k self.K)
                exp3 =  sum(((-sum(self.C.e[r,s]*(2*model.q[r, s, k] - ADMM.q[u,r,s,k]) for r in self.R for s in self.S ) + ADMM.g[u, k])/2)** 2 for k in self.K )
                exp2 = sum(ADMM.rho[u,k] * sum(self.C.e[r,s]*model.q[r,s,k] for r in self.R for s in self.S) for k in self.K)
                return self.C.b1/self.C.b3*(exp0 + exp1) - exp2 + ADMM.r_2/2.0 * exp3

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
        x1 = {}
        x2 = {}
                     
        for u in self.Scn.U:
            
            self.C.model[u].del_component('obj')
            def obj_rule_tap(model):
                exp0 = sum(self.C.tff[r,s]*(model.v[r,s]+(self.C.b[r,s]/(self.C.alpha[r,s]+1.0))*(model.v[r,s]**(self.C.alpha[r,s]+1))/(self.C.cap[r,s]**(self.C.alpha[r,s]))) for (r,s) in self.C.A)
                exp1 = 1.0/self.C.b1*( sum( model.q[r,s,k]*( model.q[r,s,k] - 1.0 -self.C.b0[k]) for k in self.K for s in self.S for r in self.R))
#                 exp2 = -sum(ADMM.rho[u,k] * sum(self.C.e[r,s]*model.q[r,s,k] for r in self.R for s in self.S) for k self.K)
                exp3 =  sum(((-sum(self.C.e[r,s]*(2*model.q[r, s, k] - ADMM.q[u,r,s,k]) for r in self.R for s in self.S ) + ADMM.g[u, k])/2)** 2 for k in self.K )
                exp2 = sum(ADMM.rho[u,k] * sum(self.C.e[r,s]*model.q[r,s,k] for r in self.R for s in self.S) for k in self.K)
                return self.C.b1/self.C.b3*(exp0 + exp1) - exp2 + ADMM.r_2/2.0 * exp3
            self.C.model[u].add_component('obj', Objective(rule=obj_rule_tap, sense=minimize))

            
            q_u, x_u, v_u, x1_u, x2_u = self.C.traffic_problem(u, ipoptObjs[u])

            for (i,j) in self.A:
                v[u,i,j] = v_u[i,j]
                self.C.model[u].v[i,j].set_value(v_u[i,j])
            for r in self.R:
                for s in self.S:
                    for k in self.K:
                        q[u,r,s,k] = q_u[r,s,k]
                        self.C.model[u].q[r,s,k].set_value(q_u[r,s,k])
                        for (i,j) in self.A:
                            x[u, i,j, r,s,k] = x_u[i,j,r,s,k]
                            x1[u, i,j, r,s,k] = x1_u[i,j,r,s,k]
                            x2[u, i,j, r,s,k] = x2_u[i,j,r,s,k]
                            self.C.model[u].x[i,j,r,s,k].set_value(x_u[i,j,r,s,k])
                            self.C.model[u].x1[i,j,r,s,k].set_value(x1_u[i,j,r,s,k])
                            self.C.model[u].x2[i,j,r,s,k].set_value(x2_u[i,j,r,s,k])
        return q, x, v, c, g, x1, x2

    def step_2(self, ADMM):
        print ('Starting Step 2:')
        print ('Updating z, rho, lambda')
        z_new = {}
        rho_new = {}
        la_new = {}
        for k in self.K:
            z_new[k] = sum(ADMM.c[u, k] * self.Scn.pr[u] for u in self.Scn.U) # modified  on Jul. 14
        for u in self.Scn.U:
            for k in self.K:
                rho_new[u,k] = ADMM.rho[u,k] + ADMM.r_2 *(sum(sum(self.C.e[r,s]*ADMM.q[u,r,s,k] for r in self.R) for s in self.S) - ADMM.g[u,k])
                la_new[u,k] = ADMM.lbda[u,k] + ADMM.r_1 * (ADMM.c[u,k]-z_new[k] )
        return z_new, rho_new, la_new

    def init_ADMM(self):
        rho_0={}
        y = 500
        for u in self.Scn.U:
            for k in self.K:
                rho_0[u,k]=y+2000*np.random.uniform(low = -1, high = 1)
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
                    r1=1,#ADMM
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
                ES[u,k] = ADMM.g[u,k] - sum( self.C.e[r,s]*ADMM.q[u,r,s,k] for r in self.R for s in self.S)
        return ES

    def scen_diff(self, ADMM):
        SD={}
        for u in self.Scn.U:
            for k in Ntw.K:
                SD[u, k] = ADMM.c[u,k] - ADMM.znu[k]
        return SD

    # all the output can be saved in this function.
    def write_evs(self,EvES,EvSD,EvR,c,g,res_path):
        f_exsu = open((res_path+'/Resulting_exsu.csv'),'w')
        f_scdi = open((res_path+'/Resulting_scendiff.csv'),'w')
        f_prices = open((res_path+'/Resulting_prices.csv'),'w')
        f_capacity = open((res_path+'/Resulting_capacity.csv'),'w')
        f_services = open((res_path+'/Resulting_services.csv'),'w')
        for u in self.Scn.U:
            for k in self.K:
                aux_p = 'Scenario%i,Node%i,'%(u,k)
                aux_s = 'Scenario%i,Node%i,'%(u,k)
                aux_sd = 'Scenario%i,Node%i,'%(u,k)
                aux_g = 'Scenario%i,Node%i,'%(u,k) 
                aux_c = 'Scenario%i,Node%i,'%(u,k)
                for nu in EvES:
                    aux_p = aux_p+('%f'%(EvR[nu][u,k]))+','
                    aux_s = aux_s+('%f'%(EvES[nu][u,k]))+','
                    aux_sd = aux_sd+('%f'%(EvSD[nu][u,k]))+','
                    aux_g += ('%f'%(g[nu][u,k]))+','
                    aux_c += ('%f'%(c[nu][u,k]))+','
                f_prices.write(aux_p+',\n')
                f_exsu.write(aux_s+',\n')
                f_scdi.write(aux_sd+',\n')
                f_capacity.write(aux_c + ',\n')
                f_services.write(aux_g + ',\n')
        f_prices.close()
        f_exsu.close()
        f_prices.close()
        f_capacity.close()
        f_services.close()

class Investors(Network):
    def __init__(self, nodes, links, origin, destination, facility, costca, costcb, costga, costgb, scenarios): #Locations,
#       self.K = Locations
        self.ca = costca
        self.cb = costcb
        self.ga = costga
        self.gb = costgb
        self.model = {}
        super().__init__(nodes, links, origin, destination, facility, scenarios)

    def centralized_problem(self, c_old):
        model = AbstractModel('System_minimization')
        model.c = Var(self.K, within=NonNegativeReals,initialize=1000.0)
        model.g = Var(self.K, within=NonNegativeReals,initialize=500.0)
        model.q = Var(self.R, self.S, self.K, within=NonNegativeReals, initialize=500.0)#=1e-1)
        model.v = Var(self.A, within=NonNegativeReals,initialize=0.0) # traffic flow
        model.x = Var(self.A, self.R, self.S, self.K, within=NonNegativeReals,initialize=0.0)   # traffic flow on link a associated with O-D pairs (r,s) and passing k
        # break up x into two pars: x1 and x2
        model.x1 = Var(self.A, self.R, self.S, self.K, within=NonNegativeReals,initialize=0.0)
        model.x2 = Var(self.A, self.R, self.S, self.K, within=NonNegativeReals,initialize=0.0)

        def obj_rule(model):
            exp0 = sum(self.ca*model.c[k]**2+self.cb*model.c[k] for k in self.K)
            exp1 = sum(self.ga*model.g[k]**2+self.gb*model.g[k] for k in self.K)
            exp20 = sum( self.tff[r,s]*(model.v[r,s]+(self.b[r,s]/(self.alpha[r,s]+1.0))*(model.v[r,s]**(self.alpha[r,s]+1))/(self.cap[r,s]**(self.alpha[r,s]))) for (r,s) in Ntw.A)
            exp21 = (sum( sum( sum( model.q[r,s,k]*( (model.q[r,s,k]) - 1.0 - self.b0[k] - self.b2*c_old[k]) for k in self.K) for s in self.S) for r in self.R))
            return exp0 + exp1 + self.b1/self.b3 * ( exp20+ 1.0/self.b1*exp21)
        model.obj = Objective(rule=obj_rule, sense=minimize)

        def market_clear_rule(model, k):
            return sum( self.e[r,s]*model.q[r, s, k] for r in self.R for s in self.S) - model.g[k] == 0
        model.market_clear = Constraint(self.K, rule = market_clear_rule)

        def capacity_rule(model, k):
            return model.g[k]-model.c[k] <= 0
        model.capacity = Constraint(self.K, rule=capacity_rule)

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
            return sum( model.q[r,s,k] for k in self.K) == Ntw.d[r,s]
        model.cond = Constraint(self.R, self.S, rule=cond_rule)
        optsolver = SolverFactory(quadsol)
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
        for k in self.K:
            rho_ans[k] = inst.dual[inst.market_clear[k]]
            print('Locational price at facility %i: %f' % (k, rho_ans[k]))
        x={}
        v={}
        for (i,j) in self.A:
            v[i,j]=value(inst.v[i,j])
            for r in self.R:
                for s in self.S:
                    for k in self.K:
                        x[i,j,r,s,k]=value(inst.x[i,j,r,s,k])
        q={}
        for r in self.R:
            for s in self.S:
                for k in self.K:
                    q[r,s,k]=value(inst.q[r,s,k])
        g={}
        c={}
        for k in self.K:
            g[k]=value(inst.g[k])
            c[k]=value(inst.c[k])
            print ('Investment capacity at facility %i: %f' % (k, c[k]))
        return g, c, q, x, v, rho_ans

    def profit_maximization_model(self, ADMM, u):
        model = AbstractModel('profit_max')
        model.c = Var(self.K, within=NonNegativeReals,initialize=1e5)
        model.g = Var(self.K, within=NonNegativeReals,initialize=1e5)
        def obj_rule(model):
            exp0 = -sum(ADMM.rho[u,k]*model.g[k] - (self.ca*model.c[k]**2 + self.cb*model.c[k]) - (self.ga*model.g[k]**2+self.gb*model.g[k]) for k in self.K)
            exp1 = ADMM.r_1/2.0 * sum( ( sum(ADMM.q[u,r,s,k] for r in self.R for s in self.S) - model.g[k])** 2 for k in self.K )
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
        model = ConcreteModel('CDA')
        model.q = Var(self.R, self.S, self.K, within=NonNegativeReals, initialize=500.0)#=1e-1)
        model.v = Var(self.A, within=NonNegativeReals,initialize=0.0) # traffic flow
        model.x = Var(self.A, self.R, self.S, self.K, within=NonNegativeReals,initialize=0.0)   # traffic flow on link a associated with O-D pairs (r,s) and passing k
        # break up x into two pars: x1 and x2
        model.x1 = Var(self.A, self.R, self.S, self.K, within=NonNegativeReals,initialize=0.0)
        model.x2 = Var(self.A, self.R, self.S, self.K, within=NonNegativeReals,initialize=0.0)

        def objrule(model):
            exp0 = sum(self.tff[r,s]*(model.v[r,s]+(self.b[r,s]/(self.alpha[r,s]+1.0))*(model.v[r,s]**(self.alpha[r,s]+1))/(self.cap[r,s]**(self.alpha[r,s]))) for (r,s) in self.A)
            exp1 = 1.0/self.b1*( sum( model.q[r,s,k]*( (model.q[r,s,k]) - 1.0 + self.b3*ADMM.rho[u,k] -self.b0[k] - self.b2*ADMM.c[u,k] ) for k in self.K for s in self.S for r in self.R))
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
        #inst = model.create_instance()

        # Ipopt bound multipliers (obtained from solution) 
        model.ipopt_zL_out = Suffix(direction=Suffix.IMPORT) 
        model.ipopt_zU_out = Suffix(direction=Suffix.IMPORT)
        # Ipopt bound multipliers (sent to solver)
        model.ipopt_zL_in = Suffix(direction=Suffix.EXPORT)
        model.ipopt_zU_in = Suffix(direction=Suffix.EXPORT)
        # Obtain dual solutions from first solve and send to warm start
        model.dual = Suffix(direction=Suffix.IMPORT_EXPORT)

        return model 

    def traffic_problem(self, u, optsolver):
        x={}
        x1={}
        x2={}
        v={}
        q={}
        warm_start =True 
        #optsolver = SolverFactory(nlsol)
        
        #optsolver = ipoptObj 

        # new problem
        if (warm_start == True):
            self.model[u].ipopt_zL_in.update(self.model[u].ipopt_zL_out) 
            self.model[u].ipopt_zU_in.update(self.model[u].ipopt_zU_out) 
            self.model[u].dual.update(self.model[u].dual) 
            optsolver.options['warm_start_init_point'] = 'yes' 
            optsolver.options['warm_start_bound_push'] = 1e-6 
            optsolver.options['warm_start_mult_bound_push'] = 1e-6 
            optsolver.options['mu_init'] = 1e-6
        results = optsolver.solve(self.model[u], tee = True)
        #self.model[u].solutions.load_from(results)
        for (i,j) in self.A:
            v[i,j]=value(self.model[u].v[i,j])
            for r in self.R:
                for s in self.S:
                    for k in self.K:
                        x[i,j,r,s,k]=value(self.model[u].x[i,j,r,s,k])
                        x1[i,j,r,s,k]=value(self.model[u].x1[i,j,r,s,k])
                        x2[i,j,r,s,k]=value(self.model[u].x2[i,j,r,s,k])
        for r in self.R:
            for s in self.S:
                for k in self.K:
                    q[r,s,k]=value(self.model[u].q[r,s,k])
        return q, x, v, x1, x2

    def logit_facility_choice(self, ADMM, u):
        q={}
        for r in self.R:
            for s in self.S:
                for k in self.K:
                    #u directly calculated q using logit model may lead to numerical issue because denominator is zero.
                    q[r,s,k]= self.d[r,s] * self.Scn.growth[u] *(np.exp(self.b0[k] - self.b1*ADMM.tt[u][r,s,k] + self.b2*ADMM.c[u,k] - self.b3*ADMM.rho[u,k]))/(sum(np.exp(self.b0[k1] - self.b1*ADMM.tt[u][r,s, k1] + self.b2*ADMM.c[u,k1] - self.b3*ADMM.rho[u,k1]) for k1 in self.K))
              #      q[r,s,k]= self.d[r,s] * self.Scn.growth[u] /(sum(np.exp(self.b0[k1] - self.b1*ADMM.tt[u][r,s, k1] + self.b2*ADMM.c[u,k1] - self.b3*ADMM.rho[u,k1] - (self.b0[k] - self.b1*ADMM.tt[u][r,s, k] + self.b2*ADMM.c[u,k] - self.b3*ADMM.rho[u,k])  ) for k1 in self.K))
    
        return q 

    def calculate_travel_time(self, u):
        tt={}
        optsolver = SolverFactory(nlsol)
        results = optsolver.solve(self.model[u], tee = True)
        for r in self.R:
            for s in self.S:
                for k in self.K:
                    tt[r,s,k] = ((self.model[u].dual[self.model[u].conc2[r,s,k, r]]-self.model[u].dual[self.model[u].conc2[r,s,k, k]]) + (self.model[u].dual[self.model[u].conc3[r,s,k, k]] - self.model[u].dual[self.model[u].conc3[r,s,k, s]]))/(self.b1/self.b3)
                   
                   # debug code: see if the calculation of travel time is correct
                    #if r == 21 and (s == 20):
                   #     print(self.tff[r,k]*(1+(self.b[r,k]*(value(self.model[u].v[r,k])**(self.alpha[r,k]))/(self.cap[r,k]**(self.alpha[r,k])))))
                   #     print(value(self.model[u].x[r,k,r,s,k]))
                   #     print(value(self.model[u].x[k,s,r,s,k]))
                   #     print(value(self.model[u].q[r,s,k]))
                   #     print(self.model[u].dual[self.model[u].conc2[r,s,k, r]])
                   #     print(self.model[u].dual[self.model[u].conc2[r,s,k, k]])
                   #     print(self.model[u].dual[self.model[u].conc3[r,s,k, k]])
                   #     print(self.model[u].dual[self.model[u].conc3[r,s,k, s]])
                   #     print(self.model[u].dual[self.model[u].cond[r,s]])
                   #     for (i,j) in self.A:
                   #         print(i, j, '  ', self.model[u].dual[self.model[u].conb2[i,j, r, s, k]])


                    print('travel_time: (%i, %i, %i) at scenario %i: %f' % (r,s,k,u, tt[r,s,k]))
                

        return tt 

class ADMM:
    def __init__(self, rho, q, r1, r2, la, z, c, g):#, growth):#pr,
        self.rho = rho
        self.q = q
        self.c = c
        self.g = g
        self.x = {}
        self.x1 = {}
        self.x2 = {}
        self.v = {}
        self.r_1 = r1
        self.r_2 = r2
        self.lbda = la
        self.znu = z
        self.tt = {}
#       self.growth = growth

class Scenarios:
    def __init__(self, growth, U, pr ):
        self.growth = growth
        self.U = U
        self.pr = pr

def Example(identical_scen, congestion):
    S = [1,7,14,20,24]
    R = [2,11,13,19,21]
    K = [3,6,12,17,22]
    U = set(range(1,6))
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
    #b1 = 1 # travel time
    b1 = 0.001 # travel time
    b2 = 0.0 # capacity
    b3 = 0.001 # price

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
        if identical_scen:
            growth[u] = 1 
        else:
            growth[u] = random.uniform(1,1.2)
        print('growth at scen ', u, ': ', growth[u])
    
    # Probablity
    pr = {}
    for u in U:
        pr[u] = random.uniform(0,1) 
    pr_sum = sum(pr[u] for u in U)
    for u in U:
        if identical_scen:
            pr[u] = 1/len(U)
        else:
            pr[u] = pr[u]/pr_sum

    Scen = {}
    Scen = Scenarios(growth, U, pr)
    # travel cost function: use BPR function

    alpha = {}
    cap = {}
    tff = {}
    b = {}
    for (r,s) in A:
        alpha[r,s] = 1.0
        if congestion:
            b[r,s] = 0.15
        else:
            b[r,s] = 0

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
    
    for (i,j) in A:
        cap[i,j] = cap[i,j]*1

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
    for (i,j) in A:
        tff[i,j] = tff[i,j] *1
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
            # need to double check the objective function has e[]

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

def Example_Anaheim(identical_scen, congestion):
    S = [2,34,21,29,41]
    R = [10,25,28,39,24]
    K = [5,9,13,12,16]
    U = set(range(1,6))
    N, A, cap, tff = transportation_network_topo('Anaheim')

    I = []

    # utility parameters for charging facility choice
    b0 = {}
    for k in K:
        b0[k] = 0.0 # locational attractiveness
    #b1 = 1 # travel time
    b1 = 0.001
    b2 = 0.0 # capacity
    b3 = 0.001 # price

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
        if identical_scen:
            growth[u] = 1
        else:
            growth[u] = random.uniform(1,1.2)
        print('growth at scen ', u, ': ', growth[u])

    # Probablity
    pr = {}
    for u in U:
        pr[u] = random.uniform(0,1)
    pr_sum = sum(pr[u] for u in U)
    for u in U:
        if identical_scen:
            pr[u] = 1/len(U)
        else:
            pr[u] = pr[u]/pr_sum

    Scen = {}
    Scen = Scenarios(growth, U, pr)
    # travel cost function: use BPR function

    alpha = {}
    b = {}

    for (r,s) in A:
        alpha[r,s] = 4.0
        if congestion:
            b[r,s] = 0.15
        else:
            b[r,s] = 0

    for (i,j) in A:
        cap[i,j] = cap[i,j] * 1

    for (i,j) in A:
        tff[i,j] = tff[i,j] *1
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
            # need to double check the objective function has e[]

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


def Example_6node(identical_scen, congestion):
    R = [1,4]
    S = [3,6]
    K = [2,5]
    U = set(range(1,2))
    N = set(range(1,7))
    A = [(1,2),(2,3),(1,4),(4,5),(5,6),(6,3)]
    I = [1]
    
    b0 = {}
    for k in K:
        b0[k] = 1 # locational attractiveness
    #b1 = 1 # travel time
    b1 = 0.01# travel time
    b2 = 1 # capacity
    b3 = 0.1 # price
    
    # income parameters for orgin and destination
    inc = {}
    for r in R:
        for s in S:
            inc[r,s]=1.0

    # travel demand between od
    d = {}
#     d[1,3] = 20
#     d[4,6] = 20
    for r in R:
        for s in S:
            d[r,s] = 20

    # EV adoption rate
    random.seed(1)
    growth = {}
    for u in U:
        if identical_scen:
            growth[u] = 1 
        else:
            growth[u] = random.uniform(1,1.2)
        print('growth at scen ', u, ': ', growth[u])
    
    # Probablity
    pr = {}
    for u in U:
        pr[u] = random.uniform(0,1) 
    pr_sum = sum(pr[u] for u in U)
    for u in U:
        if identical_scen:
            pr[u] = 1/len(U)
        else:
            pr[u] = pr[u]/pr_sum

    Scen = {}
    Scen = Scenarios(growth, U, pr)
    # travel cost function: use BPR function

    alpha = {}
    cap = {}
    tff = {}
    b = {}
    for (r,s) in A:
        alpha[r,s] = 2.0
        if congestion:
            b[r,s] = 0.15
        else:
            b[r,s] = 0
    cap={}
    tff={}
    for (i,j) in A:        
        cap[i,j]=700
        tff[i,j]=8

    
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
            # need to double check the objective function has e[]

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

def Example_3node(identical_scen, congestion):
    R = [1]
    S = [3]
    K = [2]
    U = set(range(1,2))
    N = set(range(1,4))
    A = [(1,2),(2,3),(3,1)]
    I = [1]
    
    b0 = {}
    for k in K:
        b0[k] = 0 # locational attractiveness
    #b1 = 1 # travel time
    b1 = 0.001# travel time
    b2 = 1 # capacity
    b3 = 0.1 # price
    
    # income parameters for orgin and destination
    inc = {}
    for r in R:
        for s in S:
            inc[r,s]=1.0

    # travel demand between od
    d = {}
#     d[1,3] = 20
#     d[4,6] = 20
    for r in R:
        for s in S:
            d[r,s] = 20

    # EV adoption rate
    random.seed(1)
    growth = {}
    for u in U:
        if identical_scen:
            growth[u] = 1 
        else:
            growth[u] = random.uniform(1,1.2)
        print('growth at scen ', u, ': ', growth[u])
    
    # Probablity
    pr = {}
    for u in U:
        pr[u] = random.uniform(0,1) 
    pr_sum = sum(pr[u] for u in U)
    for u in U:
        if identical_scen:
            pr[u] = 1/len(U)
        else:
            pr[u] = pr[u]/pr_sum

    Scen = {}
    Scen = Scenarios(growth, U, pr)
    # travel cost function: use BPR function

    alpha = {}
    cap = {}
    tff = {}
    b = {}
    for (r,s) in A:
        alpha[r,s] = 2.0
        if congestion:
            b[r,s] = 0.15
        else:
            b[r,s] = 0
    cap={}
    tff={}
    for (i,j) in A:        
        cap[i,j]=700
        tff[i,j]=8

    
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
            # need to double check the objective function has e[]

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
    # this is the main function
    congestion = True 
    identical_scen = False 
#     Ntw = Example(identical_scen, congestion)
    Ntw = Example_6node(identical_scen, congestion)
    Algo = Ntw.init_ADMM()
    time_bq = {}
    start = time.time()
    SD_tol = 1
    ES_tol = 1
    print ('Stopping critieria %f' % ES_tol)
    Maxit = 50
    Pre_Iter = 10
    EE = {} # excess supply
    SS = {} # scenario difference
    RR = {} # prices
    RR[0] = Algo.rho
    CC={}
    GG={}
    EM = {} # max excess supply for all the locations
    SM = {} # max scenario difference for all the locations
    ipoptObjs = {} # ipopt solver objects for each scenarios
    h_max = 4 # maximum hours for running this program
    for u in Ntw.Scn.U:
        Ntw.I.model[u] = Ntw.I.profit_maximization_model(Algo,u)
        Ntw.C.model[u] = Ntw.C.traffic_problem_model(Algo,u)
        ipoptObjs[u] = SolverFactory(nlsol)
        Algo.tt[u] = Ntw.C.calculate_travel_time(u)
    for iter in range(Maxit):
        time_bq[iter] = time.time()
        print ('Start iteration %i\t' % iter)

        if iter < Pre_Iter:
            #Algo.r_1 = min(1.05*Algo.r_1, 5)
            #Algo.r_2 = min(1.05*Algo.r_2, 5)
            Algo.q, Algo.c, Algo.g = Ntw.step_0(Algo)
        else:
            Algo.q, Algo.x, Algo.v, Algo.c, Algo.g, Algo.x1, Algo.x2 = Ntw.step_1(Algo, ipoptObjs)
        Algo.znu, Algo.rho, Algo.lbda = Ntw.step_2(Algo)
        ES = Ntw.exsu(Algo)
        SD = Ntw.scen_diff(Algo)
        maxee = max( abs(ES[u, k]) for k in Ntw.K for u in Ntw.Scn.U)
        maxsd = max( abs(SD[u, k]) for k in Ntw.K for u in Ntw.Scn.U)
        RR[iter] = Algo.rho
        CC[iter] = Algo.c
        GG[iter] = Algo.g
        EE[iter] = ES
        SS[iter] = SD
        EM[iter] = maxee
        SM[iter] = maxsd
        #Algo.r_1, Algo.r_2 = calculate_r(RR, EE, SS)
        print ('Iteration %i ES_{max} %f, SD_{max} %f'% (iter,maxee, maxsd))
        if iter > 5:
            if maxee <= ES_tol and maxsd <= SD_tol:
                if EM[iter-1] <= 1*ES_tol and EM[iter-2] <= 2*ES_tol and EM[iter-3] <= 3*ES_tol and SM[iter-1] <= 1*SD_tol and SM[iter-2] <= 2*SD_tol and SM[iter-3] <= 3*SD_tol:
                    if iter > Pre_Iter or (not congestion):
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
                    c=CC,
                    g=GG,
                    res_path=pname)
    os.system('cp -r '+pname+'/ Results') # this line creates a copy of results so that we don't need to change code in plotting. TODA, a more efficient way is to write code in python to visulize the results directly.
