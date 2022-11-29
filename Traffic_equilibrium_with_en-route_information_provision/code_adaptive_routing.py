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
from time import perf_counter, strftime,localtime
from pyomo.environ import log as pyolog
from pyomo.environ import exp as pyoexp
#from pyomo.bilevel import *
from collections import defaultdict
from pyomo.core.base.symbolic import differentiate
from pyomo.core.expr.current import identify_variables
from pyomo.core.expr.current import evaluate_expression
import networkx as nx
import math

quadsol = 'cplex'
nlpsol = 'ipopt'
minlpsol = 'scipampl'
mipsol = 'cplex'
bisol = 'bilevel_blp_global'
infinity = float('inf')

class Lib:
	def Network_topo(N, A):
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
		return NodesOut, NodesIn

	def Network_prev_nodes(P, A, P_set):
		# nodes on a path before reaching link a, including a_i
		N_ap = {}
		for p in P_set:
			for (i,j) in A:
				N_ap[i,j,p] = []
			for k in range(len(P[p])):
				for m in range(k+1):# including a_i
					N_ap[P[p][k][0], P[p][k][1],p].append(P[p][m][0])
		#print(N_ap)
		return N_ap

	def Network_incidence_link_path(A,P, P_set):
	# link-path incident scaler
		d = {}
		for p in P_set:
			for (i,j) in A:
				d[i,j,p] = 0
			for k in range(len(P[p])):
				d[P[p][k][0], P[p][k][1],p] = 1
		#print(d)
		return d

	def Network_incidence_path_od(P,R,S, P_set):
	# link-path incident scaler
		d = {} # p, r, s
		for p in P_set:
			for r in R:
				for s in S:
					d[p,r,s] = 0
			d[p,P[p][0][0],P[p][len(P[p])-1][1]]= 1
		#print(d)
		return d

	def get_derivatives(expr):
		obj = evaluate_expression(expr)
		varList = list( identify_variables(expr) )
		varList_num = []
		firstDerivs = differentiate(expr, wrt_list=varList)
		# Note this calculates d^2/dx_i^2; if you want the full Hessian matrix
		#   ( \delta^2/{\delta x_i \delta x_j} ) replace "wrt=v" with "wrt_list=varList"
		secondDerivs = [ differentiate(firstDerivs[i], wrt_list=varList) for i,v in enumerate(varList) ]
		var_i = -1
		firstDerivs_num = {}
		secondDerivs_num = defaultdict(dict)
		for var in varList:
			var_i += 1
			varList_num.append(value(var))
			firstDerivs_num[var_i] = evaluate_expression(firstDerivs[var_i])
			#print("1st Derivs: %s (%f)" % (firstDerivs[var_i],firstDerivs_num[var_i]))
			var_j = -1
			for var2 in varList:
				#print("1st Var: %s (%f), 2nd Var: %s (%f)" % (var.name, value(var), var2.name, value(var2)))
				var_j += 1
				secondDerivs_num[var_i][var_j] = evaluate_expression(secondDerivs[var_i][var_j])
				#print("2nd Derivs: %s (%f)" % (secondDerivs[var_i][var_j], secondDerivs_num[var_i][var_j]))
		return obj, varList, varList_num, firstDerivs_num, secondDerivs_num

class IO():
	v = {}
	def write_file(f1, dict1):
		for i in dict1.keys():
			# TODO need to figure out what's the format for keys...
			f1.write("%i, %f," % (i, dict1[i]))
		f1.write("\n")

	

	def Orlando():
		#S = [1, 3, 11,12,17,18]
		#R = [3, 11,17,18]
		S = [3,11,17]
		#S = [3]   # destination
		R = [3,11,17]
		#R = [17]   # origin 
		#TODO can we change to list(range(4))
		N = set(range(1,19))
		A = [(1,2),(1,3),
				(2,1),(2,4),
				(3,1),(3,4),(3,5),
				(4,3),(4,2),(4,6),(4,7),
				(5,3),(5,11),(5,10),
				(6,4),(6,10),(6,7),(6,9),
				(7,4),(7,6),(7,8),
				(8,7),(8,9),(8,16),
				(9,6),(9,8),(9,14),(9,15),
				(10,5),(10,13),(10,14),(10,6),
				(11,5),(11,12),
				(12,11),(12,13),
				(13,12),(13,10),(13,14),
				(14,13),(14,10),(14,9),(14,15),(14,17),
				(15,9),(15,14),(15,16),(15,18),
				(16,8),(16,15),
				(17,14),(17,18),
				(18,15),(18,17)]
		# utility parameters for charging facility choice

		def all_path(G, R, S, max_len):
			i = 0
			P ={}
			for r in R:
				for s in S:
					paths = nx.all_simple_paths(G, source=r, target=s, cutoff=max_len)
					for path in map(nx.utils.pairwise, paths):
						P[i] = list(path)
						print("Path ", i, " : ",  P[i])
						i = i+1
			return P

		G = nx.Graph()
		G.add_nodes_from(list(N))
		G.add_edges_from(A)
		P = all_path(G, R,S, 7)  #6.46
		Q = {} # slope of demand function: demand = demand - demand_slope*rho

		for s in S:
			for r in R:
				Q[r,s]=0
		Q[17,11]=1455
		Q[11,17]=854
		Q[17,3]=9000
		Q[3,17]=8610
		Q[11,3]=5186
		Q[3,11]=4029

		# travel cost function: use BPR function: t = tff(c+b(v/cap)**alpha)
		alpha = {}
		cap = {}
		tff = {}
		b = {}
		c = {}
		l_n = {}
		factor = 0.24

		for (r,s) in A:
				#tff[r,s] = 1
				alpha[r,s] = 1 #4.0
				b[r,s] = 0.15
				#cap[r,s] = 1
				c[r,s] = 1
				l_n[r,s] =  1

		def inti_lane():
			l_n[1,2] = 3
			l_n[2,1] = 3
			l_n[2,4] = 2
			l_n[4,2] = 2
			l_n[3,4] = 3
			l_n[4,3] = 3
			l_n[1,3] = 3
			l_n[3,1] = 3
			l_n[3,5] = 3
			l_n[5,3] = 3
			l_n[4,6] = 3
			l_n[6,4] = 3
			l_n[4,7] = 2
			l_n[7,4] = 2
			l_n[6,7] = 3
			l_n[7,6] = 3
			l_n[7,8] = 2
			l_n[8,7] = 2
			l_n[8,9] = 3
			l_n[9,8] = 3
			l_n[6,9] = 3
			l_n[9,6] = 3
			l_n[9,15] = 3
			l_n[15,9] = 3
			l_n[8,16] = 2
			l_n[16,8] = 3
			l_n[15,16] = 3
			l_n[16,15] = 3
			l_n[10,6] = 3
			l_n[6,10] = 3
			l_n[5,10] = 4
			l_n[10,5] = 4
			l_n[5,11] = 3
			l_n[11,5] = 3
			l_n[11,12] = 3
			l_n[12,11] = 3
			l_n[12,13] = 3
			l_n[13,12] = 3
			l_n[13,10] = 3
			l_n[10,13] = 3
			l_n[10,14] = 4
			l_n[14,10] = 4
			l_n[13,14] = 3
			l_n[14,13] = 3
			l_n[9,14] = 3
			l_n[14,9] = 3
			l_n[14,15] = 3
			l_n[15,14] = 3
			l_n[15,18] = 2
			l_n[18,15] = 2
			l_n[17,18] = 2
			l_n[18,17] = 2
			l_n[14,17] = 3
			l_n[17,14] = 3

		inti_lane()




		# capacity for the mixed traffic is not bounded by the real traffic calculated here.
		# this is because the CAV allows the capacity of the road to increase (smaller headway)
		def init_cap():
			cap[3,5] = 18400*factor
			cap[5,3] = 18400*factor
			cap[1,3] = 18400*factor
			cap[3,1] = 18400*factor
			cap[1,2] = 18400*factor
			cap[2,1] = 18400*factor
			cap[2,4] = 8800*factor
			cap[4,2] = 8800*factor
			cap[4,3] = 13200*factor
			cap[3,4] = 13200*factor
			cap[4,7] = 8800*factor
			cap[7,4] = 8800*factor
			cap[4,6] = 9200*factor
			cap[6,4] = 9200*factor
			cap[5,11] = 18400*factor
			cap[11,5] = 18400*factor
			cap[11,12] = 17600*factor
			cap[12,11] = 17600*factor
			cap[12,13] = 16800*factor
			cap[13,12] = 16800*factor
			cap[13,10] = 19200*factor
			cap[10,13] = 19200*factor
			cap[13,14] = 16800*factor
			cap[14,13] = 16800*factor
			cap[10,14] = 19200*factor
			cap[14,10] = 19200*factor
			cap[5,10] = 19200*factor
			cap[10,5] = 19200*factor
			cap[6,10] = 14400*factor
			cap[10,6] = 14400*factor
			cap[14,15] = 14400*factor
			cap[15,14] = 14400*factor
			cap[9,14] = 12000*factor
			cap[14,9] = 12000*factor
			cap[6,9] = 9600*factor
			cap[9,6] = 9600*factor
			cap[9,15] = 12000*factor
			cap[15,9] = 12000*factor
			cap[15,18] = 9600*factor
			cap[18,15] = 9600*factor
			cap[17,18] = 10000*factor
			cap[18,17] = 10000*factor
			cap[6,7] = 14400*factor
			cap[7,6] = 14400*factor
			cap[7,8] = 8800*factor
			cap[8,7] = 8800*factor
			cap[8,16] = 8800*factor
			cap[16,8] = 8800*factor
			cap[9,8] = 12000*factor
			cap[8,9] = 12000*factor
			cap[14,17] = 18000*factor
			cap[17,14] = 18000*factor
			cap[15,16] = 12000*factor
			cap[16,15] = 12000*factor

			for (r,s) in A:
				#cap[r,s] = cap[r,s]* 0.5
				cap[r,s] = cap[r,s]*0.5/ factor
		init_cap()


		def init_tff():
			tff[3,5] = 360
			tff[5,3] = 360
			tff[1,3] = 60
			tff[3,1] = 60
			tff[1,2] = 120
			tff[2,1] = 120
			tff[2,4] = 60
			tff[4,2] = 60
			tff[4,3] = 120
			tff[3,4] = 120
			tff[4,7] = 480
			tff[7,4] = 480
			tff[4,6] = 300
			tff[6,4] = 300
			tff[5,11] = 180
			tff[11,5] = 180
			tff[11,12] = 60
			tff[12,11] = 60
			tff[12,13] = 120
			tff[13,12] = 120
			tff[13,10] = 180
			tff[10,13] = 180
			tff[13,14] = 120
			tff[14,13] = 120
			tff[10,14] = 60
			tff[14,10] = 60
			tff[5,10] = 180
			tff[10,5] = 180
			tff[6,10] = 240
			tff[10,6] = 240
			tff[14,15] = 240
			tff[15,14] = 240
			tff[9,14] = 480
			tff[14,9] = 480
			tff[6,9] = 60
			tff[9,6] = 60
			tff[9,15] = 60
			tff[15,9] = 60
			tff[15,18] = 240
			tff[18,15] = 240
			tff[17,18] = 420
			tff[18,17] = 420
			tff[6,7] = 180
			tff[7,6] = 180
			tff[7,8] = 120
			tff[8,7] = 120
			tff[8,16] = 180
			tff[16,8] = 180
			tff[9,8] = 240
			tff[8,9] = 240
			tff[14,17] = 360
			tff[17,14] = 360
			tff[15,16] = 120
			tff[16,15] = 120
			
		init_tff()

		a0 = {}
		a1 = {}
		n = {}

		for (r,s) in A:
			a0[r,s] = -10.9543
			a1[r,s] = 0.000006
			n[r,s] = 3

		# a0
		a0[12,13] = -5.2844
		a0[13,12] = -5.2844
		a0[13,14] = -5.6073
		a0[14,13] = -5.6073
		a0[13,10] = -4.7818
		a0[10,13] = -4.7818
		a0[8,16] = -4.4609
		a0[16,8] = -4.4609
		a0[15,16] = -5.3167
		a0[16,15] = -5.3167
		a0[14,15] = -4.4551
		a0[15,14] = -4.4551
		a0[9,15] = -6.4915
		a0[15,9] = -6.4915
		a0[14,17] = -3.9633
		a0[17,14] = -3.9633
		a0[1,2] = -5.8835
		a0[2,1] = -5.8835
		a0[2,4] = -8.4944
		a0[4,2] = -8.4944
		a0[5,11] = -4.6374
		a0[11,5] = -4.6374
		a0[1,3] = -5.8791
		a0[3,1] = -5.8791
		a0[3,5] = -4.2791
		a0[5,3] = -4.2791
		a0[11,12] = -5.8288
		a0[12,11] = -5.8288
		a0[9,14] = -4.7099
		a0[14,9] = -4.7099
		a0[4,6] = -5.8051
		a0[6,4] = -5.8051
		a0[4,3] = -7.6803
		a0[3,4] = -7.6803
		a0[10,14] = -6.1302
		a0[14,10] = -6.1302
		a0[5,10] = -5.1920
		a0[10,5] = -5.1920

		# a1
		a1[12,13] = 0.0004
		a1[13,12] = 0.0004
		a1[13,14] = 0.0003
		a1[14,13] = 0.0003
		a1[13,10] = 0.0002 
		a1[10,13] = 0.0002
		a1[8,16] = 0.0006
		a1[16,8] = 0.0006
		a1[15,16] = 0.0006
		a1[16,15] = 0.0006 
		a1[14,15] = 0.0006
		a1[15,14] = 0.0006
		a1[9,15] = -0.0001
		a1[15,9] = -0.0001
		a1[14,17] = 0.0012
		a1[17,14] = 0.0012
		a1[1,2] = 0.0007
		a1[2,1] = 0.0007
		a1[2,4] = 0.0014
		a1[4,2] = 0.0014
		a1[5,11] = 0.0004
		a1[11,5] = 0.0004
		a1[1,3] = 0.0005
		a1[3,1] = 0.0005
		a1[3,5] = 0.0003
		a1[5,3] = 0.0003
		a1[11,12] = 0.0003
		a1[12,11] = 0.0003
		a1[9,14] = 0.0014
		a1[14,9] = 0.0014
		a1[4,6] = 0.0011
		a1[6,4] = 0.0011
		a1[4,3] = 0.0029
		a1[3,4] = 0.0029
		a1[10,14] = 0.0079
		a1[14,10] = 0.0079
		a1[5,10] = 0.0006
		a1[10,5] = 0.0006

		return Network_transport(N, A, R, S, P, Q, alpha,  b, cap, tff, c,a0,a1)
		

	def scen_four_node_safer_sim(ntw_1):
		Xi = [1,2]
		prob = {} # probability of each sceanrio
		link = [(2,3)]
		link_tff ={}
		link_c ={}
		link_b ={}
		link_alpha ={}
		link_cap ={}
		link_a1 = {}
		link_a0 = {}

		for xi in Xi:
			prob[xi] = 1/len(Xi)

		#prob[1] = 0.5
		#prob[2] = 0.5
		# actual
		for (i,j) in link:
			link_tff[i,j,1] = ntw_t.tff[i,j]
			link_tff[i,j,2] = ntw_t.tff[i,j]
			link_c[i,j,1] = ntw_t.c[i,j] # 1, 2 link nodes, 1 scenario id
			link_c[i,j,2] = ntw_t.c[i,j] # 1, 2 link nodes, 1 scenario id
			link_b[i,j,1] = ntw_t.b[i,j] # 1, 2 link nodes, 1 scenario id
			link_b[i,j,2] = ntw_t.b[i,j]
			link_alpha[i,j,1] = ntw_t.alpha[i,j]
			link_alpha[i,j,2] = ntw_t.alpha[i,j]
			link_cap[i,j,1] = ntw_t.cap[i,j]*0.1 
			link_cap[i,j,2] = ntw_t.cap[i,j]*0.5
			link_a0[i,j,1] = ntw_t.a0[i,j]
			link_a0[i,j,2] = ntw_t.a0[i,j]
			link_a1[i,j,1] = ntw_t.a0[i,j]
			link_a1[i,j,2] = ntw_t.a0[i,j]



		
		
		return Scen(Xi, prob,link, link_tff, link_c, link_b, link_cap, link_alpha)

	

	

	def scen_Orlando(ntw_t):
		Xi = [1,2]
		prob = {} # probability of each sceanrio
		link = [(17,14),(14,17),(3,5),(5,3)]
		link_tff ={}
		link_c ={}
		link_b ={}
		link_alpha ={}
		link_cap ={}
		link_a1 = {}
		link_a0 = {}

		for xi in Xi:
			prob[xi] = 1/len(Xi)

		#prob[1] = 0.5
		#prob[2] = 0.5
		# actual
		for (i,j) in link:
			link_tff[i,j,1] = ntw_t.tff[i,j]
			link_tff[i,j,2] = ntw_t.tff[i,j]
			link_c[i,j,1] = ntw_t.c[i,j] # 1, 2 link nodes, 1 scenario id
			link_c[i,j,2] = ntw_t.c[i,j] # 1, 2 link nodes, 1 scenario id
			link_b[i,j,1] = ntw_t.b[i,j] # 1, 2 link nodes, 1 scenario id
			link_b[i,j,2] = ntw_t.b[i,j]
			link_alpha[i,j,1] = ntw_t.alpha[i,j]
			link_alpha[i,j,2] = ntw_t.alpha[i,j]
			link_cap[i,j,1] = ntw_t.cap[i,j]*0.1 
			link_cap[i,j,2] = ntw_t.cap[i,j]*0.5
			link_a0[i,j,1] = ntw_t.a0[i,j]
			link_a0[i,j,2] = ntw_t.a0[i,j]
			link_a1[i,j,1] = ntw_t.a0[i,j]
			link_a1[i,j,2] = ntw_t.a0[i,j]
		return Scen(Xi, prob,link, link_tff, link_c, link_b, link_cap, link_alpha)

	



	def four_node_safer_sim():
		R = [1]
		S = [3]
		N = list(range(1,5))
		A= [(1,2),(2,3),(2,4),(1,4),(4,3)]
		def all_path(G, R, S, max_len):
		    i = 0
		    P ={}
		    for r in R:
		        for s in S:
		            paths = nx.all_simple_paths(G, source=r, target=s, cutoff=max_len)
		            for path in map(nx.utils.pairwise, paths):
		                P[i] = list(path)
		                #print("path:", P[i])
		                i = i+1
		    return P
		G = nx.DiGraph()
		G.add_nodes_from(list(N))
		G.add_edges_from(A)
		P = all_path(G, R,S, 6)
		# travel cost function: use BPR function: t = tff(c+b(v/cap)**alpha)
		alpha = {}
		cap = {}
		tff = {}
		b = {}
		c = {}
		Q = {} # slope of demand function: demand = demand - demand_slope*rho

		for s in S:
		    for r in R:
		        Q[r,s]=8000
		#Q[1,3] = 1000
		#Q[2,3] = 800
		#Q[4,3] = 200

		for (r,s) in A:
		        tff[r,s] = 1
		        alpha[r,s] = 4.0
		        b[r,s] = 1.0
		        cap[r,s] = 1
		        c[r,s] = 1
		# e3
		tff[1,4] = 146.84299191
		#alpha[1,4] = 3.41053692
		b[1,4] = 2.70103383
		cap[1,4] = 1950.0

		#e1
		tff[1,2] = 141.0499002
		#alpha[1,2] = 3.94580842
		b[1,2] = 0.76710326
		cap[1,2] = 3920.0

		#e2
		tff[2,3] = 165.29733326
		#alpha[2,3] = 3.87066037
		b[2,3] = 2.42170115
		cap[2,3] = 1950.0

		#e4
		tff[4,3] = 140.98809353
		#alpha[4,3] = 4.17299034
		b[4,3] = 0.78246818
		cap[4,3] = 3920.0

		#e5
		tff[2,4] = 165.07866445
		#alpha[2,4] = 3.67176982
		b[2,4] = 2.30109182
		cap[2,4] = 1950.0

		# collision risk: CR = a0 + a1 * v + a2 * (v)^n
		a0 = {}
		a1 = {}
		a2 = {}
		n = {}

		for (r,s) in A:
		    a0[r,s] = 1
		    a1[r,s] = 1
		    a2[r,s] = 1
		    n[r,s] = 2
		#e1
		a0[1,2] = 2.99992*10**(-11)
		n[1,2] = 3

		#e2
		a0[2,3] = 4.6021*10**(-10)
		n[2,3] = 3


		#e3
		a0[1,4] = 3.721*10**(-10)
		n[1,4] = 3

		#e4
		a0[4,3] = 2.833*10**(-11)
		n[4,3] = 3
		#
		#e5
		a0[2,4] = 4.628*10**(-10)
		n[2,4] = 3

		return Network_transport(N, A, R, S, P, Q, alpha,  b, cap, tff, c, a0,a1)

class Network_transport:
	#def __init__(self, nodes, links, origin, destination, path, Q, alpha,b, cap, tff, c, a0, a1, n):
	#def __init__(self, nodes, links, origin, destination, path, Q, alpha,b, cap, tff, c,l_n):
	def __init__(self, nodes, links, origin, destination, path, Q, alpha,b, cap, tff, c,a0,a1):
		#self.N = nodes
		self.N = nodes
		self.A = links
		self.R = origin
		self.S = destination
		self.P = path
		self.Q = Q
		self.alpha = alpha
		self.b = b
		self.c = c
		self.cap = cap
		self.tff = tff
		self.a0 = a0
		self.a1 = a1
		#self.n = n
		self.P_set = list(self.P.keys())
		self.delta_path_od = Lib.Network_incidence_path_od(self.P,self.R,self.S, self.P_set)
		self.delta_link_path = Lib.Network_incidence_link_path(self.A, self.P, self.P_set)
		self.NodesOut,self.NodesIn = Lib.Network_topo(self.N,self.A)
		self.N_ap = Lib.Network_prev_nodes(self.P, self.A, self.P_set)
		#self.l_n = l_n

class Scen:

	def __init__(self, Xi, prob,link, link_tff, link_c, link_b, link_cap, link_alpha):
		self.Xi = Xi
		self.prob = prob
		self.link = link
		self.link_tff = link_tff
		self.link_c = link_c
		self.link_b = link_b
		self.link_cap = link_cap
		self.link_alpha = link_alpha
		#self.link_a2 = link_a2
		#self.link_a0 = link_a0
		#self.link_a1 = link_a1

class Problem:

	def __init__(self, ntw_transport, scen):
		self.ntw_t = ntw_transport
		self.scen = scen
		self.max_info_nodes = 24
		self.opt_gap = 0.001
		self.congestion = True
		self.capacity_scaler = 1
		self.M = 10000
		self.epsilon = 0.000001
		self.approx =True
		self.model = ConcreteModel('relax')
	# function to calculate if at link (i,j) on path p, whether a driver has receive information or not. 
	
	def calc_z(self, y):
		z = {}
		for (i,j) in self.ntw_t.A:
			for p in self.ntw_t.P_set:
				temp = 1
				for n in self.ntw_t.N_ap[i,j,p]:
					temp = temp * (1-y[n])
				z[i,j,p] = round(1 -temp)
		return z

	def update_travel_time(self, xi):

		for (r,s) in self.scen.link:
			self.ntw_t.tff[r,s] = self.scen.link_tff[r,s,xi]
			self.ntw_t.c[r,s] = self.scen.link_c[r,s,xi]
			self.ntw_t.b[r,s] = self.scen.link_b[r,s,xi]
			self.ntw_t.cap[r,s] = self.scen.link_cap[r,s,xi]
			self.ntw_t.alpha[r,s] = self.scen.link_alpha[r,s,xi]
			#self.ntw_t.a2[r,s] = self.scen.link_a2[r,s,xi]
			#self.ntw_t.a0[r,s] = self.scen.link_a0[r,s,xi]
			#self.ntw_t.a1[r,s] = self.scen.link_a1[r,s,xi]
		if not self.congestion:
			for (i,j) in self.ntw_t.A:
				self.ntw_t.b[i,j] = 0
		else:
			for (i,j) in self.ntw_t.A:
				self.ntw_t.cap[i,j] = self.ntw_t.cap[i,j] * self.capacity_scaler

	def uer_fix_upper(self, z):
		print("Start uer_fix_upper...")
		model = AbstractModel('uer_fix_upper')
		# variables
		print("Variables...")
		# traffic flow on path p the connects r,s at secen
		model.x_p = Var(self.ntw_t.P_set,self.scen.Xi, within=NonNegativeReals, initialize = 1)
		#traffic flow on link a from r to s that haven't received infromation on scen
		model.x_a = Var(self.ntw_t.R, self.ntw_t.S, self.ntw_t.A, within=NonNegativeReals, initialize = 1)
		# traffic flow of link a
		model.v = Var(self.ntw_t.A, self.scen.Xi, within=NonNegativeReals, initialize=1.0)
		print("Objective...")

		def obj_rule(model):
			exp = {}
			for xi in self.scen.Xi:
				self.update_travel_time(xi)
				exp[xi] = sum(self.ntw_t.tff[r,s]*(self.ntw_t.c[r,s]*model.v[r,s,xi]+(self.ntw_t.b[r,s]/(self.ntw_t.alpha[r,s]+1.0))*(model.v[r,s,xi]**(self.ntw_t.alpha[r,s]+1))/(self.ntw_t.cap[r,s]**(self.ntw_t.alpha[r,s]))) for (r,s) in self.ntw_t.A)
				print('exp',exp[xi])
			return sum(self.scen.prob[xi]*(exp[xi]) for xi in self.scen.Xi)
		model.obj = Objective(rule=obj_rule, sense=minimize)
		print("Constraints...")
		print("- Traffic flow on link..")
		# traffic flow path constraints

		def con_link_path_rule(model,i,j, xi):
			return (0, model.v[i,j, xi] - sum(self.ntw_t.delta_link_path[i,j,p]*model.x_p[p,xi] for p in self.ntw_t.P_set), 0)
		model.con_link_path = Constraint(self.ntw_t.A, self.scen.Xi, rule = con_link_path_rule)
		print("- Traffic path OD..")
		#traffic path OD constraints

		def con_path_od_rule(model,r,s,xi):
			#if r==s or s==r:
				#return Constraint.Skip
			#else:
			return (0, sum(self.ntw_t.delta_path_od[p,r,s] * model.x_p[p,xi] for p in self.ntw_t.P_set) - self.ntw_t.Q[r,s],0)
		model.con_path_od = Constraint(self.ntw_t.R, self.ntw_t.S,self.scen.Xi,rule = con_path_od_rule)
		print("- Traffic non-anticipitivity...")
		# traffic non anti constraints; z=1 means haven't receive information when they reach link a following path p
		
		def con_non_anti_rule(model,r, s,i,j,xi):
			if sum((1-z[i,j, p])*self.ntw_t.delta_link_path[i,j,p]*self.ntw_t.delta_path_od[p,r,s] for p in self.ntw_t.P_set) > 0:
				return (0, sum((1-z[i,j,p])*self.ntw_t.delta_link_path[i,j,p]*self.ntw_t.delta_path_od[p,r,s]* model.x_p[p,xi] for p in self.ntw_t.P_set)   - model.x_a[r,s,i,j] ,0)
			else:
				return Constraint.Skip
		model.con_non_anti = Constraint(self.ntw_t.R, self.ntw_t.S, self.ntw_t.A, self.scen.Xi, rule = con_non_anti_rule)
		print("Select solver...")
		optsolver = SolverFactory(nlpsol)
		#optsolver = SolverFactory(mipsol)
		print("Create instance...")
		inst = model.create_instance(report_timing=False)
		#inst.pprint()
		#print("Dual direction...")
		inst.dual = Suffix(direction=Suffix.IMPORT)
		print("Start solving...")
		results = optsolver.solve(inst, tee=True, load_solutions=True)

		print(inst.dual.display())
		v={}
		v_cap={}
		print("Solution: print v...")
		#print("Solution: print v_cap...")
		for (i,j) in self.ntw_t.A:
			for xi in self.scen.Xi:
				v[i,j,xi]=value(inst.v[i,j,xi])
				v_cap[i,j,xi] = v[i,j,xi]/self.ntw_t.cap[i,j] 
				print ("%i, %i, %i : %f" % (i, j, xi, v[i,j,xi]))
					
		

		crl_1={}
		print("Solution: print cr for link in xi 1...")

		for (i,j) in self.ntw_t.A:
			for xi in self.scen.Xi:
				v[i,j,xi]=value(inst.v[i,j,xi])
				

		crl_2={}
		print("Solution: print cr for link in xi 2...")

		for (i,j) in self.ntw_t.A:
			for xi in self.scen.Xi:
				v[i,j,xi]=value(inst.v[i,j,xi])
				



		x_p={}
		bi_x_p={}
		print("Solution: print x_p...")
		for p in self.ntw_t.P_set:
			for xi in self.scen.Xi:
				x_p[p, xi] = value(inst.x_p[p, xi])
				if x_p[p,xi] <=self.epsilon:
					bi_x_p[p,xi] = 0
				else:
					bi_x_p[p,xi] = 1
				print ("%i,%i : %f" % (p, xi, x_p[p, xi]))
		x_a={}
		bi_x_a={}
		print("Solution: print x_a...")
		for r in self.ntw_t.R:
			for s in self.ntw_t.S:
				for (i,j) in self.ntw_t.A:
					if sum((1-z[i,j, p])*self.ntw_t.delta_link_path[i,j,p]*self.ntw_t.delta_path_od[p,r,s] for p in self.ntw_t.P_set) > 0:
						x_a[r, s, i,j] = value(inst.x_a[r,s,i,j])
						if x_a[r, s, i,j] <=self.epsilon:
							bi_x_a[r, s, i,j] = 0
						else:
							bi_x_a[r, s, i,j] = 1
						print ("%i,%i,%i,%i : %f" % (r, s, i,j, x_a[r,s,i,j]))
					else:
						x_a[r, s, i,j] =0
						bi_x_a[r,s,i,j] = 0
						print ("%i,%i,%i,%i : %f" % (r, s, i,j, x_a[r,s,i,j]))
		
		'''dual_path_od={}
								print("Solution: print dual_path_od...")
								for r in self.ntw_t.R:
									for s in self.ntw_t.S:
										for xi in self.scen.Xi:
											dual_path_od[r,s,xi]= inst.dual[inst.con_path_od[r,s,xi]]
											print ("%i, %i, %i: %f" % (r,s, xi, dual_path_od[r,s,xi]))
								dual_non_anti={}
								print("Solution: print dual_non_anti...")
								for r in self.ntw_t.R:
									for s in self.ntw_t.S:
										for (i,j) in self.ntw_t.A:
											for xi in self.scen.Xi:
												if sum((1-z[i,j, p])*self.ntw_t.delta_link_path[i,j,p]*self.ntw_t.delta_path_od[p,r,s] for p in self.ntw_t.P_set) > 0:
													dual_non_anti[r,s,i,j,xi]= inst.dual[inst.con_non_anti[r,s,i,j,xi]]
													print ("%i, %i, %i, %i, %i: %f" % (r,s,i,j, xi, dual_non_anti[r,s,i,j,xi]))
												else:
													dual_non_anti[r,s,i,j,xi] = 0'''
		return v, x_p, x_a, bi_x_p, bi_x_a,v_cap,crl_1,crl_2
		#return v,v2, x_p, x_a, bi_x_p, bi_x_a,v_cap,v_cap2,crl_1,crl_2


	def check(self, z, LOA, h1, p_r,share):
		print("Start uer_fix_upper...")
		model = AbstractModel('uer_fix_upper')
		# variables
		print("Variables...")
		# traffic flow on path p the connects r,s at secen
		#model.x_p = Var(self.ntw_t.P_set,self.scen.Xi, within=NonNegativeReals, initialize = 1)
		model.x_p_h = Var(self.ntw_t.P_set,self.scen.Xi, within=NonNegativeReals, initialize = 1)
		model.x_p_c = Var(self.ntw_t.P_set,self.scen.Xi, within=NonNegativeReals, initialize = 1)
		#traffic flow on link a from r to s that haven't received infromation on scen
		model.x_a = Var(self.ntw_t.R, self.ntw_t.S, self.ntw_t.A, within=NonNegativeReals, initialize = 1)
		# traffic flow of link a
		model.v = Var(self.ntw_t.A, self.scen.Xi, within=NonNegativeReals, initialize=1.0)
		# capacity
		#model.cap = Var(self.ntw_t.A, self.scen.Xi, within=NonNegativeReals, initialize=1.0)
		l = 0.5
		#LOA = 1
		#p_r = 1

		print("Objective...")



		def obj_rule(model):
			exp = {}
			exp1 = {}
			for xi in self.scen.Xi:
				self.update_travel_time(xi)
				exp[xi] = sum(self.ntw_t.tff[r,s]*(self.ntw_t.c[r,s]*model.v[r,s,xi]+(self.ntw_t.b[r,s]/(self.ntw_t.alpha[r,s]+1.0))*(model.v[r,s,xi]**(self.ntw_t.alpha[r,s]+1))/(3600*l*self.ntw_t.l_n[r,s]/((sum(self.ntw_t.delta_link_path[r,s,p]*model.x_p_h[p,xi] for p in self.ntw_t.P_set)*1.8/ model.v[r,s, xi])+(1.8 - h1)*(sum(self.ntw_t.delta_link_path[r,s,p]*model.x_p_c[p,xi] for p in self.ntw_t.P_set)/ model.v[r,s, xi])))**(self.ntw_t.alpha[r,s])) for (r,s) in self.ntw_t.A)
				exp1[xi] = sum(model.v[r,s,xi]*self.ntw_t.tff[r,s]*(self.ntw_t.c[r,s]+self.ntw_t.b[r,s]*(model.v[r,s,xi]**(self.ntw_t.alpha[r,s]))/(3600*l*self.ntw_t.l_n[r,s]/((sum(self.ntw_t.delta_link_path[r,s,p]*model.x_p_h[p,xi] for p in self.ntw_t.P_set)*1.8/ model.v[r,s, xi])+(1.8 - h1)*(sum(self.ntw_t.delta_link_path[r,s,p]*model.x_p_c[p,xi] for p in self.ntw_t.P_set)/ model.v[r,s, xi])))**(self.ntw_t.alpha[r,s])) for (r,s) in self.ntw_t.A)
			#return sum(self.scen.prob[xi]*(p_r*((1-LOA)* exp[xi]+LOA*exp1[xi])) for xi in self.scen.Xi)
			return sum(self.scen.prob[xi]*((1-p_r) * exp[xi] + p_r*((1-LOA)* exp[xi]+LOA*exp1[xi])) for xi in self.scen.Xi)
			#return sum(self.scen.prob[xi]*(p_r*((1-LOA)* exp[xi]+LOA*exp1[xi])) for xi in self.scen.Xi)
		model.obj = Objective(rule=obj_rule, sense=minimize)
		print("Constraints...")
		print("- Traffic flow on link..")
		# traffic flow path constraints

		

		def con_link_path_rule(model,i,j, xi):
			#return (0, model.v[i,j, xi] - - sum(self.ntw_t.delta_link_path[i,j,p]*model.x_p_c[p,xi] for p in self.ntw_t.P_set), 0)
			#return (0, model.v[i,j, xi] - sum(self.ntw_t.delta_link_path[i,j,p]*model.x_p_c[p,xi] for p in self.ntw_t.P_set), 0)
			return (0, model.v[i,j, xi] - sum(self.ntw_t.delta_link_path[i,j,p]*model.x_p_h[p,xi] for p in self.ntw_t.P_set) - sum(self.ntw_t.delta_link_path[i,j,p]*model.x_p_c[p,xi] for p in self.ntw_t.P_set), 0)
		model.con_link_path = Constraint(self.ntw_t.A, self.scen.Xi, rule = con_link_path_rule)
		print("- Traffic path OD..")
		#traffic path OD constraints
		

		def con_path_od_hdv_rule(model,r,s,xi):
			return (0, sum(self.ntw_t.delta_path_od[p,r,s] * model.x_p_h[p,xi] for p in self.ntw_t.P_set) - self.ntw_t.Q[r,s]*(1-p_r),0)
		model.con_path_od_hdv = Constraint(self.ntw_t.R, self.ntw_t.S,self.scen.Xi,rule = con_path_od_hdv_rule)

		
		def con_path_od_cav_rule(model,r,s,xi):
			return (0, sum(self.ntw_t.delta_path_od[p,r,s] * model.x_p_c[p,xi] for p in self.ntw_t.P_set) - self.ntw_t.Q[r,s]*p_r,0)
		model.con_path_od_cav = Constraint(self.ntw_t.R, self.ntw_t.S,self.scen.Xi,rule = con_path_od_cav_rule)
		print("- Traffic non-anticipitivity...")

		#(3600*l*self.ntw_t.l_n[r,s]/((sum(self.ntw_t.delta_link_path[i,j,p]*model.x_p_h[p,xi] for p in self.ntw_t.P_set)*1.8/ model.v[i,j, xi])+(1.8 - LOA)*(sum(self.ntw_t.delta_link_path[r,s,p]*model.x_p_c[p,xi] for p in self.ntw_t.P_set)/ model.v[r,s, xi])))

		# traffic non anti constraints; z=1 means haven't receive information when they reach link a following path p
		'''
		def con_capacity_rule(model,i,j,xi):
			return (0,model.cap[i,j, xi] - (3600*l*self.ntw_t.l_n[r,s]/(1.8 - LOA)*(sum(self.ntw_t.delta_link_path[r,s,p]*model.x_p_c[p,xi] for p in self.ntw_t.P_set)/ model.v[r,s, xi])), 0)
			#return (0,model.cap[i,j, xi] - 3600*l*self.ntw_t.l_n[i,j]/((sum(self.ntw_t.delta_link_path[i,j,p]*model.x_p_h[p,xi] for p in self.ntw_t.P_set)/ model.v[i,j, xi])*1.8), 0)	
			#return (0,model.cap[i,j, xi] - (3600*l*self.ntw_t.l_n[i,j]/(sum(self.ntw_t.delta_link_path[i,j,p]*model.x_p_h[p,xi] for p in self.ntw_t.P_set)*1.8/ model.v[i,j, xi])), 0)	
			#return (0,model.cap[i,j, xi] - 3600*l*self.ntw_t.l_n[i,j]/((sum(self.ntw_t.delta_link_path[i,j,p]*model.x_p_h[p,xi] for p in self.ntw_t.P_set)/ model.v[i,j, xi])*1.8+(sum(self.ntw_t.delta_link_path[i,j,p]*model.x_p_c[p,xi] for p in self.ntw_t.P_set)/ model.v[i,j, xi]) *(1.8 - LOA)), 0)			
		model.con_capacity = Constraint(self.ntw_t.A, self.scen.Xi, rule = con_capacity_rule)'''
		
		def con_non_anti_rule(model,r, s,i,j,xi):
			if sum((1-z[i,j, p])*self.ntw_t.delta_link_path[i,j,p]*self.ntw_t.delta_path_od[p,r,s] for p in self.ntw_t.P_set) > 0:
				#return (0, sum((1-z[i,j,p])*self.ntw_t.delta_link_path[i,j,p]*self.ntw_t.delta_path_od[p,r,s]* model.x_p_c[p,xi] for p in self.ntw_t.P_set)  - model.x_a[r,s,i,j] ,0)
				return (0, sum((1-share*z[i,j,p])*self.ntw_t.delta_link_path[i,j,p]*self.ntw_t.delta_path_od[p,r,s]* model.x_p_h[p,xi] for p in self.ntw_t.P_set) + sum((1-z[i,j,p])*self.ntw_t.delta_link_path[i,j,p]*self.ntw_t.delta_path_od[p,r,s]* model.x_p_c[p,xi] for p in self.ntw_t.P_set)   - model.x_a[r,s,i,j] ,0)
			else:
				return Constraint.Skip
		model.con_non_anti = Constraint(self.ntw_t.R, self.ntw_t.S, self.ntw_t.A, self.scen.Xi, rule = con_non_anti_rule)
		print("Select solver...")
		optsolver = SolverFactory(nlpsol)
		#optsolver = SolverFactory(mipsol)
		print("Create instance...")
		inst = model.create_instance(report_timing=False)
		#inst.pprint()
		#print("Dual direction...")
		inst.dual = Suffix(direction=Suffix.IMPORT)
		print("Start solving...")
		'''
		results = optsolver.solve(inst, options={'acceptable_dual_inf_tol':10e20,
			'acceptable_compl_inf_tol':1e-08,'acceptable_constr_viol_tol':1e-08,
			'linear_solver':'ma27','acceptable_tol': 1e-09,'halt_on_ampl_error':'no'}, 
			tee=True)'''
		results = optsolver.solve(inst, options={'linear_solver':'ma27','halt_on_ampl_error':'no'},tee=True)
		#results = optsolver.solve(inst, tee=True, load_solutions=True)

		v={}
		print("Solution: print v...")
		for (i,j) in self.ntw_t.A:
			for xi in self.scen.Xi:
				v[i,j,xi]=value(inst.v[i,j,xi])
				print ("%i, %i, %i : %f" % (i, j, xi, v[i,j,xi]))

		v={}
		v_cap={}
		print("Solution: print v1...")
		#print("Solution: print v_cap...")
		for (i,j) in self.ntw_t.A:
			for xi in self.scen.Xi:
				v[i,j,xi]=value(inst.v[i,j,xi])
				if xi == 1:
					print ("%i, %i, %i : %f" % (i, j, xi, v[i,j,xi]))
					#print ("%i, %i, %i : %f" % (i, j, xi, v_cap[i,j,xi]))
				else:
					pass
		print("Solution: print v2...")
		v2={}
		v_cap2={}
		#print("Solution: print v...")
		#print("Solution: print v_cap...")
		for (i,j) in self.ntw_t.A:
			for xi in self.scen.Xi:
				v[i,j,xi]=value(inst.v[i,j,xi])
				v_cap[i,j,xi] = v[i,j,xi]/self.ntw_t.cap[i,j] 
				if xi == 2:
					print ("%i, %i, %i : %f" % (i, j, xi, v[i,j,xi]))
					#print ("%i, %i, %i : %f" % (i, j, xi, v_cap[i,j,xi]))
				else:
					pass

		crl_1={}
		print("Solution: print cr for link in xi 1...")

		for (i,j) in self.ntw_t.A:
			for xi in self.scen.Xi:
				v[i,j,xi]=value(inst.v[i,j,xi])
				crl_1[i,j,xi] = 1.0/ (1.0 + math.exp(-(self.ntw_t.a0[i,j] + self.ntw_t.a1[i,j]  * v[i,j,xi])))
				if xi == 1:
					print ("%i, %i, %i : %f" % (i, j, xi, crl_1[i,j,xi]))
				else:
					pass

		crl_2={}
		print("Solution: print cr for link in xi 2...")

		for (i,j) in self.ntw_t.A:
			for xi in self.scen.Xi:
				v[i,j,xi]=value(inst.v[i,j,xi])
				crl_2[i,j,xi] = 1.0/ (1.0 + math.exp(-(self.ntw_t.a0[i,j] + self.ntw_t.a1[i,j]  * v[i,j,xi])))
				if xi == 2:
					print ("%i, %i, %i : %f" % (i, j, xi, crl_2[i,j,xi]))
				else:
					pass


		
		x_p_h={}
		bi_x_p_h={}
		print("Solution: print x_p...")
		for p in self.ntw_t.P_set:
			for xi in self.scen.Xi:
				x_p_h[p, xi] = value(inst.x_p_h[p, xi])
				if x_p_h[p,xi] <=self.epsilon:
					bi_x_p_h[p,xi] = 0
				else:
					bi_x_p_h[p,xi] = 1
				print ("%i,%i : %f" % (p, xi, x_p_h[p, xi]))
		
		x_p_c={}
		bi_x_p_c={}
		print("Solution: print x_p...")
		for p in self.ntw_t.P_set:
			for xi in self.scen.Xi:
				x_p_c[p, xi] = value(inst.x_p_c[p, xi])
				if x_p_c[p,xi] <=self.epsilon:
					bi_x_p_c[p,xi] = 0
				else:
					bi_x_p_c[p,xi] = 1
				print ("%i,%i : %f" % (p, xi, x_p_c[p, xi]))

		x_a={}
		bi_x_a={}
		print("Solution: print x_a...")
		for r in self.ntw_t.R:
			for s in self.ntw_t.S:
				for (i,j) in self.ntw_t.A:
					if sum((1-z[i,j, p])*self.ntw_t.delta_link_path[i,j,p]*self.ntw_t.delta_path_od[p,r,s] for p in self.ntw_t.P_set) > 0:
						x_a[r, s, i,j] = value(inst.x_a[r,s,i,j])
						if x_a[r, s, i,j] <=self.epsilon:
							bi_x_a[r, s, i,j] = 0
						else:
							bi_x_a[r, s, i,j] = 1
						#print ("%i,%i,%i,%i : %f" % (r, s, i,j, x_a[r,s,i,j]))
					else:
						x_a[r, s, i,j] =0
						bi_x_a[r,s,i,j] = 0
		
		dual_path_od_hdv={}
		print("Solution: print dual_path_od...")
		for r in self.ntw_t.R:
			for s in self.ntw_t.S:
				for xi in self.scen.Xi:
					dual_path_od_hdv[r,s,xi]= inst.dual[inst.con_path_od_hdv[r,s,xi]]
					#print ("%i, %i, %i: %f" % (r,s, xi, dual_path_od[r,s,xi]))
		
		dual_path_od_cav={}
		print("Solution: print dual_path_od...")
		for r in self.ntw_t.R:
			for s in self.ntw_t.S:
				for xi in self.scen.Xi:
					dual_path_od_cav[r,s,xi]= inst.dual[inst.con_path_od_cav[r,s,xi]]
					#print ("%i, %i, %i: %f" % (r,s, xi, dual_path_od[r,s,xi]))

		dual_non_anti={}
		print("Solution: print dual_non_anti...")
		for r in self.ntw_t.R:
			for s in self.ntw_t.S:
				for (i,j) in self.ntw_t.A:
					for xi in self.scen.Xi:
						if sum((1-z[i,j, p])*self.ntw_t.delta_link_path[i,j,p]*self.ntw_t.delta_path_od[p,r,s] for p in self.ntw_t.P_set) > 0:
							dual_non_anti[r,s,i,j,xi]= inst.dual[inst.con_non_anti[r,s,i,j,xi]]
							#print ("%i, %i, %i, %i, %i: %f" % (r,s,i,j, xi, dual_non_anti[r,s,i,j,xi]))
						else:
							dual_non_anti[r,s,i,j,xi] = 0
		return v, x_p_h,x_p_c, x_a, dual_path_od_hdv,dual_path_od_cav, dual_non_anti, bi_x_p_h,bi_x_p_c, bi_x_a,v_cap,v_cap2,crl_1,crl_2
		#return v,v2, x_p_h, x_a, dual_path_od_hdv, dual_non_anti, bi_x_p_h, bi_x_a,v_cap,v_cap2,crl_1,crl_2
		#return v,v2, x_p_c, x_a, dual_path_od_cav, dual_non_anti, bi_x_p_c, bi_x_a,v_cap,v_cap2,crl_1,crl_2


	def uer_mixed_traffic(self, z):
		print("Start uer_fix_upper...")
		model = AbstractModel('uer_fix_upper')
		# variables
		print("Variables...")
		# traffic flow on path p the connects r,s at secen
		#model.x_p = Var(self.ntw_t.P_set,self.scen.Xi, within=NonNegativeReals, initialize = 1)
		model.x_p_h = Var(self.ntw_t.P_set,self.scen.Xi, within=NonNegativeReals, initialize = 1)
		model.x_p_c = Var(self.ntw_t.P_set,self.scen.Xi, within=NonNegativeReals, initialize = 1)
		#traffic flow on link a from r to s that haven't received infromation on scen
		model.x_a = Var(self.ntw_t.R, self.ntw_t.S, self.ntw_t.A, within=NonNegativeReals, initialize = 1)
		# traffic flow of link a
		model.v = Var(self.ntw_t.A, self.scen.Xi, within=NonNegativeReals, initialize=1.0)
		# capacity
		model.cap = Var(self.ntw_t.A, self.scen.Xi, within=NonNegativeReals, initialize=1.0)
		l = 0.5
		LOA = 1
		p_r = 0

		print("Objective...")



		def obj_rule(model):
			exp = {}
			exp1 = {}
			for xi in self.scen.Xi:
				self.update_travel_time(xi)
				exp[xi] = sum(self.ntw_t.tff[r,s]*(self.ntw_t.c[r,s]*model.v[r,s,xi]+(self.ntw_t.b[r,s]/(self.ntw_t.alpha[r,s]+1.0))*(model.v[r,s,xi]**(self.ntw_t.alpha[r,s]+1))/(model.cap[r,s,xi]**(self.ntw_t.alpha[r,s]))) for (r,s) in self.ntw_t.A)
				exp1[xi] = sum(model.v[r,s,xi]*self.ntw_t.tff[r,s]*(self.ntw_t.c[r,s]+self.ntw_t.b[r,s]*(model.v[r,s,xi]**(self.ntw_t.alpha[r,s]))/(model.cap[r,s,xi]**(self.ntw_t.alpha[r,s]))) for (r,s) in self.ntw_t.A)
			return sum(self.scen.prob[xi]*((1-p_r) * exp[xi] + p_r*((1-LOA)* exp[xi]+LOA*exp1[xi])) for xi in self.scen.Xi)
		model.obj = Objective(rule=obj_rule, sense=minimize)
		print("Constraints...")
		print("- Traffic flow on link..")
		# traffic flow path constraints

		

		def con_link_path_rule(model,i,j, xi):
			return (0, model.v[i,j, xi] - sum(self.ntw_t.delta_link_path[i,j,p]*model.x_p_h[p,xi] for p in self.ntw_t.P_set) - sum(self.ntw_t.delta_link_path[i,j,p]*model.x_p_c[p,xi] for p in self.ntw_t.P_set), 0)
		model.con_link_path = Constraint(self.ntw_t.A, self.scen.Xi, rule = con_link_path_rule)
		print("- Traffic path OD..")
		#traffic path OD constraints

		def con_path_od_hdv_rule(model,r,s,xi):
			return (0, sum(self.ntw_t.delta_path_od[p,r,s] * model.x_p_h[p,xi] for p in self.ntw_t.P_set) - self.ntw_t.Q[r,s]*(1-p_r),0)
		model.con_path_od_hdv = Constraint(self.ntw_t.R, self.ntw_t.S,self.scen.Xi,rule = con_path_od_hdv_rule)

		def con_path_od_cav_rule(model,r,s,xi):
			return (0, sum(self.ntw_t.delta_path_od[p,r,s] * model.x_p_c[p,xi] for p in self.ntw_t.P_set) - self.ntw_t.Q[r,s]*p_r,0)
		model.con_path_od_cav = Constraint(self.ntw_t.R, self.ntw_t.S,self.scen.Xi,rule = con_path_od_cav_rule)
		print("- Traffic non-anticipitivity...")
		# traffic non anti constraints; z=1 means haven't receive information when they reach link a following path p


		def con_capacity_rule(model,i,j,xi):
			return (0,model.cap[i,j, xi] - 3600*l*self.ntw_t.l_n[i,j]/((sum(self.ntw_t.delta_link_path[i,j,p]*model.x_p_h[p,xi] for p in self.ntw_t.P_set)/ model.v[i,j, xi])*1.8+(sum(self.ntw_t.delta_link_path[i,j,p]*model.x_p_c[p,xi] for p in self.ntw_t.P_set)/ model.v[i,j, xi]) *(1.8 - LOA)), 0)			
		model.con_capacity = Constraint(self.ntw_t.A, self.scen.Xi, rule = con_capacity_rule)
		
		def con_non_anti_rule(model,r, s,i,j,xi):
			if sum((1-z[i,j, p])*self.ntw_t.delta_link_path[i,j,p]*self.ntw_t.delta_path_od[p,r,s] for p in self.ntw_t.P_set) > 0:
				return (0, sum((1-z[i,j,p])*self.ntw_t.delta_link_path[i,j,p]*self.ntw_t.delta_path_od[p,r,s]* model.x_p_h[p,xi] for p in self.ntw_t.P_set) + sum((1-z[i,j,p])*self.ntw_t.delta_link_path[i,j,p]*self.ntw_t.delta_path_od[p,r,s]* model.x_p_c[p,xi] for p in self.ntw_t.P_set)   - model.x_a[r,s,i,j] ,0)
			else:
				return Constraint.Skip
		model.con_non_anti = Constraint(self.ntw_t.R, self.ntw_t.S, self.ntw_t.A, self.scen.Xi, rule = con_non_anti_rule)
		print("Select solver...")
		optsolver = SolverFactory(nlpsol)
		#optsolver = SolverFactory(mipsol)
		print("Create instance...")
		inst = model.create_instance(report_timing=False)
		#inst.pprint()
		#print("Dual direction...")
		inst.dual = Suffix(direction=Suffix.IMPORT)
		print("Start solving...")
		'''
		results = optsolver.solve(inst, options={'acceptable_dual_inf_tol':10e20,
			'acceptable_compl_inf_tol':1e-08,'acceptable_constr_viol_tol':1e-08,
			'linear_solver':'ma27','acceptable_tol': 1e-09,'halt_on_ampl_error':'no'}, 
			tee=True)'''
		results = optsolver.solve(inst, options={'linear_solver':'ma27','halt_on_ampl_error':'no'},tee=True)
		#results = optsolver.solve(inst, tee=True, load_solutions=True)

		v={}
		print("Solution: print v...")
		for (i,j) in self.ntw_t.A:
			for xi in self.scen.Xi:
				v[i,j,xi]=value(inst.v[i,j,xi])
				print ("%i, %i, %i : %f" % (i, j, xi, v[i,j,xi]))

		v={}
		v_cap={}
		print("Solution: print v1...")
		#print("Solution: print v_cap...")
		for (i,j) in self.ntw_t.A:
			for xi in self.scen.Xi:
				v[i,j,xi]=value(inst.v[i,j,xi])
				if xi == 1:
					print ("%i, %i, %i : %f" % (i, j, xi, v[i,j,xi]))
					#print ("%i, %i, %i : %f" % (i, j, xi, v_cap[i,j,xi]))
				else:
					pass
		print("Solution: print v2...")
		v2={}
		v_cap2={}
		#print("Solution: print v...")
		#print("Solution: print v_cap...")
		for (i,j) in self.ntw_t.A:
			for xi in self.scen.Xi:
				v[i,j,xi]=value(inst.v[i,j,xi])
				v_cap[i,j,xi] = v[i,j,xi]/self.ntw_t.cap[i,j] 
				if xi == 2:
					print ("%i, %i, %i : %f" % (i, j, xi, v[i,j,xi]))
					#print ("%i, %i, %i : %f" % (i, j, xi, v_cap[i,j,xi]))
				else:
					pass

		crl_1={}
		print("Solution: print cr for link in xi 1...")

		for (i,j) in self.ntw_t.A:
			for xi in self.scen.Xi:
				v[i,j,xi]=value(inst.v[i,j,xi])
				crl_1[i,j,xi] = 1.0/ (1.0 + math.exp(-(self.ntw_t.a0[i,j] + self.ntw_t.a1[i,j]  * v[i,j,xi])))
				if xi == 1:
					print ("%i, %i, %i : %f" % (i, j, xi, crl_1[i,j,xi]))
				else:
					pass

		crl_2={}
		print("Solution: print cr for link in xi 2...")

		for (i,j) in self.ntw_t.A:
			for xi in self.scen.Xi:
				v[i,j,xi]=value(inst.v[i,j,xi])
				crl_2[i,j,xi] = 1.0/ (1.0 + math.exp(-(self.ntw_t.a0[i,j] + self.ntw_t.a1[i,j]  * v[i,j,xi])))
				if xi == 2:
					print ("%i, %i, %i : %f" % (i, j, xi, crl_2[i,j,xi]))
				else:
					pass



		x_p_h={}
		bi_x_p_h={}
		print("Solution: print x_p...")
		for p in self.ntw_t.P_set:
			for xi in self.scen.Xi:
				x_p_h[p, xi] = value(inst.x_p_h[p, xi])
				if x_p_h[p,xi] <=self.epsilon:
					bi_x_p_h[p,xi] = 0
				else:
					bi_x_p_h[p,xi] = 1
				print ("%i,%i : %f" % (p, xi, x_p_h[p, xi]))

		x_p_c={}
		bi_x_p_c={}
		print("Solution: print x_p...")
		for p in self.ntw_t.P_set:
			for xi in self.scen.Xi:
				x_p_c[p, xi] = value(inst.x_p_c[p, xi])
				if x_p_c[p,xi] <=self.epsilon:
					bi_x_p_c[p,xi] = 0
				else:
					bi_x_p_c[p,xi] = 1
				print ("%i,%i : %f" % (p, xi, x_p_c[p, xi]))

		x_a={}
		bi_x_a={}
		print("Solution: print x_a...")
		for r in self.ntw_t.R:
			for s in self.ntw_t.S:
				for (i,j) in self.ntw_t.A:
					if sum((1-z[i,j, p])*self.ntw_t.delta_link_path[i,j,p]*self.ntw_t.delta_path_od[p,r,s] for p in self.ntw_t.P_set) > 0:
						x_a[r, s, i,j] = value(inst.x_a[r,s,i,j])
						if x_a[r, s, i,j] <=self.epsilon:
							bi_x_a[r, s, i,j] = 0
						else:
							bi_x_a[r, s, i,j] = 1
						#print ("%i,%i,%i,%i : %f" % (r, s, i,j, x_a[r,s,i,j]))
					else:
						x_a[r, s, i,j] =0
						bi_x_a[r,s,i,j] = 0
		dual_path_od_hdv={}
		print("Solution: print dual_path_od...")
		for r in self.ntw_t.R:
			for s in self.ntw_t.S:
				for xi in self.scen.Xi:
					dual_path_od_hdv[r,s,xi]= inst.dual[inst.con_path_od_hdv[r,s,xi]]
					#print ("%i, %i, %i: %f" % (r,s, xi, dual_path_od[r,s,xi]))

		dual_path_od_cav={}
		print("Solution: print dual_path_od...")
		for r in self.ntw_t.R:
			for s in self.ntw_t.S:
				for xi in self.scen.Xi:
					dual_path_od_cav[r,s,xi]= inst.dual[inst.con_path_od_cav[r,s,xi]]
					#print ("%i, %i, %i: %f" % (r,s, xi, dual_path_od[r,s,xi]))

		dual_non_anti={}
		print("Solution: print dual_non_anti...")
		for r in self.ntw_t.R:
			for s in self.ntw_t.S:
				for (i,j) in self.ntw_t.A:
					for xi in self.scen.Xi:
						if sum((1-z[i,j, p])*self.ntw_t.delta_link_path[i,j,p]*self.ntw_t.delta_path_od[p,r,s] for p in self.ntw_t.P_set) > 0:
							dual_non_anti[r,s,i,j,xi]= inst.dual[inst.con_non_anti[r,s,i,j,xi]]
							#print ("%i, %i, %i, %i, %i: %f" % (r,s,i,j, xi, dual_non_anti[r,s,i,j,xi]))
						else:
							dual_non_anti[r,s,i,j,xi] = 0
		return v,v2, x_p_h,x_p_c, x_a, dual_path_od_hdv,dual_path_od_cav, dual_non_anti, bi_x_p_h,bi_x_p_c, bi_x_a,v_cap,v_cap2,crl_1,crl_2
	
	# linearized bilinear terms
	def uer_fix_upper2(self, z):
		print("Start uer_fix_upper2...")
		model = ConcreteModel('uer_fix_upper2')
		# variables
		print("Variables...")
		# traffic flow on path p the connects r,s at secen
		model.x_p = Var(self.ntw_t.P_set,self.scen.Xi, within=NonNegativeReals, initialize = 1)
		#traffic flow on link a from r to s that haven't received infromation on scen
		model.x_a = Var(self.ntw_t.R, self.ntw_t.S, self.ntw_t.A, within=NonNegativeReals, initialize = 1)
		# traffic flow of link a
		model.v = Var(self.ntw_t.A, self.scen.Xi, within=NonNegativeReals, initialize=1.0)
		# m = (1-z)x_p
		model.m = Var(self.ntw_t.A, self.ntw_t.P_set, self.scen.Xi, within=NonNegativeReals, initialize=1.0)
		print("Objective...")
		
		def obj_rule(model):
			exp = {}
			for xi in self.scen.Xi:
				self.update_travel_time(xi)
				exp[xi] = sum(self.ntw_t.tff[r,s]*(self.ntw_t.c[r,s]*model.v[r,s,xi]+(self.ntw_t.b[r,s]/(self.ntw_t.alpha[r,s]+1.0))*(model.v[r,s,xi]**(self.ntw_t.alpha[r,s]+1))/(self.ntw_t.cap[r,s]**(self.ntw_t.alpha[r,s]))) for (r,s) in self.ntw_t.A)
			return sum(self.scen.prob[xi]*(exp[xi]) for xi in self.scen.Xi)
		model.obj = Objective(rule=obj_rule, sense=minimize)
		print("Constraints...")
		print("- Traffic flow on link..")
		# traffic flow path constraints
		
		def con_link_path_rule(model,i,j, xi):
			return (0, model.v[i,j, xi] - sum(self.ntw_t.delta_link_path[i,j,p]*model.x_p[p,xi] for p in self.ntw_t.P_set), 0)
		model.con_link_path = Constraint(self.ntw_t.A, self.scen.Xi, rule = con_link_path_rule)
		print("- Traffic path OD..")
		#traffic path OD constraints
		
		def con_path_od_rule(model,r,s,xi):
			return (0, sum(self.ntw_t.delta_path_od[p,r,s] * model.x_p[p,xi] for p in self.ntw_t.P_set) - self.ntw_t.Q[r,s],0)
		model.con_path_od = Constraint(self.ntw_t.R, self.ntw_t.S,self.scen.Xi,rule = con_path_od_rule)
		print("- Traffic non-anticipitivity...")
		# traffic non anti constraints; z=1 means haven't receive information when they reach link a following path p
		
		def con_non_anti_rule(model,r, s,i,j,xi):
			if sum((1-z[i,j, p])*self.ntw_t.delta_link_path[i,j,p]*self.ntw_t.delta_path_od[p,r,s] for p in self.ntw_t.P_set) > 0:
				return (0, sum(self.ntw_t.delta_link_path[i,j,p]*self.ntw_t.delta_path_od[p,r,s]* model.m[i,j,p,xi] for p in self.ntw_t.P_set)   - model.x_a[r,s,i,j] ,0)
			else:
				return (0, sum(self.ntw_t.delta_link_path[i,j,p]*self.ntw_t.delta_path_od[p,r,s]* model.m[i,j,p,xi] for p in self.ntw_t.P_set)   - model.x_a[r,s,i,j] ,0)
				#return Constraint.Skip
		model.con_non_anti = Constraint(self.ntw_t.R, self.ntw_t.S, self.ntw_t.A, self.scen.Xi, rule = con_non_anti_rule)
		
		def con1_linear_rule(model,i,j,p,xi):
			#return Constraint.Skip
			return (-infinity, model.m[i,j,p,xi] - model.x_p[p, xi], 0)
		model.con1_linear = Constraint(self.ntw_t.A, self.ntw_t.P_set,self.scen.Xi,rule = con1_linear_rule)
		
		def con2_linear_rule(model,i,j,p,xi):
			return (-infinity, -model.m[i,j,p,xi] + model.x_p[p, xi], self.M*z[i,j,p])
		model.con2_linear = Constraint(self.ntw_t.A, self.ntw_t.P_set,self.scen.Xi,rule = con2_linear_rule)
		
		def con3_linear_rule(model,i,j,p,xi):
			print(z[i,j,p])
			return (-infinity, model.m[i,j,p,xi], self.M*(1-z[i,j,p]))
		model.con3_linear = Constraint(self.ntw_t.A, self.ntw_t.P_set,self.scen.Xi,rule = con3_linear_rule)
		print("Select solver...")
		optsolver = SolverFactory(nlpsol)
		#optsolver = SolverFactory(mipsol)
		print("Create instance...")
		#print("Dual direction...")
		model.dual = Suffix(direction=Suffix.IMPORT_EXPORT)
		print("Start solving...")
		results = optsolver.solve(model, tee=True, load_solutions=True)
		v={}
		print("Solution: print v...")
		for (i,j) in self.ntw_t.A:
			for xi in self.scen.Xi:
				v[i,j,xi]=value(model.v[i,j,xi])
				if v[i,j,xi] <= self.epsilon:
					v[i,j,xi] = 0
				print ("%i, %i, %i : %f" % (i, j, xi, v[i,j,xi]))
		x_p={}
		bi_x_p={}
		print("Solution: print x_p...")
		for p in self.ntw_t.P_set:
			for xi in self.scen.Xi:
				x_p[p, xi] = value(model.x_p[p, xi])
				if x_p[p,xi] <=self.epsilon:
					bi_x_p[p,xi] = 0
					x_p[p,xi] = 0
				else:
					bi_x_p[p,xi] = 1
				print ("%i,%i : %f" % (p, xi, x_p[p, xi]))
		x_a={}
		bi_x_a={}
		print("Solution: print x_a...")
		for r in self.ntw_t.R:
			for s in self.ntw_t.S:
				for (i,j) in self.ntw_t.A:
					if sum((1-z[i,j, p])*self.ntw_t.delta_link_path[i,j,p]*self.ntw_t.delta_path_od[p,r,s] for p in self.ntw_t.P_set) > 0:
						x_a[r, s, i,j] = value(model.x_a[r,s,i,j])
						if x_a[r, s, i,j] <=self.epsilon:
							bi_x_a[r, s, i,j] = 0
							x_a[r,s,i,j] = 0
						else:
							bi_x_a[r, s, i,j] = 1
						print ("%i,%i,%i,%i : %f" % (r, s, i,j, x_a[r,s,i,j]))
					else:
						x_a[r, s, i,j] =0
						bi_x_a[r,s,i,j] = 0
		m = {}
		dual1_linear={}
		dual2_linear={}
		dual3_linear={}
		print("Solution: print dual1_linear...")
		for (i,j) in self.ntw_t.A:
			for p in self.ntw_t.P_set:
				for xi in self.scen.Xi:
					m[i,j,p,xi] = value(model.m[i,j,p,xi])
					if (m[i,j,p,xi] <= self.epsilon):
						m[i,j,p,xi] = 0
					dual1_linear[i,j,p,xi]= -model.dual[model.con1_linear[i,j,p,xi]]
					#dual1_linear[i,j,p,xi]= 0
					if (dual1_linear[i,j,p,xi] <= self.epsilon and dual1_linear[i,j,p,xi] >= -self.epsilon):
						dual1_linear[i,j,p,xi] = 0
					dual2_linear[i,j,p,xi]= -model.dual[model.con2_linear[i,j,p,xi]]
					if (dual2_linear[i,j,p,xi] <= self.epsilon and dual2_linear[i,j,p,xi] >= -self.epsilon):
						dual2_linear[i,j,p,xi] = 0
					dual3_linear[i,j,p,xi]= -model.dual[model.con3_linear[i,j,p,xi]]
					if (dual3_linear[i,j,p,xi] <= self.epsilon and dual3_linear[i,j,p,xi] >= -self.epsilon):
						dual3_linear[i,j,p,xi] = 0
					print ("dual1_linear: %i, %i, %i, %i: %f" % (i,j,p, xi, dual1_linear[i,j,p,xi]))
					print ("zmdual2_linear: %i, %i, %i, %i: %f" % (i,j,p, xi, dual2_linear[i,j,p,xi]))
					print ("dual3_linear: %i, %i, %i, %i: %f" % (i,j,p, xi, dual3_linear[i,j,p,xi]))
		return v, x_p, x_a, m, dual1_linear, dual2_linear, dual3_linear
	
	'''def calc_obj_upper(self, v):
		print("Start calc_obj_upper...")
		exp = {}
		for xi in self.scen.Xi:
			self.update_travel_time(xi)
			exp[xi] = sum(v[r,s,xi]*self.ntw_t.tff[r,s]*(self.ntw_t.c[r,s]+self.ntw_t.b[r,s]*(v[r,s,xi]**(self.ntw_t.alpha[r,s]))/(self.ntw_t.cap[r,s]**(self.ntw_t.alpha[r,s]))) for (r,s) in self.ntw_t.A)
		return sum(self.scen.prob[xi]*(exp[xi]) for xi in self.scen.Xi)'''

	def calc_obj_upper(self, v):
		print("Start calc_obj_upper...")
		exp = {}
		for xi in self.scen.Xi:
			self.update_travel_time(xi)
			#exp[xi] = sum(self.ntw_t.a0[r,s]  * v[r,s,xi]** self.ntw_t.n[r,s] for (r,s) in self.ntw_t.A)
			exp[xi] = sum(v[r,s,xi]*self.ntw_t.tff[r,s]*(self.ntw_t.c[r,s]+self.ntw_t.b[r,s]*(v[r,s,xi]**(self.ntw_t.alpha[r,s]))/(self.ntw_t.cap[r,s]**(self.ntw_t.alpha[r,s]))) for (r,s) in self.ntw_t.A)
			#exp[xi] = sum((math.exp(self.ntw_t.a0[r,s] + self.ntw_t.a1[r,s]  * v[r,s,xi]) / (1 + math.exp(self.ntw_t.a0[r,s] + self.ntw_t.a1[r,s]  * v[r,s,xi]))) * v[r,s,xi] for (r,s) in self.ntw_t.A)
			#exp[xi] = sum(1.0/ (1.0 + math.exp(-(self.ntw_t.a0[r,s] + self.ntw_t.a1[r,s]  * v[r,s,xi]))) * v[r,s,xi] for (r,s) in self.ntw_t.A)
		return sum(self.scen.prob[xi]*(exp[xi]) for xi in self.scen.Xi)
	
	def build_relax_problem(self, y_ini, z_ini, mipgap):
		print("Start build_relax_problem...")
		#self.model = ConcreteModel('relax')
		model = self.model
		# variables
		print("Variables...")
		
		# information node
		def y_init_rule(model, i):
			return y_ini[i]
		model.y = Var(self.ntw_t.N, domain=Binary, initialize = y_init_rule)
		
		# z_ap: = 1 if a on p has info;
		def z_init_rule(model, i,j,p):
			return z_ini[i,j,p]
		model.z = Var(self.ntw_t.A, self.ntw_t.P_set, domain=Binary, initialize = z_init_rule)
		
		# upper level obj
		model.obj_fn = Var(within = NonNegativeReals,initialize = 0)
		def obj_rule(model):
			exp = model.obj_fn
			return exp
		model.obj = Objective(rule=obj_rule, sense=minimize)
		print("model cuts")
		model.cuts = ConstraintList()
		
		# has info if any node before link a is infomation node
		def has_info_rule(model,i,j, p, n):
			if n in self.ntw_t.N_ap[i,j,p]:
				return (0,model.z[i,j,p] - model.y[n], infinity)
			else:
				return Constraint.Skip
		model.has_info = Constraint(self.ntw_t.A, self.ntw_t.P_set, self.ntw_t.N, rule = has_info_rule)
		
		# no information: when all the nodes before link a are not information nodes
		def no_info_rule(model,i,j, p):
			return (-infinity,model.z[i,j,p] - sum(model.y[n] for n in self.ntw_t.N_ap[i,j,p]), 0)
		model.no_info = Constraint(self.ntw_t.A, self.ntw_t.P_set, rule = no_info_rule)
		
		# total number of information nodes
		def con_max_info_nodes_rule(model):
			return (-infinity, sum(model.y[n] for n in self.ntw_t.N)-self.max_info_nodes,0)
		model.con_max_info_nodes = Constraint(rule=con_max_info_nodes_rule)
	
	def solve_relax_problem(self,model,mipgap):
		print("Start solve_relax_problem...")
		print("Select solver...")
		optsolver = SolverFactory(nlpsol)
		#optsolver = SolverFactory(mipsol)
		#optsolver.options['CPXPARAM_Simplex_Tolerances_Feasibility']= 0.000001
		#optsolver.options['mipgap']= mipgap
		#print("Print Model...")
		#model.pprint()
		print("Start solving...")
		results = optsolver.solve(model, tee=True, load_solutions=True)
		if (results.solver.status == SolverStatus.ok) and (results.solver.termination_condition == TerminationCondition.optimal):
			# Do something when the solution in optimal and feasible
			print("Optimal")
		elif (results.solver.termination_condition == TerminationCondition.infeasible):
			# Do something when model in infeasible
			print("infeasible")
		else:
			# Something else is wrong
			print("Solver Status: ",  results.solver.status)
		obj_upper = value(model.obj)
		print("Solution: print y...")
		y = {}
		for n in self.ntw_t.N:
			y[n] = round(value(model.y[n]))
			print ("%i : %i" % (n, y[n]))
		print("Solution: print z...")
		z = {}
		for (i,j) in self.ntw_t.A:
			for p in self.ntw_t.P_set:
				z[i,j,p] = round(value(model.z[i,j,p]))
				print ("%i, %i, %i : %i" % (i, j, p, z[i,j,p]))
		print("zm2 Solution: print obj_fn...")
		obj_upper = value(model.obj_fn)
		print ("%f" % obj_upper)
		return obj_upper, y, z
	
	def calc_support_lower(self, model,z, v, x_p, x_a, dual_non_anti):
		exp = {}
		exp1 = {}
		for xi in self.scen.Xi:
			self.update_travel_time(xi)
			exp[xi] = sum(self.ntw_t.tff[r,s]*(self.ntw_t.c[r,s]*v[r,s,xi]+(self.ntw_t.b[r,s]/(self.ntw_t.alpha[r,s]+1.0))*(v[r,s,xi]**(self.ntw_t.alpha[r,s]+1))/(self.ntw_t.cap[r,s]**(self.ntw_t.alpha[r,s]))) for (r,s) in self.ntw_t.A)
			exp1[xi] = 0
			for (i,j) in self.ntw_t.A:
				if sum((1-z[i,j, p])*self.ntw_t.delta_link_path[i,j,p] for p in self.ntw_t.P_set) > 0:
					exp1[xi] = exp1[xi]+ sum(dual_non_anti[self.ntw_t.P[p][0][0], self.ntw_t.P[p][len(self.ntw_t.P[p])-1][1], i, j, xi] *( self.ntw_t.delta_link_path[i,j,p]*(1-model.z[i,j,p])*x_p[p,xi] - x_a[self.ntw_t.P[p][0][0], self.ntw_t.P[p][len(self.ntw_t.P[p])-1][1], i,j]) for p in self.ntw_t.P_set)
		return sum(self.scen.prob[xi]*(exp[xi]) + exp1[xi] for xi in self.scen.Xi)
	
	def calc_support_lower2(self, model,z, v, x_p, x_a,m, dual2_linear, dual3_linear):
		exp = {}
		exp1 = {}
		for xi in self.scen.Xi:
			self.update_travel_time(xi)
			exp[xi] = sum(self.ntw_t.tff[r,s]*(self.ntw_t.c[r,s]*v[r,s,xi]+(self.ntw_t.b[r,s]/(self.ntw_t.alpha[r,s]+1.0))*(v[r,s,xi]**(self.ntw_t.alpha[r,s]+1))/(self.ntw_t.cap[r,s]**(self.ntw_t.alpha[r,s]))) for (r,s) in self.ntw_t.A)
			#exp1[xi] = sum((dual1_linear[i,j,p,xi]-dual2_linear[i,j,p,xi])*(m[i,j, p, xi]-x_p[p,xi] - model.z[i,j,p]*self.M) + dual2_linear[i,j,p,xi]*(-m[i,j,p,xi] + x_p[p,xi]-model.z[i,j,p]*self.M) + dual3_linear[i,j,p,xi]*(m[i,j,p,xi] - (1-model.z[i,j,p])*self.M) for (i,j) in self.ntw_t.A for p in self.ntw_t.P_set)
			#exp1[xi] = sum((dual1_linear[i,j,p,xi]-dual2_linear[i,j,p,xi])*(m[i,j, p, xi]-x_p[p,xi] - model.z[i,j,p]*self.M)+dual3_linear[i,j,p,xi]*(m[i,j,p,xi] - (1-model.z[i,j,p])*self.M) for (i,j) in self.ntw_t.A for p in self.ntw_t.P_set)
			exp1[xi] = sum((-dual2_linear[i,j,p,xi])*(m[i,j, p, xi]-x_p[p,xi] - model.z[i,j,p]*self.M)+dual3_linear[i,j,p,xi]*(m[i,j,p,xi] - (1-model.z[i,j,p])*self.M) for (i,j) in self.ntw_t.A for p in self.ntw_t.P_set)
		return sum(self.scen.prob[xi]*(exp[xi]) + exp1[xi] for xi in self.scen.Xi)

# relax second term: lower bound may be too low
	
	def calc_support_upper(self, support_lower, z):
		model = ConcreteModel('support_upper')
		# variables
		print("Variables...")
		# traffic flow on path p the connects r,s at secen
		model.x_p = Var(self.ntw_t.P_set,self.scen.Xi, within=NonNegativeReals, initialize = 1)
		#traffic flow on link a from r to s that haven't received infromation on scen
		model.x_a = Var(self.ntw_t.R, self.ntw_t.S, self.ntw_t.A, within=NonNegativeReals, initialize = 1)
		# traffic flow of link a
		model.v = Var(self.ntw_t.A, self.scen.Xi, within=NonNegativeReals, initialize=1.0)
		print("Objective...")
		
		def obj_rule(model):
			exp = {}
			for xi in self.scen.Xi:
				self.update_travel_time(xi)
				exp[xi] = sum(self.ntw_t.tff[r,s]*self.ntw_t.b[r,s]*self.ntw_t.alpha[r,s]/(self.ntw_t.alpha[r,s]+1.0)*model.v[r,s,xi]**(self.ntw_t.alpha[r,s]+1)/(self.ntw_t.cap[r,s]**(self.ntw_t.alpha[r,s])) for (r,s) in self.ntw_t.A)
			return sum(self.scen.prob[xi]*(exp[xi]) for xi in self.scen.Xi)
		model.obj = Objective(rule=obj_rule, sense=minimize)
		print("Constraints...")
		print("- Traffic flow on link..")
		
		# traffic flow path constraints
		def con_link_path_rule(model,i,j, xi):
			return (0, model.v[i,j, xi] - sum(self.ntw_t.delta_link_path[i,j,p]*model.x_p[p,xi] for p in self.ntw_t.P_set), 0)
		model.con_link_path = Constraint(self.ntw_t.A, self.scen.Xi, rule = con_link_path_rule)
		print("- Traffic path OD..")
		
		#traffic path OD constraints
		def con_path_od_rule(model,r,s,xi):
			return (0, sum(self.ntw_t.delta_path_od[p,r,s] * model.x_p[p,xi] for p in self.ntw_t.P_set) - self.ntw_t.Q[r,s],0)
		model.con_path_od = Constraint(self.ntw_t.R, self.ntw_t.S,self.scen.Xi,rule = con_path_od_rule)
		print("- Traffic non-anticipitivity...")
		
		# traffic non anti constraints; z=1 means haven't receive information when they reach link a following path p
		def con_non_anti_rule(model,r, s,i,j,xi):
			if sum((1-z[i,j, p])*self.ntw_t.delta_link_path[i,j,p]*self.ntw_t.delta_path_od[p,r,s] for p in self.ntw_t.P_set) > 0:
				return (0, sum((1-z[i,j,p])*self.ntw_t.delta_link_path[i,j,p]*self.ntw_t.delta_path_od[p,r,s]* model.x_p[p,xi] for p in self.ntw_t.P_set)   - model.x_a[r,s,i,j] ,0)
			else:
				return Constraint.Skip
		model.con_non_anti = Constraint(self.ntw_t.R, self.ntw_t.S, self.ntw_t.A, self.scen.Xi, rule = con_non_anti_rule)
		print("Select solver...")
		#optsolver = SolverFactory(nlpsol)
		optsolver = SolverFactory(mipsol)
		print("Create instance...")
		model.pprint()
		print("Start solving...")
		results = optsolver.solve(model, tee=True, load_solutions=True)
		v={}
		print("Solution: print v...")
		for (i,j) in self.ntw_t.A:
			for xi in self.scen.Xi:
				v[i,j,xi]=value(model.v[i,j,xi])
				print ("%i, %i, %i : %f" % (i, j, xi, v[i,j,xi]))
		x_p={}
		bi_x_p={}
		print("Solution: print x_p...")
		for p in self.ntw_t.P_set:
			for xi in self.scen.Xi:
				x_p[p, xi] = value(model.x_p[p, xi])
				if x_p[p,xi] <=self.epsilon:
					bi_x_p[p,xi] = 0
				else:
					bi_x_p[p,xi] = 1
				print ("%i,%i : %f" % (p, xi, x_p[p, xi]))
		obj_val = value(model.obj)
		return support_lower + obj_val
	
	# combine first and second terms: lower bound may be too high, higher than upper bound. This seems wrong, but the paper by Gao is using this approach. Don't know why yet.
	def calc_support_upper2(self, z, dual_non_anti):
		model = ConcreteModel('support_upper2')
		# variables
		print("Variables...")
		# traffic flow on path p the connects r,s at secen
		model.x_p = Var(self.ntw_t.P_set,self.scen.Xi, within=NonNegativeReals, initialize = 1)
		#traffic flow on link a from r to s that haven't received infromation on scen
		model.x_a = Var(self.ntw_t.R, self.ntw_t.S, self.ntw_t.A, within=NonNegativeReals, initialize = 1)
		# traffic flow of link a
		model.v = Var(self.ntw_t.A, self.scen.Xi, within=NonNegativeReals, initialize=1.0)
		print("Objective...")
		
		def obj_rule(model):
			exp = {}
			exp1 = {}
			exp2 = {}
			for xi in self.scen.Xi:
				self.update_travel_time(xi)
				exp[xi] = sum(self.ntw_t.tff[r,s]*(self.ntw_t.c[r,s]*model.v[r,s,xi]+(self.ntw_t.b[r,s]/(self.ntw_t.alpha[r,s]+1.0))*(model.v[r,s,xi]**(self.ntw_t.alpha[r,s]+1))/(self.ntw_t.cap[r,s]**(self.ntw_t.alpha[r,s]))) for (r,s) in self.ntw_t.A)
				exp1[xi] = sum(self.ntw_t.tff[r,s]*self.ntw_t.b[r,s]*self.ntw_t.alpha[r,s]/(self.ntw_t.alpha[r,s]+1.0)*model.v[r,s,xi]**(self.ntw_t.alpha[r,s]+1)/(self.ntw_t.cap[r,s]**(self.ntw_t.alpha[r,s])) for (r,s) in self.ntw_t.A)
				exp2[xi] = 0
				for (i,j) in self.ntw_t.A:
					if sum((1-z[i,j, p])*self.ntw_t.delta_link_path[i,j,p] for p in self.ntw_t.P_set) > 0:
						exp2[xi] = exp2[xi]+ sum(dual_non_anti[self.ntw_t.P[p][0][0], self.ntw_t.P[p][len(self.ntw_t.P[p])-1][1], i, j, xi] *( self.ntw_t.delta_link_path[i,j,p]*(1-z[i,j,p])*model.x_p[p,xi] - model.x_a[self.ntw_t.P[p][0][0], self.ntw_t.P[p][len(self.ntw_t.P[p])-1][1], i,j]) for p in self.ntw_t.P_set)
			return sum(self.scen.prob[xi]*(exp[xi] + exp1[xi]) + exp2[xi] for xi in self.scen.Xi)
		model.obj = Objective(rule=obj_rule, sense=minimize)
		print("Constraints...")
		print("- Traffic flow on link..")
		
		# traffic flow path constraints
		def con_link_path_rule(model,i,j, xi):
			return (0, model.v[i,j, xi] - sum(self.ntw_t.delta_link_path[i,j,p]*model.x_p[p,xi] for p in self.ntw_t.P_set), 0)
		model.con_link_path = Constraint(self.ntw_t.A, self.scen.Xi, rule = con_link_path_rule)
		print("- Traffic path OD..")
		
		#traffic path OD constraints
		def con_path_od_rule(model,r,s,xi):
			return (0, sum(self.ntw_t.delta_path_od[p,r,s] * model.x_p[p,xi] for p in self.ntw_t.P_set) - self.ntw_t.Q[r,s],0)
		model.con_path_od = Constraint(self.ntw_t.R, self.ntw_t.S,self.scen.Xi,rule = con_path_od_rule)
		print("- Traffic non-anticipitivity...")
		
		# traffic non anti constraints; z=1 means haven't receive information when they reach link a following path p
		def con_non_anti_rule(model,r, s,i,j,xi):
			if sum((1-z[i,j, p])*self.ntw_t.delta_link_path[i,j,p]*self.ntw_t.delta_path_od[p,r,s] for p in self.ntw_t.P_set) > 0:
				return (0, sum((1-z[i,j,p])*self.ntw_t.delta_link_path[i,j,p]*self.ntw_t.delta_path_od[p,r,s]* model.x_p[p,xi] for p in self.ntw_t.P_set)   - model.x_a[r,s,i,j] ,0)
			else:
				return Constraint.Skip
		model.con_non_anti = Constraint(self.ntw_t.R, self.ntw_t.S, self.ntw_t.A, self.scen.Xi, rule = con_non_anti_rule)
		print("Select solver...")
		optsolver = SolverFactory(nlpsol)
		#optsolver = SolverFactory(mipsol)
		print("Create instance...")
		model.pprint()
		##print("Dual direction...")
		#inst.dual = Suffix(direction=Suffix.IMPORT)
		print("Start solving...")
		results = optsolver.solve(model, tee=True, load_solutions=True)
		v={}
		print("Solution: print v...")
		for (i,j) in self.ntw_t.A:
			for xi in self.scen.Xi:
				v[i,j,xi]=value(model.v[i,j,xi])
				print ("%i, %i, %i : %f" % (i, j, xi, v[i,j,xi]))
		x_p={}
		print("Solution: print x_p...")
		for p in self.ntw_t.P_set:
			for xi in self.scen.Xi:
				x_p[p, xi] = value(model.x_p[p, xi])
				print ("%i,%i : %f" % (p, xi, x_p[p, xi]))
		x_a={}
		print("Solution: print x_a...")
		for r in self.ntw_t.R:
			for s in self.ntw_t.S:
				for (i,j) in self.ntw_t.A:
					if sum((1-z[i,j, p])*self.ntw_t.delta_link_path[i,j,p]*self.ntw_t.delta_path_od[p,r,s] for p in self.ntw_t.P_set) > 0:
						x_a[r, s, i,j] = value(model.x_a[r,s,i,j])
						print ("%i,%i,%i,%i : %f" % (r, s, i,j, x_a[r,s,i,j]))
					else:
						x_a[r,s,i,j] = 0
		obj_val = value(model.obj)
		exp = {}
		exp1 = {}
		exp2 = {}
		for xi in self.scen.Xi:
			self.update_travel_time(xi)
			exp[xi] = sum(self.ntw_t.tff[r,s]*(self.ntw_t.c[r,s]*v[r,s,xi]+(self.ntw_t.b[r,s]/(self.ntw_t.alpha[r,s]+1.0))*(v[r,s,xi]**(self.ntw_t.alpha[r,s]+1))/(self.ntw_t.cap[r,s]**(self.ntw_t.alpha[r,s]))) for (r,s) in self.ntw_t.A)
			exp1[xi] = sum(self.ntw_t.tff[r,s]*self.ntw_t.b[r,s]*self.ntw_t.alpha[r,s]/(self.ntw_t.alpha[r,s]+1.0)*v[r,s,xi]**(self.ntw_t.alpha[r,s]+1)/(self.ntw_t.cap[r,s]**(self.ntw_t.alpha[r,s])) for (r,s) in self.ntw_t.A)
			exp2[xi] = 0
			for (i,j) in self.ntw_t.A:
				#if sum((1-z[i,j, p])*self.ntw_t.delta_link_path[i,j,p] for p in self.ntw_t.P_set) > 0:
				exp2[xi] = exp2[xi]+ sum(dual_non_anti[self.ntw_t.P[p][0][0], self.ntw_t.P[p][len(self.ntw_t.P[p])-1][1], i, j, xi] *( self.ntw_t.delta_link_path[i,j,p]*(1-self.model.z[i,j,p])*x_p[p,xi] - x_a[self.ntw_t.P[p][0][0], self.ntw_t.P[p][len(self.ntw_t.P[p])-1][1], i,j]) for p in self.ntw_t.P_set)
		return sum(self.scen.prob[xi]*(exp[xi] + exp1[xi]) + exp2[xi] for xi in self.scen.Xi)
	
	# relax second term, but use optimal solution from lower level: lower bound may be low, but should be better
	def calc_support_upper3(self, support_lower, v):
		def obj_rule():
			exp = {}
			for xi in self.scen.Xi:
				self.update_travel_time(xi)
				exp[xi] = sum(self.ntw_t.tff[r,s]*self.ntw_t.b[r,s]*self.ntw_t.alpha[r,s]/(self.ntw_t.alpha[r,s]+1.0)*v[r,s,xi]**(self.ntw_t.alpha[r,s]+1)/(self.ntw_t.cap[r,s]**(self.ntw_t.alpha[r,s])) for (r,s) in self.ntw_t.A)
			return sum(self.scen.prob[xi]*(exp[xi]) for xi in self.scen.Xi)
		obj_val = obj_rule()
		return support_lower + obj_val

# linearized bilinear term and follow the paper by Gao.
	def calc_support_upper4(self, z, dual2_linear, dual3_linear):
		print("Start calc_support_upper4...")
		model = ConcreteModel('support_upper4')
		# variables
		print("Variables...")
		# traffic flow on path p the connects r,s at secen
		model.x_p = Var(self.ntw_t.P_set,self.scen.Xi, within=NonNegativeReals, initialize = 1)
		#traffic flow on link a from r to s that haven't received infromation on scen
		model.x_a = Var(self.ntw_t.R, self.ntw_t.S, self.ntw_t.A, within=NonNegativeReals, initialize = 1)
		# traffic flow of link a
		model.v = Var(self.ntw_t.A, self.scen.Xi, within=NonNegativeReals, initialize=1.0)
		# m = (1-z)x_p
		model.m = Var(self.ntw_t.A, self.ntw_t.P_set, self.scen.Xi, within=NonNegativeReals, initialize=1.0)
		print("Objective...")
		
		def obj_rule(model):
			exp = {}
			exp1 = {}
			exp2 = {}
			for xi in self.scen.Xi:
				self.update_travel_time(xi)
				exp[xi] = sum(self.ntw_t.tff[r,s]*(self.ntw_t.c[r,s]*model.v[r,s,xi]+(self.ntw_t.b[r,s]/(self.ntw_t.alpha[r,s]+1.0))*(model.v[r,s,xi]**(self.ntw_t.alpha[r,s]+1))/(self.ntw_t.cap[r,s]**(self.ntw_t.alpha[r,s]))) for (r,s) in self.ntw_t.A)
				exp1[xi] = sum(self.ntw_t.tff[r,s]*self.ntw_t.b[r,s]*self.ntw_t.alpha[r,s]/(self.ntw_t.alpha[r,s]+1.0)*model.v[r,s,xi]**(self.ntw_t.alpha[r,s]+1)/(self.ntw_t.cap[r,s]**(self.ntw_t.alpha[r,s])) for (r,s) in self.ntw_t.A)
				exp2[xi] = sum(dual2_linear[i,j,p,xi]*(-model.m[i,j,p,xi] + model.x_p[p,xi]-z[i,j,p]*self.M) + dual3_linear[i,j,p,xi]*(model.m[i,j,p,xi] - (1-z[i,j,p])*self.M) for (i,j) in self.ntw_t.A for p in self.ntw_t.P_set)
			return sum(self.scen.prob[xi]*(exp[xi] + exp1[xi]) + exp2[xi] for xi in self.scen.Xi)
		model.obj = Objective(rule=obj_rule, sense=minimize)
		print("Constraints...")
		print("- Traffic flow on link..")
		
		# traffic flow path constraints
		def con_link_path_rule(model,i,j, xi):
			return (0, model.v[i,j, xi] - sum(self.ntw_t.delta_link_path[i,j,p]*model.x_p[p,xi] for p in self.ntw_t.P_set), 0)
		model.con_link_path = Constraint(self.ntw_t.A, self.scen.Xi, rule = con_link_path_rule)
		print("- Traffic path OD..")
		
		#traffic path OD constraints
		def con_path_od_rule(model,r,s,xi):
			return (0, sum(self.ntw_t.delta_path_od[p,r,s] * model.x_p[p,xi] for p in self.ntw_t.P_set) - self.ntw_t.Q[r,s],0)
		model.con_path_od = Constraint(self.ntw_t.R, self.ntw_t.S,self.scen.Xi,rule = con_path_od_rule)
		print("- Traffic non-anticipitivity...")
		
		# traffic non anti constraints; z=1 means haven't receive information when they reach link a following path p
		def con_non_anti_rule(model,r, s,i,j,xi):
			if sum((1-z[i,j, p])*self.ntw_t.delta_link_path[i,j,p]*self.ntw_t.delta_path_od[p,r,s] for p in self.ntw_t.P_set) > 0:
				return (0, sum(self.ntw_t.delta_link_path[i,j,p]*self.ntw_t.delta_path_od[p,r,s]* model.m[i,j,p,xi] for p in self.ntw_t.P_set)   - model.x_a[r,s,i,j] ,0)
			else:
				return Constraint.Skip
		model.con_non_anti = Constraint(self.ntw_t.R, self.ntw_t.S, self.ntw_t.A, self.scen.Xi, rule = con_non_anti_rule)
		
		def con1_linear_rule(model,i,j,p,xi):
			#return Constraint.Skip
			return (-infinity, model.m[i,j,p,xi] - model.x_p[p, xi], 0)
		model.con1_linear = Constraint(self.ntw_t.A, self.ntw_t.P_set,self.scen.Xi,rule = con1_linear_rule)
		print("Select solver...")
		optsolver = SolverFactory(nlpsol)
		#optsolver = SolverFactory(mipsol)
		print("Create instance...")
		#print("Dual direction...")
		#model.dual = Suffix(direction=Suffix.IMPORT)
		print("Start solving...")
		results = optsolver.solve(model, tee=True, load_solutions=True)
		v={}
		print("Solution: print v...")
		for (i,j) in self.ntw_t.A:
			for xi in self.scen.Xi:
				v[i,j,xi]=value(model.v[i,j,xi])
				print ("%i, %i, %i : %f" % (i, j, xi, v[i,j,xi]))
		x_p={}
		bi_x_p={}
		print("Solution: print x_p...")
		for p in self.ntw_t.P_set:
			for xi in self.scen.Xi:
				x_p[p, xi] = value(model.x_p[p, xi])
				if x_p[p,xi] <=self.epsilon:
					bi_x_p[p,xi] = 0
				else:
					bi_x_p[p,xi] = 1
				print ("%i,%i : %f" % (p, xi, x_p[p, xi]))
		x_a={}
		bi_x_a={}
		print("Solution: print x_a...")
		for r in self.ntw_t.R:
			for s in self.ntw_t.S:
				for (i,j) in self.ntw_t.A:
					if sum((1-z[i,j, p])*self.ntw_t.delta_link_path[i,j,p]*self.ntw_t.delta_path_od[p,r,s] for p in self.ntw_t.P_set) > 0:
						x_a[r, s, i,j] = value(model.x_a[r,s,i,j])
						if x_a[r, s, i,j] <=self.epsilon:
							bi_x_a[r, s, i,j] = 0
						else:
							bi_x_a[r, s, i,j] = 1
						print ("%i,%i,%i,%i : %f" % (r, s, i,j, x_a[r,s,i,j]))
					else:
						x_a[r, s, i,j] =0
						bi_x_a[r,s,i,j] = 0
		m = {}
		print("Solution: print dual1_linear...")
		for (i,j) in self.ntw_t.A:
			for p in self.ntw_t.P_set:
				for xi in self.scen.Xi:
					m[i,j,p,xi] = value(model.m[i,j,p,xi])
		exp = {}
		exp1 = {}
		exp2 = {}
		for xi in self.scen.Xi:
			self.update_travel_time(xi)
			exp[xi] = sum(self.ntw_t.tff[r,s]*(self.ntw_t.c[r,s]*v[r,s,xi]+(self.ntw_t.b[r,s]/(self.ntw_t.alpha[r,s]+1.0))*(v[r,s,xi]**(self.ntw_t.alpha[r,s]+1))/(self.ntw_t.cap[r,s]**(self.ntw_t.alpha[r,s]))) for (r,s) in self.ntw_t.A)
			exp1[xi] = sum(self.ntw_t.tff[r,s]*self.ntw_t.b[r,s]*self.ntw_t.alpha[r,s]/(self.ntw_t.alpha[r,s]+1.0)*v[r,s,xi]**(self.ntw_t.alpha[r,s]+1)/(self.ntw_t.cap[r,s]**(self.ntw_t.alpha[r,s])) for (r,s) in self.ntw_t.A)
			exp2[xi] = sum(dual2_linear[i,j,p,xi]*(-m[i,j,p,xi] + x_p[p,xi]-self.model.z[i,j,p]*self.M) + dual3_linear[i,j,p,xi]*(m[i,j,p,xi] - (1-self.model.z[i,j,p]*self.M)) for (i,j) in self.ntw_t.A for p in self.ntw_t.P_set)
		return sum(self.scen.prob[xi]*(exp[xi] + exp1[xi]) + exp2[xi] for xi in self.scen.Xi)

	def calc_level_headway(self, level):
		if level ==0:
			LOA = 0
			h1 = 0
		if level ==1:
			LOA = 0.2
			h1 = 0.4
		if level ==3:
			LOA = 0.4
			h1 = 0.5
		if level ==4:
			LOA = 0.6
			h1 = 0.6
		if level ==5:
			LOA = 0.8
			h1 = 0.7
		if level ==6:
			LOA = 1
			h1 = 0.8
		return LOA, h1

	
	def algorithm(self):
		# solve for the initialization
		y = {}
		z = {}
		v = {}
		v_rs = {}
		x_p = {}
		x_a = {}
		dual_path_od = {}
		dual_non_anti = {}
		bi_v_p = {}
		bi_v_a = {}
		obj_approx = 0
		c1 = 0.5
		c2 = 0.7
		c3 = 0.7
		c4 = 1.3
		gap = 0.005
		mipgap = 0.005
		# Step 1: initialize y, z upper level decisions
		y = defaultdict(dict)
		z = defaultdict(dict)
		obj_lb = defaultdict(dict)
		obj_ub = defaultdict(dict)
		temp = 0
		# TODO change to randomly select n from N
		# y: if information is shared at node i
		'''for i in self.ntw_t.N:
									temp = temp+1
									if temp <= self.max_info_nodes:
										y[0][i] = 1
									else:
										y[0][i] = 0'''
		for i in self.ntw_t.N:
			y[0][i] = 0
		
		'''
		cases = {'c0':[],'c1': [1],'c2': [2],'c3': [3],'c4': [4],'c5': [5],'c6': [6],'c7': [7],'c8': [8],
		'c9': [9],'c10': [10],'c11': [11],'c12': [12],'c13': [13],'c14': [14],'c15': [15],'c16': [16],
		'c17': [17],'c18': [18],'c19':[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18]}'''
		#cases = {'c0':[],'c1': [1],'c2': [2],'c3': [3],'c4': [4]}
		cases = {'c0':[],'c1': [1],'c2': [2],'c3': [3],'c4': [4],'c14': [14],'c15': [15],'c16': [16],
		'c17': [17],'c18': [18],'c19':[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18]}

		v_cap = {}
		v_cap2={}

		#form a relationship between LOA and headway and p_r
		LOA, h1= self.calc_level_headway(0)

		# half automated vehicles
		#LOA = 0.5 # level 4
		#h1 = 0.5

		#if LOA is '0' then p_r=0
		if LOA == 0:
			p_r = 0



		
		for key,value in cases.items():
			for i in self.ntw_t.N:
				y[0][i] = 0
			for n in value:
				y[0][n]=1
			z[0] = self.calc_z(y[0])
			print("y[0]", y[0])
			#print("z[0]", z[0])
			#v,v2, x_p, x_a, m, dual1_linear, dual2_linear, dual3_linear,v_cap,v_cap2, crl_1,crl_2= self.uer_fix_upper(z[0])
			#v,v2, x_p, x_a, bi_x_p, bi_x_a,v_cap,v_cap2,crl_1,crl_2
			#v,v2, x_p, x_a, m, dual1_linear, v_cap,v_cap2, crl_1,crl_2= self.uer_fix_upper(z[0])
			v, x_p, x_a, bi_x_p, bi_x_a,v_cap,crl_1,crl_2= self.uer_fix_upper(z[0])
			#v,v2, x_p_h,x_p_c, x_a, dual_path_od_hdv,dual_path_od_cav, dual_non_anti, bi_x_p_h,bi_x_p_c, bi_x_a,v_cap,v_cap2,crl_1,crl_2 = self.check(z[0], LOA, h1,p_r,share)
			obj_ub[key] = self.calc_obj_upper(v)
			print(obj_ub[key])
		return  obj_ub, v,x_p,x_a

if __name__ == "__main__":
	#ntw_t = IO.four_node_transport_braess()
	#scen = IO.scen_four_braess()
	#ntw_t = IO.four_node_safer_sim()
	#scen = IO.scen_four_node_safer_sim(ntw_t)
	ntw_t = IO.Orlando()
	scen = IO.scen_Orlando(ntw_t)
	#ntw_t = IO.Sioux_Fall()
	#scen = IO.scen_Sioux_Fall()
	prob = Problem(ntw_t,scen)
	time_bq = {}
	start = time.time()
	tol = 0.1
	print ('Stopping critieria %f' % tol)
	Maxit = 1
	aux_str=strftime("%Y%m%d_%H%M", localtime())
	res_path=''+aux_str
	os.system('mkdir '+res_path)
	file_name = res_path+'/rho.txt'
	file1 = open(file_name,"w")
	for iter in range(Maxit):
		time_bq[iter] = time.time()
		print ('Start iteration %i\t' % iter)
		IO.obj, IO.v,IO.x_p,IO.x_a = prob.algorithm()
		print(IO.obj)
		#print(IO.v)
		#IO.write_file(file1, IO.rho)
	file1.close() #to change file access modes
	end = time.time()
	el_time = end - start
	print ('Elapsed time %f' % el_time)
