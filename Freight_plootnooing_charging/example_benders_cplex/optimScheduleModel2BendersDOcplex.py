#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 12 21:46:50 2022

@author: rakibulalam
"""



# write in docplex to call benders in CPPLEX new version
# import module
import docplex.mp
from docplex.mp import model_reader

#readl .lp file
m = model_reader.ModelReader.read('testModel.lp', ignore_names=True)

# select benders strategy, 3= Full (integer in master problem, all continuous in one subproblem)
m.parameters.benders.strategy = 3

# print model information
m.print_information()

# solve and see progres of convergence
msol = m.solve(clean_before_solve=True, log_output =True)
assert msol is not None, "model can't solve"
m.report()
