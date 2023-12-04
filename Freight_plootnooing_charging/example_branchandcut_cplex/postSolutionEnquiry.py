#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 12 21:46:50 2022

@author: rakibulalam
"""



# import module
import docplex.mp
from docplex.mp import model_reader

#readl .lp file
m = model_reader.ModelReader.read('testModel150_10_10_10.lp', ignore_names=True)

# strategy.search =1 for conventional branch and cut; 2: for dynamic
# ref: http://www-eio.upc.edu/lceio/manuals/cplex-11/html/refparameterscplex/refparameterscplex68.html

m.parameters.mip.strategy.search = 1
m.parameters.mip.display = 2

# print model information
m.print_information()
# solve and see progres of convergence
msol = m.solve(clean_before_solve=True, log_output =True)
assert msol is not None, "model can't solve"
m.report()

# print optimal values of variables
import pandas as pd
outputs = {}
print("status:", msol.solve_details.status)                                    
if msol:                                                                         
    print("obj val:", msol.get_objective_value())
    data = []
    for var,value in msol.iter_var_values():
        data.append([var,value])
    df_sol = pd.DataFrame(data=data, columns=['var', 'value'])
    outputs['solution'] = df_sol

# find objective value: this value is avaialble in solution
msol.get_objective_value()





