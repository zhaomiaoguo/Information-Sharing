import os
import sys
import numpy as np
import pandas as pd
import openmatrix as omx

def import_matrix(matfile):
    # Function to import OMX matrices
    f = open(matfile, 'r')
    all_rows = f.read()
    blocks = all_rows.split('Origin')[1:]
    matrix = {}
    for k in range(len(blocks)):
        orig = blocks[k].split('\n')
        dests = orig[1:]
        orig=int(orig[0])

        d = [eval('{'+a.replace(';',',').replace(' ','') +'}') for a in dests]
        destinations = {}
        for i in d:
            destinations = {**destinations, **i}
        matrix[orig] = destinations
    zones = max(matrix.keys())
    mat = np.zeros((zones, zones))
    for i in range(zones):
        for j in range(zones):
            # We map values to a index i-1, as Numpy is base 0
            mat[i, j] = matrix.get(i+1,{}).get(j+1,0)

    index = np.arange(zones) + 1

    myfile = omx.open_file('demand.omx','w')
    myfile['matrix'] = mat
    myfile.create_mapping('taz', index)
    myfile.close()
    
def transportation_network_topo(network_name):
    # Assume you have not changed the relative paths after cloning the repository
    root = os.path.dirname(os.getcwd())
    
    
    
    # Importing the networks into a Pandas dataframe consists of a single line of code
    # but we can also make sure all headers are lower case and without trailing spaces
    
    netfile = os.path.join(root,'TransportationNetworks', network_name, network_name+'_net.tntp')
    net = pd.read_csv(netfile, skiprows=8, sep='\t')

    trimmed= [s.strip().lower() for s in net.columns]
    net.columns = trimmed

    # And drop the silly first andlast columns
    net.drop(['~', ';'], axis=1, inplace=True)
    n = max(max(net.init_node),max(net.term_node))
    N = list(range(1,n+1))
    A=[]
    cap={}
    tff={}
    for i in range(net.shape[0]):
        tmp = (int(net.iloc[i].init_node),int(net.iloc[i].term_node))
        A.append(tmp)
        cap[int(net.iloc[i].init_node),int(net.iloc[i].term_node)] = net.iloc[i].capacity
        tff[int(net.iloc[i].init_node),int(net.iloc[i].term_node)] = net.iloc[i].free_flow_time
        
    return(N,A,cap,tff)
#N, A, capp, tff = transportation_network_topo('Anaheim')
#print('Success')

