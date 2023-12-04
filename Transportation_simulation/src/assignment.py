'''
Network assignment

methods to carry out traffic assignment and user equilibrium
'''


from time import time
import path4gmns as pg
from classes import NetworkDescriptor
from consts import  *

__all__ = ['user_equilibrium', 'user_equilibrium_DTA', 'DTALite_UE_DTA']


def user_equilibrium(descriptor: NetworkDescriptor):
    '''
    Perform user equilibrium

    agent and link performance data will be output in GMNS format
    '''
    st = time()

    network = pg.read_network(length_unit=descriptor.length_unit,
                              speed_unit=descriptor.speed_unit)
    
    pg.read_zones(network)
    pg.load_demand(network)

    pg.perform_column_generation(column_gen_num=descriptor.column_gen_num,
                                 column_update_num=descriptor.column_update_num,
                                 ui=network)
    
    pg.output_columns(network)
    pg.output_link_performance(network)

    if time_log:
        print(f'Path-based user equilibrium processing time: {time()-st:.2f} seconds')


def user_equilibrium_DTA(descriptor: NetworkDescriptor):
    '''
    Perform user equilibrium and dynamic traffic assignment

    trajectory data will be output in GMNS format
    '''
    st = time()

    network = pg.read_network(length_unit=descriptor.length_unit,
                              speed_unit=descriptor.speed_unit)
    
    pg.read_zones(network)
    pg.load_demand(network)

    pg.perform_column_generation(column_gen_num=descriptor.column_gen_num,
                                 column_update_num=descriptor.column_update_num,
                                 ui=network)
    pg.perform_simple_simulation(network)
    print('complete dynamic simulation')

    print('writing agent trajectories')
    pg.output_agent_trajectory(network)

    if time_log:
        print(f'User equilibrium and DTA processing time: {time()-st:.2f} seconds')


def DTALite_UE_DTA(descriptor: NetworkDescriptor, mode=1):
    ''''
    Perform traffic assignment using Path4gmns' DTALite API

    A settings.csv file is required.
    agent and link performance data will be output in GMNS format
    '''
    st = time()

    pg.perform_network_assignment_DTALite(mode,
                                          column_gen_num=descriptor.column_gen_num,
                                          column_update_num=descriptor.column_update_num)

    if time_log:
        print(f'User equilibrium and DTA through DTALite processing time: {time()-st:.2f} seconds')