import osm2gmns as og
import path4gmns as pg
import subprocess
from time import time
import os

osm05 = 'map05.osm'
osm06 = 'map06.osm'
target = 'map.osm'

osmconvert = 'osmconvert'
osmosis = 'osmosis'

time_log = True
og.og_settings.verbose = False

networks = [
    {
        'name': 'UCF',
        'id': 1159228,
        # zone and demand generation parameters
        'grid_dimension': 8,
        'total_demand': 10_000,
        'time_budget': 120
    },

    {
        'name': 'Ocala',
        'id': 119060,
        # zone and demand generation parameters
        'grid_dimension': 16,
        'total_demand': 40_000,
        'time_budget': 360
    },

    {
        'name': 'Orlando',
        'id': 1128379,
        # zone and demand generation parameters
        'grid_dimension': 36,
        'total_demand': 60_000,
        'time_budget': 360
    },
]

def prepare_OSM_data(id):
    st = time()

    og.downloadOSMData(id, output_filename=osm05)

    try:
        subprocess.run(
            [
                osmconvert,
                osm05,
                '--fake-author',
                f'-o={osm06}'
            ],
            timeout=10
        )

        print(f'map data has been converted from OSM 0.5 to OSM 0.6 and written to {osm06}')
    except FileNotFoundError as e:
        print(f'Process failed because the executable {osmconvert} could not be found.\n{e}')
    except subprocess.CalledProcessError as e:
        print(
            f'Process failed, return code was not successful'
            f'Returned {e.returncode}\n{e}'
        )
    except subprocess.TimeoutExpired as e:
        print(f'Process timed out\n{e}')

    try:
        subprocess.run(
            [
                osmosis,
                '--read-xml',
                osm06,
                '--tf',
                'accept-ways',
                'highway=*',
                '--used-node',
                '--write-xml',
                target
            ]
        )

        print(f'OSM 0.6 map data has been filtered to include highways only and written to {target}')
    except FileNotFoundError as e:
        print(f'Process failed because the executable {osmosis} could not be found.\n{e}')
    except subprocess.CalledProcessError as e:
        print(
            f'Process failed, return code was not successful'
            f'Returned {e.returncode}\n{e}'
        )

    if time_log:
        print(f'Prepare OSM Data...\tProcessing Time:\t{time()-st:.2f} seconds')


def prepare_POI_data(network_name='gmns'):    
    network = og.getNetFromFile(filename=osm05, POI=True)
    og.connectPOIWithNet(network)
    og.outputNetToCSV(network, output_folder=network_name)


def prepare_filtered_data(network_name='gmns'):
    network = og.getNetFromFile(filename=target)
    og.consolidateComplexIntersections(network, auto_identify=True)
    og.outputNetToCSV(network, output_folder=network_name)


def prepare_demand_data(network_name='gmns', grid_dimension=8, total_demand=10_000, time_budget=120):

    print(f'''

netork name:    {network_name}
grid dimensiom: {grid_dimension}
total demand:   {total_demand}
time budget:    {time_budget}

''')

    pwd = os.getcwd()
    # print(network_name)
    os.chdir(network_name)

    network = pg.read_network()

    pg.network_to_zones(network, 
                        grid_dimension=grid_dimension, 
                        total_demand=total_demand, 
                        time_budget=time_budget
                        )

    pg.output_zones(network)
    pg.output_synthesized_demand(network)

    os.chdir(pwd)

    demand = network.get_ODMatrix()
    print(len(demand))


def DTALite_perform_UE_DTA(network_name='gmns', mode=1):
    st = time()
    pwd = os.getcwd()
    os.chdir(network_name)

    # mode = 1
    column_gen_num = 10
    column_update_num = 10

    pg.perform_network_assignment_DTALite(mode, 
                                          column_gen_num, 
                                          column_update_num
                                          )

    os.chdir(pwd)

    if time_log:
        print(f'Perform UE and DTA...\tProcessing Time:\t{time()-st:.2f} seconds')


def perform_UE_DTA(network_name='gmns'):
    st = time()
    pwd = os.getcwd()
    os.chdir(network_name)

    network = pg.read_network()

    pg.read_zones(network)
    pg.load_demand(network)

    # UE + DTA
    column_gen_num = 10
    column_update_num = 10
    pg.perform_column_generation(column_gen_num, column_update_num, network)
    pg.perform_simple_simulation(network)
    print('complete dynamic simulation')

    print('writing agent trajectories')
    pg.output_agent_trajectory(network)

    if time_log:
        print(f'Perform UE and DTA...\tProcessing Time:\t{time()-st:.2f} seconds')


def perform_user_equilibrium(network_name='gmns'):
    st = time()
    pwd = os.getcwd()
    os.chdir(network_name)

    network = pg.read_network()
    pg.read_zones(network)
    pg.load_demand(network)
    
    column_gen_num = 10
    column_update_num = 10

    pg.perform_column_generation(column_gen_num, column_update_num, network)

    pg.output_columns(network)
    pg.output_link_performance(network)

    if time_log:
        print(f'Finished path-based user equilibrium\n'
              f'Processing time: {time()-st:.2f}')

def main():

    network = networks[2]
    network['name'] += '_test'

    # prepare_directory(network['name'])
    # prepare_OSM_data(network['id'])

    st = time()
    # prepare_POI_data(network['name'])
    # prepare_demand_data(network['name'], 
    #                     network['grid_dimension'], 
    #                     network['total_demand'], 
    #                     network['time_budget']
    #                     )
    print(f'Convert to GMNS...\tProcessing Time:\t{time()-st:.2f} seconds')

    # prepare_filtered_data(network['name'])
    perform_UE_DTA(network['name'])


if __name__ == '__main__':
    main()