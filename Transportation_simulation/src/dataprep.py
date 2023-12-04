'''
Data preparation

Methods to obtain and synthesize network data
'''


from time import time
import subprocess
import osm2gmns as og
import path4gmns as pg
from classes import NetworkDescriptor
from consts import  *


__all__ = ['prep_OSM_network', 'filter_OSM_network', 'extract_POI_data', 'filtered_data_to_csv', 'synthesize_zones_demand']


def prep_OSM_network(descriptor: NetworkDescriptor):
    ''' 
    Download network data from OpenStreetMap in OSM format and add dummy data to it

    the dummy data will be added using osmconvert, a command line tool.  
    osmconvert installation: https://wiki.openstreetmap.org/wiki/Osmconvert#Download

    After installing, make sure to update the osmconvert global variable in the consts 
    module to match the executable's name (and path if necessary)
    '''
    st = time()

    og.downloadOSMData(descriptor.network_id, osm05)

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
            f'Process failed, return code was not successful\n'
            f'Returned {e.returncode}\n{e}'
        )
    except subprocess.TimeoutExpired as e:
        print(f'Process timed out\n{e}')

    descriptor.set_osm_target(osm06)

    if time_log:
        print(f'OSM download processing time: {time()-st:.2f} seconds')


def filter_OSM_network(descriptor: NetworkDescriptor):
    '''
    Filter network data in OSM network to include highway data only

    The filtering will be done with osmosis, a command line tool.
    osmosis installation: https://wiki.openstreetmap.org/wiki/Osmosis#How_to_install

    After installing, make sure to update the osmosis global variable in the consts
    module to match the executable's name (and path if necessary)

    Note: calling osmosis from Python's subprocess module results in a FileNotFoundError
    in Windows. In this case, you can skip using this function, although you will be 
    working with an unfiltered network.
    '''
    st = time()
    
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
                osm_filtered
            ]
        )

        print(f'OSM 0.6 map data has been filtered to include highways only and written to {osm_filtered}')
    except FileNotFoundError as e:
        print(f'Process failed because the executable {osmosis} could not be found.\n{e}')
    except subprocess.CalledProcessError as e:
        print(
            f'Process failed, return code was not successful'
            f'Returned {e.returncode}\n{e}'
        )

    descriptor.set_osm_target(osm_filtered)

    if time_log:
        print(f'OSM data filtering processing time: {time()-st:.2f} seconds')


def extract_POI_data(descriptor: NetworkDescriptor):
    '''
    Extract Point of Interest data from the OSM network

    The POIs will be connected with the rest of the network.
    The link, node and poi data will be output in GMNS format.
    '''
    st = time()

    network = og.getNetFromFile(descriptor.osm_target, POI=True)
    og.connectPOIWithNet(network)
    og.outputNetToCSV(network)

    if time_log:
        print(f'POI data extraction processing time: {time()-st:.2f} seconds')


def filtered_data_to_csv(descriptor: NetworkDescriptor):
    '''
    Consolidate the filtered OSM network and output it in GMNS format

    The network's link and node data will be output.
    '''
    st = time()

    network = og.getNetFromFile(descriptor.osm_target, combine=True)
    og.consolidateComplexIntersections(network, auto_identify=True)
    og.outputNetToCSV(network)

    if time_log:
        print(f'filtered data conversion to gmns processing time: {time()-st:.2f} seconds')


def synthesize_zones_demand(descriptor: NetworkDescriptor):
    '''
    Synthesize zone and demand data for a network, and output them in GMNS format
    '''
    st = time()

    network = pg.read_network(length_unit=descriptor.length_unit,
                              speed_unit=descriptor.speed_unit)
    
    pg.network_to_zones(network,
                        grid_dimension=descriptor.grid_dimension,
                        total_demand=descriptor.total_demand,
                        time_budget=descriptor.time_budget)
    
    pg.output_zones(network)
    pg.output_synthesized_demand(network)

    demand = network.get_ODMatrix()
    if len(demand) == 0:
        print(f'Warning! demand data was not successfuly generated.\n'
              f'Consider increasing the values in grid_dimension, total_demand, and/or time_budget')

    if time_log:
        print(f'demand and zone data synthesis processing time: {time()-st:.2f} seconds')
