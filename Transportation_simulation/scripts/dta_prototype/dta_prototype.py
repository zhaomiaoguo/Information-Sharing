'''
Small script that prepares network data and performs dynamic traffic assign assignment

Command line arguments:


Required data:

'''

import osm2gmns as og
import path4gmns as pg
import argparse
import os
from os.path import exists

import subprocess

# Note: filtering data with osmosis clashes with the POI data generated (links and POIs)
def prepare_OSM_data(network_id: int, gmns_target_dir: str):
    '''
    Download data from OSM and convert to GMNS format 
    : param network_id : a network id from Open Street Map
    : param gmns_target_dir : a path to a directory to store the GMNS data
    : return : None
    '''

    osm5_target = 'map_05.osm'
    osm6_target = 'map_06.osm'
    filtered_target = 'map_filtered.osm'

    osmconvert = 'osmconvert'
    osmosis = 'osmosis'

    og.downloadOSMData(area_id=network_id, output_filename=osm5_target)

    # call osmconvert to add dummy data to osm5_target
    # this step is necessary as osmosis relies on that data to be present
    try:
        subprocess.run(
            [
                osmconvert,
                osm5_target,
                '--fake-author',
                '-o='+osm6_target
            ],
            timeout=10,
            check=True
        )
        print(f'map data has been converted from OSM 0.5 to OSM 0.6 and written to {osm6_target}')
    except FileNotFoundError as e:
        print(f'Process failed because the executable {osmconvert} could not be found.\n{e}')
    except subprocess.CalledProcessError as e:
        print(
            f'Process failed, return code was not successful'
            f'Returned {e.returncode}\n{e}'
        )
    except subprocess.TimeoutExpired as e:
        print(f'Process timed out\n{e}')        

    # call osmosis to filter data in Linux and MacOS
    if os.name != 'nt':
        try:
            subprocess.run(
                [
                    osmosis,
                    '--read-xml',
                    osm6_target,
                    '--tf',
                    'accept-ways',
                    'highway=*',
                    '--used-node',
                    '--write-xml',
                    filtered_target
                ],
                timeout=10,
                check=True
            )
            print(f'OSM 0.6 map data has been filtered to include highways only and written to {filtered_target}')
        except FileNotFoundError as e:
            print(f'Process failed because the executable {osmosis} could not be found.\n{e}')
        except subprocess.CalledProcessError as e:
            print(
                f'Process failed, return code was not successful'
                f'Returned {e.returncode}\n{e}'
            )
        except subprocess.TimeoutExpired as e:
            print(f'Process timed out\n{e}')

    # load filtered data
    # consolidate intesections
    # output network to GMNS format
    if os.name == 'nt': network = og.getNetFromFile(filename=osm6_target, POI=True)
    else: network = og.getNetFromFile(filename=filtered_target)
    # network = og.getNetFromFile(filename=osm6_target, POI=True)

    # POIs need to be connected to generate significant demand data
    # og.connectPOIWithNet(network)
    og.consolidateComplexIntersections(network, auto_identify=True)
    og.outputNetToCSV(network, output_folder=gmns_target_dir)


def perform_UE_DTA(gmns_dir: str):
    '''
    Perform User Equilibrium and Dynamic Traffic Assignment using DTALite
    : param gmns_dir : a path to the directory containing the GMNS data
    : return : None
    '''

    pwd = os.getcwd()
    os.chdir(gmns_dir)

    network = pg.read_network()

    # autogenerate zones and demand data if they are not provided
    if not exists('zone.csv') or not exists('demand.csv'):
        pg.network_to_zones(network)
        pg.output_zones(network)
        pg.output_synthesized_demand(network)

    mode = 1
    column_gen_num = 10
    column_update_num = 10

    pg.perform_network_assignment_DTALite(mode, column_gen_num, column_update_num)

    os.chdir(pwd)


def init_parser() -> argparse.ArgumentParser:
    '''
    Initialize a parser object and set up command line arguments
    : return : parser object
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--download', action='store_true')
    parser.add_argument('--network_id', default=1159228)
    parser.add_argument('--gmns_dir', default='gmns')

    return parser


def main():
    parser = init_parser()
    args = parser.parse_args()

    download = args.download
    network_id = int(args.network_id)
    gmns_dir = args.gmns_dir

    if download: prepare_OSM_data(network_id, gmns_dir)
    perform_UE_DTA(gmns_dir)


if __name__ == '__main__':
    main()