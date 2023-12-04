'''
Small script that downloads an OSM network and filters it to include only highway data,
using the osmconvert and osmosis command line tools.

Notes: osmconvert doesn't seem to have a MacOS binary version, 
       though the source code can be compiled in MacOS without issues
       [https://wiki.openstreetmap.org/wiki/Osmconvert#Source]
       [better compilation command: $ cc osmconvert.c -O3 -o osmconvert]
       osmosis can't be called from the subprocess module from Windows

Compatibility Status:
Windows:    osmconvert runs without issues, but calling osmosis returns a FileNotFoundError
MacOS:      Script runs, but you need to compile osmconvert in your system first
Linux:      Script runs without issues
'''

import osm2gmns as og
import subprocess
from time import time

networks = [
    {'id': 1159228, 'name': 'UCF'},
    {'id': 119060, 'name': 'Ocala'},
    {'id': 1128379, 'name': 'Orlando'}
]

def prepare_data(area_id):
    '''
    : param area_id : The network's OSM id 
    : return : A string detailing the processing time in seconds
    '''

    target_05 = 'map05.osm'
    target_06 = 'map06.osm'
    highways_target = 'map.osm'

    # change this to the names of the binaries in your system
    # including the absolute path if needed
    osmconvert = 'osmconvert'
    osmosis = 'osmosis'

    st = time()
    og.downloadOSMData(area_id=area_id, output_filename=highways_target)

    # osmconvert: add dummy data to network data downloaded with the overpass API
    try:
        subprocess.run(
            [
                osmconvert,
                target_05,
                '--fake-author',
                '-o='+target_06
            ],
            timeout=10,
            check=True
        )

        print(f'map data has been converted from OSM 0.5 to OSM 0.6 and written to {target_06}')
    except FileNotFoundError as e:
        print(f'Process failed because the executable {osmconvert} could not be found.\n{e}')
    except subprocess.CalledProcessError as e:
        print(
            f'Process failed, return code was not successful'
            f'Returned {e.returncode}\n{e}'
        )
    except subprocess.TimeoutExpired as e:
        print(f'Process timed out\n{e}')

    # osmosis: filter data to contain highways info only
    try:
        subprocess.run(
            [
                osmosis,
                '--read-xml',
                target_06,
                '--tf',
                'accept-ways',
                'highway=*',
                '--used-node',
                '--write-xml',
                highways_target
            ],
            # 10 seconds is not enough for larger networks, such as Orlando
            # timeout=10,
            check=True
        )
        print(f'OSM 0.6 map data has been filtered to include highways only and written to {highways_target}')
    except FileNotFoundError as e:
        print(f'Process failed because the executable {osmosis} could not be found.\n{e}')
    except subprocess.CalledProcessError as e:
        print(
            f'Process failed, return code was not successful'
            f'Returned {e.returncode}\n{e}'
        )
    except subprocess.TimeoutExpired as e:
        print(f'Process timed out\n{e}')

    network = og.getNetFromFile()
    og.outputNetToCSV(network, output_folder='gmns')

    return time()-st

def convert_to_gmns():
    '''
    '''

    st = time()
    # network = og.getNetFromFile(filename='map_highways.osm')
    network = og.getNetFromFile()

    # og.connectPOIWithNet(network)
    # og.show(network)

    # og.consolidateComplexIntersections(network, auto_identify=True)
    og.outputNetToCSV(network, output_folder='gmns')

    return time()-st


def time_all_networks(iterations):
    records = {}
    for net in networks:
        records[net['name']] = []

    print(records)
    for _ in range(iterations):
        # records.append([])

        for net in networks:
            name = net['name']
            id = net['id']

            prepare_time = prepare_data(id)
            # convert_time = convert_to_gmns()
            records[name].append(prepare_time)

    print(records)

# def main():
#     with open('log.txt', 'w') as file:

#         # iterations = 10
#         for net in networks:
#             name = net['name']
#             id = net['id']

#             # iterations = 10
#             prepare_time_list = []
#             # convert_time_list = []
            
#             file.write(f'network: {name}\t\tid: {id}\n\n')

#             # for _ in range(iterations):
#             prepare_time = prepare_data(id)
#             # convert_time = convert_to_gmns()

#             prepare_time_list.append(prepare_time)
#             # convert_time_list.append(convert_time)
        
#             file.write(f'data filtering time: {prepare_time} seconds\n')
#                         # f'data conversion time: {convert_time} seconds\n\n')
                
#             # file.write(f'data filtering average time: {sum(prepare_time_list) / iterations}\n'
#             #           f'data conversion average time: {sum(convert_time_list) / iterations}')
#     # time_all_networks(2)

def main():
    # net_id = 159228
    # prepare_data(net_id)
    time_all_networks(2)
            
if __name__ == '__main__':
    main()