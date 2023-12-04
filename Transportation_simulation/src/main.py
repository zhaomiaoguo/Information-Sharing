'''
Description Pending
'''

import os
from consts import *
from classes import *
from dataprep import *
from assignment import *


networks = {
    'UCF'       : NetworkDescriptor(network_id=1159228,
                                    input_directory='test_data/UCF',
                                    name='UCF',
                                    length_unit='mile',
                                    speed_unit='mph'),

    'Ocala'     : NetworkDescriptor(network_id=119060,
                                    input_directory='test_data/Ocala',
                                    name='Ocala',
                                    length_unit='mile',
                                    speed_unit='mph',
                                    grid_dimension=16,
                                    total_demand=40_000,
                                    time_budget=360),

    'Orlando'   : NetworkDescriptor(network_id=1128379,
                                    input_directory='test_data/Orlando',
                                    name='Orlando',
                                    length_unit='mile',
                                    speed_unit='mph',
                                    grid_dimension=36,
                                    total_demand=60_000,
                                    time_budget=360,
                                    column_gen_num=20)
}


def main():
    descriptor = networks['UCF']
    print(descriptor.name)

    pwd = os.getcwd()
    os.chdir(descriptor.input_directory)

    prep_OSM_network(descriptor)
    extract_POI_data(descriptor)

    # synthesize_zones_demand(descriptor)

    # filtered_data_to_csv(descriptor)

    # user_equilibrium(descriptor)

    os.chdir(pwd)


if __name__ == '__main__':
    main()