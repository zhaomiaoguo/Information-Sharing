'''
Class definitions
'''

__all__ = ['NetworkDescriptor']


class NetworkDescriptor():
    ''' An individual network's info and parameters to carry out traffic assignment '''
    
    def __init__(self, 
                 network_id: int, 
                 input_directory: str,
                 name: str = 'gmns', 
                 length_unit: str = 'mile', 
                 speed_unit: str = 'mph', 
                 grid_dimension: int = 8, 
                 total_demand: int = 10_000, 
                 time_budget: int = 120,
                 column_gen_num: int = 10,
                 column_update_num: int = 10):
        ''' Attributes of the Network '''

        self.network_id = network_id
        self.name = name
        self.input_directory = input_directory
        self.length_unit = length_unit
        self.speed_unit = speed_unit
        self.grid_dimension = grid_dimension
        self.total_demand = total_demand
        self.time_budget = time_budget
        self.column_gen_num = column_gen_num
        self.column_update_num = column_update_num


    def set_osm_target(self, target: str):
        ''' Specify the file to read OSM data from '''
        self.osm_target = target