# DTALite IO Specification

## Inputs

### Required 

#### link.csv

| Field | Description | Sample Values | Optional |
| --- | --- | --- | --- |
| name | Name of link (optional for visualization)| Centaurus Drive | Yes |
| link_id | Unique id of link | 0 | Yes |
| osm_way_id | | 11176652 | Yes |
| from_node_id | id of link's upstream node | 0 | No | 
| to_node_id | id of link's downstream node | 1 | No |
| dir_flag | indicaiton of directions of the link (0 $\rightarrow$ bi-direction, 1 $\rightarrow$ single direction) | 1 | No |
| length | Length of link (measured in miles or km) | 58.7 | No |
| lanes | the number of lanes on the link | 1 | No |
| free_speed | free-flow speed (measured on mph or kmph)| 48 | No |
| capacity | Link capacity (measured on vehicles/hr/lane) | 4000 | Yes |
| link_type_name | Name of link type (optional for visualization) | tertiary | Yes |
| link_type | index of link type name | 1 | No |
| geometry | Shape of the link | "LINESTRING (-81.2050301 28.6024228, -81.2054316 28.6023837, -81.2055453 28.6023780, -81.2056278 28.6023743)" | Yes | 
| allowed_uses | | auto | Yes |
| from_biway | | 1 | Yes |
| is_link | | 0 | Yes |
| VDF_fftt1 | free-flow travel time used in the volume-delay function| | Yes |
| VDF_cap1 | capacity used in the volume-delay function | | Yes |
| VDF_alpha1 | coefficient used in the volume-delay function | 0.15 | |
| VDF_beta1 | coefficient used in the volume-delay function | 4 | |
| toll | optional toll cost of the link, which could also be the cost of fuel | 0 | Yes |

#### node.csv

| Field | Description | Sample Values | Optional |
| --- | --- | --- | --- |
| name | Name of node (for visualization only) | | Yes | 
| node_id | Unique id of node | 0 | No |
| osm_node_id | | 892740184 | Yes |
| osm_highway | | traffic_signals | Yes | 
| zone_id | indication of node's physical location | 1 | Yes | 
| ctrl_type | | signal | Yes |
| node_type | optional text label for visualization | 1 | Yes |
| activity_type | | | Yes |
| is_boundary | | | Yes |
| x_coord | longitude or horizontal coordinate in any arbitrary coordinate system | -81.2050301 | No |
| y_coord | latitude or vertical coordinate in any arbitrary coordinate system | 28.6024228 | No |
| intersection_id | | | Yes | 
| poi_id | | | Yes |
| notes | | | Yes |

#### zone.csv

Note: zone data would be required as its own .csv file if the zone data is not included in node.csv already.

| Field | Description | Sample Values | Optional |
| --- | --- | --- | --- |
| zone_id | unique id of zone | 1 | No | 
| bin_index | | 0 | | 
| activity_nodes | | 10;20 | |
| x_coord | | -81.20626505000001 | No |
| y_coord | | 28.602461400000003 | No |
| geometry | | "LINESTRING (-81.21000000000001 28.605,-81.20500000000001 28.605,-81.20500000000001 28.6,-81.21000000000001 28.6,-81.21000000000001 28.605)" | No |
| production | | 833 | |
| attraction | | 833 | |

#### demand.csv

| Field | Description | Sample Values | Optional |
| --- | --- | --- | --- |
| o_zone_id | origin zone number of the link (must already be defined in node.csv or zone.csv) | 1 | No |
| d_zone_id | destination zone number of the link (must already be defined in node.csv or zone.csv) | 2 | No |
| volume | travel demand | 185.09 | No |

#### settings.csv

[Settings.csv description](https://github.com/asu-trans-ai-lab/DTALite/tree/main/user_guide#43-assignment-and-simulation-configuration-file)

### Optional

## Outputs

#### static_link_performance.csv

As described by the DTALite [repo](https://github.com/asu-trans-ai-lab/DTALite/tree/main/user_guide)

| Field | Description | Sample Values |
| --- | --- | --- |
| link_id | link id number of the road | 1 |
| from_node_id | upstream node of the link | 1 |
| to_node_id | downstream node of the link | 2 |
| time_period | simulation period of the agent in HHMM format | 0700_0800 |
| volume | link based flow volume for the defined period | 5600 |
| travel_time | link travel time in minutes | 15 |
| speed | average travel speed on the link | 38 | 
| VOC | volume / capacity ratio | 0.4 |
| notes | explanatory text | period-based |

#### path_flow.csv

As described by the DTALite [repo](https://github.com/asu-trans-ai-lab/DTALite/tree/main/user_guide)

| Field | Description | Sample Values |
| --- | --- | --- |
| agent_id | agent id number | 1 | 
| o_zone_id | origin zone number of the agent | 1 |
| d_zone_id | destination zone number of the agent | 2 |
| path_id | path identification number | 0 |
| o_node_id | origin node number of agent | 1 | 
| d_node_id | destination node number of agent | 2 |
| agent_type | type of the agent | a |
| demand_period | name of the demand period | AM |
| volume | flow volume assigned on the agent | 5600 |
| toll | the amount of money/time that agent pays (in dollars) | 360 |
| travel_time | total time from origin to destination of the agent | 31.5 | 
| distance | total travel distance from origin to destination of the agent (in miles or km) | 20 |
| node_sequence | the nodes the agent passes | 1;3;2; |
| link_sequence | the links the agent passes | 1003;1002; |
| time sequence | time points through which agent passes | 0730:00;0745:45;0801:31; |
| time_decimal_sequence | the decimal times through which the agent passes | 450.00;465.76;481.52; |

#### link_performance.csv

As obtained from testing

| Field | Description | Sample Values |
| --- | --- | --- |
| link_id | link id number of the road | 1 |
| from_node_id | upstream node of the link | 1 |
| to_node_id | downstream node of the link | 2 |
| time_period | simulation period of the agent in HHMM format | 0700_0800 |
| volume | link based flow volume for the defined period | 5600 |
| CA | | |
| CD | | | 
| density | | 0 | 
| queue | | 0 | 
| travel_time | link travel time in minutes | 15 | 
| waiting_time_in_sec | | |
| speed | average travel speed on the link | 38 |
| voc | volume / capacity ratio | 0.4 |
| geometry | | "LINESTRING (-81.2050301 28.6024228, -81.2054316 28.6023837, -81.2055453 28.6023780, -81.2056278 28.6023743)" |
| notes | explanatory text | simulation-based |

#### agent.csv 

As obtained from testing 

| Field | Description | Sample Values |
| --- | --- | --- |
| agent_id | agent unique id | 1 |
| o_zone_id | | 5 |
| d_zone_id | | 1 |
| path_id | | 0 |
| agent_type | | a |
| demand_period | | AM |
| volume | | 6109.71 |
| toll | | 0 |
| travel_time | | 1308.89 |
| distance | | 1308.89 | 
| node sequence | | 3058;1461;1462;889;1463;1897;1898;2393;727;1577;1994;973 |
| link_sequence | | 12235;2404;2406;2408;3329;3331;4532;4534;4536;4538;3555 |
| time_sequence | | |
| time_decimal_sequence | | |
| link_travel_time_sequence | | |
| geometry | | "LINESTRING (-82.1360235 29.1854731, -82.1369651 29.1854898, -82.1378378 29.1854898, -82.1386639 29.1855017, -82.139483 29.1855064, -82.1394987 29.1847693, -82.1395026 29.1840455, -82.1404216 29.1840534, -82.1423263 29.1840502, -82.1447711 29.1840125, -82.147196 29.1839332, -82.147191 29.1833717)" |

#### trajectory.csv 

As obtained from testing 

| Field | Description | Sample Values |
| --- | --- | --- |
| agent_id | | 1 |
| o_zone_id | | 5 |
| d_zone_id | | 1 |
| departure_time_in_min | | 420 |
| arrival_time_in_min | | 507.8 |
| complete_trip | | n |
| travel_time_in_min | total travel time in minutes | 87.80000000000001 |
| PCE | | 1 |
| distance | total distance traveled | 1517.26 |
| number of nodes | total amount of nodes traversed | |
| node_sequence | | 3058;1404;1923;5177;2801;2939;3212;4494;1310;2609;2610;160;3662;2909;2910;2911;2412;2413;1;0 |
| geometry | | "LINESTRING (-82.1360235 29.1854731, -82.1360155 29.1862657, -82.1360131 29.1868095, -82.136012 29.1870216, -82.136015 29.187854, -82.1359997 29.1885662, -82.1369328 29.1885792, -82.1378235 29.1885899, -82.1378164 29.1892993, -82.1386345 29.1893169, -82.1394692 29.1893211, -82.1401545 29.1893369, -82.1403444 29.1893427, -82.1407076 29.191194, -82.1409565 29.1919894, -82.1413772 29.1934315, -82.1415897 29.1943188, -82.1417171 29.194323, -82.1419064 29.1943291, -82.1417649 29.1935266)" |
| time_sequence | | 420.0;507.8 |
| time_sequence_hhmmss | time sequence in hours:minutes:seconds format | 7:00:00;8:27:48 |