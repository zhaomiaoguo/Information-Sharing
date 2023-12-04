# Path4GMNS Data IO Table

TODO: check which functions require settings.yml

| function | requires | outputs |
| :------ | :------ | :----- |
| read_network | link.csv, node.csv, poi.csv (optional) | - | 
| output_zones | - | zones.csv |
| output_synthesized_demand | - | demand.csv |
| read_zones | zones.csv | - |
| load_demand | demand.csv | - |
| output_columns | - | agent.csv |
| output_link_performance | - | link_perfomance.csv | 
| evaluate_accessibility | - | zone_accessibility.csv | 
| evaluate_equity | - | equity_str.csv |
| download_sample_setting_file | - | data/ |
| output_agent_trajectory | - | trajectory.csv |
| perform_network_assignment_DTALite(mode=0) | link.csv, node.csv, (check for more) | link_perfomance.csv |
| perform_network_assignment_DTALite(mode=1) | link.csv, node.csv, demand.csv, (check for more) | agent.csv, link_perfomance.csv |
| perform_network_assignment_DTALite(mode=2) | link.csv, node.csv, (check for more) | agent.csv, link_perfomance.csv |
| perform_network_assignment_DTALite(mode=3) | link.csv, node.csv, (check for more) | (check) |