import os
import sys
import optparse

# we need to import some python modules from the $SUMO_HOME/tools directory
if 'SUMO_HOME' in os.environ:
    tools = os.path.join(os.environ['SUMO_HOME'], 'tools')
    sys.path.append(tools)
else:
    sys.exit("please declare environment variable 'SUMO_HOME'")


from sumolib import checkBinary  # Checks for the binary in environ vars
import traci


def get_options():
    opt_parser = optparse.OptionParser()
    opt_parser.add_option("--nogui", action="store_true",
                         default=False, help="run the commandline version of sumo")
    options, args = opt_parser.parse_args()
    return options


#def update_real_life_info():
    # Simulate an accident at link 1-2-3 (edge e23) at simulation time 100
    #if traci.simulation.getTime() >= 15:
        #for vehicle_id in traci.vehicle.getIDList():
            #if traci.vehicle.getRoadID(vehicle_id) == "route1":
                # Reroute the vehicle to link 1-2-4-3
                #traci.vehicle.changeTarget(vehicle_id, "route3")


def run():
    step = 0
    while traci.simulation.getMinExpectedNumber() > 0:
        traci.simulationStep()
        print("Simulation time:", traci.simulation.getTime())

        if step == 15:
            traci.vehicle.changeTarget("carflow.1", "e24")
            traci.vehicle.changeTarget("carflow.3", "e24")



        # Update real-life information
        #update_real_life_info()

        step += 1

    traci.close()


# main entry point
if __name__ == "__main__":
    options = get_options()

    # check binary
    if options.nogui:
        sumoBinary = checkBinary('sumo')
    else:
        sumoBinary = checkBinary('sumo-gui')

    # traci starts sumo as a subprocess and then this script connects and runs
    traci.start([sumoBinary, "-c", "fnn.sumocfg",
                             "--tripinfo-output", "tripinfo.xml"])
    run()