# MPI-Wireless-Sensor-Network

#Introduction
#Selected Charging and Base Station Simulation Architecture

 Figure 1: Architecture for Charging and Base Station Simulation
 As shown in Figure 1, the Architecture is as follows, the simulation consists of an M X N
 cartesian grid for the charging stations, these can communicate with adjacent nodes, nodes
 to the right, left, up and down from their current position as shown using blue arrows.
 Each charging station can send a report to the base station and receive data from the base
 station as denoted using the red arrows, the red arrow also denotes the termination signal
 that will be sent to all the nodes in the grid from the base station once a timer that is
 selected by the user is complete, this causes all the nodes to finish up and terminate.

 
 Figure 2: Demonstration of a simulated run of the architecture.
 Architecture Functionality Using Demonstration
 The algorithm starts as follows, first, the grid for the simulation is set up using the grid size
 provided by the user, m x n, each node in the grid is an MPI thread and the Base station is
 an extra MPI thread. As spoken previously the nodes in the grid can only communicate with
 their adjacent nodes and the base station.
 Step 0- At each Iteration, the nodes will update each of their ports, the number of ports at
 each node is chosen as 10, and these nodes will update each port number in parallel using
 OpenMP threads to randomly generate a value between 0 and 1 at every iteration and
 update an array shared by each thread of the node. This array is then used to calculate the
 total available ports for the node.
 
 Step 1- Each node will check using a condition if they are considered to be full, as a
 default the condition is set as if there are no available ports. For all nodes that this condition
 is true for (Node 0, 1, 2, 3, 4) they will send a message to all adjacent nodes (the node to
 their top, bottom, left and right in the grid) requesting their port availability concurrently.
 
 Step 2- Each Node that has been requested for its port availability will then concurrently
 send its port availability to the nodes that requested for it.
 
 Step 3- all nodes will then check if all their neighbour nodes have at least a port available,
 if it is such as Node 3 it will print all the adjacent nodes available, here it would say Nodes 6
 and 4 are available. If a Node finds that all adjacent charging station nodes are also not
 available it will send a report with required data to the Base Station

  Step 4- Once the base station gets a report it will save the report in a report list for the
 current iteration, it will then wait for 3 seconds to allow all the reports for the iteration to
 arrive, once the wait time is over it will process and calculate the ranks of all the nearby
 charging station nodes for the reporting node, it will then process which of these nearby
 ranks have also sent a report within the same iteration.
 
 Step 5- The base station will then check for each nearby node have they sent a report to
 the base station within the same iteration as well, if they have they are no longer
 considered as available nearby nodes. If all said nearby nodes for a reporting node are not
 available the base station will send a message back to the reporting node stating “no
 nearby nodes are available”, if there are available nearby nodes for the reporting node that
 have not sent a report to the base station, the base station will send back to the reporting
 node all the ranks of the nearby nodes that are available back to the reporting node, the
 reporting node will then print out the available nearby nodes to the user.

 Step 6- Once all the reports for the iteration are processed, the base station will log all this
 information to the log.txt file. this will include all required information for a report sent in a
 specific iteration.
