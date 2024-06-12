# MPI-Wireless-Sensor-Network

## Introduction

### Selected Charging and Base Station Simulation Architecture

![Architecture](<link-to-figure-1>)

As shown in Figure 1, the Architecture is as follows: the simulation consists of an M x N cartesian grid for the charging stations. These can communicate with adjacent nodes - nodes to the right, left, up, and down from their current position as shown using blue arrows. Each charging station can send a report to the base station and receive data from the base station as denoted using the red arrows. The red arrow also denotes the termination signal that will be sent to all the nodes in the grid from the base station once a timer selected by the user is complete. This causes all the nodes to finish up and terminate.

![Demonstration](<link-to-figure-2>)

### Architecture Functionality Using Demonstration

The algorithm starts as follows:

**Step 0:** At each Iteration, the nodes will update each of their ports. The number of ports at each node is chosen as 10, and these nodes will update each port number in parallel using OpenMP threads to randomly generate a value between 0 and 1 at every iteration and update an array shared by each thread of the node. This array is then used to calculate the total available ports for the node.

**Step 1:** Each node will check using a condition if they are considered to be full, as a default, the condition is set as if there are no available ports. For all nodes that this condition is true for (Node 0, 1, 2, 3, 4), they will send a message to all adjacent nodes (the node to their top, bottom, left, and right in the grid) requesting their port availability concurrently.

**Step 2:** Each Node that has been requested for its port availability will then concurrently send its port availability to the nodes that requested it.

**Step 3:** All nodes will then check if all their neighbour nodes have at least a port available. If it is such as Node 3, it will print all the adjacent nodes available. Here it would say Nodes 6 and 4 are available. If a Node finds that all adjacent charging station nodes are also not available, it will send a report with the required data to the Base Station.

**Step 4:** Once the base station gets a report, it will save the report in a report list for the current iteration. It will then wait for 3 seconds to allow all the reports for the iteration to arrive. Once the wait time is over, it will process and calculate the ranks of all the nearby charging station nodes for the reporting node. It will then process which of these nearby ranks have also sent a report within the same iteration.

**Step 5:** The base station will then check for each nearby node if they have sent a report to the base station within the same iteration as well. If they have, they are no longer considered as available nearby nodes. If all said nearby nodes for a reporting node are not available, the base station will send a message back to the reporting node stating “no nearby nodes are available”. If there are available nearby nodes for the reporting node that have not sent a report to the base station, the base station will send back to the reporting node all the ranks of the nearby nodes that are available back to the reporting node. The reporting node will then print out the available nearby nodes to the user.

**Step 6:** Once all the reports for the iteration are processed, the base station will log all this information to the log.txt file. This will include all information necessary for a report sent in a specific iteration.

## Roadmap

- [x] Implement MPI communication for grid nodes
- [x] Define data structures for node and base station communication
- [x] Develop algorithms for port availability updates and report generation
- [x] Implement timer functionality for termination signal
- [x] Test and debug MPI communication and algorithmic logic
- [x] Optimize code for performance and scalability
- [x] Document code with comments and explanations

## License

Licensed under MIT
