#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <mpi.h>
#include <time.h>
#include <unistd.h>
#include <omp.h>
#include <math.h>

// Define the maximum grid dimensions.
#define MAX_GRID_X 100
#define MAX_GRID_Y 100
// Define the maximum number of ports per ChargingNode.
#define MAX_PORTS 4
#define REPORT_INTERVAL 1
#define SHIFT_ROW 0
#define SHIFT_COL 1
#define DISP 1
// Define the value for the number of ports to consider a node full.
#define FULL_PORT 3
// Define the time limit for termination.
#define TERMINATE_TIME 10

/**
 * Structure representing a charging port.
 * @param port_id The port's unique ID.
 * @param availability The availability status of the port.
 */

typedef struct {
    int port_id;
    int availability;
} ChargingPort;

/**
 * Structure representing a ChargingNode.
 * @param x The X-coordinate of the ChargingNode in the grid.
 * @param y The Y-coordinate of the ChargingNode in the grid.
 * @param ports An array of ChargingPort structures representing the charging ports.
 */
typedef struct {
    int x;
    int y;
    ChargingPort ports[MAX_PORTS];
} ChargingNode;

/**
 * Structure for storing data for reporting.
 * @param YYYY Year.
 * @param MM Month.
 * @param DD Day.
 * @param HH Hour.
 * @param M Minute.
 * @param SS Second.
 * @param Availability Availability status.
 */
struct Data {
    int YYYY, MM, DD, HH, M, SS, Availability;
};

/**
 * Structure for reporting information.
 * @param adjacent An array of adjacent nodes.
 * @param time_taken Time taken for communication.
 * @param sent Number of reports sent.
 * @param received Number of reports received.
 * @param iteration Current iteration.
 * @param num_adj Number of adjacent nodes.
 * @param rank Rank of the ChargingNode.
 * @param ports Number of ports.
 * @param adj_val An array of adjacency values.
 * @param alert_time Alert time information.
 */
typedef struct {
    int adjacent[4];
    double time_taken;
    int sent;
    int received;
    int iteration;
    int num_adj;
    int rank;
    int ports;
    int adj_val[4];
    char alert_time[];
} Report;

/**
 * Simulates the behavior of a ChargingNode.
 * @param node A pointer to the ChargingNode to simulate.
 * @param rank The rank of the ChargingNode.
 * @param grid_size_x The X dimension of the grid.
 * @param grid_size_y The Y dimension of the grid.
 * @param master_comm The master MPI communicator.
 * @param grid_comm The grid-specific MPI communicator.
 * @param size The total number of processes.
 */
void simulateChargingNode(ChargingNode *node, int rank, int grid_size_x, int grid_size_y, MPI_Comm master_comm, MPI_Comm grid_comm, int size);

/**
 * Simulates the behavior of the base station.
 * @param rank The rank of the base station.
 * @param num_nodes The number of ChargingNodes.
 */
void simulateBaseStation(int rank, int num_nodes);

/**
 * Handles the logic for the base station.
 * @param master_comm The master MPI communicator.
 * @param comm The MPI communicator for the base station.
 * @param grid_size_x The X dimension of the grid.
 * @param grid_size_y The Y dimension of the grid.
 */
void base_station(MPI_Comm master_comm, MPI_Comm comm, int grid_size_x, int grid_size_y);

/**
 * Handles the logic for the charging station (ChargingNode).
 * @param master_comm The master MPI communicator.
 * @param comm The MPI communicator for the ChargingNode.
 * @param grid_size_x The X dimension of the grid.
 * @param grid_size_y The Y dimension of the grid.
 */
void charging_station(MPI_Comm master_comm, MPI_Comm comm, int grid_size_x, int grid_size_y);

// Global termination signal flag.
int termination_signal = 0;

/**
 * The main function for the simulation.
 * @param argc The number of command-line arguments.
 * @param argv An array of command-line arguments.
 * @return 0 on successful execution.
 */
int main(int argc, char *argv[]) {

    int grid_size_x, grid_size_y;
    int size, rank;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    
    if (argc != 3) {
    if (rank == 0) {
        printf("Usage: %s <grid_size_x> <grid_size_y>\n", argv[0]);
    }
    MPI_Finalize();
    return 1;
    }

    grid_size_x = atoi(argv[1]);
    grid_size_y = atoi(argv[2]);

    if (size != (grid_size_x * grid_size_y + 1)) {
        if (rank == 0) {
            printf("Number of processes should be %d\n", (grid_size_x * grid_size_y + 1));
        }
        MPI_Finalize();
        return 1;
    }
    
    
    MPI_Comm new_comm;
    MPI_Comm_split(MPI_COMM_WORLD, rank == size - 1, 0, &new_comm);
    if (rank == size - 1){
        base_station(MPI_COMM_WORLD, new_comm, grid_size_x, grid_size_y);
        printf("Base Station closed\n");
    } else {
        charging_station(MPI_COMM_WORLD, new_comm, grid_size_x, grid_size_y);
        
    }

    MPI_Finalize();
    return 0;
}

void charging_station(MPI_Comm master_comm, MPI_Comm comm, int grid_size_x, int grid_size_y){

    int ndims = 2, size, rank, reorder, my_cart_rank, ierr;
    int wrap_around[ndims];

    MPI_Comm grid_comm;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    
    int dims[2] = {grid_size_x, grid_size_y};
    int periods[2] = {0, 0};
    int coords[2];
    MPI_Dims_create(size, ndims, dims);
    if (rank == 0){
        printf("Root Rank: %d, Comm Size: %d: Grid Dimensions = [%d x %d] \n", rank, size, dims[0], dims[1]);
    }
    
    /* creating the cartesion mapping*/
    wrap_around[0] = wrap_around[1] = 0;
    reorder = 1;
    ierr = 0;
    
    /* the MPI_Cart_create to create the cartesian communicator*/
    ierr = MPI_Cart_create(comm, ndims, dims, wrap_around, reorder, &grid_comm);
    
    /* the MPI_Cart_coords to pass in the communicator and the current rank, to get the coordinates for the node*/
    MPI_Cart_coords(grid_comm, rank, ndims, coords);
    
    /*use the current coordinates of the node to get its rank*/
    MPI_Cart_rank(grid_comm, coords, &my_cart_rank);
   
    
    //if (rank != size - 1) {
    ChargingNode node;
    node.x = coords[0];
    node.y = coords[1];
    simulateChargingNode(&node, rank, grid_size_x, grid_size_y, master_comm, grid_comm, size);

}

// code for the base station that will be contacted if a quadrant is fully utilized
void base_station(MPI_Comm master_comm, MPI_Comm comm, int grid_size_x, int grid_size_y) {
    int reports_received = 0;
    MPI_Status status;
    MPI_Request request;
    FILE *log_file;

    int size, rank;

    MPI_Comm_rank(master_comm, &rank);
    MPI_Comm_size(master_comm, &size);

    // Open a log file for writing
    FILE *logFile = fopen("log.txt", "w+");

    if (logFile == NULL) {
        printf("Error opening the log file.\n");
        return;
    }

    fprintf(logFile, "Base station started.\n");

    double report_message[14];
    Report reportArray[100];
    int array_size = 0;



    #pragma omp parallel num_threads(2) shared(report_message, reportArray, array_size)
    {
        int thread_id = omp_get_thread_num();
        
        if (thread_id == 1) {
            // Receiving thread
            while (!termination_signal) {
                if (logFile == NULL) {
                    printf("Error opening the log file.\n");
                }

                int numAdjacentNodes = 0;
                int availabilityFull = 0;
                int iteration = 0;
                int messages_sent = 0;
                int messages_received = 0;
                double sent_time, recv_time;
                int report_rank;
                Report report;
                int nearby_ranks[8];  
                int num_adjacent = 0;
                
                int sender_rank;

                // Set the timeout (in seconds)
                double timeout = 4.0;
                double start_time = MPI_Wtime();

                // Start non-blocking receive
                MPI_Irecv(report_message, 14, MPI_DOUBLE, MPI_ANY_SOURCE, 0, master_comm, &request);
                time_t now;
                struct tm *timeinfo;

                time(&now);
                timeinfo = localtime(&now);

                // Format the current time in the desired format
                strftime(report.alert_time, 80, "%a %Y-%m-%d %H:%M:%S", timeinfo);
                
                
                

                int message_received = 0;
                while (1) {
                    // Check if the message has been received
                    int flag;
                    MPI_Test(&request, &flag, &status);

                    if (flag) {
                        message_received = 1;
                        break;  // Message received, exit the loop
                    }

                    // Check for the timeout
                    double current_time = MPI_Wtime();
                    if (current_time - start_time > timeout) {
                        break;  // Timeout exceeded, exit the loop
                    }
                }
                
                if (!message_received) {
                    // Process the received message
                    
                    for (int i = 0; i < array_size; i++){
                        report = reportArray[i];
                        report_rank = report.rank;
                        
                        int current_rank = report_rank;  
                        int max_distance = 2;

                        // Calculate the row and column values for the current rank
                        int current_row = current_rank / grid_size_x;
                        int current_col = current_rank % grid_size_x;
                        
                        num_adjacent = 0;




                        // Calculate and validate ranks of adjacent nodes at 2-unit distance
                        for (int row = current_row - max_distance; row <= current_row + max_distance; row++) {
                            for (int col = current_col - max_distance; col <= current_col + max_distance; col++) {
                                if (row >= 0 && row < grid_size_y && col >= 0 && col < grid_size_x) {
                                    int rank = row * grid_size_x + col;
                                    int row_diff = abs(row - current_row);
                                    int col_diff = abs(col - current_col);

                                    if (rank != current_rank && (row_diff + col_diff) == max_distance) {
                                        nearby_ranks[num_adjacent++] = rank;
                                    }
                                }
                            }
                        }
                        
                        // Create arrays

                        // Fill in your reportArray and nearby_ranks with data

                        // Initialize a new array to store unmatched nearby_ranks
                        int unmatchedArray[8];
                        int unmatchedCount = 0;

                        // Loop through nearby_ranks
                        for (int i = 0; i < num_adjacent; i++) {
                            int nearbyRank = nearby_ranks[i];
                            bool found = false;

                            // Loop through reportArray to check if nearbyRank is found
                            for (int j = 0; j < array_size; j++) {
                                if (reportArray[j].rank == nearbyRank) {
                                    found = true;
                                    break;
                                }
                            }

                            // If nearbyRank is not found in reportArray, add it to unmatchedArray
                            if (!found) {
                                unmatchedArray[unmatchedCount] = nearbyRank;
                                unmatchedCount++;
                            }
                        }
                        
                        // Print the header
                        fprintf(logFile, "------------------------------------------------------------------------------------------------------------\n");
                        fprintf(logFile, "Iteration : %d\n", report.iteration);
                        
                        time_t now;
                        struct tm *timeinfo;
                        char buffer[80];

                        time(&now);
                        timeinfo = localtime(&now);

                        // Format the current time in the desired format
                        strftime(buffer, 80, "%a %Y-%m-%d %H:%M:%S", timeinfo);
                        
                        fprintf(logFile, "Logged time : %s\n", buffer);
                        fprintf(logFile, "Alert reported time : %s\n", report.alert_time);
                        fprintf(logFile, "Number of adjacent nodes: %d\n", report.num_adj);

                        
                        fprintf(logFile, "\n");
                        
                        fprintf(logFile, "Availability to be considered full: %d\n\n", FULL_PORT);
                        
                        // Print reporting node data
                        fprintf(logFile, "Reporting Node \tCoord\tPort Value\tAvailable Port\n");
                        // Add your reporting node data here
                        
                        int x_coord = report.rank % grid_size_x;
                        int y_coord = report.rank / grid_size_x;
                        fprintf(logFile, "\t\t%d\t\t(%d, %d)\t\t%d\t\t\t\t%d\n", report.rank, x_coord, y_coord, MAX_PORTS, report.ports);
                        
                        // Print adjacent nodes data
                        fprintf(logFile, "\nAdjacent Nodes\tCoord\tPort Value\tAvailable Port\n");
                        for (int i = 0; i < 4; i++) {
                            if (report.adjacent[i] >= 0) {
                                int x_coord = report.adjacent[i] % grid_size_x;
                                int y_coord = report.adjacent[i] / grid_size_x;
                                fprintf(logFile, "\t\t%d\t\t(%d,%d)\t\t%d\t\t\t\t%d\n", report.adjacent[i], x_coord, y_coord, MAX_PORTS, report.adj_val[i]);
                            }
                        }
                        // Add your adjacent nodes data here
                        
                        // Print nearby nodes data
                        fprintf(logFile, "\nNearby Nodes \tCoord\n");
                        
                        // Add your nearby nodes data here
                        for (int i = 0; i < num_adjacent; i++) {
                            int x_coord = nearby_ranks[i] % grid_size_x;
                            int y_coord = nearby_ranks[i] / grid_size_x;
                            fprintf(logFile, "\t\t%d\t\t(%d,%d)\n", nearby_ranks[i], x_coord, y_coord);
                        }
                        
                        // Available station nearby data
                        fprintf(logFile, "\nAvailable station nearby (no report received in the iteration):");
                        
                        if (unmatchedCount == 0){
                            fprintf(logFile, " N/A\n");
                        }
                        
                        for (int i = 0; i < unmatchedCount; i++) {
                            int unmatchedValue = unmatchedArray[i];
                            // Process or print the unmatchedValue
                            if (i == unmatchedCount - 1){
                                fprintf(logFile, " %d\n", unmatchedValue);
                                break;
                            }
                            fprintf(logFile, " %d,", unmatchedValue);
                        }

                        // Add your available station nearby data here
                        fprintf(logFile, "Communication Time (seconds) : %f\n", report.time_taken/100);
                        fprintf(logFile, "Total Messages sent between reporting node and base station: 2\n\n");
                        // Print the header
                        fprintf(logFile, "------------------------------------------------------------------------------------------------------------\n");


                    }                    
                    
                    array_size = 0;

                    continue;
                    // Process report_message as needed
                }
                
                sent_time = report_message[7];
                recv_time = MPI_Wtime();
                double time_taken = recv_time - sent_time;
                
                report.adjacent[0] = (int)report_message[0];
                report.adjacent[1] = (int)report_message[1];
                report.adjacent[2] = (int)report_message[2];
                report.adjacent[3] = (int)report_message[3];
                
                report.sent = (int)report_message[4];
                report.received = (int)report_message[5];
                report.iteration = (int)report_message[6];
                report.rank = (int)report_message[8];
                
                report.ports = (int)report_message[9];
                
                report.adj_val[0] = (int)report_message[10];
                report.adj_val[1] = (int)report_message[11];
                report.adj_val[2] = (int)report_message[12];
                report.adj_val[3] = (int)report_message[13];
                
                report.time_taken = fabs(time_taken);
 
                for (int i = 0; i < 4; i++){
                    if (report.adjacent[i] >= 0){
                        numAdjacentNodes ++;
                    }
                }
                
                report.num_adj = numAdjacentNodes;
                
                if (array_size < 100) { // Check if the array is not full
                    reportArray[array_size] = report;
                    array_size++; // Increment the array size
                } else {
                    printf("Report array is full. Cannot append more reports.\n");
                }
                
                
                // Handle the data
            }
        }
        else if (thread_id == 0) {
            // Sending thread        
            sleep(TERMINATE_TIME);
            termination_signal = 1;
            
            MPI_Request send_request[size - 1];

            for (int i = 0; i < size - 1; i++){
                 MPI_Isend(&termination_signal, 1, MPI_INT, i, 0, master_comm, &send_request[i]);
            }
            
            // Wait for all non-blocking sends to complete
             MPI_Waitall(size - 1, send_request, MPI_STATUSES_IGNORE);
            //MPI_Bcast(&termination_signal, 1, MPI_INT, size - 1, master_comm);
            
        }
    }

    fclose(logFile);
}



void simulateChargingNode(ChargingNode *node, int rank, int grid_size_x, int grid_size_y, MPI_Comm master_comm, MPI_Comm grid_comm, int size) {
    int neighbor_data[MAX_PORTS];
    int neighbor_avail[MAX_PORTS];
    int neighbor_rank;
    int total_messages_sent = 1;
    int total_messages_received = 1;
    MPI_Status status;
    int iteration = 1;
    

    // Initialize charging node data
    srand(time(NULL) + rank);

    // Main simulation loop
    int reports_sent = 0;
   

    // Non-blocking receive
    MPI_Request recv_request;
    // Non-blocking receive to check for termination signal
    MPI_Irecv(&termination_signal, 1, MPI_INT, size, 0, master_comm, &recv_request);

    int recv_complete = 0;
    
    while (!termination_signal) {
        int total_available_ports = 0;
        
        

        // Check if the receive operation is complete
        MPI_Test(&recv_request, &recv_complete, MPI_STATUS_IGNORE);
        
        if (termination_signal){
            return;
        }

        
        struct Data dataBuffer[15] = {0};
        
        int currentRow = 0; // Index for the current row
        int numRows = 0;    // Number of rows currently in the buffer
        
        time_t now = time(NULL);
        struct tm* timeInfo = localtime(&now);
        
        // Simulate periodic port updates
        #pragma omp parallel num_threads(MAX_PORTS) reduction(+:total_available_ports) shared(dataBuffer, timeInfo)
        { 
            #pragma omp for
            for (int i = 0; i < MAX_PORTS; i++) {
                // Simulate port availability updates
                
                int availability_change = rand() % 2; // Randomly change to 0 or 1
                node->ports[i].availability = availability_change;
                total_available_ports += availability_change; // Update total available ports
            }
            
            dataBuffer[currentRow].YYYY = timeInfo->tm_year + 1900; // Year
            dataBuffer[currentRow].MM = timeInfo->tm_mon + 1;        // Month
            dataBuffer[currentRow].DD = timeInfo->tm_mday;           // Day
            dataBuffer[currentRow].HH = timeInfo->tm_hour;           // Hour
            dataBuffer[currentRow].M = timeInfo->tm_min;             // Minute
            dataBuffer[currentRow].SS = timeInfo->tm_sec;            // Second
            dataBuffer[currentRow].Availability = total_available_ports;      

        }
        
        int coords[2];
        MPI_Cart_coords(grid_comm, rank, 2, coords);
        int adjacent_rank[4] = {-1, -1, -1, -1};

        MPI_Cart_shift(grid_comm, SHIFT_ROW, DISP, &adjacent_rank[0], &adjacent_rank[1]);
        MPI_Cart_shift(grid_comm, SHIFT_COL, DISP, &adjacent_rank[2], &adjacent_rank[3]);

        MPI_Request reqs[8]; // 8 requests for 4 sends and 4 receives
        
        // if total available ports is 0 send a message to 4 adjacent ranks
        if (total_available_ports <= FULL_PORT){
            for (int i = 0; i < 4; i++){
                MPI_Isend(&total_available_ports, 1, MPI_INT, adjacent_rank[i], 0, grid_comm, &reqs[0]);
            }
            total_messages_sent += 4;
        }
        
        
        int adjacent_value[4] = {-1, -1, -1, -1};
        MPI_Request recvReqs[4];
        
        // revc messages from the 4 adjacent ranks whichever are zer0
        for (int i = 0; i < 4; i++){
            MPI_Irecv(&adjacent_value[i], 1, MPI_DOUBLE, adjacent_rank[i], 0, grid_comm, &recvReqs[i]);
        }

        // Wait for the receive operations to complete or timeout 
        int completed[4] = {0, 0, 0, 0};
        double start_time = MPI_Wtime();
        double timeout = 0.2;  // Set your desired timeout value in seconds

        while (completed[0] + completed[1] + completed[2] + completed[3] < 4) {
            for (int i = 0; i < 4; i++) {
                if (!completed[i]) {
                    MPI_Test(&recvReqs[i], &completed[i], MPI_STATUS_IGNORE);
                }
            }

            double elapsed_time = MPI_Wtime() - start_time;
            if (elapsed_time >= timeout) {
                // Handle the timeout (e.g., break from the loop or take appropriate action)
                break;
            }
        }
        
        for (int i = 0; i < 4; i++) {
            total_messages_received += completed[i];
        }
        
        // send my current ports to a rank that requesting for my nodes
        for (int i = 0; i < 4; i++){
            if (adjacent_value[i] <= FULL_PORT) {
                
                MPI_Isend(&total_available_ports, 1, MPI_INT, adjacent_rank[i], 1, grid_comm, &reqs[0]);
                total_messages_sent += 1;
            } 
        }
        
        // recv back the port numbers from the ranks node sent requests to
        for (int i = 0; i < 4; i++){
            MPI_Irecv(&adjacent_value[i], 1, MPI_INT, adjacent_rank[i], 1, grid_comm, &recvReqs[i]);
        }


        int completed2[4] = {0, 0, 0, 0};
        start_time = MPI_Wtime();
        
        while (completed2[0] + completed2[1] + completed2[2] + completed2[3] < 4) {
            for (int i = 0; i < 4; i++) {
                if (!completed2[i]) {
                    MPI_Test(&recvReqs[i], &completed2[i], MPI_STATUS_IGNORE);
                }
            }

            double elapsed_time = MPI_Wtime() - start_time;
            if (elapsed_time >= timeout) {
                // Handle the timeout (e.g., break from the loop or take appropriate action)
                break;
            }
        }
        
        for (int i = 0; i < 4; i++) {
            total_messages_received += completed2[i];
        }
   
        if (total_available_ports <= FULL_PORT){
          
      
            if (adjacent_value[2] <= FULL_PORT && adjacent_value[3] <= FULL_PORT && adjacent_value[0] <= FULL_PORT && adjacent_value[1] <= FULL_PORT) {
                
                double available_ranks[14] = {adjacent_rank[0],adjacent_rank[1],adjacent_rank[2],adjacent_rank[3], total_messages_sent, total_messages_received, iteration, MPI_Wtime(), rank, total_available_ports, adjacent_value[0],adjacent_value[1],adjacent_value[2],adjacent_value[3]};
                MPI_Send(&available_ranks, 14, MPI_DOUBLE, size, 0, master_comm);
                total_messages_sent = total_messages_received = 0;
                
            }else {
                if (adjacent_value[2] > 0){
                    printf("For Rank: %d Left Node with rank %d is available. \n", rank, adjacent_rank[2]);
                }
                if (adjacent_value[3] > 0){
                    printf("For Rank: %d Right Node with rank %d is available. \n", rank, adjacent_rank[3]);
                }
                if (adjacent_value[0] > 0){
                    printf("For Rank: %d Upper Node with rank %d is available. \n", rank, adjacent_rank[0]);
                } 
                if (adjacent_value[1] > 0){
                    printf("For Rank: %d Lower Node with rank %d is available. \n", rank, adjacent_rank[1]);
                }
            }
        
        } 
       
        iteration++;
        
        int flag;
        MPI_Test(&recv_request, &flag, MPI_STATUS_IGNORE);
        if (flag) {
            // Terminate the process when the termination signal is received
            break;
        }
        
        sleep(5);
    }
}


