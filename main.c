#define LAB4_EXTEND

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include "Lab4_IO.h"
#include "timer.h"

#define EPSILON 0.00001
#define DAMPING_FACTOR 0.85

int main (int argc, char* argv[]){
    // Variables to hold data_input information
    struct node *nodehead;
    int nodecount;
    int *num_in_links, *num_out_links;

    // Pointers to data that hold iterative data of PageRank Algorithm
    double *r, *r_pre;

    // Looping variables, chunk_size, and starting and ending index of split
    int i, j, chunk_size, start_index, end_index;

    // Set inital run flag
    int array_flag = 0;

    // Damping Constant
    double damp_const;

    // Variables to hold start time and end time
    double start, end;

    // Variables for size of processes, and rank of current process
    int comm_sz, my_rank;

    // Initialize MPI
    MPI_Init(&argc, &argv);

    // Find Comm size with pass-by-reference
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);

    // Find Comm rank with pass-by-reference
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    // Read the data_input file and get the number of nodes, list of in links,
    // and list of out links using pass-by-reference
    if (get_node_stat(&nodecount, &num_in_links, &num_out_links)) return 254;

    // Create the structures to hold the information from data_input file
    // using pass-by-reference
    if (node_init(&nodehead, num_in_links, num_out_links, 0, nodecount)) return 254;
    
    // Calculate chunk_size, start and end index
    chunk_size = nodecount/comm_sz;
    start_index = chunk_size * my_rank;
    end_index = start_index + chunk_size;   

	// Receive buffer	
	double recv_buffer[nodecount];	

    // Allocate memory for the iterative steps of PageRank
    r = malloc(nodecount * sizeof(double));
    r_pre = malloc(nodecount * sizeof(double));

    // Create first iteration of PageRank values
    for ( i = 0; i < nodecount; ++i)
        r[i] = 1.0 / nodecount;

    // Calculate damp_const from Damping factor and node count
    damp_const = (1.0 - DAMPING_FACTOR) / nodecount;

    // Get start time
    GET_TIME(start);
	
	// CORE CALCULATION
    // Do PageRank Algorithm calculations
    // Cycle through only until error of r and r_pre are less than EPSILON
	do{

		// Check if this is the first run, if so set array_flag to 1 and copy r into r_pre
		// In other iterations copy recv_buffer into r_pre
		if (array_flag == 1) {
			vec_cp(recv_buffer, r_pre, nodecount);
		}
		else {
			array_flag = 1;
			vec_cp(r, r_pre, nodecount);
		}

        // Cycle through all of nodecount
        for ( i = 0; i < nodecount; ++i){

            // Clear r[i] first
            r[i] = 0;

            // Cycle through all of the inlinks and sum the probability into r[i]
            // PageRank Algorithm step
            for ( j = 0; j < nodehead[i].num_in_links; ++j)
                r[i] += r_pre[nodehead[i].inlinks[j]] / num_out_links[nodehead[i].inlinks[j]];

            // Adjust with damping factor and constant
            r[i] *= DAMPING_FACTOR;
            r[i] += damp_const;
        }

        // MPI_Allgather to gather all of the data and distribute so everyone has the data
		MPI_Allgather(&r[start_index], chunk_size, MPI_DOUBLE, &recv_buffer, chunk_size, MPI_DOUBLE, MPI_COMM_WORLD);
		
    }while(rel_error(recv_buffer, r_pre, nodecount) >= EPSILON);

	// Checks if you are the master process
	if (my_rank == 0) {

		// Get end time
		GET_TIME(end);    

	    // Save result to the output file
		Lab4_saveoutput(r, nodecount, end-start);
	}

    // Post processing
    // Freeing up pointers
    node_destroy(nodehead, nodecount);
    free(num_in_links); free(num_out_links);
    free(r); free(r_pre);

    // Finalize MPI
	MPI_Finalize();
	
    return 0;
}
