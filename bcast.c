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

    // Looping variables, and chunk_size of split
    int i, j, chunk_size;

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
    
    // Calculate chunk_size
    chunk_size = nodecount/comm_sz;   

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
        
        // Make copy of previous r data into r_pre
        vec_cp(r, r_pre, nodecount);

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

        // Cycle through all of the processes and use MPI_Bcast to send and receive result
		for (i = 0; i < comm_sz; i++) {
			MPI_Bcast(&r[i * chunk_size], chunk_size, MPI_DOUBLE, i, MPI_COMM_WORLD);
		}
    }while(rel_error(r, r_pre, nodecount) >= EPSILON);

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
