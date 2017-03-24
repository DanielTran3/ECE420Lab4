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
    struct node *nodehead;
    int nodecount;
    int *num_in_links, *num_out_links;
    double *r, *r_pre;
    int i, j, chunk_size, start_index, end_index;
	int array_flag = 0;    
	double damp_const;
    double start, end;
	int comm_sz, my_rank;

	MPI_Init(&argc, &argv);

	MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);

	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

	// Read the data_input file and get the number of nodes, list of in links,
	// and list of out links using pass-by-reference
    if (get_node_stat(&nodecount, &num_in_links, &num_out_links)) return 254;

    // Calculate the result
    if (node_init(&nodehead, num_in_links, num_out_links, 0, nodecount)) return 254;
    
	chunk_size = nodecount/comm_sz;
	start_index = chunk_size * my_rank;
	end_index = start_index + chunk_size;

	//double send_buffer[chunk_size];	
	double recv_buffer[nodecount];	

    r = malloc(nodecount * sizeof(double));
    r_pre = malloc(nodecount * sizeof(double));
    for ( i = 0; i < nodecount; ++i)
        r[i] = 1.0 / nodecount;
    damp_const = (1.0 - DAMPING_FACTOR) / nodecount;

    // CORE CALCULATION
	GET_TIME(start);
	int loop = 0;
	do{   
		if (array_flag == 1) {
			vec_cp(recv_buffer, r_pre, nodecount);
		}
		else {
			array_flag = 1;
			vec_cp(r, r_pre, nodecount);
		}
        for ( i = start_index; i < end_index; ++i){
            r[i] = 0;
            for ( j = 0; j < nodehead[i].num_in_links; ++j)
                r[i] += r_pre[nodehead[i].inlinks[j]] / num_out_links[nodehead[i].inlinks[j]];
            r[i] *= DAMPING_FACTOR;
            r[i] += damp_const;
        }
		MPI_Allgather(&r[start_index], chunk_size, MPI_DOUBLE, &recv_buffer, chunk_size, MPI_DOUBLE, MPI_COMM_WORLD);
		//vec_cp(recv_buffer, r, nodecount);
		loop++;
    }while(rel_error(recv_buffer, r_pre, nodecount) >= EPSILON);

	if (my_rank == 0) {
		GET_TIME(end);
		//if (array_flag == 1) {
			//Lab4_saveoutput(recv_buffer, nodecount, end-start);
		//}
		//else {
			//Lab4_saveoutput(r, nodecount, end-start);
		//}	
		printf("%f\n", end-start);
    }
    // post processing
	node_destroy(nodehead, nodecount);
    free(num_in_links); free(num_out_links);
    free(r); free(r_pre);

	MPI_Finalize();
    return 0;
}
