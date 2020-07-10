#include "QuadGrid.h"

int main(int argc, char **argv) {
	int N = 20;
	if (argc > 1) {
		N = std::stoi(argv[1]);
	}
	long long timeMPI, timeNoMPI;
	QuadGrid withMPI = QuadGrid(N, 0.001);
	withMPI.initMPI(argc, argv);
	timeMPI = withMPI.start(1);
	std::cout << "time taken with MPI: " << timeMPI << std::endl;
	withMPI.writeToCSV("withMPI_" + std::to_string(N), timeMPI);

	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (rank == 0) {
		QuadGrid noMPI = QuadGrid(N, 0.001);
		timeNoMPI = noMPI.start(1);
		noMPI.writeToCSV("noMPI_" + std::to_string(N), timeNoMPI);
		std::cout << "time taken with no MPI: " << timeNoMPI << std::endl;
	}

	if (timeMPI < timeNoMPI) {
		std::cout << "MPI is faster at N: " << N << std::endl;
	}
	else {
		std::cout << "no MPI is faster at N: " << N << std::endl;
	}


	return 0;
}
//#define _CRT_SECURE_NO_DEPRECATE
//#define _CRT_SECURE_NO_WARNINGS
//
//#include <stdio.h>
//#include <stdlib.h>
//#include <math.h>
//#include <iostream>
//#include <mpi.h>
//
//#define MAX_ITER 100
//
//// Maximum value of the matrix element
//#define MAX 100
//#define TOL 0.000001
//
//
//void printMatrix(float *matrix, int rows, int cols) {
//	for (int row = 0; row < rows; row++) {
//		for (int col = 0; col < cols; col++) {
//			printf(" %.2f ", matrix[row * cols + col]);
//		}
//		printf("\n");
//	}
//}
//
//
//
//// Generate a random float number with the maximum value of max
//float rand_float(int max) {
//	return ((float)rand() / (float)(RAND_MAX)) * max;
//}
//
//
//
//
//// Calculates how many rows are given, as maximum, to each node
//int get_max_rows(int num_nodes, int n) {
//	return (int)(ceil((n - 2) / num_nodes) + 2);
//}
//
//
//
//
//// Gets the position from which elements are gonna be sent / received
//int get_node_offset(int node_id, int n, int max_rows) {
//	return node_id * n * (max_rows - 2);
//}
//
//
//
//
//// Calculates how many elements are going to a given node
//int get_node_elems(int node_id, int n, int max_rows) {
//
//	int node_offset = get_node_offset(node_id, n, max_rows);
//	int node_elems = max_rows * n;
//
//	// Case in which the node receive the full set of elements
//	if (node_offset + node_elems <= (n*n)) {
//		return node_elems;
//	}
//
//	// Case of the last node, which could get less elements
//	else {
//		return (n*n) - node_offset;
//	}
//}
//
//
//
//
//// Allocate 2D matrix in the master node
//void allocate_root_matrix(float **mat, int n, int m) {
//
//	*mat = (float *)malloc(n * m * sizeof(float));
//	//for (int i = 0; i < (n*m); i++) {
//	//	(*mat)[i] = rand_float(MAX);
//	//}
//
//	for (int i = 0; i < (n*m); i++) {
//		(*mat)[i] = 0;
//	}
//
//
//	int dim = n;
//	float distance = 1.0 / n;
//	(*mat)[0] = 1.0;
//	(*mat)[n * m - 1] = 1.0;
//	for (int i = 1; i < n; i++) {
//		(*mat)[i] = (*mat)[i - 1] - distance;
//		(*mat)[(n - 1) * dim + (dim - 1) - i] = (*mat)[(dim - 1) * dim + dim - i] - distance;
//		(*mat)[i * dim] = (*mat)[(i - 1) * dim] - distance;
//		int sec = (dim - i) * dim + dim - 1;
//		int first = (dim - 1 - i) * dim + dim - 1;
//		(*mat)[first] =(*mat)[sec] - distance;
//	}
//}
//
//
//
//
//// Allocate 2D matrix in the slaves nodes
//void allocate_node_matrix(float **mat, int num_elems) {
//	*mat = (float *)malloc(num_elems * sizeof(float));
//}
//
//
//
//
//// Solves as many elements as specified in "num_elems"
//void solver(float **mat, int n, int num_elems) {
////
////	float diff = 0, temp;
////	int done = 0, cnt_iter = 0, myrank;
////
////	while (!done && (cnt_iter < MAX_ITER)) {
////		diff = 0;
////
////		// Neither the first row nor the last row are solved
////		// (that's why it starts at "n" and it goes up to "num_elems - 2n")
////		for (int i = n; i < num_elems - (2 * n); i++) {
////
////			// Additionally, neither the first nor last column are solved
////			// (that's why the first and last positions of "rows" are skipped)
////			if ((i % n == 0) || (i + 1 % n == 0)) {
////				continue;
////			}
////
////			int pos_up = i - n;
////			int pos_do = i + n;
////			int pos_le = i - 1;
////			int pos_ri = i + 1;
////
////			temp = (*mat)[i];
////			(*mat)[i] = 0.2 * ((*mat)[i] + (*mat)[pos_le] + (*mat)[pos_up] + (*mat)[pos_ri] + (*mat)[pos_do]);
////			diff += abs((int)((*mat)[i] - temp));
////		}
////
////		if (diff / n / n < TOL) {
////			done = 1;
////		}
////		cnt_iter++;
////	}
////
////
////	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
////	if (done) {
////		printf("Node %d: Solver converged after %d iterations\n", myrank, cnt_iter);
////	}
////	else {
////		printf("Node %d: Solver not converged after %d iterations\n", myrank, cnt_iter);
////	}
//	while (true) {
//		for (int i = n; i < num_elems - (2 * n); i++) {
//			residuum(i, n, mat, )
//		}
//	}
//
//}
//
//float faultFunction(int i, int n, float **mat, float distance) {
//	int pos_up = i - n;
//	int pos_do = i + n;
//	int pos_le = i - 1;
//	int pos_ri = i + 1;
//	return (1 / sqrt(distance)) * (4 * (*mat)[i] - (*mat)[pos_up] - (*mat)[pos_do] - (*mat)[pos_le] - (*mat)[pos_ri]);
//}
//
//double residuum(int i, int n, float **mat, float distance) {
//	int pos_up = i - n;
//	int pos_do = i + n;
//	int pos_le = i - 1;
//	int pos_ri = i + 1;
//	return faultFunction(i, n, mat, distance) * (distance * distance) - (4 * (*mat)[i] - (*mat)[pos_up] - (*mat)[pos_do] - (*mat)[pos_le] - (*mat)[pos_ri]);
//}
//
//
//int main(int argc, char *argv[]) {
//
//	int np, myrank, n, communication, i;
//	float *a = nullptr, *b = nullptr;
//
//	MPI_Init(&argc, &argv);
//	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
//	MPI_Comm_size(MPI_COMM_WORLD, &np);
//
//	if (argc < 3) {
//		if (myrank == 0) {
//			printf("Call this program with two parameters: matrix_size communication \n");
//			printf("\t matrix_size: Add 2 to a power of 2 (e.g. : 18, 1026)\n");
//			printf("\t communication:\n");
//			printf("\t\t 0: initial and final using point-to-point communication\n");
//			printf("\t\t 1: initial and final using collective communication\n");
//		}
//
//		MPI_Finalize();
//		exit(1);
//	}
//
//
//	n = atoi(argv[1]);
//	communication = atoi(argv[2]);
//
//	if (myrank == 0) {
//		printf("Matrix size = %d communication = %d\n", n, communication);
//	}
//
//
//	// Calculate common relevant values for each node
//	int max_rows = get_max_rows(np, n);
//
//	// Array of containing the offset and number of elements per node
//	int *nodes_offsets = new int[np];
//	int *nodes_elems = new int[np];
//	for (i = 0; i < np; i++) {
//		nodes_offsets[i] = get_node_offset(i, n, max_rows);
//		nodes_elems[i] = get_node_elems(i, n, max_rows);
//	}
//
//	// Variable to store the node local number of elements
//	int num_elems = nodes_elems[myrank];
//
//
//	double tscom1 = MPI_Wtime();
//
//
//	switch (communication) {
//	case 0: {
//		if (myrank == 0) {
//
//			// Allocating memory for the whole matrix
//			allocate_root_matrix(&a, n, n);
//
//			std::cout << "before calculating matrix:" << std::endl << std::endl;
//			printMatrix(a, n, n);
//			std::cout << std::endl << std::endl;
//
//
//
//			// Master sends chuncks to every other node
//			for (i = 1; i < np; i++) {
//				int i_offset = nodes_offsets[i];
//				int i_elems = nodes_elems[i];
//				MPI_Send(&a[i_offset], i_elems, MPI_FLOAT, i, 0, MPI_COMM_WORLD);
//			}
//		}
//		else {
//
//			// Allocating the exact memory to the rows receiving
//			allocate_node_matrix(&a, num_elems);
//			MPI_Status status;
//
//			// Receiving the data from the master node
//			MPI_Recv(a, num_elems, MPI_FLOAT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
//		}
//		break;
//	}
//
//	case 1: {
//		if (myrank == 0) {
//			// Allocating memory for the whole matrix
//			allocate_root_matrix(&a, n, n);
//		}
//		// Allocating the exact memory where the receiving rows are computed
//		allocate_node_matrix(&b, num_elems);
//
//		// Collective communication for scattering the matrix
//		// Info: https://www.mpich.org/static/docs/v3.1/www3/MPI_Scatterv.html
//		MPI_Scatterv(a, nodes_elems, nodes_offsets, MPI_FLOAT, b, num_elems, MPI_FLOAT, 0, MPI_COMM_WORLD);
//		break;
//	}
//	}
//
//
//	double tfcom1 = MPI_Wtime();
//	double tsop = MPI_Wtime();
//
//	// --------- SOLVER ---------
//	if (communication == 0) {
//		solver(&a, n, num_elems);
//	}
//	else {
//		solver(&b, n, num_elems);
//	}
//
//	double tfop = MPI_Wtime();
//	double tscom2 = MPI_Wtime();
//
//
//	switch (communication) {
//	case 0: {
//		if (myrank == 0) {
//
//			MPI_Status status;
//
//			// Master sends chuncks to every other node
//			for (i = 1; i < np; i++) {
//				int i_offset = nodes_offsets[i] + n; 		// +n to skip cortex values
//				int i_elems = nodes_elems[i] - (2 * n);		// -2n to skip cortex values
//				MPI_Recv(&a[i_offset], i_elems, MPI_FLOAT, i, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
//			}
//		}
//		else {
//			int solved_offset = n;							// Start at n to skip cortex values
//			int solved_elems = num_elems - (2 * n);			// Reach num_elems-2n to skip cortex values
//
//			// Compute num_elems sin la corteza
//			MPI_Send(&a[solved_offset], solved_elems, MPI_FLOAT, 0, 0, MPI_COMM_WORLD);
//		}
//
//		break;
//	}
//
//	case 1: {
//
//		// Collective communication for gathering the matrix
//		// Info: http://www.mpich.org/static/docs/v3.2.1/www/www3/MPI_Gatherv.html
//		MPI_Gatherv(b, num_elems, MPI_FLOAT, a, nodes_elems, nodes_offsets, MPI_FLOAT, 0, MPI_COMM_WORLD);
//		break;
//	}
//	}
//
//	double tfcom2 = MPI_Wtime();
//
//
//	if (myrank == 0) {
//		float com_time = (tfcom1 - tscom1) + (tfcom2 - tscom2);
//		float ops_time = tfop - tsop;
//		float total_time = com_time + ops_time;
//
//		printf("Communication time: %f\n", com_time);
//		printf("Operations time: %f\n", ops_time);
//		printf("Total time: %f\n", total_time);
//
//		FILE *f;
//		f = fopen("results.csv", "a");
//		fprintf(f, "Communication;Nodes;Processes;Size;Communication-time;Operations-time;Total-time;\n");
//
//
//		fprintf(f, "%d;%d;;%d;%f;%f;%f;\n", communication, np, n, com_time, ops_time, total_time);
//		fclose(f);
//	}
//
//	if (myrank == 0) {
//		std::cout << "finished calculating matrix:" << std::endl << std::endl;
//		printMatrix(a, n, n);
//		std::cout << std::endl;
//	}
//
//
//	if (a != nullptr)
//		free(a);
//	if (b != nullptr)
//		free(b);
//
//	MPI_Finalize();
//	return 0;
//}
