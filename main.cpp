#include "QuadGrid.h"

int main(int argc, char **argv) {
	int N = 20;
	int method = SCATTER;
	const double PRECISION = 0.0001;

	// first command line argument size of Matrix N second the method used for computing
	// for example: mpiexec -n 8 Relaxation 999 1 (calcualtes a 1000x1000 matrix with scatter method)
	if (argc > 1) {
		N = std::stoi(argv[1]);
		if (argc >= 3) {
			method = std::stoi(argv[2]);
		}
	}
	QuadGrid grid = QuadGrid(N, PRECISION);
	long long time;
	switch (method) {
	case SEQ: 
		std::cout << "Sequential used" << std::endl;
		time = grid.start(1);
		grid.writeToCSV("noMPI_" + std::to_string(N), time);
		std::cout << "time taken with no MPI: " << time << std::endl;
		break;

	case SCATTER:
		grid.initMPI(argc, argv);
		time = grid.start(1, SCATTER);
		std::cout << "time taken with MPI: " << time << std::endl;
		grid.writeToCSV("withMPI_scatter_" + std::to_string(N), time);
		break;

	case SEND_RECV:
		grid.initMPI(argc, argv);
		time = grid.start(1, SEND_RECV);
		std::cout << "time taken with MPI: " << time << std::endl;
		grid.writeToCSV("withMPI_send_recv_" + std::to_string(N), time);
		break;

	default:
		grid.initMPI(argc, argv);
		time = grid.start(1, SCATTER);
		std::cout << "time taken with MPI: " << time << std::endl;
		grid.writeToCSV("withMPI_scatter_" + std::to_string(N), time);
		break;
	}
	
	return 0;
}
