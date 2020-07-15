#include "RelaxationMethod.h"

int main(int argc, char **argv) {
	int N = 20;
	int method = SCATTER;
	double PRECISION = 0.0001;
	int edgeCase = 0;

	// first command line argument size of Matrix N second the method used for computing
	// for example: mpiexec -n 8 Relaxation 999 1 0.0001 0 (calcualtes a 1000x1000 matrix with scatter method, for case 0)
	if (argc > 1) {
		N = std::stoi(argv[1]);
		if (argc >= 3) {
			method = std::stoi(argv[2]);
		}
		if (argc >= 4) {
			PRECISION = std::atof(argv[3]);
		}
		if (argc >= 5) {
			edgeCase = std::stoi(argv[4]);
		}
	}
	auto programm_start = std::chrono::system_clock::now();
	RelaxationMethod grid = RelaxationMethod(N, PRECISION);
	long long time;
	switch (method) {
	case SEQ: 
		std::cout << "Sequential used" << std::endl;
		time = grid.start(edgeCase);
		//grid.writeToCSV("noMPI_edgeCase" + std::to_string(edgeCase) + "_"  + std::to_string(N), time);
		std::cout << "time taken with no MPI: " << time << std::endl;
		break;

	case SCATTER:
		grid.initMPI(argc, argv);
		time = grid.start(edgeCase, SCATTER);
		if (time != -1) {
			std::cout << "time taken with MPI: " << time << std::endl;
		}
		grid.writeToCSV("withMPI_scatter_edgeCase" + std::to_string(edgeCase) + "_" + std::to_string(N), time);
		break;

	case SEND_RECV:
		grid.initMPI(argc, argv);
		time = grid.start(edgeCase, SEND_RECV);
		if (time != -1) {
			std::cout << "time taken with MPI: " << time << std::endl;
		}
		//grid.writeToCSV("withMPI_send_recv_edgeCase" + std::to_string(edgeCase) + "_" + std::to_string(N), time);
		break;

	default:
		grid.initMPI(argc, argv);
		time = grid.start(edgeCase, SCATTER);
		if (time != -1) {
			std::cout << "time taken with MPI: " << time << std::endl;
		}
		//grid.writeToCSV("withMPI_scatter_edgeCase" + std::to_string(edgeCase) + "_" + std::to_string(N), time);
		break;
	}
	auto programm_end = std::chrono::system_clock::now();
	std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(programm_end - programm_start).count() << " is the complete time of the programm" << std::endl;
	
	return 0;
}
