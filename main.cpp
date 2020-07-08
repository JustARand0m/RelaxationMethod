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