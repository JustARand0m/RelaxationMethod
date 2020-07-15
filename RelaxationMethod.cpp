#include "RelaxationMethod.h"

RelaxationMethod::RelaxationMethod(int n, double _limit): dim(n+1), matrix(dim, dim) {
	world_rank = -1;
	world_size = -1;
	distance = 1.0 / n;
	limit = _limit;
}

long long RelaxationMethod::start(int _edgeCase, int method) {
	long long timeTaken = -1;
	double r = 0;
	Matrix m(5, 5);
	edgeCase = _edgeCase;
	if (world_rank == 0) {
		std::cout << "N was set to " << dim - 1 << std::endl;
		std::cout << "method was set to " << method << std::endl;
		std::cout << "precition was set to " << limit << std::endl;
		std::cout << "edge case was set to " << edgeCase << std::endl;
		std::cout << "processes used " << world_size << std::endl;
	}
	switch (edgeCase) {
	case 0:
		// initialize Grid
		matrix[0] = 1.0;
		matrix[dim * dim - 1] = 1.0;
		for (int i = 1; i < dim; i++) {
			matrix[i] = matrix[i-1] - distance;
			matrix[(dim - 1) * dim + (dim-1) - i] = matrix[(dim - 1) * dim + dim - i] - distance;
			matrix[i * dim] = matrix[(i - 1) * dim] - distance;
			int sec = (dim - i) * dim + dim - 1;
			int first = (dim - 1 - i) * dim + dim - 1;
			matrix[first] = matrix[sec] - distance;
		}
		if (init) {
			std::chrono::system_clock::time_point start;
			if (world_rank == MASTER_NODE) {
				start = std::chrono::system_clock::now();
			}
			if (method == SCATTER) {
				if(world_rank == 0)
					std::cout << "scatter used" << std::endl;
				partitionMPIScatter();
			}
			else if(method == SEND_RECV) {
				if(world_rank == 0)
					std::cout << "send and recv used" << std::endl;
				partitionMPI();
			}
			else {
				if(world_rank == 0)
					std::cout << "scatter used" << std::endl;
				partitionMPIScatter();
			}
			if (world_rank == MASTER_NODE) {
				auto end = std::chrono::system_clock::now();
				timeTaken = (long long)std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
			}
		}
		else {
			auto start = std::chrono::system_clock::now();
			int loopcount = 0;
			do {
				loopcount++;
				r = calcualteCell(1, dim - 2, matrix);
			} while (r >= limit);
			std::cout << "looped for " << loopcount << std::endl;

			auto end = std::chrono::system_clock::now();
			timeTaken = (long long)std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
		}
		break;
	case 1:
		matrix[0] = 1.0;
		matrix[dim * dim - 1] = 1.0;
		for (int i = 1; i < dim; i++) {
			matrix[i] = matrix[i-1] - distance;
			matrix[(dim - 1) * dim + (dim-1) - i] = matrix[(dim - 1) * dim + dim - i] - distance;
			matrix[i * dim] = matrix[(i - 1) * dim] - distance;
			int sec = (dim - i) * dim + dim - 1;
			int first = (dim - 1 - i) * dim + dim - 1;
			matrix[first] = matrix[sec] - distance;
		}

		if (init) {
			std::chrono::system_clock::time_point start;
			if (world_rank == MASTER_NODE) {
				start = std::chrono::system_clock::now();
			}
			if (method == SCATTER) {
				if(world_rank == 0)
					std::cout << "scatter used" << std::endl;
				partitionMPIScatter();
			}
			else if(method == SEND_RECV) {
				if(world_rank == 0)
					std::cout << "send and recv used" << std::endl;
				partitionMPI();
			}
			else {
				if(world_rank == 0)
					std::cout << "scatter used" << std::endl;
				partitionMPIScatter();
			}
			if (world_rank == MASTER_NODE) {
				auto end = std::chrono::system_clock::now();
				timeTaken = (long long)std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
			}
		}
		else {
			auto start = std::chrono::system_clock::now();
			int loopcount = 0;
			do {
				loopcount++;
				r = calcualteCell(1, dim - 2, matrix);
			} while (r >= limit);
			std::cout << "looped for " << loopcount << std::endl;

			auto end = std::chrono::system_clock::now();
			timeTaken = (long long)std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
		}

		break;
	}
	return timeTaken;
}

double RelaxationMethod::calcualteCell(int startRow, int lastRow, Matrix &_matrix) {
	double r = std::numeric_limits<double>::min();

	for (int row = startRow; row <= lastRow; row++) {
		for (int col = 1; col < dim - 1; col++) {
			double rnew = residuum(row, col, _matrix);
			_matrix.set(row, col, _matrix.get(row, col)  + (1.0 / 4.0) * rnew);
			r = std::max(r, rnew);
		}
	}
	return r;
}

void RelaxationMethod::partitionMPI() {
	int size = world_size - 1;
	int rows = dim / (size);
	rows = rows <= 0 ? 1 : rows;

	// master that partitions data and sends to slave processes
	if (world_rank == MASTER_NODE) {
		std::cout << "starting ..." << std::endl;
		while (true) {
			double r = std::numeric_limits<double>::min();

			// send submatrices to all slaves
			for (int i = 1; i < world_size; i++) {
				int startRow, endRow;
				startEndRows(startRow, endRow, rows, i);
				Matrix subMatrix = matrix.getSubMatrix(startRow, endRow);
				MPI_Send(subMatrix.getData(), subMatrix.size(), MPI_DOUBLE, i, MATRIX_BEFORE_CALC_TAG, MPI_COMM_WORLD);
			}

			// receive the results form the slaves
			std::vector<Matrix> subMatrices;
			std::vector<double> residuums(size, 0.0);
			MPI_Status status;
			for (int i = 1; i < world_size; i++) {
				int startRow, endRow;
				startEndRows(startRow, endRow, rows, i);
				Matrix subMatrix = matrix.getSubMatrix(startRow, endRow);
				subMatrices.push_back(subMatrix);
				MPI_Recv(subMatrices[i - 1].getData(), subMatrix.size(), MPI_DOUBLE, i, MATRIX_AFTER_CALC_TAG, MPI_COMM_WORLD, &status);
				MPI_Recv(&(residuums[i-1]), 1, MPI_DOUBLE, i, RESIDUUM_TAG, MPI_COMM_WORLD, &status);
			}

			// insert the result matrices back into the original one for the next iteration
			for (int i = 1; i < world_size; i++) {
				int startRow, endRow;
				startEndRows(startRow, endRow, rows, i);
				Matrix importantContent = subMatrices[i - 1].getSubMatrix(1, endRow - startRow - 1);
				matrix.insertSubmatrix(importantContent, startRow + 1);
				r = std::max(r, residuums[i - 1]);
			}
			if (r <= limit && r != std::numeric_limits<double>::min()) {
				//matrix.printMatrix();
				return;
			}
		}
	}
	else {
		MPI_Status status;
		int startRow, endRow;
		startEndRows(startRow, endRow, rows, world_rank);
		Matrix subMatrix = matrix.getSubMatrix(startRow, endRow);
		double r = std::numeric_limits<double>::min();

		while (true) {
			r = std::numeric_limits<double>::min(); 

			//wait for Master to send submatrix
			MPI_Recv(subMatrix.getData(), subMatrix.size(), MPI_DOUBLE, MASTER_NODE, MATRIX_BEFORE_CALC_TAG, MPI_COMM_WORLD, &status);
			int start = 1, end = endRow - startRow - 1;

			double rNew = calcualteCell(start, end, subMatrix);
			r = std::max(rNew, r);

			// send results back to master
			MPI_Send(subMatrix.getData(), subMatrix.size(), MPI_DOUBLE, MASTER_NODE, MATRIX_AFTER_CALC_TAG, MPI_COMM_WORLD);
			MPI_Send(&r, 1, MPI_DOUBLE, MASTER_NODE, RESIDUUM_TAG, MPI_COMM_WORLD);
		}
	}
}


void RelaxationMethod::partitionMPIScatter() {
	int rows = dim / world_size;
	rows = rows <= 0 ? 1 : rows;
	int startRow, endRow;
	startEndRowsScatter(startRow, endRow, rows, world_rank);

	std::vector<int> sendCounts = calculateSendCounts(rows);
	std::vector<int> sendDispls = calculateSendDispls(rows);
	std::vector<int> recvCounts = calculateRecvCounts(rows);
	std::vector<int> recvDispls = calculateRecvDispls(rows);
	std::vector<double> rValues(world_size);

	int start = 1, end = endRow - startRow - 1;

	int done = 0;
	int loopcount = 0;
	while (true) {
		Matrix subMatrix = matrix.getSubMatrix(startRow, endRow);

		// send the partitioned matrix to all processes
		MPI_Scatterv(matrix.getData(), sendCounts.data(), sendDispls.data(), MPI_DOUBLE, subMatrix.getData(), sendCounts[world_rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);

		double r = calcualteCell(start, end, subMatrix);

		double *matrixStart = subMatrix.getData() + dim;
		matrixStart = world_rank == 0 ? subMatrix.getData() : subMatrix.getData() + dim;
		int matrixSize = subMatrix.size() - (2 * dim);
		matrixSize = world_rank == 0 || world_rank == world_size - 1 ? subMatrix.size() - dim : subMatrix.size() - (2 * dim);

		// get the submatrices and merge them into the old matrix, also get all calcualted r values
		MPI_Gatherv(matrixStart, matrixSize, MPI_DOUBLE, matrix.getData(), recvCounts.data(), recvDispls.data(), MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
		MPI_Gather(&r, 1, MPI_DOUBLE, rValues.data(), 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);

		// search the max out of all r values
		double rMax = std::numeric_limits<double>::min();

	
		if (world_rank == 0) {
			for (auto value : rValues) {
				rMax = std::max(rMax, value);
			}
		}

		loopcount++;

		if (world_rank == 0 && rMax <= limit && rMax != std::numeric_limits<double>::min()) {
			std::cout << "loopcount: " << loopcount << std::endl;
			std::cout << "limit reached" << std::endl;
			//matrix.printMatrix();
			done = 1;
		}
		MPI_Bcast(&done, 1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD);
		if (done == 1) {
			std::cout << "Process: " << world_rank << " is done!"<< std::endl;
			return;
		}
	}
}


std::vector<int> RelaxationMethod::calculateSendCounts(int rows) {
	std::vector<int> buffer(world_size);
	int startRow, endRow;

	for (int i = 0; i < world_size; i++) {
		startEndRowsScatter(startRow, endRow, rows, i);
		buffer[i] = (endRow - startRow + 1) * dim;
	}
	return buffer;
}


std::vector<int> RelaxationMethod::calculateSendDispls(int rows) {
	std::vector<int> buffer(world_size);
	int startRow, endRow;
	for (int i = 0; i < world_size; i++) {
		startEndRowsScatter(startRow, endRow, rows, i);
		buffer[i] = dim * startRow;
	}
	return buffer;
}

std::vector<int> RelaxationMethod::calculateRecvCounts(int rows) {
	std::vector<int> buffer(world_size);
	int startRow, endRow;

	for (int i = 0; i < world_size; i++) {
		startEndRowsGather(startRow, endRow, rows, i);
		buffer[i] = (endRow - startRow + 1) * dim;
	}
	return buffer;
}

std::vector<int> RelaxationMethod::calculateRecvDispls(int rows) {
	std::vector<int> buffer(world_size);
	int startRow, endRow;
	for (int i = 0; i < world_size; i++) {
		startEndRowsGather(startRow, endRow, rows, i);
		buffer[i] = dim * startRow;
	}
	return buffer;
}


void RelaxationMethod::initMPI(int argc, char **argv) {
	if (init == false) {
		MPI_Init(&argc, &argv);
		MPI_Comm_size(MPI_COMM_WORLD, &world_size);
		MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	}
	init = true;
}

RelaxationMethod::~RelaxationMethod()
{
	if (init) {
		MPI_Finalize();
	}
}


void RelaxationMethod::writeToCSV(std::string name, long long time) {
	std::ofstream csvFile;
	csvFile.open(name + ".csv");
	if (init) {
		csvFile << "with MPI,time taken:," << time << ",dimensions:," << dim - 1 << "x" << dim - 1 << ",\n";
	}
	else {
		csvFile << "no MPI,time taken:," << time << ",dimensions:," << dim - 1 << "x" << dim - 1 << ",\n";
	}
	for (int row = 0; row < dim; row++) {
		for (int col = 0; col < dim; col++) {
			csvFile << matrix.get(row, col) << ",";
		}
		csvFile << "\n";
	}
	csvFile.close();
}


void RelaxationMethod::startEndRows(int &startRow, int &endRow, int rows, int _world_rank) {
	startRow = ((_world_rank - 1) * rows) - 1;
	if (startRow < 0) {
		startRow = 0;
	}
	endRow = (_world_rank * rows);
	if (endRow >= dim) {
		endRow = dim - 1;
	}
	// if uneven numbers the last process rows until the end
	if (_world_rank == world_size - 1) {
		endRow = dim - 1;
	}
}

void RelaxationMethod::startEndRowsScatter(int &startRow, int &endRow, int rows, int _world_rank) {
	startRow = ((_world_rank) * rows) - 1;
	if (startRow < 0) {
		startRow = 0;
	}
	endRow = ((_world_rank + 1) * rows);
	if (endRow >= dim) {
		endRow = dim - 1;
	}
	// if uneven numbers the last process rows until the end
	if (_world_rank == world_size - 1) {
		endRow = dim - 1;
	}
}


void RelaxationMethod::startEndRowsGather(int &startRow, int &endRow, int rows, int _world_rank) {
	startRow = ((_world_rank) * rows);
	if (startRow < 0) {
		startRow = 0;
	}
	endRow = ((_world_rank + 1) * rows) - 1;
	if (endRow >= dim) {
		endRow = dim - 1;
	}
	// if uneven numbers the last process rows until the end
	if (_world_rank == world_size - 1) {
		endRow = dim - 1;
	}
}

double RelaxationMethod::faultFunction(int row, int col, Matrix &_matrix) {
	if (edgeCase == 1) {
		const double PI = 3.1415926535897932384626433832795028841971693993751058209;
		return 2 * (PI * PI) * std::sin(PI * row) * std::sin(PI * col);
	} 
	return 0;
}

double RelaxationMethod::residuum(int row, int col, Matrix &_matrix) {
	return faultFunction(row, col, _matrix) * (distance * distance) - (4 * _matrix.get(row, col) - _matrix.get(row - 1, col) - _matrix.get(row + 1, col) - _matrix.get(row, col - 1) - _matrix.get(row, col + 1));
}
