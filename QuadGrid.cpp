#include "QuadGrid.h"

QuadGrid::QuadGrid(int n, double _limit): dim(n+1), matrix(dim, dim) {
	world_rank = -1;
	world_size = -1;
	distance = 1.0 / n;
	limit = _limit;
}

int QuadGrid::start(int edgeCase) {
	long long timeTaken = -1;
	double r = 0;
	Matrix m(5, 5);
	switch (edgeCase) {
	case 1:
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
			partitionMPI();
			if (world_rank == MASTER_NODE) {
				auto end = std::chrono::system_clock::now();
				timeTaken = (long long)(end - start).count();
			}
		}
		else {
			auto start = std::chrono::system_clock::now();
			do {
				r = calcualteCell(1, dim - 2, matrix, NULL);
			} while (r >= limit);

			auto end = std::chrono::system_clock::now();
			timeTaken = (long long)(end - start).count();
		}
		break;
	case 2:
		break;
	}
	return timeTaken;
}

double QuadGrid::calcualteCell(int startRow, int lastRow, Matrix &_matrix, MPI_Request* request) {
	MPI_Status status;
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

void QuadGrid::partitionMPI() {
	int size = world_size - 1;
	int rows = dim / (size);
	rows = rows <= 0 ? 1 : rows;

	// master that partitions data and sends to slave processes
	if (world_rank == MASTER_NODE) {
		while (true) {
			double r = std::numeric_limits<double>::min();
			MPI_Request *requests_send = new MPI_Request[size];
			MPI_Status *statuses_send = new MPI_Status[size];

			// send submatrices to all slaves
			for (int i = 1; i < world_size; i++) {
				int startRow, endRow;
				startEndRows(startRow, endRow, rows, i);
				Matrix subMatrix = matrix.getSubMatrix(startRow, endRow);
				MPI_Isend(subMatrix.getData(), subMatrix.size(), MPI_DOUBLE, i, MATRIX_BEFORE_CALC_TAG, MPI_COMM_WORLD, &requests_send[i - 1]);
			}

			// wait until all messages are sent to all slaves
			MPI_Waitall(size, requests_send, statuses_send);

			delete[] requests_send;
			delete[] statuses_send;

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

			// check if r reached the limit
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

void QuadGrid::initMPI(int argc, char **argv) {
	if (init == false) {
		MPI_Init(&argc, &argv);
		MPI_Comm_size(MPI_COMM_WORLD, &world_size);
		MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	}
	init = true;
}

QuadGrid::~QuadGrid()
{
	if (init) {
		MPI_Abort(MPI_COMM_WORLD, 0);
		MPI_Finalize();
	}
}


void QuadGrid::writeToCSV(std::string name, long long time) {
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


void QuadGrid::startEndRows(int &startRow, int &endRow, int rows, int _world_rank) {
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

double QuadGrid::faultFunction(int row, int col, Matrix &_matrix) {
	return (1 / sqrt(distance)) * (4 * _matrix.get(row, col) - _matrix.get(row - 1, col) - _matrix.get(row + 1, col) - _matrix.get(row, col - 1) - _matrix.get(row, col + 1));
}

double QuadGrid::residuum(int row, int col, Matrix &_matrix) {
	return faultFunction(row, col, _matrix) * (distance * distance) - (4 * _matrix.get(row, col) - _matrix.get(row - 1, col) - _matrix.get(row + 1, col) - _matrix.get(row, col - 1) - _matrix.get(row, col + 1));
}
