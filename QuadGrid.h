#pragma once
#include <vector>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <mpi.h>
#include <chrono>
#include <fstream>
#include <string>

#include "Matrix.h"

const int MASTER_NODE = 0;

const int MATRIX_BEFORE_CALC_TAG = 0;
const int MATRIX_AFTER_CALC_TAG = 1;
const int RESIDUUM_TAG = 2;

#define SEQ 0
#define SCATTER 1 
#define SEND_RECV 2

class QuadGrid
{
public:
	QuadGrid(int n, double limit);
	int start(int edgeCase, int method = 0);
	void initMPI(int argc, char **argv);
	void writeToCSV(std::string name, long long time);
	~QuadGrid();

private:
	int world_size;
	int world_rank;
	bool init = false;
	double limit;
	int dim;
	double distance;
	Matrix matrix;

	void partitionMPI();
	void partitionMPIScatter();
	std::vector<int> calculateSendCounts(int rows);
	std::vector<int> calculateSendDispls(int rows);
	std::vector<int> calculateRecvCounts(int rows);
	std::vector<int> calculateRecvDispls(int rows);
	double faultFunction(int row, int col, Matrix &matrix);
	double residuum(int row, int col, Matrix &matrix);
	double calcualteCell(int startRow, int lastRow, Matrix &_matrix, MPI_Request* request = NULL);
	void startEndRows(int &startRow, int &endRow, int rows, int _world_rank);
	void startEndRowsScatter(int &startRow, int &endRow, int rows, int _world_rank);
	void startEndRowsGather(int &startRow, int &endRow, int rows, int _world_rank);
};

