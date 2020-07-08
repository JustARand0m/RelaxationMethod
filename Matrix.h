#pragma once
#include <vector>
#include <iostream>
#include <iomanip>

class Matrix {
public: 
	Matrix(const Matrix &m);
	Matrix(int rows, int cols);
	Matrix& operator=(const Matrix &other);
	Matrix getSubMatrix(int startRow, int endRow);
	bool insertSubmatrix(Matrix m, int startRow);
	bool set(int row, int col, double value);
	double get(int row, int col);
	double& operator [] (int i);
	void printMatrix();
	double *getData();
	int maxRows();
	int maxCols();
	int size();
private:
	void setVector(std::vector<double> vec);

	std::vector<double> matrix;
	int rows, cols;
};
