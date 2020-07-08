#include "Matrix.h"


Matrix::Matrix(const Matrix &m) {
	rows = m.rows;
	cols = m.cols;
	matrix = m.matrix;
}

Matrix::Matrix(int _rows, int _cols) {
	rows = _rows;
	cols = _cols;
	matrix = std::vector<double>(rows * cols, 0.0);
}

Matrix& Matrix::operator=(const Matrix &other) {
	if (this != &other) {
		matrix = other.matrix;
		cols = other.cols;
		rows = other.rows;
	}
	return *this;
}

bool Matrix::set(int row, int col, double value) {
	if (row > rows || col > cols || row < 0 || col < 0) {
		return false;
	}
	matrix[row * cols + col] = value;
}

double Matrix::get(int row, int col) {
	return matrix[row * cols + col];
}

void Matrix::printMatrix() {
	for (unsigned int i = 0; i < matrix.size(); i++) {
		std::cout << std::fixed << std::setprecision(3) << matrix[i] << " ";
		if ((i+1) % cols == 0) {
			std::cout << std::endl;
		}
	}
	std::cout << std::endl;
}

Matrix Matrix::getSubMatrix(int startRow, int endRow) {
	Matrix newMatrix(endRow + 1 - startRow, cols);
	if (startRow > rows || endRow > rows) {
		newMatrix;
	}
	std::vector<double> newVecMatrix((matrix.begin() + (startRow * cols)), (matrix.begin() + ((endRow + 1) * cols)));
	newMatrix.setVector(newVecMatrix);
	return newMatrix;
}

bool Matrix::insertSubmatrix(Matrix m, int startRow) {
	if (startRow > rows || startRow < 0) {
		return false;
	}
	for (int row = 0; row < m.rows; row++) {
		for (int col = 0; col < cols; col++) {
			set(row + startRow, col, m.get(row, col));
		}
	}
	return true;
}


void Matrix::setVector(std::vector<double> vec) {
	matrix = vec;
}

double& Matrix::operator[] (int i) { 
	return matrix[i]; 
}

double *Matrix::getData() {
	return matrix.data();
}

int Matrix::maxCols() {
	return cols;
}

int Matrix::maxRows() {
	return rows;
}

int Matrix::size() {
	return cols * rows;
}
