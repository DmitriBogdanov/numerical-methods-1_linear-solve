#pragma once

#include "cmatrix.hpp"
#include "math_helpers.hpp"


 
template<typename T>
class GaussianEliminationSolver {
	CMatrix<T> _matrix;

public:
	// Creates solver for passed linear system
	GaussianEliminationSolver(const CMatrix<T> &matrix) {
		// Ensure correct matrix dimensions
		if (matrix.cols != matrix.rows + 1)
			throw std::runtime_error("ERROR: Cannot create linear system with dimensions (" +
				std::to_string(matrix.rows) + ", " + std::to_string(matrix.cols) + ").");

		_matrix = matrix;
	}

	CMatrix<T> solve() {
		this->gauss_forward_elimination();

		this->gauss_backward_elimination();

		// Return solution as a column
		CMatrix<T> solution(_matrix.rows, 1);
		for (size_t i = 0; i < _matrix.rows; ++i) solution[i][0] = _matrix[i][_matrix.cols - 1];

		return solution;
	}

	// When we need to solve multiple systems it's better to reset matrix without reallocating the entire object
	void reset_matrix(const CMatrix<T> &matrix) {
		_matrix = matrix;
	}

	void add_disturbance(T magnitude) {
		// Add rand[-magnitude, +magnitude] elementwise to the last column 
		for (size_t i = 0; i < _matrix.rows; ++i) {
			
			const T randCoef = static_cast<T>(rand()) / RAND_MAX * 2 - 1; // random value in [-1, 1] range
			_matrix[i][_matrix.cols - 1] += magnitude * randCoef;
		}
	}

	CMatrix<T> get_b() const {
		// Returns last column as its own object
		CMatrix<T> b(_matrix.rows, 1);
		for (size_t i = 0; i < _matrix.rows; ++i) b[i][0] = _matrix[i][_matrix.cols - 1];

		return b;
	}

private:
	void gauss_forward_elimination() {
		for (size_t k = 0; k < _matrix.rows; ++k) {
			// Select leading row and swap it with current row if needed
			const auto leadingRow = this->find_leading_row(k);
			if (leadingRow != k) _matrix.swap_rows(leadingRow, k);

			// If the diagonal element of a leading row is zero, matrix is singular => throw
			if (isZero(_matrix[k][k]))
				throw std::runtime_error("ERROR: Could not solve the system with singular matrix");

			// Set diagonal element to 1 by multiplying its row
			const T factor = static_cast<T>(1) / _matrix[k][k];
			for (size_t j = k; j < _matrix.cols; ++j) _matrix[k][j] *= factor;

			// Substract current row from the following ones completing the step
			for (size_t i = k + 1; i < _matrix.rows; ++i) {
				const T firstElement = _matrix[i][k];

				for (size_t j = k; j < _matrix.cols; ++j) _matrix[i][j] += _matrix[k][j] * firstElement;
			}
		}
	}

	void gauss_backward_elimination() {
		// Go backwards until our matrix is diagonal
		for (size_t k = _matrix.rows - 1; k > 0; --k)
			for (size_t i = 1; i <= k; ++i) {
				const T factor = -_matrix[k - i][k];
				for (size_t j = k; j <= _matrix.cols; ++j) _matrix[k - i][j] += _matrix[k][j] * factor;
			}
	}

	void add_to_row(size_t destRow, size_t sourceRow, T factor) {
		// Add rows in vector-like fashion: 'destRow += sourceRow * factor'
		for (size_t j = 0; j < _matrix.cols; ++j) _matrix[destRow][j] += _matrix[sourceRow][j] * factor;
	}

	size_t find_leading_row(size_t k) {
		// Select row with the largest 'first' element
		size_t indexOfMax = k; 

		for (size_t i = k; i < _matrix.rows; ++i)
			if (std::abs(_matrix[i][k]) > std::abs(_matrix[indexOfMax][k]))
				indexOfMax = i;

		return indexOfMax;
	}
};