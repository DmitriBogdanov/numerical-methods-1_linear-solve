#pragma once

#include "cmatrix.hpp"



// # LinearSystem #
// Class that handles general functionality of a linear system, parses, prints, calculates solution residuals and etc
template<typename T>
class LinearSystem {
	CMatrix<T> _matrix; // square
	CMatrix<T> _column;
		// square matrix and column are separated for the purpose of easier substitution when
		// checking solutions and calculating residues

public:
	// Initializes system based on (n, n + 1) matrix
	LinearSystem(const CMatrix<T> &matrix) {
		// Ensure correct matrix dimensions
		if (matrix.cols != matrix.rows + 1)
			throw std::runtime_error("ERROR: Cannot create linear system with dimensions (" +
				std::to_string(matrix.rows) + ", " + std::to_string(matrix.cols) + ").");

		// Fill square matrix and column
		const size_t N = matrix.rows;
		_matrix.allocate_for_size(N, N);
		_column.allocate_for_size(N, 1);

		
		for (size_t i = 0; i < N; ++i) {
			memcpy(_matrix[i], matrix[i], sizeof(T) * N); // copy N elements at once to the square matrix
			_column[i][0] = matrix[i][matrix.cols - 1]; // copy last column to _column
		}
	}

	// Initializes system based on (n, n) matrix and a (n, 1) column
	LinearSystem(const CMatrix<T> &A, const CMatrix<T> &b) {
		// Ensure correct matrix dimensions
		if (A.cols != A.rows)
			throw std::runtime_error("ERROR: Cannot create linear system with dimensions (" +
				std::to_string(A.rows) + ", " + std::to_string(A.cols) + ").");

		_matrix = A;
		_column = b;
	}

	// Print matrix to console
	void print() const {
		// Last column is visually separated 
		for (size_t i = 0; i < _matrix.rows; ++i) {
			std::cout << "[ ";
			for (size_t j = 0; j < _matrix.cols; ++j) std::cout << std::setw(10) << _matrix[i][j] << " ";
			std::cout << "|" << std::setw(10) << _column[i][0] << " ]\n";
		}
	}

	T residual_cubic(const CMatrix<T> &solution) const {
		// Substitute given solution and find resulting column
		CMatrix<T> column = _matrix * solution;

		// Find difference between calculated column and original column
		for (size_t i = 0; i < _matrix.rows; ++i)
			column[i][0] -= _column[i][0];

		// Return octahedric norm
		return column.norm_cubic();
	}
};

