#pragma once

#include "linear_system.hpp"


 
template<typename T>
class GaussianEliminationSolver : public LinearSystem<T> {
	using LinearSystem<T>::_matrix;
		// allows us to use '_matrix' insted of writing 'LinearSystem<T>::_matrix' or 'this->_matrix' every time

public:
	// Fills matrix based on given file
	explicit GaussianEliminationSolver(const std::string &filepath) :
		LinearSystem<T>(filepath)
	{}

	CMatrix<T> solve_through_gaussian_elimination() {
		this->gauss_forward_elimination();

		this->gauss_backward_elimination();

		// Return solution as a column
		CMatrix<T> solution(_matrix.rows, 1);
		for (size_t i = 0; i < _matrix.rows; ++i) solution[i][0] = _matrix[i][_matrix.cols - 1];

		return solution;
	}

private:
	void gauss_forward_elimination() {
		for (size_t k = 0; k < _matrix.rows; ++k) {
			// Select leading row and swap it with current row if needed
			const auto leadingRow = this->find_leading_row(k);
			if (leadingRow != k) _matrix.swap_rows(leadingRow, k);

			// If the diagonal element of a leading row is zero, matrix is singular => throw
			if (_matrix[k][k] == static_cast<T>(0))
				throw std::runtime_error("ERROR: Could not solve the system with singular matrix");

			// Set diagonal element to 1 by multiplying its row
			const T factor = static_cast<T>(1) / _matrix[k][k];
			for (size_t j = k; j < _matrix.cols; ++j) _matrix[k][j] *= factor;

			// Substract current row from the following ones completing the step
			for (size_t i = k + 1; i < _matrix.rows; ++i) {
				const T firstElement = _matrix[i][k];
				this->add_to_row(i, k, -firstElement);
			}
		}
	}

	void gauss_backward_elimination() {
		// Go backwards untill our matrix is diagonal
		for (size_t k = _matrix.rows - 1; k > 0; --k)
			for (size_t i = 1; i <= k; ++i)
				this->add_to_row(k - i, k, - _matrix[k - i][k]);
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