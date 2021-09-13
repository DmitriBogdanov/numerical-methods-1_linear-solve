#pragma once

#include <tuple>

#include "cmatrix.hpp"
#include "math_helpers.hpp"

#include "linear_system.hpp"



template<typename T>
class QRDecompositionSolver {
	CMatrix<T> _matrix;
	CMatrix<T> _b;
	///CMatrix<T> _Q;
	///CMatrix<T> _R;

	size_t _N;

public:
	QRDecompositionSolver(const CMatrix<T> &matrix) {
		// Ensure correct matrix dimensions
		if (matrix.cols != matrix.rows + 1)
			throw std::runtime_error("ERROR: Cannot create linear system with dimensions (" +
				std::to_string(matrix.rows) + ", " + std::to_string(matrix.cols) + ").");

		_N = matrix.rows;

		// Allocate matrices
		_matrix.allocate_for_size(_N, _N);
		_b.allocate_for_size(_N, 1);
		///_Q.allocate_for_size(_N, _N);
		///_R.allocate_for_size(_N, _N);

		// Fill square matrix and column
		for (size_t i = 0; i < _N; ++i) {
			memcpy(_matrix[i], matrix[i], sizeof(T) * _N); // copy N elements at once to the square matrix
			_b[i][0] = matrix[i][_N]; // copy last column to _column
		}

		
	}

	// @return 1 => solution
	// @return 2 => Q
	// @return 3 => R
	std::tuple<CMatrix<T>, CMatrix<T>, CMatrix<T>> solve() {
		// Decompose matrix into QR
		// Q - orthogonal matrix => Q^-1 == Q^T
		// R - upper triangular matrix
		auto [Q, R] = this->QR_decompose();

		// Go from system
		// Ax = b
		// To an equivalent system
		// Rx = Q^T b
		_matrix = R;
		_b = Q.get_transposed() * _b;

		// Check that system matrix isn't singular
		const T det = this->get_det();
		if (isZero(det))
			throw std::runtime_error("ERROR: Could not solve the system with singular matrix");

		// Solve through Gauss backwards elimination
		for (size_t k = 0; k < _N; ++k) {
			const T factor = static_cast<T>(1) / _matrix[k][k];
			for (size_t j = k; j < _N; ++j) _matrix[k][j] *= factor;
			_b[k][0] *= factor;
		}

		for (size_t k = _N - 1; k > 0; --k) {
			for (size_t i = 1; i <= k; ++i) {
				const T factor = _matrix[k - i][k];
				_matrix[k - i][k] -= _matrix[k][k] * factor;
				_b[k - i][0] -= _b[k][0] * factor;
			}
		}

		return { _b, Q, R };
	}

private:
	// @return 1 => Q
	// @return 1 => R
	std::tuple<CMatrix<T>, CMatrix<T>> QR_decompose() {
		// Set Q to identity matrix
		CMatrix<T> Q(_N, _N);
		Q.fill_with_zeroes();

		for (size_t i = 0; i < _N; ++i)
			Q[i][i] = 1;

		// Copy original matrix into R
		CMatrix R(_matrix);

		// Givens rotations go brrr
		T c(0);
		T s(0);
		for (size_t i = 0; i < _N - 1; ++i)
			for (size_t j = i + 1; j < _N; ++j)
				if (!isZero(R[j][i])) {
					const T inverseNorm = static_cast<T>(1) / std::sqrt(sqr(R[i][i]) + sqr(R[j][i]));
					c = R[i][i] * inverseNorm;
					s = R[j][i] * inverseNorm;

					for (size_t j_ = 0; j_ < _N; ++j_) {
						const T tempR = c * R[i][j_] + s * R[j][j_];
						R[j][j_] = - s * R[i][j_] + c * R[j][j_];
						R[i][j_] = tempR;

						const T tempQ = c * Q[i][j_] + s * Q[j][j_];
						Q[j][j_] = - s * Q[i][j_] + c * Q[j][j_];
						Q[i][j_] = tempQ;
					}
				}

		Q.transpose_square();

		return { Q, R };
	}

	T get_det() {
		T det(1);
		for (size_t i = 0; i < _N; ++i) det *= _matrix[i][i];
		return det;
	}
};