#pragma once

#include "cmatrix.hpp"



// # LinearSystem #
// Class that handles general functionality of a linear system, parses, prints, calculates solution residuals and etc
template<typename T>
class LinearSystem {
protected:
	CMatrix<T> _matrix;

public:
	// Fills matrix based on given file
	explicit LinearSystem(const std::string& filepath) {
		_matrix.parse_from_file(filepath);
	}

	// Print matrix to console
	void print() const {
		// Unlike _matrix.print() last column is visually separated 
		for (size_t i = 0; i < _matrix.rows; ++i) {
			std::cout << "[ ";
			for (size_t j = 0; j < _matrix.cols - 1; ++j) std::cout << std::setw(10) << _matrix[i][j] << " ";
			std::cout << "|" << std::setw(10) << _matrix[i][_matrix.cols - 1] << " ]\n";
		}
	}
};

