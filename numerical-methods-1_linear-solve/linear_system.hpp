#pragma once

#include "cmatrix.hpp"



// Класс обобщающий функционал СЛАУ
// Содержит методы инициализации, печати, проверки невязки решений и т.д.
template<typename T>
class LinearSystem {
protected:
	CMatrix<T> _matrix;

public:
	// Конструктор заполняющий матрицу на основе файла
	explicit LinearSystem(const std::string& filepath) {
		_matrix.parse_from_file(filepath);
	}

	// Вывод системы в консоль
	void print() const {
		// Отличие от вывода непосредственно матрицы заключается в визуальном отделении последнего столбца
		for (size_t i = 0; i < _matrix.rows; ++i) {
			std::cout << "[ ";
			for (size_t j = 0; j < _matrix.cols - 1; ++j) std::cout << std::setw(10) << _matrix[i][j] << " ";
			std::cout << "|" << std::setw(10) << _matrix[i][_matrix.cols - 1] << " ]\n";
		}

		std::cout << std::endl;
	}
};

/*
template<typename T>
class Solution {
	CMatrix<T> _column;

public:
	Solution(size_t rows) {
		_column.allocate_for_size(rows, 1);
	}

	T& operator[] (size_t i) { // позволяет использовать индексацию вида matrix[i][j]
		return _column[i][0];
	}

	const T& operator[] (size_t i) const { // позволяет использовать индексацию вида matrix[i][j]
		return _column[i][0];
	}

	void print() const {
		_column.print();
	}
};
*/

