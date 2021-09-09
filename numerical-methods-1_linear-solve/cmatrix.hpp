#pragma once

#include <type_traits> // std::is_arithmetic<T> is used to ensure matrix elements are of a numeric type
#include <stdexcept> // std::runtime_error
#include <cstring> // memcpy() and memmove() for copying-moving entire matrices at once

#include <fstream> // required for parsing matrices from files
#include <string>

#include <iostream> // required for printing matrices to console
#include <iomanip>



// # CMatrix #
// Matrix class
template<typename T>
class CMatrix {
	// Template compiles only if T is numeric
	static_assert(std::is_arithmetic<T>::value, "T must be numeric");

	size_t _rows;
	size_t _cols;
	T** _data;

public:
	// Default constructor
	CMatrix() :
		_rows(0),
		_cols(0),
		_data(nullptr),
		rows(_rows),
		cols(_cols)
	{}

	// Copy constructor
	CMatrix(const CMatrix<T> &other) :
		rows(_rows),
		cols(_cols)
	{
		this->allocate_for_size(other._rows, other._cols);

		memcpy(_data[0], other._data[0], sizeof(T) * rows * cols);
	}

	// Move constructor
	CMatrix(CMatrix<T> &&other) :
		rows(_rows),
		cols(_cols)
	{
		this->allocate_for_size(other._rows, other._cols);

		auto t = sizeof(other._data[0]);

		memmove(_data[0], other._data[0], sizeof(T) * rows * cols);
	}

	CMatrix(size_t rows, size_t cols) :
		rows(_rows),
		cols(_cols)
	{
		this->allocate_for_size(rows, cols);
	}

	// Allocates memory for given dimensions
	void allocate_for_size(size_t rows, size_t cols) {
		_rows = rows;
		_cols = cols;

		_data = new T*[rows];

		// Allocate memory in a contiguous chunk for the sake of performance
		_data[0] = new T[rows * cols];
		for (size_t i = 1; i < rows; ++i) _data[i] = _data[i - 1] + cols;
	}

	~CMatrix() {
		if (!_data) return;

		delete[] _data[0];
		delete[] _data;
	}

	T* operator[] (size_t i) { // allows us to use matrix[i][j] style of indexation
		return _data[i];
	}

	const T* operator[] (size_t i) const {
		return _data[i];
	}

	// Const references are used as 'getters' for the purpose of beauty
	const size_t &rows;
	const size_t &cols;

	void print() const {
		for (size_t i = 0; i < _rows; ++i) {
			std::cout << "[ ";
			for (size_t j = 0; j < _cols; ++j) std::cout << std::setw(10) << _data[i][j] << " ";
			std::cout << " ]\n";
		}

		std::cout << std::endl;
	}

	void parse_from_file(const std::string & filepath) {
		// Open file, throw if unsuccesfull
		std::ifstream inFile(filepath);

		if (!inFile.is_open())
			throw std::runtime_error("ERROR: Could not open file <" + filepath + ">");

		// Read matrix size and allocate memory
		size_t n;
		inFile >> n;
		this->allocate_for_size(n, n + 1);

		// Read matrix elements, since it's stored contiguously in memory we only need 1 loop
		for (size_t i = 0; i < _rows * _cols; ++i) inFile >> _data[0][i];
	}

	void swap_rows(size_t row1, size_t row2) {
		for (size_t j = 0; j < _cols; ++j) std::swap(_data[row1][j], _data[row2][j]);
			/// Perhaps swapping pointers would be more efficient, but it complicates deallocation
			/// as _data[0] is no longer guaranteed to point to the beginning of allocated memory,
			/// also effect of such swaps on cache utilization should probably be bechmarked for safety
	}	
};