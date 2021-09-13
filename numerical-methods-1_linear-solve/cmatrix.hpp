#pragma once

#include <type_traits> // std::is_arithmetic<T> is used to ensure matrix elements are of a numeric type
#include <stdexcept> // std::runtime_error
#include <cstring> // memcpy() and memmove() for copying-moving entire matrices at once

#include <fstream> // required for parsing matrices from files
#include <string>

#include <iostream> // required for printing matrices to console
#include <iomanip>



// # CMatrix #
// Simplistic matrix implementation
// - uses C-style arrays internally
// - stores data in a contiguous fashion
// - [i][j] style of indexation 
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
		CMatrix()
	{
		this->allocate_for_size(other._rows, other._cols);

		// Copy data a single contiguous chunk
		memcpy(_data[0], other._data[0], sizeof(T) * rows * cols);
	}

	// Move constructor
	CMatrix(CMatrix<T> &&other) noexcept :
		CMatrix()
	{
		// No need co copy/move data itself, we can just rewire pointers
		this->_rows = other._rows;
		this->_cols = other._cols;
		this->_data = other._data;

		// Reset pointer for other object so it doesn't destroy our data when going out of scope
		other._data = nullptr;
	}

	CMatrix(size_t rows, size_t cols) :
		CMatrix()
	{
		this->allocate_for_size(rows, cols);
	}

	// Allocates memory for given dimensions, deallocates old data if necessary
	void allocate_for_size(size_t rows, size_t cols) {
		// Deallocate old data if present
		if (_data) {
			delete[] _data[0];
			delete[] _data;
		}

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

	// Const references are used as 'getters' for the purpose of beauty
	const size_t& rows;
	const size_t& cols;

	CMatrix<T>& operator= (const CMatrix<T>& other) {;
		this->allocate_for_size(other._rows, other._cols);

		// Copy data a single contiguous chunk
		memcpy(_data[0], other._data[0], sizeof(T) * rows * cols);

		return *this;
	}

	T* operator[] (size_t i) { // allows us to use matrix[i][j] style of indexation
		return _data[i];
	}

	const T* operator[] (size_t i) const {
		return _data[i];
	}

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

	CMatrix<T> operator* (const CMatrix<T>& other) const {
		if (_cols != other._rows)
			throw std::runtime_error("ERROR: Matrix multipliation impossible for given dimensions");

		CMatrix<T> result(_rows, other._cols);
		result.fill_with_zeroes();

		for (size_t i = 0; i < _rows; ++i)
			for (size_t k = 0; k < _cols; ++k)
				for (size_t j = 0; j < other._cols; ++j)
					result[i][j] += _data[i][k] * other._data[k][j];
		// note that naive loop order would be [i]->[j]->[k], swapping [k] and [j]
		// loops reduces the number of cache misses since we access contiguously
		// stored elements in the inner-most loop

		return result;
	}

	CMatrix<T> get_transposed() const { // NORE: constructs new matrix, does not modify the original one
		CMatrix<T> transposed(_cols, _rows);

		for (size_t i = 0; i < _rows; ++i)
			for (size_t j = 0; j < _cols; ++j)
				transposed[i][j] = _data[j][i];

		return transposed;
	}

	CMatrix<T>& transpose_square() {
		for (size_t i = 0; i < _rows; ++i)
			for (size_t j = i; j < _cols; ++j)
				std::swap(_data[i][j], _data[j][i]);

		return *this;
	}

	void swap_rows(size_t row1, size_t row2) {
		for (size_t j = 0; j < _cols; ++j) std::swap(_data[row1][j], _data[row2][j]);
			/// Perhaps swapping pointers would be more efficient, but it complicates deallocation
			/// as _data[0] is no longer guaranteed to point to the beginning of allocated memory,
			/// also effect of such swaps on cache utilization should probably be bechmarked for safety
	}	

	void fill_with_zeroes() {
		// Using memset allows us to fill matrix with zeroes much faster than element-wise operations
		memset(_data[0], 0, sizeof(T) * _rows * _cols);
	}

	T norm_cubic() const {
		// Cubic norm is max_i { SUM_j |matrix[i][j]|}
		T maxSum(0);

		// Go over each row calculating SUM_j |matrix[i][j]|, select max_i
		for (size_t i = 0; i < _rows; ++i) {
			T sum(0);

			for (size_t j = 0; j < _cols; ++j) sum += std::abs(_data[i][j]);

			if (sum > maxSum) maxSum = sum;
		}

		return maxSum;
	}
};