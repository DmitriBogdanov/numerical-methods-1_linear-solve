#pragma once

#include <type_traits>
#include <stdexcept>
#include <cstring> // memcpy � memmove

#include <fstream>
#include <string>

#include <iostream>
#include <iomanip>



template<typename T>
class CMatrix {
	// ������ ������������� ������ ��� ��������� �����
	static_assert(std::is_arithmetic<T>::value, "T must be numeric");

	size_t _rows;
	size_t _cols;
	T** _data;

public:
	CMatrix() :
		_rows(0),
		_cols(0),
		_data(nullptr),
		rows(_rows),
		cols(_cols)
	{}

	// Copy-�����������
	CMatrix(const CMatrix<T> &other) :
		rows(_rows),
		cols(_cols)
	{
		this->allocate_for_size(other._rows, other._cols);

		memcpy(_data[0], other._data[0], sizeof(T) * rows * cols);
	}

	// Move-�����������
	CMatrix(CMatrix<T> &&other) :
		rows(_rows),
		cols(_cols)
	{
		this->allocate_for_size(other._rows, other._cols);

		auto t = sizeof(other._data[0]);

		memmove(_data[0], other._data[0], sizeof(T) * rows * cols);
	}

	explicit CMatrix(size_t rows, size_t cols) :
		rows(_rows),
		cols(_cols)
	{
		this->allocate_for_size(rows, cols);
	}

	void allocate_for_size(size_t rows, size_t cols) {
		_rows = rows;
		_cols = cols;

		_data = new T*[rows];

		// �� ����������� ������������� �������� ������ ��� ������ ����� ����������� ������
		_data[0] = new T[rows * cols];
		for (size_t i = 1; i < rows; ++i) _data[i] = _data[i - 1] + cols;
	}

	~CMatrix() {
		if (!_data) return;

		delete[] _data[0];
		delete[] _data;
	}

	T* operator[] (size_t i) { // ��������� ������������ ���������� ���� matrix[i][j]
		return _data[i];
	}

	const T* operator[] (size_t i) const { // ��������� ������������ ���������� ���� matrix[i][j]
		return _data[i];
	}

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
		// ��������� ����, ������� ������ � ������ �������
		std::ifstream inFile(filepath);

		if (!inFile.is_open())
			throw std::runtime_error("ERROR: Could not open file <" + filepath + ">");

		// ��������� ������ �������, �������� ������
		size_t n;
		inFile >> n;
		this->allocate_for_size(n, n + 1);

		// ��������� �������� �������
		for (size_t i = 0; i < _rows * _cols; ++i) inFile >> _data[0][i];
	}

	void swap_rows(size_t row1, size_t row2) { /// ������� ����� ���������
		for (size_t j = 0; j < _cols; ++j) std::swap(_data[row1][j], _data[row2][j]);
	}	
};