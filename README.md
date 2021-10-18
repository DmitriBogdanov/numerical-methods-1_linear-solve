# Numerical methods 1 / Linear Solver

Contains implementations of following methods for solving linear systems:
* Gaussian elimination
* QR-Decomposition using Givens rotations

Note that present implementations are intended for study and analyzing properties of aforementioned methods, as such they are not meant to be used in any sort of high-performance production code.

## Compilation

* Recommended compiler: MSVC v142
* Requires C++17 support

## Usage

Для управления программой используется config-файл следующего формата:

Input is a .dat file, containing floating-point matrix that represents any linear system. To configure input file, output path and other parameters, place config file of the following format into the same folder as executable:

* Line 1: [input relative path]
* Line 2: [output relative path]
* Line 3: [number of iterations used when estimating condition number]

## Version history

* 00.07
    * Implemented computation of conditional number through its definition
    * Optimized matrix inverse and condition number estimate methods though the use of QR-decomposition with various right–hand sides 

* 00.06
    * Matrix norm function (octahedrical norm)

* 00.05
    * Implemented condition number estimate thorough the use of residuals
    * Added matrix inverse function as a testing ground for implemeted methods
* 00.04
    * Implemented QR-decomposition function (modification using Givens rotations)
    * Implemented QR-decomposition method for solving linear systems

* 00.03
    * Memory leak checking in debug mode

* 00.02
    * Matrix norm function (cubic norm)
    * Implemented matrix arithmetic operators
    * Added functionality of calculating solution residuals
    * All results are now saved as .dat files
    * Structural changes
    * Optimizations

* 00.01
    * Implemented Gaussian elimination