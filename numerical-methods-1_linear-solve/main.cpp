#include <chrono>

#include "static_timer.hpp"
#include "linear_system.hpp"
#include "gaussian_elimination_solver.hpp"
#include "qr_decomposition_solver.hpp"

#ifdef _DEBUG
#define CHECK_MEMORY
#endif

#ifdef CHECK_MEMORY
#include <crtdbg.h>
#define _CRTDBG_MAP_ALLOC
#define new new (_NORMAL_BLOCK, __FILE__, __LINE__)
#endif



// Parse config and return input/output filepaths
// @return 1 => input filepath
// @return 2 => output filepath
std::tuple<std::string, std::string> parse_config() {
	// Parse config file
	std::string inputFilepath;
	std::string outputFilepath;

	std::ifstream inConfig("config.ini");

	if (!inConfig.is_open())
		throw std::runtime_error("ERROR: Could not oper config file.");

	std::getline(inConfig, inputFilepath);
	std::getline(inConfig, outputFilepath);

	return { inputFilepath , outputFilepath };
}


// Parse double/float matrices from given file
// @return 1 => double matrix
// @return 2 => float matrix
std::tuple<CMatrix<double>, CMatrix<float>> parse_matrices(const std::string &inputFilepath) {
	CMatrix<double> doubleMatrix;
	CMatrix<float> floatMatrix;

	doubleMatrix.parse_from_file(inputFilepath);
	floatMatrix.parse_from_file(inputFilepath);

	return { doubleMatrix, floatMatrix };
}


template<typename T>
CMatrix<T> solve_throught_gaussian_elimination(const CMatrix<T>& matrix) {
	std::cout << "\n##### Method -> Gaussian elimination\n##### Type   -> " <<
		typeid(T).name() << "\n>>> Solving...\n";

	GaussianEliminationSolver solver(matrix);

	StaticTimer::start();
	const auto solution = solver.solve();
	const auto elapsed = StaticTimer::elapsed();

	std::cout << ">>> Solved in " << elapsed << " ms\n>>> Solution:\n";
	solution.print();

	return solution;
}


template<typename T>
CMatrix<T> solve_through_QR_decomposition(const CMatrix<T>& matrix) {
	std::cout << "\n##### Method -> QR decomposition\n##### Type   -> " <<
		typeid(T).name() << "\n>>> Solving...\n";

	QRDecompositionSolver solver(matrix);

	StaticTimer::start();
	const auto [ solution, Q, R ] = solver.solve();
	const auto elapsed = StaticTimer::elapsed();

	std::cout << ">>> Solved in " << elapsed << " ms\n>>> Solution:\n";
	solution.print();

	std::cout << ">>> Q:\n";
	Q.print();

	std::cout << ">>> R:\n";
	R.print();

	return solution;
}


template<typename T>
void save_solution_to_file(const LinearSystem<T> &system, const std::string &filepath, const std::string &method, const CMatrix<T> &solution) {
	const std::string pathWithoutExtension = filepath.substr(0, filepath.find_last_of("."));
	const std::string extension = filepath.substr(filepath.find_last_of("."));

	std::ofstream outFile(pathWithoutExtension + "[" + method + "]" + "{" + typeid(T).name() + "}" + extension);

	outFile << "Solution:\n";
	for (size_t i = 0; i < solution.rows; ++i) outFile << solution[i][0] << '\n';
	outFile << "\nResidual (cubic norm):\n" << system.residual_cubic(solution);
	outFile.close();
}


template<typename T>
T calculate_ñondition_number(const LinearSystem<T> &system, const CMatrix<T> &matrix, size_t iterations, T disturbanceMagnitude) {
	std::cout << "\n##### Calculating condition number\n##### Systems tested        -> " << iterations
		<< "\n##### Disturbance magnitude -> " << disturbanceMagnitude << "\n##### Precision             -> "
		<< typeid(T).name() << "\n>>> Calculating...\n";

	GaussianEliminationSolver<T> solver(matrix);

	StaticTimer::start();

	const auto b = solver.get_b();
	const auto x = solver.solve();

	const auto inverseNormB = static_cast<T>(1) / b.norm_cubic();
	const auto inverseNormX = static_cast<T>(1) / x.norm_cubic();

	// Solve <iterations> disturbed systems and select max error
	T ñonditionNumber(0);
	for (size_t k = 0; k < iterations; ++k) {
		solver.reset_matrix(matrix);

		solver.add_disturbance(disturbanceMagnitude);
		const auto db = solver.get_b() - b;
		const auto dx = solver.solve() - x;

		const auto delta_b = db.norm_cubic() * inverseNormB;
		const auto delta_x = dx.norm_cubic() * inverseNormX;

		ñonditionNumber = std::max(ñonditionNumber, delta_x / delta_b);
	}

	const auto elapsed = StaticTimer::elapsed();

	std::cout << ">>> Calculated in " << elapsed << " ms\n>>> Condition number:\n" << ñonditionNumber << '\n';

	return ñonditionNumber;
}

template<typename T>
CMatrix<T> calculate_matrix_inverse(const CMatrix<T> &matrix) {
	std::cout << "\n##### Calculating matrix inverse\n##### Precision -> "
		<< typeid(T).name() << "\n>>> Calculating...\n";

	// Ensure we don't get a non-square matrix
	// det != 0 will be ensured by solver
	if (matrix.rows != matrix.cols)
		throw std::runtime_error("ERROR: Cannot inverse non-square matrix");

	const auto N = matrix.rows;
	CMatrix<T> inverseMatrix(N, N);

	// Matrix (N, N + 1) that will hold systems used for calculating inverse
	CMatrix<T> systemMatrix(N, N + 1); 

	for (size_t i = 0; i < N; ++i) {
		memcpy(systemMatrix[i], matrix[i], sizeof(T) * N); // copy N elements at once to the square matrix
	}

	GaussianEliminationSolver<T> solver(systemMatrix);

	StaticTimer::start();

	for (size_t k = 0; k < N; ++k) {
		// Fill last column with 0's and a single 1 is a correct position
		for (size_t i = 0; i < N; ++i)
			systemMatrix[i][N] = (i == k) ? 1 : 0;

		// Solver system
		solver.reset_matrix(systemMatrix);
		const auto solution = solver.solve();

		// Fill a column of the inverse matrix
		for (size_t i = 0; i < N; ++i)
			inverseMatrix[i][k] = solution[i][0];
	}

	const auto elapsed = StaticTimer::elapsed();

	std::cout << ">>> Calculated in " << elapsed << " ms\n>>> Inverse matrix:\n";
	inverseMatrix.print();
	std::cout << ">>> A * A^-1:\n";
	(matrix * inverseMatrix).print();


	return inverseMatrix;
}


int main(int argc, char** argv) {
	#ifdef CHECK_MEMORY
	_CrtMemState _ms;
	_CrtMemCheckpoint(&_ms);
	#endif

	try {
		// Parse config file
		const auto [ inputFilepath, outputFilepath ] = parse_config(); // legal since C++17

		const std::string pathWithoutExtension = outputFilepath.substr(0, outputFilepath.find_last_of("."));
		const std::string extension = outputFilepath.substr(outputFilepath.find_last_of("."));

		// Fill matrices from input file
		std::cout << ">>> Parsing matrix...\n";

		const auto [doubleMatrix, floatMatrix] = parse_matrices(inputFilepath);
		
		// Print parsed matrix, if it's too big supress matrix printing
		std::cout << ">>> Parsed:\n";
		floatMatrix.print(); 

		// Initialize linear systems
		LinearSystem doubleSystem(doubleMatrix);
		LinearSystem floatSystem(floatMatrix);

		// Solve through various methods
		const auto solutionGaussianDouble = solve_throught_gaussian_elimination(doubleMatrix);
		const auto solutionGaussianFloat = solve_throught_gaussian_elimination(floatMatrix);
		const auto solutionQRDouble = solve_through_QR_decomposition(doubleMatrix);
		const auto solutionQRFloat = solve_through_QR_decomposition(floatMatrix);

		// Save solutions to files
		save_solution_to_file(doubleSystem, outputFilepath, "Gaussian_elimination", solutionGaussianDouble);
		save_solution_to_file(floatSystem, outputFilepath, "Gaussian_elimination", solutionGaussianFloat);
		save_solution_to_file(doubleSystem, outputFilepath, "QR_decomposition", solutionQRDouble);
		save_solution_to_file(floatSystem, outputFilepath, "QR_decomposition", solutionQRFloat);

		// Find system condition number
		constexpr size_t ITERATIONS = 1000;
		constexpr double DISTURBANCE_MAGNITUDE = 1e-3;
		const auto conditionNumber = calculate_ñondition_number(doubleSystem, doubleMatrix, ITERATIONS, DISTURBANCE_MAGNITUDE);

		// Save it as file
		constexpr double UNSTABLE_IF_UNDER = 1e2;

		std::ofstream outFile(pathWithoutExtension + "[condition_number]" + extension);
		outFile << "Condition number:\n" << conditionNumber;
		if (conditionNumber < UNSTABLE_IF_UNDER) outFile << "\n\nSystem is stable";
		else outFile << "\n\nSystem is unstable";
		outFile.close();

		// Find inverse matrix and confirm its validity
		const auto N = doubleMatrix.rows;
		CMatrix<double> squareMatrix(N, N);
		for (size_t i = 0; i < N; ++i)
			memcpy(squareMatrix[i], doubleMatrix[i], sizeof(double) * N); // copy N elements at once to the square matrix

		const auto matrixInverse = calculate_matrix_inverse(squareMatrix);

		// Save it as file
		outFile.open(pathWithoutExtension + "[inverse_matrix]" + extension);
		for (size_t i = 0; i < matrixInverse.rows; ++i) {
			for (size_t j = 0; j < matrixInverse.cols; ++j) outFile << std::setw(14) << matrixInverse[i][j];
			outFile << '\n';
		}
		outFile.close();
	}
	// If caught any errors, show error message
	catch (const std::runtime_error& err) {
		std::cerr << err.what() << std::endl;
	}

	#ifdef CHECK_MEMORY
	_CrtMemDumpAllObjectsSince(&_ms);
	#endif

	return 0;
}