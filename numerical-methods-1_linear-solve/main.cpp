#include <chrono>
#include <tuple>

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
// @return 3 => condition number estimate iterations
std::tuple<std::string, std::string, size_t> parse_config() {
	// Parse config file
	std::string inputFilepath;
	std::string outputFilepath;
	size_t iterations;

	std::ifstream inConfig("config.txt");

	if (!inConfig.is_open())
		throw std::runtime_error("ERROR: Could not oper config file.");

	std::getline(inConfig, inputFilepath);
	std::getline(inConfig, outputFilepath);
	inConfig >> iterations;

	return { inputFilepath , outputFilepath, iterations };
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


// @return 1 => solution
// @return 2 => residual (cubic)
// @return 3 => residual (coctahedral)
template<typename T>
std::tuple<CMatrix<T>, T, T> solve_throught_gaussian_elimination(const CMatrix<T>& matrix, const LinearSystem<T> &system) {
	std::cout << "\n##### Method -> Gaussian elimination\n##### Type   -> " <<
		typeid(T).name() << "\n>>> Solving...\n";

	GaussianEliminationSolver solver(matrix);

	StaticTimer::start();
	const auto solution = solver.solve();
	const auto elapsed = StaticTimer::elapsed();

	std::cout << ">>> Solved in " << elapsed << " ms\n>>> Solution:\n";
	solution.print();

	const auto residualCubic = system.residual_cubic(solution);
	const auto residualOctahedral = system.residual_octahedral(solution);
	
	std::cout << ">>> Residual (cubic norm):\n" << residualCubic << "\n>>> Residual (octahedral norm):\n" << residualOctahedral << '\n';

	return { solution, residualCubic, residualOctahedral };
}


// @return 1 => solution
// @return 2 => residual (cubic)
// @return 3 => residual (coctahedral)
template<typename T>
std::tuple<CMatrix<T>, T, T> solve_through_QR_decomposition(const CMatrix<T>& matrix, const LinearSystem<T> &system) {
	std::cout << "\n##### Method -> QR decomposition\n##### Type   -> " <<
		typeid(T).name() << "\n>>> Solving...\n";

	QRDecompositionSolver solver(matrix);

	StaticTimer::start();
	const auto solution = solver.solve();
	const auto [ Q, R ] = solver.get_QR();
	const auto elapsed = StaticTimer::elapsed();

	std::cout << ">>> Solved in " << elapsed << " ms\n>>> Solution:\n";
	solution.print();

	std::cout << ">>> Q:\n";
	Q.print();

	std::cout << ">>> R:\n";
	R.print();

	const auto residualCubic = system.residual_cubic(solution);
	const auto residualOctahedral = system.residual_octahedral(solution);

	std::cout << ">>> Residual (cubic norm):\n" << residualCubic << "\n>>> Residual (octahedral norm):\n" << residualOctahedral << '\n';

	return { solution, residualCubic, residualOctahedral };
}


template<typename T>
void save_solution_to_file(const std::string &filepath, const std::string &method, const CMatrix<T> &solution, T residualCubic, T residualOctahedral) {
	const std::string pathWithoutExtension = filepath.substr(0, filepath.find_last_of("."));
	const std::string extension = filepath.substr(filepath.find_last_of("."));

	std::ofstream outFile(pathWithoutExtension + "[" + method + "]" + "{" + typeid(T).name() + "}" + extension);

	outFile << "Solution:\n";
	for (size_t i = 0; i < solution.rows; ++i) outFile << solution[i][0] << '\n';
	outFile << "\nResidual (cubic norm):\n" << residualCubic;
	outFile << "\n\nResidual (octahedral norm):\n" << residualOctahedral;
	outFile.close();
}


template<typename T>
CMatrix<T> calculate_matrix_inverse(const CMatrix<T>& matrix) {
	std::cout << "\n##### Calculating matrix inverse\n##### Precision -> "
		<< typeid(T).name() << "\n>>> Calculating...\n";

	// Ensure we don't get a non-square matrix
	// det != 0 will be ensured by solver
	if (matrix.rows != matrix.cols)
		throw std::runtime_error("ERROR: Cannot inverse non-square matrix");

	const auto N = matrix.rows;
	CMatrix<T> inverseMatrix(N, N);

	// Right side
	CMatrix<T> column(N, 1);

	QRDecompositionSolver<T> solver(matrix, column);

	StaticTimer::start();

	for (size_t k = 0; k < N; ++k) {
		// Fill last column with 0's and a single 1 is a correct position
		column.fill_with_zeroes();
		column[k][0] = static_cast<T>(1);

		// Solve system
		solver.reset_right_side(column);
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


// @return 1 => ñondition number estimate (cubic norm)
// @return 2 => ñondition number estimate (octahedral norm)
template<typename T>
std::tuple<T, T> estimate_ñondition_number(const LinearSystem<T> &system, const CMatrix<T> &matrix, size_t iterations, T disturbanceMagnitude) {
	std::cout << "\n##### Condition number estimate\n##### Systems tested        -> " << iterations
		<< "\n##### Disturbance magnitude -> " << disturbanceMagnitude << "\n##### Precision             -> "
		<< typeid(T).name() << "\n>>> Calculating...\n";

	QRDecompositionSolver<T> solver(matrix);

	StaticTimer::start();

	const CMatrix<T> b = solver.get_b();
	const CMatrix<T> x = solver.solve();

	const auto inverseNormB = static_cast<T>(1) / b.norm_cubic();
	const auto inverseNormX = static_cast<T>(1) / x.norm_cubic();

	// Solve <iterations> disturbed systems and select max error
	T ñonditionNumber_cubic(0);
	T ñonditionNumber_octahedral(0);
	for (size_t k = 0; k < iterations; ++k) {
		solver.reset_right_side(b);

		solver.add_disturbance(disturbanceMagnitude);
		const auto db = solver.get_b() - b;
		const auto dx = solver.solve() - x;

		const auto delta_b_cubic = db.norm_cubic() * inverseNormB;
		const auto delta_x_cubic = dx.norm_cubic() * inverseNormX;

		const auto delta_b_octahedral = db.norm_octahedral() * inverseNormB;
		const auto delta_x_octahedral = dx.norm_octahedral() * inverseNormX;

		ñonditionNumber_cubic = std::max(ñonditionNumber_cubic, delta_x_cubic / delta_b_cubic);
		ñonditionNumber_octahedral = std::max(ñonditionNumber_octahedral, delta_x_octahedral / delta_b_octahedral);
	}

	const auto elapsed = StaticTimer::elapsed();

	std::cout << ">>> Calculated in " << elapsed << " ms\n>>> Condition number estimate (cubic norm):\n" << ñonditionNumber_cubic
		<< "\n>>> Condition number estimate (octahedral norm) :\n"  << ñonditionNumber_octahedral  << '\n';

	return { ñonditionNumber_cubic, ñonditionNumber_octahedral };
}


// @return 1 => ñondition number (cubic norm)
// @return 2 => ñondition number (octahedral norm)
template<typename T>
std::tuple<T, T> calculate_ñondition_number(const CMatrix<T>& matrix, const CMatrix<T>& matrixInverse) {
	std::cout << "\n##### Condition number\n##### Precision -> "
		<< typeid(T).name() << "\n>>> Calculating...\n";

	StaticTimer::start();

	const auto cubic = matrix.norm_cubic() * matrixInverse.norm_cubic();
	const auto octahedral = matrix.norm_octahedral() * matrixInverse.norm_octahedral();

	const auto elapsed = StaticTimer::elapsed();

	std::cout << ">>> Calculated in " << elapsed << " ms\n>>> Condition number (cubic norm):\n" << cubic
		<< "\n>>> Condition number (octahedral norm):\n" << octahedral << '\n';

	return { cubic, octahedral };
}


int main(int argc, char** argv) {
	#ifdef CHECK_MEMORY // only checks in DEBUG, not defined in release
	_CrtMemState _ms;
	_CrtMemCheckpoint(&_ms);
	#endif
	 
	try {
		// Parse config file
		const auto [ inputFilepath, outputFilepath, iterations] = parse_config(); // legal since C++17

		const std::string pathWithoutExtension = outputFilepath.substr(0, outputFilepath.find_last_of("."));
		const std::string extension = outputFilepath.substr(outputFilepath.find_last_of("."));

		// Fill matrices from input file
		std::cout << ">>> Parsing matrix...\n";

		const auto [ doubleMatrix, floatMatrix ] = parse_matrices(inputFilepath);
		
		// Print parsed matrix, if it's too big supress matrix printing
		std::cout << ">>> Size:\n(" << floatMatrix.rows << ", " << floatMatrix.cols << ")\n>>> Parsed:\n";
		floatMatrix.print(); 

		// Initialize linear systems
		LinearSystem doubleSystem(doubleMatrix);
		LinearSystem floatSystem(floatMatrix);

		// Solve through various methods
		const auto [ solutionGD, resGD_cubic, resGD_oct ] = solve_throught_gaussian_elimination(doubleMatrix, doubleSystem);
		const auto [ solutionGF, resGF_cubic, resGF_oct ] = solve_throught_gaussian_elimination(floatMatrix, floatSystem);
		const auto [ solutionQRD, resQRD_cubic, resQRD_oct ] = solve_through_QR_decomposition(doubleMatrix, doubleSystem);
		const auto [ solutionQRF, resQRF_cubic, resQRF_oct ] = solve_through_QR_decomposition(floatMatrix, floatSystem);

		// Save solutions to files
		save_solution_to_file(outputFilepath, "Gaussian_elimination", solutionGD, resGD_cubic, resGD_oct);
		save_solution_to_file(outputFilepath, "Gaussian_elimination", solutionGF, resGF_cubic, resGF_oct);
		save_solution_to_file(outputFilepath, "QR_decomposition", solutionQRD, resQRD_cubic, resQRD_oct);
		save_solution_to_file(outputFilepath, "QR_decomposition", solutionQRF, resQRF_cubic, resQRF_oct);

		// Find inverse matrix
		const auto N = doubleMatrix.rows;
		CMatrix<double> squareMatrix(N, N);
		for (size_t i = 0; i < N; ++i)
			memcpy(squareMatrix[i], doubleMatrix[i], sizeof(double) * N); // copy N elements at once to the square matrix

		const auto matrixInverse = calculate_matrix_inverse(squareMatrix);

		// Save it as file
		std::ofstream outFile(pathWithoutExtension + "[inverse_matrix]" + extension);
		for (size_t i = 0; i < matrixInverse.rows; ++i) {
			for (size_t j = 0; j < matrixInverse.cols; ++j) outFile << std::setw(14) << matrixInverse[i][j];
			outFile << '\n';
		}
		outFile.close();

		// Estimate condition number
		constexpr double DISTURBANCE_MAGNITUDE = 1e-2;
		const auto [ estimateCubic, estimateOctahedral ] = estimate_ñondition_number(doubleSystem, doubleMatrix, iterations, DISTURBANCE_MAGNITUDE);

		// Find exact condition number
		const auto [exactCubic, exactOctahedral] = calculate_ñondition_number(squareMatrix, matrixInverse);

		// Save it to file
		constexpr double STABLE_IF_UNDER = 1e5;

		outFile.open(pathWithoutExtension + "[condition_number]" + extension);

		outFile << "Condition number estimate (cubic norm):\n" << estimateCubic
			<< "\n\nCondition number estimate (octahedral norm):\n" << estimateOctahedral
			<< "\n\nCondition number (cubic norm):\n" << exactCubic
			<< "\n\nCondition number  (octahedral norm):\n" << exactOctahedral;

		if (estimateCubic < STABLE_IF_UNDER) outFile << "\n\nSystem is stable";
		else outFile << "\n\nSystem is unstable";

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