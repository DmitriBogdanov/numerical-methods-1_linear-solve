#include <chrono>

#include "static_timer.hpp"
#include "linear_system.hpp"
#include "gaussian_elimination_solver.hpp"

#ifdef _DEBUG
#define CHECK_MEMORY
#endif

#ifdef CHECK_MEMORY
#include <crtdbg.h>
#define _CRTDBG_MAP_ALLOC
#define new new (_NORMAL_BLOCK, __FILE__, __LINE__)
#endif



// Parse config and return input/output filepaths
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
std::tuple<CMatrix<double>, CMatrix<float>> parse_matrices(const std::string &inputFilepath) {
	CMatrix<double> doubleMatrix;
	CMatrix<float> floatMatrix;

	doubleMatrix.parse_from_file(inputFilepath);
	floatMatrix.parse_from_file(inputFilepath);

	return { doubleMatrix, floatMatrix };
}

// Solve given system through gauss elimination
template<typename T>
CMatrix<T> solve_throught_gaussian_elimination(const CMatrix<T>& matrix) {
	std::cout << "Method -> Gaussian elimination\n" << "Type   -> " <<
		typeid(T).name() << "\n>>> Solving...\n";

	GaussianEliminationSolver solver(matrix);

	StaticTimer::start();
	const auto solution = solver.solve();
	const auto elapsed = StaticTimer::elapsed();

	std::cout << ">>> Solved in " << elapsed << " ms\n\n";

	return solution;
}

/// Solve given system through QR-decomposition

int main(int argc, char** argv) {
	#ifdef CHECK_MEMORY
	_CrtMemState _ms;
	_CrtMemCheckpoint(&_ms);
	#endif

	try {
		// Parse config file
		const auto [ inputFilepath, outputFilepath ] = parse_config(); // legal since C++17

		// Fill matrices from input file
		std::cout << ">>> Parsing matrices...\n";

		const auto [doubleMatrix, floatMatrix] = parse_matrices(inputFilepath);

		// Print parsed matrix, if it's too big supress matrix printing
		constexpr size_t MAX_DISPLAYED_SIZE = 6;
		const bool supressMatrixPrinting = floatMatrix.cols > MAX_DISPLAYED_SIZE;

		if (supressMatrixPrinting) std::cout << "[ Matrix output supressed due to large size ]\n";
		else floatMatrix.print(); 

		// Initialize linear systems
		LinearSystem doubleSystem(doubleMatrix);
		LinearSystem floatSystem(floatMatrix);

		// Solve through various methods
		const auto solutionGaussianDouble = solve_throught_gaussian_elimination(doubleMatrix);
		const auto solutionGaussianFloat = solve_throught_gaussian_elimination(floatMatrix);
		
		// Print solutions to console
		if (supressMatrixPrinting) std::cout << "[ Matrix output supressed due to large size ]\n";
		else {
			std::cout << "Solution 1 // Gaussian elimination // double\n";
			solutionGaussianDouble.print();

			std::cout << "Solution 2 // Gaussian elimination // float\n";
			solutionGaussianFloat.print();
		}

		// Save solutions to files
		std::cout << ">>> Saving solutions to files, calculating residuals...\n";
		const auto dotPos = outputFilepath.find_last_of(".");
		const std::string pathWithoutExtension = outputFilepath.substr(0, dotPos);
		const std::string extension = outputFilepath.substr(dotPos);
		std::ofstream outFile;
		
		outFile.open(pathWithoutExtension + "[GaussElimination]{double}" + extension);
		outFile << "Solution:\n";
		for (size_t i = 0; i < solutionGaussianDouble.rows; ++i) outFile << solutionGaussianDouble[i][0] << '\n';
		outFile << "\nResidual (cubic norm):\n" << doubleSystem.residual_cubic(solutionGaussianDouble);
		outFile.close();

		outFile.open(pathWithoutExtension + "[GaussElimination]{float}" + extension);
		outFile << "Solution:\n";
		for (size_t i = 0; i < solutionGaussianFloat.rows; ++i) outFile << solutionGaussianFloat[i][0] << '\n';
		outFile << "\nResidual (cubic norm):\n" << floatSystem.residual_cubic(solutionGaussianFloat);
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