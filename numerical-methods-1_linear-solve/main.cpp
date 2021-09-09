#include <chrono>

#include "static_timer.hpp"

#include "gaussian_elimination_solver.hpp"



int main(int argc, char** argv) {
	// Method: Gaussian elimination
	// Precision: double
	try {
		std::cout << "### GAUSSIAN ELIMINATION /// DOUBLE ###\n";
		std::cout << ">>> Parsing...\n";

		GaussianEliminationSolver<double> solver("[test]/DATA2.txt");

		std::cout << ">>> Parsed system:\n";
		solver.print();

		std::cout << "\n>>> Solving...\n";

		StaticTimer::start();
		const auto solution = solver.solve_through_gaussian_elimination();
		const auto timeElapsed = StaticTimer::elapsed();

		std::cout << ">>> Solved in " << timeElapsed << " ms:\n";
		solution.print();
	}
	catch (const std::runtime_error& err) {
		// If caught any errors, show error message
		std::cerr << err.what() << std::endl;
	}

	// Method: Gaussian elimination
	// Precision: float
	try {
		std::cout << "### GAUSSIAN ELIMINATION /// FLOAT ###\n";
		std::cout << ">>> Parsing...\n";

		GaussianEliminationSolver<double> solver("[test]/DATA2.txt");

		std::cout << ">>> Parsed system:\n";
		solver.print();

		std::cout << "\n>>> Solving...\n";

		StaticTimer::start();
		const auto solution = solver.solve_through_gaussian_elimination();
		const auto timeElapsed = StaticTimer::elapsed();

		std::cout << ">>> Solved in " << timeElapsed << " ms:\n";
		solution.print();
	}
	catch (const std::runtime_error& err) {
		// If caught any errors, show error message
		std::cerr << err.what() << std::endl;
	}

	return 0;
}