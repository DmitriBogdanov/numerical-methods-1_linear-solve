#include <chrono>

#include "static_timer.hpp"

#include "gaussian_elimination_solver.hpp"



int main(int argc, char** argv) {
	// Вычисление в точности double
	try {
		std::cout << "### GAUSSIAN ELIMINATION /// DOUBLE ###\n";

		GaussianEliminationSolver<double> solver("[test]/DATA2.txt");

		std::cout << ">>> Parsed linear system:\n";
		solver.print();

		StaticTimer::start();
		const auto solution = solver.solve_through_gaussian_elimination();
		const auto timeElapsed = StaticTimer::elapsed();

		std::cout << ">>> System solved in " << timeElapsed << " ms:\n";
		solution.print();
	}
	catch (const std::runtime_error& err) {
		std::cerr << err.what() << std::endl;
	}

	// Вычисление в точности float
	try {
		std::cout << "### GAUSSIAN ELIMINATION /// FLOAT ###\n";

		GaussianEliminationSolver<float> solver("[test]/DATA2.txt");

		std::cout << ">>> Parsed linear system:\n";
		solver.print();

		StaticTimer::start();
		solver.solve_through_gaussian_elimination();
		const auto timeElapsed = StaticTimer::elapsed();

		std::cout << ">>> System solved in " << timeElapsed << " ms:\n";
		solver.print();
	}
	catch (const std::runtime_error& err) {
		std::cerr << err.what() << std::endl;
	}

	return 0;
}