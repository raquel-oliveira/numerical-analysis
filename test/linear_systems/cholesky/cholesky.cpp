#include "NaiveLinearSystemSolver.h"
#include "LinearSystemSolver.h"
#include "Matrix.h"
#include <iostream>
#include <chrono>

int main(void) {

    numerical_analysis::Matrix<double> A = {
        {4, -1, 1},
		{-1, 4.25, 2.75},
		{1, 2.75, 3.5}
    };

    numerical_analysis::Matrix<double> b = {
        {4}, 
        {6},
		{7.25}
    };

    std::cout << A;
    std::cout << b;

	auto start = std::chrono::steady_clock::now();
    numerical_analysis::Matrix<double> x {b.rows, 1, 0};
    numerical_analysis::NaiveLinearSystemSolver<double>::solve_by_cholesky(A, b, x);
    //numerical_analysis::LinearSystemSolver<double>::solve_by_cholesky(A, b, x);
	auto duration = std::chrono::duration_cast<std::chrono::duration<double>> 
						(std::chrono::steady_clock::now() - start);
    std::cout << x;
	std::cout << "Elapsed time:" << duration.count() << std::endl;

    return 0;
};
