#include "LinearSystemSolver.h"
#include "NaiveLinearSystemSolver.h"
#include "Matrix.h"
#include <iostream>

int main(void) {

    numerical_analysis::Matrix<double> A = {
        {0.0030, 30.0000},
        {1.0000,  4.0000}
    };

    numerical_analysis::Matrix<double> b = {
        {5.0010}, 
        {1.0000}, 
    };

    std::cout << A;
    std::cout << b;

    numerical_analysis::Matrix<double> x {1, b.rows, 0};
    //numerical_analysis::NaiveLinearSystemSolver<double>::solve_by_lu(A, b, x);

	numerical_analysis::LinearSystemSolver<double>::solve_by_lu(A, b, x, false);

    std::cout << x;

    return 0;
}
