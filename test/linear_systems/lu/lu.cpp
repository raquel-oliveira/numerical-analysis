#include "LinearSystemSolver.h"
#include "NaiveLinearSystemSolver.h"
#include "Matrix.h"
#include <iostream>

int main(void) {

	numerical_analysis::Matrix<long double> A = {
        {7.154, 1.075, .638, 1.146, .482},
        {1.075, 7.751, 1.521, .348, 1.616},
        {.638, 1.521, 8.467, 1.257, 1.540},
        {1.146, .348, 1.257, 5.528, 1.639},
        {.482, 1.616, 1.540, 1.639, 7.138}
    };

    numerical_analysis::Matrix<long double> b = {
        {10495}, 
        {12311},
        {13423},
        {9918},
        {12415} 
    };

    std::cout << A;
    std::cout << b;

    numerical_analysis::Matrix<long double> x {b.rows, 1, 0};
    //numerical_analysis::NaiveLinearSystemSolver<double>::solve_by_lu(A, b, x);

	numerical_analysis::LinearSystemSolver<long double>::solve_by_lu(A, b, x, true);

    std::cout << x;

    return 0;
}
