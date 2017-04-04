#include "NaiveLinearSystemSolver.h"
#include "LinearSystemSolver.h"
#include "Matrix.h"
#include <iostream>

int main(void) {

    numerical_analysis::Matrix<double> A = {
        {2, 1},
        {-1,  3}
    };

    numerical_analysis::Matrix<double> b = {
        {3}, 
        {2}, 
    };

    numerical_analysis::Matrix<double> x = { {3}, {3}};

    double e = 0.01;

    std::cout << "A:" << std::endl;
    std::cout << A;
    std::cout << "b:" << std::endl; 
    std::cout << b;
    std::cout << "correctness: " << e << std::endl;
    std::cout << "start x:" << std::endl; 
    std::cout << x;

    numerical_analysis::NaiveLinearSystemSolver<double>::solve_by_jacobi(A, b, e, x);

    std::cout << "Answer" << std::endl;
    std::cout << x;

    return 0;
}
