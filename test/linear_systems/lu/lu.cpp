#include "LinearSystemsMethods.h"
#include "Matrix.h"
#include <iostream>

int main(void) {

    numerical_analysis::Matrix<double> A {3, 3, 1};
    numerical_analysis::Matrix<double> L {3, 3, 1};
    numerical_analysis::Matrix<double> U {3, 3, 1};

    
    numerical_analysis::Matrix<double> T = {
        {1, 0, 3},
        {3, 4, 4},
        {1, 2, 3}
    };

    numerical_analysis::Matrix<double> b = {{1}, {2}, {3}};

    std::cout << T;
    std::cout << b;

    numerical_analysis::LinearSystemsMethods<double>::getLU(A, L, U);

    std::cout << A;
    std::cout << L;
    std::cout << U;

    return 0;
}
