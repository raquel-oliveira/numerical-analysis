#include "Matrix.h"
#include <iostream>

int main(void) {
    
    numerical_analysis::Matrix<double> m1 {6, 6, 1.0}; 
    std::cout << m1; 
	std::cout << -m1;

    numerical_analysis::Matrix<double> m2 {6, 6, 1.0};
    std::cout << m2;

    numerical_analysis::Matrix<double> m3 = m1 * m2;
    std::cout << m3;

    m2 += m1;
    m2 += m1;
    std::cout << m2;
    std::cout << m1;

    m2 -= m1;
    std::cout << m2;

    numerical_analysis::Matrix<double> m4 {1, 2, 3};
    numerical_analysis::Matrix<double> m5 {2, 1, 3};
    numerical_analysis::Matrix<double> m6 = m4 * m5;
    std::cout << m6;

    numerical_analysis::Matrix<double> m7 = 0.5 * m1;
    numerical_analysis::Matrix<double> m8 = m1 * 0.5;

    std::cout << m7;
    std::cout << m8;

    m7 = 0.5 * m4;

    std::cout << m7;

    numerical_analysis::Matrix<double> m9 {5, 1, 0};
    std::cout << m9;

    numerical_analysis::Matrix<double> m10 {{1,2,3},{4,5,6},{7,8,9}};
    std::cout << m10;
    std::cout << m10.isSymmetric() << std::endl;

    // Testing [] operator
	std::cout << m9[1][0] << std::endl;
	m9[4][0] = 4;
	std::cout << m9[1][0] << std::endl;
	
    return 0;
}
