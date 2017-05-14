#include "Matrix.h"
#include <functional>
#include <vector>

using namespace numerical_analysis;

int main(void) {

	Matrix<std::function<double(const Matrix<double> &)>> J = {
		{[](const Matrix<double> p){return p[0][0]*p[0][0];},
		[](const Matrix<double> p){return p[0][0]*p[1][0];}},
		{[](const Matrix<double> p){return p[0][0]*p[1][0];},
		[](const Matrix<double> p){return p[0][0]*p[1][0];}},
	};

	Matrix<double> M = {{1, 2}, {2, 1}};

	Matrix<double> v {
		{20},{2}
	};

	std::cout << eval<double>(J,v) << std::endl;

	return 0;
}
