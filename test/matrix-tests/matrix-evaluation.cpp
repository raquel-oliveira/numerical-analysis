#include "Matrix.h"
#include <functional>
#include <vector>

using namespace numerical_analysis;

int main(void) {

	Matrix<std::function<double(std::vector<double>)>> J = {
		{[](std::vector<double> p){return p[0]*p[0];},
		[](std::vector<double> p){return p[0]*p[1];}},
		{[](std::vector<double> p){return p[0]*p[1];},
		[](std::vector<double> p){return p[0]*p[1];}},
	};

	Matrix<double> M = {{1, 2}, {2, 1}};

	std::vector<double> v {
		{20, 2}
	};

	std::cout << eval<double>(J,v) << std::endl;

	return 0;
}
