#include "NonlinearSystemSolver.h"
#include "Matrix.h"

int main(void) {

	// A vector field
	numerical_analysis::Matrix<std::function<double (double, double)>> F = {
		{[](double x, double y){return x + y;}},
		{[](double x, double y){return x*x + y;}},
		{[](double x, double y){return x*y + y;}}
	};

	// It's Jacobian
	numerical_analysis::Matrix<std::function<double (double, double)>> J = {
		{ [](double x, double y){return 1;},	[](double x, double y){return 1;}   },
		{ [](double x, double y){return 2*x;},	[](double x, double y){return 1;}	},
		{ [](double x, double y){return y;},	[](double x, double y){return x+1;} }
	};

	// The root
	numerical_analysis::Matrix<long double> root {3, 1, 0};

	// Solve by Newton's method
	numerical_analysis::NonlinearSystemSolver<double(double,double)>::newton(
		F, J,
		1, 0.001,
		root
	);

	return 0;
}
