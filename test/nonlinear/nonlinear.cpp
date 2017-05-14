#include "NonlinearSystemSolver.h"
#include "Matrix.h"

using namespace numerical_analysis;

int main(void) {

	// A vector field
	numerical_analysis::Matrix<std::function<double (const Matrix<double>&)>> F = {
		{[](const Matrix<double> & p){return p[0][0] + p[1][0];}},
		{[](const Matrix<double> & p){return p[0][0]*p[0][0] + p[1][0];}},
		{[](const Matrix<double> & p){return p[0][0]*p[1][0] + p[1][0];}}
	};

	// It's Jacobian
	numerical_analysis::Matrix<std::function<double (const Matrix<double>&) >> J = {
		{ [](const Matrix<double> & p){return 1;},	[](const Matrix<double> & p){return 1;}   },
		{ [](const Matrix<double> & p){return 2*p[0][0];},	[](const Matrix<double> & p){return 1;}	},
		{ [](const Matrix<double> & p){return p[1][0];},	[](const Matrix<double> & p){return p[0][0]+1;} }
	};

	// The root
	numerical_analysis::Matrix<double> root {3, 1, 0};

	numerical_analysis::Matrix<double> initial {0, 0, 0};
	std::cout << numerical_analysis::eval<double>(J,initial) << std::endl;

	// Solve by Newton's method
	numerical_analysis::NonlinearSystemSolver<double>::newton(
		F, J, initial, root,
		1, 0.001, 1000
	);

	return 0;
}
