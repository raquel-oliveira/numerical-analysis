#include "NonlinearSystemSolver.h"
#include "Matrix.h"

int main(void) {

	// A vector field
	numerical_analysis::Matrix<std::function<double (const std::vector<double>&)>> F = {
		{[](const std::vector<double> & p){return p[0] + p[1];}},
		{[](const std::vector<double> & p){return p[0]*p[0] + p[1];}},
		{[](const std::vector<double> & p){return p[0]*p[1] + p[1];}}
	};

	// It's Jacobian
	numerical_analysis::Matrix<std::function<double (const std::vector<double>&) >> J = {
		{ [](const std::vector<double> & p){return 1;},	[](const std::vector<double> & p){return 1;}   },
		{ [](const std::vector<double> & p){return 2*p[0];},	[](const std::vector<double> & p){return 1;}	},
		{ [](const std::vector<double> & p){return p[1];},	[](const std::vector<double> & p){return p[0]+1;} }
	};

	// The root
	numerical_analysis::Matrix<double> root {3, 1, 0};
	numerical_analysis::Matrix<double> initial {0, 0, 0};

	// Solve by Newton's method
	numerical_analysis::NonlinearSystemSolver<double>::newton(
		F, J, initial, root,
		1, 0.001, 1000
	);

	return 0;
}
