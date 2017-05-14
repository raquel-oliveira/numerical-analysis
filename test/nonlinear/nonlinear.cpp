#include "NonlinearSystemSolver.h"
#include "Matrix.h"

using namespace numerical_analysis;

int main(void) {

	// A vector field
	numerical_analysis::Matrix<std::function<double (const Matrix<double>&)>> F = {
		{
			[](const Matrix<double> & p){return 3*p[0][0] - std::cos(p[1][0]*p[2][0]) - 0.5;}
		},
		{
			[](const Matrix<double> & p){
					return p[0][0]*p[0][0] - 81*std::pow(p[1][0] + 0.1,2) + std::sin(p[2][0]) + 1.06;}
		},
		{
			[](const Matrix<double> & p){return std::exp(-p[0][0]*p[1][0]) + 20*p[2][0] + (10*M_PI - 3)/3;}
		}
	};

	// It's Jacobian
	numerical_analysis::Matrix<std::function<double (const Matrix<double>&) >> J = {
		{ [](const Matrix<double> & p){return 3;},
		  [](const Matrix<double> & p){return p[2][0]*std::sin(p[1][0]*p[2][0]);},
		  [](const Matrix<double> & p){return p[1][0]*std::sin(p[1][0]*p[2][0]);}   
		},
		{ [](const Matrix<double> & p){return 2*p[0][0];},
		  [](const Matrix<double> & p){return -162*(p[1][0] + 0.1);},
		  [](const Matrix<double> & p){return std::cos(p[2][0]);}   
		},
		{ [](const Matrix<double> & p){return -p[1][0]*std::exp(-p[0][0]*p[1][0]);},
		  [](const Matrix<double> & p){return -p[0][0]*std::exp(-p[0][0]*p[1][0]);},   
		  [](const Matrix<double> & p){return 20;}   
		}
	};

	// The root
	numerical_analysis::Matrix<double> root {3, 1, 0};

	numerical_analysis::Matrix<double> initial {{0.1}, {0.1}, {-0.1}};

	// Solve by Newton's method
	numerical_analysis::NonlinearSystemSolver<double>::newton(
		F, J, initial, root,
		1, 0.000000001, 1000
	);

	std::cout << root << std::endl;

	return 0;
}
