#include "NonlinearSystemSolver.h"
#include "Matrix.h"
#include <chrono>

using namespace numerical_analysis;

int main(void) {

	// F1
	numerical_analysis::Matrix<std::function<long double (const Matrix<long double>&)>> F = {
		{
			[](const Matrix<long double> & p){return 3*p[0][0] - std::cos(p[1][0]*p[2][0]) - 0.5;}
		},
		{
			[](const Matrix<long double> & p){
					return p[0][0]*p[0][0] - 81*std::pow(p[1][0] + 0.1,2) + std::sin(p[2][0]) + 1.06;}
		},
		{
			[](const Matrix<long double> & p){return std::exp(-p[0][0]*p[1][0]) + 20*p[2][0] + (10*M_PI - 3)/3;}
		}
	};

	numerical_analysis::Matrix<std::function<long double (const Matrix<long double>&) >> J = {
		{ [](const Matrix<long double> & p){return 3;},
		  [](const Matrix<long double> & p){return p[2][0]*std::sin(p[1][0]*p[2][0]);},
		  [](const Matrix<long double> & p){return p[1][0]*std::sin(p[1][0]*p[2][0]);}   
		},
		{ [](const Matrix<long double> & p){return 2*p[0][0];},
		  [](const Matrix<long double> & p){return -162*(p[1][0] + 0.1);},
		  [](const Matrix<long double> & p){return std::cos(p[2][0]);}   
		},
		{ [](const Matrix<long double> & p){return -p[1][0]*std::exp(-p[0][0]*p[1][0]);},
		  [](const Matrix<long double> & p){return -p[0][0]*std::exp(-p[0][0]*p[1][0]);},   
		  [](const Matrix<long double> & p){return 20;}   
		}
	};

	// F2
	numerical_analysis::Matrix<std::function<long double (const Matrix<long double> &)>> F2 = {
		{
			[](const Matrix<long double> & p){return p[0][0]*p[0][0] + p[1][0]*p[1][0] - 1;}
		},
		{
			[](const Matrix<long double> & p){return p[0][0] + p[1][0] - 1;}
		}
	};


	numerical_analysis::Matrix<std::function<long double (const Matrix<long double> &)>> J2 = {
		{[](const Matrix<long double> & p){return 2*p[0][0];} , [](const Matrix<long double> & p){return 2*p[1][0] ;} },
		{[](const Matrix<long double> & p){return 1;} , [](const Matrix<long double> & p){return 1;} }
	};


	// The roots
	numerical_analysis::Matrix<long double> rootnewton1 {3, 1, 0};
	numerical_analysis::Matrix<long double> rootbroyden1 {3, 1, 0};

	numerical_analysis::Matrix<long double> rootnewton2 {2, 1, 0};
	numerical_analysis::Matrix<long double> rootbroyden2 {2, 1, 0};

	numerical_analysis::Matrix<long double> initial1 {{0.1}, {0.1}, {-0.1}};
	//numerical_analysis::Matrix<long double> initial1 {{0.5}, {0.5}, {-0.5}};
	//numerical_analysis::Matrix<long double> initial1 {{1}, {1}, {1}};
	//numerical_analysis::Matrix<long double> initial1 {{5}, {5}, {-5}};
	//numerical_analysis::Matrix<long double> initial1 {{50}, {50}, {50}};

	numerical_analysis::Matrix<long double> initial2 {{1},{2}};
	//numerical_analysis::Matrix<long double> initial2 {{1},{1}};
	//numerical_analysis::Matrix<long double> initial2 {{-15000},{-10000}};

	std::cout << "//----- Primeiro sistema -----//" << std::endl;

	std::cout << "Inicial: " << std::endl;
	std::cout << initial1 << std::endl;

	std::cout << "---- Newton results ----" << std::endl;
	auto start = std::chrono::steady_clock::now();
	// Solve by Newton's method
	numerical_analysis::NonlinearSystemSolver<long double>::newton(
		F, J, initial1, rootnewton1,
		1, 0.0001, 10000
	);
	auto end = std::chrono::steady_clock::now();
	std::chrono::duration<double, std::milli> elapsed = (end - start);
	std::cout << "Root: \n" << rootnewton1 << "Time: " << elapsed.count() << "ms" << std::endl << std::endl;
	std::cout << "F(root) = \n" << eval<long double>(F,rootnewton1) << std::endl;

	std::cout << "---- Broyden results ----" << std::endl;
	start = std::chrono::steady_clock::now();
	// Solve by Broyden's method
	numerical_analysis::NonlinearSystemSolver<long double>::broyden(
		F, J, initial1, rootbroyden1,
		1, 0.0001, 10000
	);
	end = std::chrono::steady_clock::now();
	elapsed = (end - start);
	std::cout << "Root: \n" << rootbroyden1 << "Time: " << elapsed.count() << "ms" << std::endl << std::endl;
	std::cout << "F(root) = \n" << eval<long double>(F,rootbroyden1) << std::endl;

	std::cout << "//----- Segundo sistema -----//" << std::endl;

	std::cout << "Inicial: " << std::endl;
	std::cout << initial2 << std::endl;

	std::cout << "---- Newton results ----" << std::endl;
	start = std::chrono::steady_clock::now();
	// Solve by Newton's method
	numerical_analysis::NonlinearSystemSolver<long double>::newton(
		F2, J2, initial2, rootnewton2,
		1, 0.0001, 1000
	);
	end = std::chrono::steady_clock::now();
	elapsed = (end - start);
	std::cout << "Root: \n" << rootnewton2 << "Time: " << elapsed.count() << "ms" << std::endl << std::endl;
	std::cout << "F(root) = \n" << eval<long double>(F2,rootnewton2) << std::endl;

	std::cout << "---- Broyden results ----" << std::endl;
	start = std::chrono::steady_clock::now();
	// Solve by Broyden's method
	numerical_analysis::NonlinearSystemSolver<long double>::broyden(
		F2, J2, initial2, rootbroyden2,
		1, 0.0001, 1000
	);
	end = std::chrono::steady_clock::now();
	elapsed = (end - start);
	std::cout << "Root: \n" << rootbroyden2 << "Time: " << elapsed.count() << "ms" << std::endl << std::endl;
	std::cout << "F(root) = \n" << eval<long double>(F2,rootbroyden2) << std::endl;
	return 0;
}
