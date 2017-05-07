#include "FunctionZeroesFinder.h"
#include <functional>
#include <iostream>
#include "Matrix.h"

int main(void) {

	// Your function
	std::function<double (const double &)> f = [](const double & x){return x*x - 2;};

	// The root we want
	double root;

	// Domain interval
	std::pair<double, double> interval;
	// TODO Apply methods for interval restriction here

	// Bisection Method
	numerical_analysis::FunctionZeroesFinder<double>::bisection(
		f, 
		interval, 
		numerical_analysis::FunctionZeroesFinder<double>::StopCriteria::IMAGE | numerical_analysis::FunctionZeroesFinder<double>::StopCriteria::DELTA_IMAGE,
		0.001, root
	);	

	//numerical_analysis::Matrix<double (double, double, double)> m {5, [](double x, double y, double z){return x;}};
	numerical_analysis::Matrix<std::function<double (double, double, double)>> m {5, 5, [](double x, double y, double z){return 1;},
	[](double x, double y, double z){return 0;}};

	std::cout << m << std::endl;

	return 0;
}
