#include "FunctionZeroesFinder.h"
#include <functional>
#include <iostream>

int main(void) {

	// Your function
	std::function<double (const double &)> f = [](const double & x){return x*x - 2;};

	// The root we want
	double root;

	// Domain interval
	std::pair<double, double> interval;
	// TODO Apply methods for interval restriction here

	// Max iterations
	int max_ite = 10;

	// Bisection Method
	numerical_analysis::FunctionZeroesFinder<double>::bisection(
		f, 
		interval, 
		numerical_analysis::FunctionZeroesFinder<double>::StopCriteria::IMAGE | numerical_analysis::FunctionZeroesFinder<double>::StopCriteria::DELTA_IMAGE,
		0.001, root, max_ite
	);	

	return 0;
}
