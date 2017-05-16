#include "FunctionZeroesFinder.h"
#include <functional>
#include <iostream>
#include <cmath>

int main(void) {

	// Your function
	std::function<double (const double &)> f = [](const double & x){return cos(x) - x;};
	// Your derivative function
	std::function<double (const double &)> df = [](const double & x){return -sin(x) - 1;};

	// The root we want
	double root;

	// Domain interval
	//std::pair<double, double> interval;
	// TODO Apply methods for interval restriction here

	double aproximation = 3.14159265/4;

	// Max iterations
	int max_ite = 10;

	// Bisection Method
	numerical_analysis::FunctionZeroesFinder<double>::newton(
		f,
		df, 
		aproximation, 
		numerical_analysis::FunctionZeroesFinder<double>::StopCriteria::IMAGE | numerical_analysis::FunctionZeroesFinder<double>::StopCriteria::DELTA_IMAGE,
		0.001, 
		root, 
		max_ite
	);	

	std::cout << "Root: " << root << std::endl;

	return 0;
}
