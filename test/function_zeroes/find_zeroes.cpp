#include "FunctionZeroesFinder.h"
#include <functional>
#include <iostream>
#include "Matrix.h"
#include <cmath>
#include <exception>

#define PI 3.14159265

int main(void) {

	// Your function
	std::function<double (const double &)> f = [](const double & x){return x*x*x + 4*(x*x) - 10;};
	// Your derivative function
	std::function<double (const double &)> df = [](const double & x){return 3*(x*x) + 8*x;};

	// The root we want
	double root;

	// Domain interval
	//std::pair<double, double> interval;
	// TODO Apply methods for interval restriction here

	double aproximation = 1.5;

	// Max iterations
	int max_ite = 100;

	try{
		// Newton Method
		numerical_analysis::FunctionZeroesFinder<double>::newton(
			f,
			df,
			aproximation,
			numerical_analysis::FunctionZeroesFinder<double>::StopCriteria::IMAGE | numerical_analysis::FunctionZeroesFinder<double>::StopCriteria::DELTA_IMAGE,
			0.0001,
			root,
			max_ite
		);

		std::cout << "Newton's root: " << root << std::endl;
	}catch (std::exception& e){
    	std::cout << e.what();
    }

    try{
    	// Your derivative function
		std::function<double (const double &)> g = [](const double & x){return x - ((x*x*x + 4*(x*x) - 10)/(3*(x*x) + 8*x));};

		// Fixed_point Method
		numerical_analysis::FunctionZeroesFinder<double>::fixed_point(
			g,
			aproximation,
			numerical_analysis::FunctionZeroesFinder<double>::StopCriteria::IMAGE |
			numerical_analysis::FunctionZeroesFinder<double>::StopCriteria::DELTA_IMAGE ,
			0.0001,
			root,
			max_ite
		);

		std::cout << "Fixed point's root: " << root << std::endl;
	}catch (std::exception& e){
    	std::cout << e.what();
    }

		try{
			// Bisection Method
			numerical_analysis::FunctionZeroesFinder<double>::bisection(
				f,
				std::make_pair(-1, 1.5),
				numerical_analysis::FunctionZeroesFinder<double>::StopCriteria::IMAGE,
				0.0001,
				root,
				max_ite
			);

			std::cout << "Bisection's root: " << root << std::endl;
		}catch (std::exception& e){
				std::cout << e.what();
			}
	return 0;
}
