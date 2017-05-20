#include <complex>
#include "FunctionZeroesFinder.h"
#include <functional>
#include <iostream>
#include "Matrix.h"
#include <cmath>
#include <exception>
#include <iomanip>


using namespace numerical_analysis;
using namespace std::complex_literals;

struct CLess {
	bool operator()(std::complex<double> z, double e){
		return std::norm<double>(z) < e;
	};
};

int main(void) {

	std::function<std::complex<double> (const std::complex<double> &)> f = 
		[](const std::complex<double> & z){return std::pow(z,3) - (1.);};

	std::function<std::complex<double> (const std::complex<double> &)> df = 
		[](const std::complex<double> & z){return 3. * std::pow(z,2);};
	
	std::complex<double> root;
	std::complex<double> initial = 1. + 1i;
	int N0 = 1000;

	try{
		// Newton Method
		numerical_analysis::FunctionZeroesFinder<std::complex<double>, CLess>::newton(
			f,
			df, 
			initial, 
			numerical_analysis::FunctionZeroesFinder<std::complex<double>>::StopCriteria::IMAGE | numerical_analysis::FunctionZeroesFinder<std::complex<double>>::StopCriteria::DELTA_IMAGE,
			0.000001, 
			root, 
			N0 
		);

		std::cout << "Newton's root: " << root << std::endl;

		std::cout << "Check: " << f(root) << std::endl;
	}catch (std::exception& e){
    	std::cout << e.what();
    }

	return 0;
}
