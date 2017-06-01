#include "FunctionZeroesFinder.h"
#include <functional>
#include <iostream>
#include "Matrix.h"
#include <cmath>
#include <exception>
#include <vector>
#include <string>

#define PI 3.14159265
using namespace std;

// Polynomial function
// Your function
std::function<double (const double &)> f = [](const double & x){return pow(x,8) - 18*pow(x,6) + 105*pow(x,4) - 232*pow(x,2) + 144;};
// Your derivative function
std::function<double (const double &)> df = [](const double & x){return 8*pow(x,7) - 108*pow(x,5) + 420*pow(x,3) - 464*x;};
// Function coefficients
std::vector<double> coeff{144, 0, -232, 0, 105, 0, -18, 0, 1};
// Positive and negative intervals to be obtained with lagrange
std::tuple<double,double,double,double> interval;	
// Positive interval obtained with signal change after lagrange
std::pair<double,double> intervalP;
// Negative interval obtained with signal change after lagrange
std::pair<double,double> intervalN;

// Non-Polynomial function
std::function<double (const double &)> h = [](const double & x){return cos(x) - x;};
// Your derivative function
std::function<double (const double &)> dh = [](const double & x){return -1 * sin(x) - 1;};
// Interval of non-polynomal function obtained after signal-change
std::pair<double,double> intervalNonPoly;
// The root we want
double root;
// Aproximation
double aproximation = 1.5;
// Max iterations
int max_ite = 200;

void polinomial_lagrange(){
	try{
		// Lagrange method
		numerical_analysis::FunctionZeroesFinder<double>::lagrange_restriction(coeff, interval);
		std::cout << "Lagrange Interval: P[" << std::get<2>(interval) << "," << std::get<3>(interval) << "] | ";
		std::cout << "N[" << std::get<0>(interval) << "," << std::get<1>(interval) << "]\n";
	}catch(std::exception& e){
		std::cout << e.what();
	}
}
void polinomial_signal_change(){
	try{
		// Signal change method
		intervalN = std::make_pair(std::get<0>(interval), std::get<1>(interval));
		intervalP = std::make_pair(std::get<2>(interval), std::get<3>(interval));
		numerical_analysis::FunctionZeroesFinder<double>::signal_change_restriction(f, intervalP, 1, 100);
		numerical_analysis::FunctionZeroesFinder<double>::signal_change_restriction(f, intervalN, 1, 100);
		std::cout << "Signal-Change Interval: P[" << intervalP.first << "," << intervalP.second << "] | ";
		std::cout << "N[" << intervalN.first << "," << intervalN.second << "]\n";
	}catch(std::exception& e){
		std::cout << e.what();
	}
}
void polinomial_newton_positive(){
	try{
		aproximation = (intervalP.first + intervalP.second) / 2;
		// Newton Method 
		numerical_analysis::FunctionZeroesFinder<double>::newton(
			f,
			df,
			aproximation,
			numerical_analysis::FunctionZeroesFinder<double>::StopCriteria::IMAGE | 
			numerical_analysis::FunctionZeroesFinder<double>::StopCriteria::DELTA_IMAGE |
			numerical_analysis::FunctionZeroesFinder<double>::StopCriteria::DELTA_DOMAIN,
			0.0001,
			root,
			max_ite
		);

		std::cout << "Root's aproximation: " << root << std::endl;
	}catch (std::exception& e){
    	std::cout << e.what();
    }
}
void polinomial_fixed_positive(){
	try{
    	// Your derivative function
		std::function<double (const double &)> g = [](const double & x){return x - ( 2*(pow(x,8) - 18*pow(x,6) + 105*pow(x,4) - 232*pow(x,2) + 144)*(8*pow(x,7) - 108*pow(x,5) + 420*pow(x,3) - 464*x) ) / 
																			   ( 2*pow(8*pow(x,7) - 108*pow(x,5) + 420*pow(x,3) - 464*x,2) - (pow(x,8) - 18*pow(x,6) + 105*pow(x,4) - 232*pow(x,2) + 144)*(56*pow(x,6)-540*pow(x,4)+1260*pow(x,2)-464) );};
		aproximation = (intervalP.first + intervalP.second) / 2;
		// Fixed_point Method
		numerical_analysis::FunctionZeroesFinder<double>::fixed_point(
			g,
			aproximation,
			numerical_analysis::FunctionZeroesFinder<double>::StopCriteria::IMAGE |
			numerical_analysis::FunctionZeroesFinder<double>::StopCriteria::DELTA_IMAGE |
			numerical_analysis::FunctionZeroesFinder<double>::StopCriteria::DELTA_DOMAIN ,
			0.0001,
			root,
			max_ite
		);

		std::cout << "Root's aproximation: " << root << std::endl;
	}catch (std::exception& e){
    	std::cout << e.what();
    }
}
void polinomial_bisection_positive(){
	try{
		// Bisection Method
		numerical_analysis::FunctionZeroesFinder<double>::bisection(
			f,
			intervalP,
			numerical_analysis::FunctionZeroesFinder<double>::StopCriteria::IMAGE | 
			numerical_analysis::FunctionZeroesFinder<double>::StopCriteria::DELTA_IMAGE |
			numerical_analysis::FunctionZeroesFinder<double>::StopCriteria::DELTA_DOMAIN,
			0.0001,
			root,
			max_ite
		);

		std::cout << "Root's aproximation: " << root << std::endl;
	}catch (std::exception& e){
			std::cout << e.what();
	}
}
void polinomial_regular_positive(){
	try{
		// Regula Falsi Method
		numerical_analysis::FunctionZeroesFinder<double>::regulaFalsi(
			f,
			intervalP,
			numerical_analysis::FunctionZeroesFinder<double>::StopCriteria::IMAGE | 
			numerical_analysis::FunctionZeroesFinder<double>::StopCriteria::DELTA_IMAGE |
			numerical_analysis::FunctionZeroesFinder<double>::StopCriteria::DELTA_DOMAIN,
			0.0001,
			root,
			max_ite
		);

		std::cout << "Root's aproximation: " << root << std::endl;
	}catch (std::exception& e){
			std::cout << e.what();
	}
}
void polinomial_hibrid_positive(){
	try{
		// Bisection-Newton Method
		numerical_analysis::FunctionZeroesFinder<double>::bissection_newton(
			f,
			df,
			intervalP,
			numerical_analysis::FunctionZeroesFinder<double>::StopCriteria::IMAGE | 
			numerical_analysis::FunctionZeroesFinder<double>::StopCriteria::DELTA_IMAGE |
			numerical_analysis::FunctionZeroesFinder<double>::StopCriteria::DELTA_DOMAIN,
			0.0001,
			root,
			7,
			max_ite
		);

		std::cout << "Root's aproximation: " << root << std::endl;
	}catch (std::exception& e){
			std::cout << e.what();
	}
}

void polinomial_newton_negative(){
	try{
		aproximation = (intervalN.first + intervalN.second) / 2;
		// Newton Method 
		numerical_analysis::FunctionZeroesFinder<double>::newton(
			f,
			df,
			aproximation,
			numerical_analysis::FunctionZeroesFinder<double>::StopCriteria::IMAGE | 
			numerical_analysis::FunctionZeroesFinder<double>::StopCriteria::DELTA_IMAGE |
			numerical_analysis::FunctionZeroesFinder<double>::StopCriteria::DELTA_DOMAIN,
			0.0001,
			root,
			max_ite
		);

		std::cout << "Root's aproximation: " << root << std::endl;
	}catch (std::exception& e){
    	std::cout << e.what();
    }
}
void polinomial_fixed_negative(){
	try{
    	// Your derivative function
		std::function<double (const double &)> g = [](const double & x){return x - ( 2*(pow(x,8) - 18*pow(x,6) + 105*pow(x,4) - 232*pow(x,2) + 144)*(8*pow(x,7) - 108*pow(x,5) + 420*pow(x,3) - 464*x) ) / 
																			   ( 2*pow(8*pow(x,7) - 108*pow(x,5) + 420*pow(x,3) - 464*x,2) - (pow(x,8) - 18*pow(x,6) + 105*pow(x,4) - 232*pow(x,2) + 144)*(56*pow(x,6)-540*pow(x,4)+1260*pow(x,2)-464) );};
		aproximation = (intervalN.first + intervalN.second) / 2;
		// Fixed_point Method
		numerical_analysis::FunctionZeroesFinder<double>::fixed_point(
			g,
			aproximation,
			numerical_analysis::FunctionZeroesFinder<double>::StopCriteria::IMAGE |
			numerical_analysis::FunctionZeroesFinder<double>::StopCriteria::DELTA_IMAGE |
			numerical_analysis::FunctionZeroesFinder<double>::StopCriteria::DELTA_DOMAIN ,
			0.0001,
			root,
			max_ite
		);

		std::cout << "Root's aproximation: " << root << std::endl;
	}catch (std::exception& e){
    	std::cout << e.what();
    }
}
void polinomial_bisection_negative(){
	try{
		// Bisection Method
		numerical_analysis::FunctionZeroesFinder<double>::bisection(
			f,
			intervalN,
			numerical_analysis::FunctionZeroesFinder<double>::StopCriteria::IMAGE | 
			numerical_analysis::FunctionZeroesFinder<double>::StopCriteria::DELTA_IMAGE |
			numerical_analysis::FunctionZeroesFinder<double>::StopCriteria::DELTA_DOMAIN,
			0.0001,
			root,
			max_ite
		);

		std::cout << "Root's aproximation: " << root << std::endl;
	}catch (std::exception& e){
			std::cout << e.what();
	}
}
void polinomial_regular_negative(){
	try{
		// Regula Falsi Method
		numerical_analysis::FunctionZeroesFinder<double>::regulaFalsi(
			f,
			intervalN,
			numerical_analysis::FunctionZeroesFinder<double>::StopCriteria::IMAGE | 
			numerical_analysis::FunctionZeroesFinder<double>::StopCriteria::DELTA_IMAGE |
			numerical_analysis::FunctionZeroesFinder<double>::StopCriteria::DELTA_DOMAIN,
			0.0001,
			root,
			max_ite
		);

		std::cout << "Root's aproximation: " << root << std::endl;
	}catch (std::exception& e){
			std::cout << e.what();
	}
}
void polinomial_hibrid_negative(){
	try{
		// Bisection-Newton Method
		numerical_analysis::FunctionZeroesFinder<double>::bissection_newton(
			f,
			df,
			intervalN,
			numerical_analysis::FunctionZeroesFinder<double>::StopCriteria::IMAGE | 
			numerical_analysis::FunctionZeroesFinder<double>::StopCriteria::DELTA_IMAGE |
			numerical_analysis::FunctionZeroesFinder<double>::StopCriteria::DELTA_DOMAIN,
			0.0001,
			root,
			7,
			max_ite
		);

		std::cout << "Root's aproximation: " << root << std::endl;
	}catch (std::exception& e){
			std::cout << e.what();
	}
}

void nonpolinomial_signal_change(){
	try{
		// Signal change method
		intervalNonPoly = std::make_pair(-100,100);
		numerical_analysis::FunctionZeroesFinder<double>::signal_change_restriction(h, intervalNonPoly, 2, 1000);
		std::cout << "Signal-Change Interval: P[" << intervalNonPoly.first << "," << intervalNonPoly.second << "] \n";
	}catch(std::exception& e){
		std::cout << e.what();
	}
}
void nonpolinomial_newton(){
	try{
		aproximation = (intervalNonPoly.first + intervalNonPoly.second) / 2;
		// Newton Method 
		numerical_analysis::FunctionZeroesFinder<double>::newton(
			h,
			dh,
			aproximation,
			numerical_analysis::FunctionZeroesFinder<double>::StopCriteria::IMAGE | 
			numerical_analysis::FunctionZeroesFinder<double>::StopCriteria::DELTA_IMAGE |
			numerical_analysis::FunctionZeroesFinder<double>::StopCriteria::DELTA_DOMAIN,
			0.0001,
			root,
			max_ite
		);

		std::cout << "Root's aproximation: " << root << std::endl;
	}catch (std::exception& e){
    	std::cout << e.what();
    }
}
void nonpolinomial_fixed(){
	try{
    	// Your derivative function
		std::function<double (const double &)> g = [](const double & x){return x - ( 2 * (cos(x)-x) * (-sin(x)-1)) / 
																			   (2 * pow(-sin(x)-1,2) - (cos(x)-x) * (-cos(x)));};
		aproximation = (intervalNonPoly.first + intervalNonPoly.second) / 2;
		// Fixed_point Method
		numerical_analysis::FunctionZeroesFinder<double>::fixed_point(
			g,
			aproximation,
			numerical_analysis::FunctionZeroesFinder<double>::StopCriteria::IMAGE |
			numerical_analysis::FunctionZeroesFinder<double>::StopCriteria::DELTA_IMAGE |
			numerical_analysis::FunctionZeroesFinder<double>::StopCriteria::DELTA_DOMAIN ,
			0.0001,
			root,
			max_ite
		);

		std::cout << "Root's aproximation: " << root << std::endl;
	}catch (std::exception& e){
    	std::cout << e.what();
    }
}
void nonpolinomial_bisection(){
	try{
		// Bisection Method
		numerical_analysis::FunctionZeroesFinder<double>::bisection(
			h,
			intervalNonPoly,
			numerical_analysis::FunctionZeroesFinder<double>::StopCriteria::IMAGE | 
			numerical_analysis::FunctionZeroesFinder<double>::StopCriteria::DELTA_IMAGE |
			numerical_analysis::FunctionZeroesFinder<double>::StopCriteria::DELTA_DOMAIN,
			0.0001,
			root,
			max_ite
		);

		std::cout << "Root's aproximation: " << root << std::endl;
	}catch (std::exception& e){
			std::cout << e.what();
	}
}
void nonpolinomial_regular(){
	try{
		// Regula Falsi Method
		numerical_analysis::FunctionZeroesFinder<double>::regulaFalsi(
			h,
			intervalNonPoly,
			numerical_analysis::FunctionZeroesFinder<double>::StopCriteria::IMAGE | 
			numerical_analysis::FunctionZeroesFinder<double>::StopCriteria::DELTA_IMAGE |
			numerical_analysis::FunctionZeroesFinder<double>::StopCriteria::DELTA_DOMAIN,
			0.0001,
			root,
			max_ite
		);

		std::cout << "Root's aproximation: " << root << std::endl;
	}catch (std::exception& e){
			std::cout << e.what();
	}
}
void nonpolinomial_hibrid(){
	try{
		// Bisection-Newton Method
		numerical_analysis::FunctionZeroesFinder<double>::bissection_newton(
			h,
			dh,
			intervalNonPoly,
			numerical_analysis::FunctionZeroesFinder<double>::StopCriteria::IMAGE | 
			numerical_analysis::FunctionZeroesFinder<double>::StopCriteria::DELTA_IMAGE |
			numerical_analysis::FunctionZeroesFinder<double>::StopCriteria::DELTA_DOMAIN,
			0.0001,
			root,
			7,
			max_ite
		);

		std::cout << "Root's aproximation: " << root << std::endl;
	}catch (std::exception& e){
			std::cout << e.what();
	}
}

int main(int argn, char ** argv) {

	if(argn > 1){
		string arg = argv[1];
		if(arg[0] == 'a'){
			std::cout << "------------ Polynomial tests ------------\n";
			std::cout << "f(x) = x^8 - 18x^6 + 105x^4 - 232x^2 + 144\n";
			std::cout << "------------------------------------------\n\n";
			// Polynomial
			polinomial_lagrange();
			polinomial_signal_change();
			std::cout << "-----> Positive interval <-----\n\n";	
			polinomial_bisection_positive();
			polinomial_regular_positive();
			polinomial_newton_positive();
			polinomial_fixed_positive();
			polinomial_hibrid_positive();
			std::cout << "-----> Negative interval <-----\n\n";
			polinomial_bisection_negative();
			polinomial_regular_negative();
			polinomial_newton_negative();
			polinomial_fixed_negative();
			polinomial_hibrid_negative();
			std::cout << "\n-- Non-Polynomial tests --\n";
			std::cout << "    f(x) = cos(x) - x\n";
			std::cout << "---------------------------\n\n";
			nonpolinomial_signal_change();
			nonpolinomial_bisection();
			nonpolinomial_regular();
			nonpolinomial_newton();
			nonpolinomial_fixed();
			nonpolinomial_hibrid();
		}else if(arg[0] == 'p'){
			std::cout << "------------ Polynomial tests ------------\n";
			std::cout << "f(x) = x^8 - 18x^6 + 105x^4 - 232x^2 + 144\n";
			std::cout << "------------------------------------------\n\n";
			// Polynomial
			polinomial_lagrange();
			polinomial_signal_change();
			if(argn > 2) arg = argv[2];
			else arg = "a";
			if(arg[0] == 'p'){
				std::cout << "-----> Positive interval <-----\n\n";
				polinomial_bisection_positive();
				polinomial_regular_positive();
				polinomial_newton_positive();
				polinomial_fixed_positive();
				polinomial_hibrid_positive();
			}else if(arg[0] == 'n'){
				std::cout << "-----> Negative interval <-----\n\n";
				polinomial_bisection_negative();
				polinomial_regular_negative();
				polinomial_newton_negative();
				polinomial_fixed_negative();
				polinomial_hibrid_negative();
			}else{
				std::cout << "-----> Positive interval <-----\n\n";	
				polinomial_bisection_positive();
				polinomial_regular_positive();
				polinomial_newton_positive();
				polinomial_fixed_positive();
				polinomial_hibrid_positive();
				std::cout << "-----> Negative interval <-----\n\n";
				polinomial_bisection_negative();
				polinomial_regular_negative();
				polinomial_newton_negative();
				polinomial_fixed_negative();
				polinomial_hibrid_negative();
			}
			
		}else if(arg[0] == 'n'){
			std::cout << "-- Non-Polynomial tests --\n";
			std::cout << "    f(x) = cos(x) - x\n";
			std::cout << "---------------------------\n\n";
			nonpolinomial_signal_change();
			nonpolinomial_bisection();
			nonpolinomial_regular();
			nonpolinomial_newton();
			nonpolinomial_fixed();
			nonpolinomial_hibrid();
		}
	}	

	return 0;
}
