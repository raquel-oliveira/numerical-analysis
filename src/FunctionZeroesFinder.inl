#include "FunctionZeroesFinder.h"
#include <bitset>
#include <climits>
#include <iostream>

template<typename TField>
void numerical_analysis::FunctionZeroesFinder<TField>::lagrange_restriction(std::vector<TField> coeff, 
	std::tuple<TField, TField, TField, TField> & bounds) {

}

template<typename TField>
void numerical_analysis::FunctionZeroesFinder<TField>::signal_change_restriction(std::function<TField (const TField &)> f,
		std::pair<TField, TField> & interval) {

}

template<typename TField>
void numerical_analysis::FunctionZeroesFinder<TField>::bisection(std::function<TField (const TField &)> f,
		const std::pair<TField, TField> & interval, 
		int criteria, const double & error,
		TField & root,
		const int iterations) {

}

template<typename TField>
void numerical_analysis::FunctionZeroesFinder<TField>::secant(std::function<TField (const TField &)> f,
		const std::pair<TField, TField> & interval, 
		int criteria, const double & error,
		TField & root,
		const int iterations) {

}

template<typename TField>
void numerical_analysis::FunctionZeroesFinder<TField>::newton(std::function<TField (const TField &)> f,
		std::function<TField (const TField &)> df,
		TField & aproximation, 
		int criteria, const double & error,
		TField & root,
		const int iterations) {

	std::bitset<4> bits(criteria); 

	TField p, p0;
	p0 = aproximation;

	for(auto i = 0; i < iterations; ++i){
		p = p0 - f(p0) / df(p0);
		
		if(bits[0]){
			if(f(p) < error){
				root = p;
				return;
			}
		}
		if(bits[1]){
			if(abs(f(p)-f(p0)) < error){
				root = p;
				return;
			}
		}
		if(bits[2]){
			if(abs(p-p0) < error){
				root = p;
				return;
			}
		}

		p0 = p;
	}

	throw std::logic_error("The Newtons's method exceeds the maximum iterations\n");
}

template<typename TField>
void numerical_analysis::FunctionZeroesFinder<TField>::fixed_point(std::function<TField (const TField &)> g,
		TField & aproximation, 
		int criteria, const double & error,
		TField & root,
		const int iterations) {

	std::bitset<4> bits(criteria); 

	TField p, p0;
	p0 = aproximation;

	for(auto i = 0; i < iterations; ++i){
		p = g(p0);
		
		if(bits[0]){
			if(g(p) < error){
				root = p;
				return;
			}
		}
		if(bits[1]){
			if(abs(g(p)-g(p0)) < error){
				root = p;
				return;
			}
		}
		if(bits[2]){
			if(abs(p-p0) < error){
				root = p;
				return;
			}
		}

		p0 = p;
	}

	throw std::logic_error("The Fixed-Point's method exceeds the maximum iterations\n");

}


