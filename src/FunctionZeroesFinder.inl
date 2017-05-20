#include "FunctionZeroesFinder.h"
#include <cmath>
#include <algorithm>
#include <bitset>
#include <climits>
#include <iostream>

template<typename TField>
void numerical_analysis::FunctionZeroesFinder<TField>::lagrange_formula(std::vector<TField> coeff,
		TField & bound) {

	char invsignal = coeff[coeff.size()-1] < 0 ? -1 : 1;

	int k = -1, b = -1;
	for (unsigned int i = 0; i < coeff.size(); ++i) {
		coeff[i] *= invsignal;
		if (coeff[i] < 0)
			k = i;
		if (coeff[i] < 0 && (b == -1 || abs(coeff[b]) < abs(coeff[i]))) 
			b = i;
	}

	if (b == -1 || k == -1) {
		std::cout << "No negative coefficients!" << std::endl;
		bound = std::numeric_limits<TField>::max();
	} else {
		bound = 1 + std::pow(std::abs(coeff[b])/coeff[coeff.size()-1], 1/(coeff.size() - 1.0 - k));
	}
}

template<typename TField, typename TLess>
void numerical_analysis::FunctionZeroesFinder<TField, TLess>::lagrange_restriction(std::vector<TField> coeff, 
	std::tuple<TField, TField, TField, TField> & bounds) {

	TField lower_pos, upper_pos, lower_neg, upper_neg;
	// upper pos
	lagrange_formula(coeff, upper_pos);
	// lower pos
	std::vector<TField> coeff1x {coeff};
	std::reverse(coeff1x.begin(), coeff1x.end());
	lagrange_formula(coeff1x, lower_pos);
	// lower neg
	std::vector<TField> coeffmx {coeff};
	for (unsigned int i = 1; i < coeff.size(); i += 2)
		coeffmx[i] *= -1;
	lagrange_formula(coeffmx, lower_neg);
	// upper neg
	std::vector<TField> coeffmx1x {coeffmx};
	std::reverse(coeffmx1x.begin(), coeffmx1x.end());
	lagrange_formula(coeffmx, upper_neg);

	bounds = std::make_tuple(-lower_neg, -1/upper_neg, 1/lower_pos, upper_pos);
}

template<typename TField, typename TLess>
void numerical_analysis::FunctionZeroesFinder<TField, TLess>::signal_change_restriction(std::function<TField (const TField &)> f,
		std::pair<TField, TField> & interval) {

}

template<typename TField, typename TLess>
void numerical_analysis::FunctionZeroesFinder<TField, TLess>::bisection(std::function<TField (const TField &)> f,
		const std::pair<TField, TField> & interval, 
		int criteria, const double & error,
		TField & root,
		const int iterations) {

}

template<typename TField, typename TLess>
void numerical_analysis::FunctionZeroesFinder<TField, TLess>::secant(std::function<TField (const TField &)> f,
		const std::pair<TField, TField> & interval, 
		int criteria, const double & error,
		TField & root,
		const int iterations) {

}

template<typename TField, typename TLess>
void numerical_analysis::FunctionZeroesFinder<TField, TLess>::newton(std::function<TField (const TField &)> f,
		std::function<TField (const TField &)> df,
		TField & aproximation, 
		int criteria, const double & error,
		TField & root,
		const int iterations) {

	TLess comp_less;
	std::bitset<4> bits(criteria); 

	TField p, p0;
	p0 = aproximation;

	for(auto i = 0; i < iterations; ++i){
		p = p0 - f(p0) / df(p0);
		
		if(bits[0]){
			if(comp_less(f(p), error)){
				root = p;
				return;
			}
		}
		if(bits[1]){
			if(comp_less(abs(f(p)-f(p0)), error)){
				root = p;
				return;
			}
		}
		if(bits[2]){
			if(comp_less(abs(p-p0), error)){
				root = p;
				return;
			}
		}

		p0 = p;
	}

	throw std::logic_error("The Newtons's method exceeds the maximum iterations\n");
}

template<typename TField, typename TLess>
void numerical_analysis::FunctionZeroesFinder<TField, TLess>::fixed_point(std::function<TField (const TField &)> g,
		TField & aproximation, 
		int criteria, const double & error,
		TField & root,
		const int iterations) {

	TLess comp_less;
	std::bitset<4> bits(criteria); 

	TField p, p0;
	p0 = aproximation;

	for(auto i = 0; i < iterations; ++i){
		p = g(p0);
		
		if(bits[0]){
			if(comp_less(g(p), error)){
				root = p;
				return;
			}
		}
		if(bits[1]){
			if(comp_less(abs(g(p)-g(p0)), error)){
				root = p;
				return;
			}
		}
		if(bits[2]){
			if(comp_less(abs(p-p0), error)){
				root = p;
				return;
			}
		}

		p0 = p;
	}

	throw std::logic_error("The Fixed-Point's method exceeds the maximum iterations\n");

}


