#include "FunctionZeroesFinder.h"
#include <cmath>
#include <algorithm>
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

template<typename TField>
void numerical_analysis::FunctionZeroesFinder<TField>::lagrange_restriction(std::vector<TField> coeff, 
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

template<typename TField>
void numerical_analysis::FunctionZeroesFinder<TField>::signal_change_restriction(std::function<TField (const TField &)> f,
		std::pair<TField, TField> & interval) {

}

template<typename TField>
void numerical_analysis::FunctionZeroesFinder<TField>::bisection(std::function<TField (const TField &)> f,
		const std::pair<TField, TField> & interval, 
		int criteria, const double & error,
		TField & root) {

}

template<typename TField>
void numerical_analysis::FunctionZeroesFinder<TField>::secant(std::function<TField (const TField &)> f,
		const std::pair<TField, TField> & interval, 
		int criteria, const double & error,
		TField & root) {

}

template<typename TField>
void numerical_analysis::FunctionZeroesFinder<TField>::newton(std::function<TField (const TField &)> f,
		std::function<TField (const TField &)> df,
		const std::pair<TField, TField> & interval, 
		int criteria, const double & error,
		TField & root) {

}

template<typename TField>
void numerical_analysis::FunctionZeroesFinder<TField>::fixed_point(std::function<TField (const TField &)> g,
		std::function<TField (const TField &)> dg,
		const std::pair<TField, TField> & interval, 
		int criteria, const double & error,
		TField & root) {

}


