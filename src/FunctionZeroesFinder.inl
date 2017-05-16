#include "FunctionZeroesFinder.h"

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

}

template<typename TField>
void numerical_analysis::FunctionZeroesFinder<TField>::fixed_point(std::function<TField (const TField &)> g,
		std::function<TField (const TField &)> dg,
		TField & aproximation, 
		int criteria, const double & error,
		TField & root,
		const int iterations) {

}


