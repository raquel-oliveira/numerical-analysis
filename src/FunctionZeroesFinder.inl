#include "FunctionZeroesFinder.h"
#include <cmath>
#include <algorithm>
#include <bitset>
#include <climits>
#include <iostream>
#include <exception>

template<typename TField, typename TLess>
void numerical_analysis::FunctionZeroesFinder<TField, TLess>::lagrange_formula(std::vector<TField> coeff,
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
	std::pair<TField, TField> & interval,
	const TField & p, const int iterations) {

	TField x;
	x = interval.first;

	for(int i = 0 ; i < iterations; ++i){
		if(f(x) * f(x+p) > 0){
			x += p;
			if(x > interval.second) std::logic_error("Step p is not good!\n");
		}else{
			interval.first = x;
			interval.second = x+p;
			return;
		}
	}
	throw std::logic_error("The Signal Change exceeds the maximum iterations\n");
}

template<typename TField, typename TLess>
bool check_interval(std::function<TField (const TField &)> f, const std::pair<TField, TField> & interval){
	TField f_first, f_second;
	f_first = f(interval.first);
	f_second = f(interval.second);
	if (f_first * f_second < 0){
		return true;
	}
	return false;
}

template<typename TField>
bool check_interval(const TField a, const TField b){
	if (a * b < 0){
		return true;
	}
	return false;
}

template<typename TField, typename TLess>
void numerical_analysis::FunctionZeroesFinder<TField, TLess>::bisection(std::function<TField (const TField &)> f,
		const std::pair<TField, TField> & interval,
		int criteria, const double & error,
		TField & root,
		const int iterations) {
			TField f_first, f_second;
			TField l, f_l, u, f_u;

			f_first = f(interval.first);
			f_second = f(interval.second);

			if (f_first * f_second > 0){
				throw std::logic_error("Interval problem in bisection method\n");
			} else {
				l = (f_first < 0) ? interval.first : interval.second;
				f_l = (f_first < 0) ? f_first : f_second;
				u = (f_first < 0) ? interval.second : interval.first;
				f_u = (f_first < 0) ? f_second : f_first;
			}

			int it = 0;
			TField f_m, m;
			TLess comp_less;
			std::bitset<4> bits(criteria);

			std::cout << "\n########################\n";
			std::cout << "### Bisection method ###\n";
			std::cout << "########################\n";

			while(it < iterations){
				//std::cout << "Iteracao " << it <<" de bisection \n";
				it++;
				m = (l+u)/2;
				f_m = f(m);
				std::cout << "Iteration: " << it << " | m = " << m << " | f(m) = " << f_m << std::endl;
				if(check_interval(f_l, f_m)){
					u = m;
					f_u = f_m;
				} else {
					if(f_m == 0){
						root = m;
						return;
					} else {
						l = m;
						f_l = f_m;
					}
				}
				/*std::cout << "A: " <<l << " f= " << f_l << "\n";
				std::cout << "B: " <<u << " f= " << f_u << "\n";
				std::cout << "P: " <<m << " f= " << f_m << "\n";*/
				if(bits[0]){
					//std::cout << "IMAGE("<<it<<"): f(p)=" <<f_m <<" abs(p)="<<std::abs(f_m)<<"   error: "<< error << "\n";
					if(comp_less(std::abs(f_m), error)){
						root = m;
						return;
					}
				}
				if(bits[1]){
					//std::cout << "DELTA IMAGE("<<it<<"):  f(a)=" <<f_l <<"  f(b)="<<f_u<<"  abs(f(a)-f(b))="<<std::abs(f_u-f_l)<<"   error: "<< error << "\n";
					if(comp_less(std::abs(f_u-f_l), error)){
						root = m;
						return;
					}
				}
				if(bits[2]){
					//std::cout << "DELTA DOMAIN("<<it<<"):  A=" <<l <<"  B="<<u<<"  abs(a-b)="<<std::abs(u-l)<<"   error: "<< error << "\n";
					if(comp_less(std::abs(u-l), error)){
						root = m;
						return;
					}
				}
			}
			root = m;
			throw std::logic_error("The Bisection's method exceeds the maximum iterations\n");
}

template<typename TField, typename TLess>
void numerical_analysis::FunctionZeroesFinder<TField, TLess>::regulaFalsi(std::function<TField (const TField &)> f,
		const std::pair<TField, TField> & interval,
		int criteria, const double & error,
		TField & root,
		const int iterations) {
			//Check interval
			TField f_first, f_second;
			f_first = f(interval.first);
			f_second = f(interval.second);
			TField l, f_l, u, f_u;
			if (f_first * f_second > 0){
				throw std::logic_error("Interval problem\n");
			} else{
				// Redirect upper and lower bound
				l = (f_first < 0) ? interval.first : interval.second;
				f_l = (f_first < 0) ? f_first : f_second;
				u = (f_first < 0) ? interval.second : interval.first;
				f_u = (f_first < 0) ? f_second : f_first;
			}

			int it = 0;
			TField f_m, m;
			TLess comp_less;
			std::bitset<4> bits(criteria);

			std::cout << "\n###########################\n";
			std::cout << "### Regula Falsi method ###\n";
			std::cout << "###########################\n";

			while(it < iterations){
				it++;
				m = l - (((u-l)/(f_u - f_l)) * f(l));
				//m = u - (((u-l)/(f_u - f_l)) * f(u));
				f_m = f(m);
				std::cout << "Iteration: " << it << " | m = " << m << " | f(m) = " << f_m << std::endl;
				if(check_interval(f_m,f_l)){
					u = m;
					f_u = f_m;
				} else{
					if(f_m == 0){
						root = m;
						return;
					} else{
						l = m;
						f_l = f_m;
					}
				}

				if(bits[0]){
					//std::cout << "IMAGE("<<it<<"): f(p)=" <<f_m <<" abs(p)="<<std::abs(f_m)<<"   error: "<< error << "\n";
					if(comp_less(std::abs(f_m), error)){
						root = m;
						return;
					}
				}

			 if(bits[1]){
				 //std::cout << "DELTA IMAGE("<<it<<"):  f(a)=" <<f_l <<"  f(b)="<<f_u<<"  abs(f(a)-f(b))="<<std::abs(f_u-f_l)<<"   error: "<< error << "\n";
					if(comp_less(std::abs(f_u-f_l), error)){
						root = m;
						return;
					}
				}

				if(bits[2]){
					//std::cout << "DELTA DOMAIN("<<it<<"):  A=" <<l <<"  B="<<u<<"  abs(a-b)="<<std::abs(u-l)<<"   error: "<< error << "\n";
					if(comp_less(std::abs(u-l), error)){
						root = m;
						return;
					}
				}
			}
			throw std::logic_error("The Regula Falsi's method (método das cordas) exceeds the maximum iterations\n");
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

	std::cout << "\n#####################\n";
	std::cout << "### Newton method ###\n";
	std::cout << "#####################\n";

	for(auto i = 0; i < iterations; ++i){
		p = p0 - f(p0) / df(p0);

		std::cout << "Iteration: " << i+1 << " | p = " << p << std::endl;

		if(bits[0]){
			if(comp_less(std::abs(f(p)), error)){
				root = p;
				return;
			}
		}
		if(bits[1]){
			if(comp_less(std::abs(f(p)-f(p0)), error)){
				root = p;
				return;
			}
		}
		if(bits[2]){
			if(comp_less(std::abs(p-p0), error)){
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

	std::cout << "\n##########################\n";
	std::cout << "### Fixed-point method ###\n";
	std::cout << "##########################\n";

	for(auto i = 0; i < iterations; ++i){
		p = g(p0);
		std::cout << "Iteration: " << i+1 << " | p = " << p << std::endl;
		if(bits[0]){
			if(comp_less(std::abs(g(p)), error)){
				root = p;
				return;
			}
		}
		if(bits[1]){
			if(comp_less(std::abs(g(p)-g(p0)), error)){
				root = p;
				return;
			}
		}
		if(bits[2]){
			if(comp_less(std::abs(p-p0), error)){
				root = p;
				return;
			}
		}

		p0 = p;
	}

	throw std::logic_error("The Fixed-Point's method exceeds the maximum iterations\n");

}

template<typename TField, typename TLess>
void numerical_analysis::FunctionZeroesFinder<TField, TLess>::bissection_newton(std::function<TField (const TField & )> f,
	std::function<TField (const TField &)> df,
	const std::pair<TField, TField> & interval,
	int criteria, const double & error,
	TField & root, const int iterationsB,
	const int iterationsN){

	std::cout << "\n###############################\n";
	std::cout << "### Bisection-newton method ###\n";
	std::cout << "###############################\n";

	try{
		bisection(f, interval, criteria, error, root, iterationsB);
	}catch(std::exception& e){
		try{
			newton(f, df, root, criteria, error, root, iterationsN);
		}catch(std::exception& e){
			throw std::logic_error("The Bisection-Newton's method exceeds the maximum iterations\n");
		}
	}

}
