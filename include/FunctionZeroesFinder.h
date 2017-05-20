#ifndef _ZEROES_FINDER_
#define _ZEROES_FINDER_

#include <functional>
#include <vector>

namespace numerical_analysis {

	/*!
	 * Contains methods for finding zeroes of functions.
	 * 
	 * @author Vitor Greati, Vinicius Campos, Raquel Oliveira
	 * @date 2017-03-23
	 * @version 1.0
	 * */
	template<typename TField, typename TLess = std::less<TField>>
	class FunctionZeroesFinder {

		public:

			/*!
			 * Enumeration that represent possible
			 * stop criteria for the methods.
			 * */
			enum StopCriteria {
				IMAGE = 1,				/*!< |f(a)| < err  */
				DELTA_IMAGE = 2,		/*!< |f(b) - f(a)| < err */
				DELTA_DOMAIN = 4		/*!< |b - a| < err */	
			};

			/*!
			 * Gives Lagrange bounds for a polynomial function.
			 * 
			 * @param coeff			List of coefficients a0,...,an
			 * @param bounds		Tuple of bounds, in the form
			 *						[lowerneg upperneg lowerpos upper pos].
			 * */
			static void lagrange_restriction(std::vector<TField> coeff, 
					std::tuple<TField, TField, TField, TField> & bounds);

			/*!
			 * Gives a restriction of the given interval based on the function's
			 * signal change.
			 * 
			 * @param f				The function.
			 * @param interval		Interval of the domain to be restricted.
			 * */
			static void signal_change_restriction(std::function<TField (const TField &)> f,
					std::pair<TField, TField> & interval);
			
			/*!
			 * Computes a root (zero) of the function given some domain
			 * interval, using the Bisection Method.
			 * 
			 * @param f				The function.
			 * @param interval		Domain interval.
			 * @param criteria		Stop criteria. Can be inserted as disjunction.
			 * @param error			Acceptable error (epsilon).
			 * @param root			The root found.
			 * */
			static void bisection(std::function<TField (const TField &)> f,
					const std::pair<TField, TField> & interval, 
					int criteria, const double & error,
					TField & root,
					const int iterations); 
		
			/*!
			 * Computes a root (zero) of the function given some domain
			 * interval, using the Secant Method.
			 * 
			 * @param f				The function.
			 * @param interval		Domain interval.
			 * @param criteria		Stop criteria. Can be inserted as disjunction.
			 * @param error			Acceptable error (epsilon).
			 * @param root			The root found.
			 * */
			static void secant(std::function<TField (const TField &)> f,
					const std::pair<TField, TField> & interval, 
					int criteria, const double & error,
					TField & root,
					const int iterations); 
	
			/*!
			 * Computes a root (zero) of the function given some domain
			 * interval, using the Newton Method.
			 * 
			 * @param f				The function.
			 * @param df			Function's first derivative.
			 * @param interval		Domain interval.
			 * @param criteria		Stop criteria. Can be inserted as disjunction.
			 * @param error			Acceptable error (epsilon).
			 * @param root			The root found.
			 * */
			static void newton(std::function<TField (const TField & )> f,
					std::function<TField (const TField &)> df,
					TField & aproximation, 
					int criteria, const double & error,
					TField & root,
					const int iterations); 

			/*!
			 * Computes a root (zero) of the function given some domain
			 * interval, using the Fixed Point Method.
			 * 
			 * @param g				The function.
			 * @param dg			Function's first derivative.
			 * @param interval		Domain interval.
			 * @param criteria		Stop criteria. Can be inserted as disjunction.
			 * @param error			Acceptable error (epsilon).
			 * @param root			The root found.
			 * */
			static void fixed_point(std::function<TField (const TField &)> g,
					TField & aproximation, 
					int criteria, const double & error,
					TField & root,
					const int iterations); 
	
	};

};

#include "../src/FunctionZeroesFinder.inl"
#endif
