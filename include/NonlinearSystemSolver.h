#ifndef _NONLINEAR_SYSTEMS_
#define _NONLINEAR_SYSTEMS_

#include "Matrix.h"
#include <functional>

namespace numerical_analysis {

	/*!
	 * Contains methods for solving nonlinear systems.
	 *
	 * @author		Vitor Greati, Vinicius Campos, Raquel Oliveira
	 * @date		2017-03-23
	 * @version		1.0
	 * */
	template<typename TFunctionSig>
	class NonlinearSystemSolver {

		public:

			/*!
			 * Solves a nonlinear system by Newton's method.
			 * 
			 * @param F				Functions.
			 * @param J				Jacobian of F.
			 * @param criteria		Stop criteria. Can be inserted as disjunction.
			 * @param error			Acceptable error (epsilon).
			 * @param root			The root.
			 * */
			static void newton(const Matrix<std::function<TFunctionSig>> & F,
					const Matrix<std::function<TFunctionSig>> & J,
					int criteria, const double & error,
					Matrix<long double> & root);
		
	};
};

#include "../src/NonlinearSystemSolver.inl"

#endif
