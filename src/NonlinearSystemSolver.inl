#include "NonlinearSystemSolver.h"

template<typename TFunctionSig>
void numerical_analysis::NonlinearSystemSolver<TFunctionSig>::newton(
				const Matrix<std::function<TFunctionSig>> & F,
				const Matrix<std::function<TFunctionSig>> & J,
				int criteria, const double & error,
				Matrix<long double> & root) {

}

