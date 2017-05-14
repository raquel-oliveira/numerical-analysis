#include "NonlinearSystemSolver.h"

template<typename TField>
void numerical_analysis::NonlinearSystemSolver<TField>::newton(
		const Matrix<std::function<TField(const std::vector<TField> &)>> & F,
		const Matrix<std::function<TField(const std::vector<TField> &)>> & J,
		Matrix<TField> & initial,
		Matrix<TField> & root,
		int criteria, double error, int iterations) {

	int k = 1;

	while (k <= iterations) {
	
	
		++k;	
	}
}

