#include "NonlinearSystemSolver.h"
#include "LinearSystemSolver.h"

template<typename TField>
void numerical_analysis::NonlinearSystemSolver<TField>::newton(
		const Matrix<std::function<TField(const Matrix<TField> &)>> & F,
		const Matrix<std::function<TField(const Matrix<TField> &)>> & J,
		Matrix<TField> & initial,
		Matrix<TField> & root,
		int criteria, double error, int iterations) {

	int k = 1;

	Matrix<TField> x {initial};
	while (k <= iterations) {
		// Calculate F(x) and J(x)
		Matrix<TField> Fx = eval<TField>(F,x);
		Matrix<TField> Jx = eval<TField>(J,x);
		// Solve the linear system J(x)y = -F(x)
		for (int i = 0; i < Fx.rows; ++i)
			for (int j = 0; j < Fx.cols; ++j)
				Fx[i][j] *= -1;
		Matrix<TField> y {F.rows, 1, 0};
		LinearSystemSolver<TField>::solve_by_lu(Jx, Fx, y);
		// x = x + y
		x += y;
		// Check tolerance
		if (y.norm_infinity() <= error) {
			root = x;
			break;
		}	
		++k;	
	}

	if (k > iterations)
		std::cout << "newton: Max iterations reached, no convergence!" << std::endl;
}

