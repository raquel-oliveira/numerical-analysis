#ifndef _NAIVE_METHODS_
#define _NAIVE_METHODS_

#include "Matrix.h"
#include "LinearSystemSolver.h"

namespace numerical_analysis {
	/**
	 * Contains various numerical methods for
	 * linear systems.
	 *
	 * @author      Vitor Greati, Vinicius Campos, Raquel Oliveira
	 * @date        2017-03-07
	 * @version     1.0
	 * */
	template<typename TField>
	class NaiveLinearSystemSolver : public LinearSystemSolver<TField> {

		public:

			/**
			 * Performs a LU decomposition with partial
			 * pivoting.
			 *
			 * @param source        Matrix to be decomposed.
			 * @param b             In case of solving a linear system, it's the b on Ax=b.
			 * @param Linv          Inverse of lower triangular matrix after decomposition.
			 * @param U             Upper triangular matrix after decomposition.
			 * */
			static void get_linvu(const Matrix<TField> & source,
					   Matrix<TField> & b,
					   Matrix<TField> & Linv,
					   Matrix<TField> & U);

			/**
			 * Solve a linear system by LU decomposition only.
			 *
			 * @param A             Matrix of coefficients.
			 * @param b             Vector.
			 * @param x             Solution of the linear system.
			 * */
			static void solve_by_lu(const Matrix<TField> A,
								  Matrix<TField> b,
								  Matrix<TField> & x);

			/**
			 * Solve a linear system by Cholesky decomposition.
			 *
			 * @param A             Matrix of coefficients.
			 * @param b             Vector.
			 * @param x             Solution of the linear system.
			 * */
			static void solve_by_cholesky(const Matrix<TField> A,
								  Matrix<TField> b,
								  Matrix<TField> &x);

	};

};

#include "../src/NaiveLinearSystemSolver.inl"

#endif
