#ifndef _METHODS_
#define _METHODS_

#include "Matrix.h"
#include "LinearSystemSolver.h"
#include <stdbool.h>

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

		private:
			/**
			* @return the euclidean norm on a n-dimensional euclidean space R^n 
			*/
			static double get_norm(Matrix<TField> c);

			/*
			* @return if it's worth find the answer using Jacobi method
			*/
			static bool check_jacobi(Matrix<TField> m);

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

			
			/**
             * Solve a linear system by Jacobi.
             *
             * @param A             Matrix of coefficients.
             * @param b             Vector.
             * @param c             Number of approximation to correctness
             * @param x             Solution of the linear system.
             * */
            static void solve_by_jacobi(Matrix<TField> A,
                          Matrix<TField> b,
                          double c,
                          Matrix<TField> & x);
	};

}

#include "../src/NaiveLinearSystemSolver.inl"

#endif
