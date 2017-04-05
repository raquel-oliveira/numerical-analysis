#ifndef _NAIVE_METHODS_
#define _NAIVE_METHODS_

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
	        static double norm_euclidean(Matrix<TField> m);

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
             * @param p             Number of precision
             * @param x             Solution of the linear system.

             * @return              number of iterations
             * */
            static long solve_by_jacobi(Matrix<TField> A,
                          Matrix<TField> b,
                          double p,
                          Matrix<TField> & x);

            /**
             * 
             * Solve a linear system by Gauss-Seidel.
             * Solution update in @param x.
             *
             * @param A             Matrix of coefficients.
             * @param b             Vector.
             * @param p             Number of precision
             * @param x             Initial guess of solution of the linear system
             * 
             * @return              number of iterations
             * */
            static long solve_by_seidel(Matrix<TField> A,
                          Matrix<TField> b,
                          double p,
                          Matrix<TField> & xe);
	};

};

#include "../src/NaiveLinearSystemSolver.inl"

#endif
