#ifndef _MATRIX_FACTORIZATION_
#define _MATRIX_FACTORIZATION_

#include "Matrix.h"
#include <cmath>

namespace numerical_analysis {

	/*!
	 * Contains methods for factoring matrices.
	 * 
	 * @author      Vitor Greati, Vinicius Campos, Raquel Oliveira
	 * @date        2017-03-23
	 * @version     1.0
	 * */
	template <typename TField>
	class MatrixDecomposer {

		/*!
		 * Performs LU factorization with partial pivoting.
		 * 
		 * Method based in the one presented in the book 'Numerical Analysis', Burden et al.
		 * Basically, returns two matrices, L and U, with L being
		 * lower triangular and U, upper triangular, such that
		 * source = L*U. Besides the basic factorization, performs also
		 * a partial pivoting to avoid round-off problems.
		 * 
		 * @param source		Matrix to be factored.
		 * @param L				Lower triangular matrix.
		 * @param U				Upper triangular matrix.
		 * */
		static void lu (const Matrix<TField> & source, 
				Matrix<TField> & L,
				Matrix<TField> & U);

		/*!
		 * Performs LDLt factorization.
		 * 
		 * Method based in the one presented in the book 'Numerical Analysis',
		 * Burden et al. It computes the matrices L and D, such that source = L*D*Lt.
		 * 
		 * @param source		Matrix to be factored.
		 * @param L				Lower triangular matrix.
		 * @param D				Diagonal matrix.
		 * */
		static void ldlt (const Matrix<TField> & source,
				Matrix<TField> & L,
				Matrix<TField> & D);

		/*!
		 * Performs Cholesky (LLt) factorization.
		 * 
		 * Method based in the one presented in the book 'Numerical Analysis',
		 * Burden et al. It computes the matrix L, such that source = L*Lt.
		 * 
		 * @param source		Matrix to be factored.
		 * @param L				Lower triangular matrix such that source = L*Lt.
		 * */
		static void cholesky (const Matrix<TField> & source,
				Matrix<TField> & L);
	};

};
#include "../src/MatrixDecomposer.inl"
#endif
