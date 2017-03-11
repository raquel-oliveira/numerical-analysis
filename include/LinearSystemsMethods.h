#ifndef _METHODS_
#define _METHODS_

#include "Matrix.h"

/**
 * Contains various numerical methods for
 * linear systems.
 *
 * @author      Vitor Greati, Vinicius Campos, Raquel Oliveira
 * @date        2017-03-07
 * @version     1.0
 * */
namespace numerical_analysis {

template<typename TField>
class LinearSystemsMethods {

    public:

        /**
         * Performs a LU decomposition with partial
         * pivoting.
         *
         * @param source        Matrix to be decomposed.
         * @param L             Lower triangular matrix after decomposition.
         * @param U             Upper triangular matrix after decomposition.
         * */
        static void getLU(const Matrix<TField> & source,
                   Matrix<TField> & L,
                   Matrix<TField> & U);

        /**
         * Performs backward substitution when receives
         * an upper triangular matrix. Solves the system
         * Ax = b.
         *
         * @param A             Matrix of coefficients.
         * @param b             Vector.
         * @param x             Solution of the linear system.
         * */
        static void backSubstitution(const Matrix<TField> A,
                              const Matrix<TField> b,
                              Matrix<TField> x);

        /**
         * Solve a linear system by LU decomposition only.
         *
         * @param A             Matrix of coefficients.
         * @param b             Vector.
         * @param x             Solution of the linear system.
         * */
        static void solveByLU(const Matrix<TField> A,
                              const Matrix<TField> b,
                              Matrix<TField> x);

        /**
         * Solve a linear system by Cholesky decomposition.
         *
         * @param A             Matrix of coefficients.
         * @param b             Vector.
         * @param x             Solution of the linear system.
         * */
        static void solveByCholesky(const Matrix<TField> A,
                              const Matrix<TField> b,
                              Matrix<TField> x);

};

}

#include "../src/LinearSystemsMethods.inl"

#endif
