#ifndef _MATRIX_
#define _MATRIX_

/**
 * Represents an m x n matrix, with its data and operations.
 *
 * It is, actually, a template class, which receives, as template
 * argument, the field of the vector space to which the represented
 * matrix belongs. 
 *
 * Users of this template class must pay attention
 * to this field, since it will only make sense when it is 
 * a subset of real numbers (like int, float, double, long double) or
 * a class representing complex numbers. That's because 
 * this class applies the definition of multiplication
 * for those kinds of fields, and it uses addition and
 * multiplication for numbers. Of course, a generalization
 * for other objects with * and + operations defined will
 * work, but maybe without any meaning.
 *
 * @author      Vitor Greati, Vinicius Campos, Raquel Oliveira
 * @date        2017-03-07
 * @version     1.0
 * */

namespace numerical_analysis {

template<typename TField>
class Matrix {

    private:

        TField * data;          /*< Matrix data. */

    public:
        
        int cols;               /*< Number of columns. */
        int rows;               /*< Number of rows. */

        /**
         * Constructor for an m x n matrix with a defined value
         * in each cell.
         *
         * @param _m            Number of lines.
         * @param _n            Number of columns.
         * @param _initial      Initial value for each cell.
         * */
        Matrix (const int & _m, const int & _n, const TField & _initial);

        /**
         * Constructor for an m x m matrix with an initial value for
         * each cell.
         *
         * @param _m            Number of lines and columns.
         * @param _initial      Initial value for each cell.
         * */
        Matrix (const int & _m, const TField & _initial);

        /**
         * Destructor for deleting the matrix data.
         *
         * */
        ~Matrix();

        /**
         * Operator for matrix addition.
         *
         * @param _rhs  The matrix to be added to this matrix.
         * */
        Matrix & operator+(const Matrix & _rhs);

        /**
         * Operator for matrix addition and assignment.
         *
         * @param _rhs  The matrix to be added to this matrix.
         * */
        Matrix & operator+=(const Matrix & _rhs);

        /**
         * Operator for matrix subtraction.
         *
         * @param _rhs  The matrix to be subtracted from this matrix.
         * */
        Matrix & operator-(const Matrix & _rhs);

        /**
         * Operator for matrix subtraction and assignment.
         *
         * @param _rhs  The matrix to be subtracted from this matrix.
         * */
        Matrix & operator-=(const Matrix & _rhs);

        /**
         * Operator for multiplication of matrices.
         *
         * @param _rhs  The matrix to right-multiply this matrix.
         * */
        Matrix & operator*(const Matrix & _rhs);

        /**
         * Operator for multiplication and assignment.
         *
         * @param _rhs  The matrix to right-multiply this matrix.
         * */
        Matrix & operator*=(const Matrix & _rhs);

        /**
         * Operator for multiplication by scalar in the form scalar * matrix.
         *
         * @param _rhs  The scalar to left-multiply this matrix.
         * */
        friend Matrix & operator*(const TField & _scalar, Matrix & _rhs);

        /**
         * Operator for multiplication by scalar in the form matrix * scalar.
         *
         * @param _rhs  The scalar to right-multiply this matrix.
         * */
        friend Matrix & operator*(Matrix & _rhs, const TField & _scalar);

        /**
         * Operator for multiplication by scalar and assignment.
         *
         * @param _rhs  The scalar to right-multiply this matrix.
         * */
        friend Matrix & operator*=(const TField & _scalar);

};

}
#endif
