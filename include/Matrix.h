#ifndef _MATRIX_
#define _MATRIX_

#include <iostream>

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

        TField ** data;          /*< Matrix data. */

    public:
        
        int rows;               /*< Number of rows. */
        int cols;               /*< Number of columns. */

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
         * Constructor for an m x m identity matrix.
         *
         * @param _m            Dimension of the identity matrix.
         * */
        Matrix(const int & _m);

        /**
         * Constructor which takes an initializer list.
         *
         * @param l             Initializer list with matrix elements.
         * */
        Matrix(const std::initializer_list<std::initializer_list<TField>> & l);

        /**
         * Copy constructor.
         *
         * */
        Matrix(const Matrix<TField> & from);

        /**
         * Destructor for deleting the matrix data.
         *
         * */
        ~Matrix();

        /**
         * Access function: takes the element data[i][j].
         *
         * @param i     Element row.
         * @param j     Element column.
         * */
        const TField & at(const int & i, const int & j) const; 

        /**
         * Set function: set data[i][j] to a value.
         *
         * @param i     Element row.
         * @param j     Element column.
         * @param value New value.
         * */
        void set(const int & i, const int & j, const TField & value);

        /**
         * Swap two lines of the matrix.
         *
         * @param i     One line.
         * @param j     Another line.
         * */
        void swap_lines(const int & i, const int & j);

        /**
        * Operator [] for accessing rows of a matrix. This
        * returns a reference.
        *
        * Since arrays have [] access defined, 
        * this overload allows using [][] for accessing and
        * modifying matrix elements.
        *
        * @param j		Row index.
        * @return 		Pointer to the first element of the row.
        **/
        TField * & operator[](const int & i);

        /**
        * Operator [] for accessing rows of a matrix. This
        * returns a copy.
        *
        * Since arrays have [] access defined, 
        * this overload allows using [][] for accessing
        * matrix elements.
        *
        * @param j		Row index.
        * @return 		Pointer to the first element of the row.
        **/
        TField * operator[](const int & i) const;

        /**
         * Operator for matrix addition.
         *
         * @param _rhs  The matrix to be added to this matrix.
         * */
        Matrix<TField> operator+(const Matrix<TField> & _rhs);

        /**
         * Operator for matrix addition and assignment.
         *
         * @param _rhs  The matrix to be added to this matrix.
         * */
        Matrix<TField> & operator+=(const Matrix<TField> & _rhs);

        /**
         * Operator for matrix subtraction.
         *
         * @param _rhs  The matrix to be subtracted from this matrix.
         * */
        Matrix<TField> operator-(const Matrix<TField> & _rhs);

        /**
         * Operator for matrix subtraction and assignment.
         *
         * @param _rhs  The matrix to be subtracted from this matrix.
         * */
        Matrix<TField> & operator-=(const Matrix<TField> & _rhs);

        /**
         * Operator for multiplication of matrices.
         *
         * @param _rhs  The matrix to right-multiply this matrix.
         * */
        Matrix<TField> operator*(const Matrix<TField> & _rhs);

        /**
         * Operator for assignment.
         *
         * @param m     The current object will be equal to m.
         * */
        Matrix<TField> & operator=(Matrix<TField> m);

        /**
         * Operator for multiplication by scalar in the form scalar * matrix.
         *
         * @param _rhs  The scalar to left-multiply this matrix.
         * */
        template<typename TFielda>
        friend Matrix<TField> operator*(const TField & _scalar, Matrix<TField> & _rhs);

        /**
         * Operator for multiplication by scalar in the form matrix * scalar.
         *
         * @param _rhs  The scalar to right-multiply this matrix.
         * */
        template<typename TFielda>
        friend Matrix<TField> operator*(Matrix<TField> & _rhs, const TField & _scalar);

        /**
         * Method to transpose a matrix and return it.
         *
         * @return Matrix<TField> Matrix transposed.  
         * */
        Matrix<TField> transpose() const;


         /**
         * Method to get diagonal of a matrix and return it.
         *
         * @return Matrix<TField> Matrix diagonal.  
         * */
        Matrix<TField> diagonal() const;

};

	/**
	 * Allows printing the matrix by stream.
	 *
	 * @param os            Output stream.
	 * @param matrix        Matrix to be printed.
	 * */
	template<typename TField>
	std::ostream& operator<<(std::ostream& os, const numerical_analysis::Matrix<TField>& matrix);

}

#include "../src/Matrix.inl"
#endif
