#include "Matrix.h"

template<typename TField>
numerical_analysis::Matrix (const int & _m, const int & _n, const TField & _initial) : rows {_m}, cols {_n} {

    // Test nullity of number of cols and rows
    if (!_m or !_n) 
        throw std::logic_error("A matrix cannot have zero rows or columns.");

    // Create data and populate it
    data = new TField * [_m];
    for (int i = 0; i < _m; ++i) {
        *(data + i) = new TField[_n];
        for (int j = 0; j < _n; ++j) {
            data[i][j] = _initial;
        }
    }
}

template<typename TField>
numerical_analysis::Matrix (const int & _m, const TField & _initial) : Matrix(m, m, _initial) { }

template<typename TField>
numerical_analysis::Matrix::~Matrix() {
    for (int i = 0; i < _m; ++i)
        delete [] (data + i);
}

void numerical_analysis::Matrix::set(const int & i, const int & j, const int & value) {
    data[i][j] = value;
}

template<typename TField>
const TField & numerical_analysis::Matrix<TField>::at(const int & i, const int & j) const {
    return data[i][j];
}

template<typename TField>
numerical_analysis::Matrix<TField> numerical_analysis::Matrix<TField>::operator+(const Matrix & _rhs) {
    Matrix<TField> sum {rows, cols, 0};
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            sum.set(i, j, this->data[i][j] + _rhs.get(i, j));  
        }
    }
    return sum;
}

template<typename TField>
numerical_analysis::Matrix<TField> & numerical_analysis::Matrix<TField>::operator+=(const Matrix<TField> & _rhs) {
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            this->data[i][j] += _rhs.get(i, j);  
        }
    }   
    return *this;
}

template<typename TField>
numerical_analysis::Matrix<TField> numerical_analysis::Matrix<TField>::operator-(const Matrix<TField> & _rhs) {
    Matrix<TField> diff {rows, cols, 0};
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            diff.set(i, j, this->data[i][j] - _rhs.get(i, j));  
        }
    }
    return diff;
}

template<typename TField>
numerical_analysis::Matrix<TField> & numerical_analysis::Matrix<TField>::operator-=(const Matrix<TField> & _rhs) {
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            this->data[i][j] -= _rhs.get(i, j);  
        }
    }   
    return *this;
}

template<typename TField>
numerical_analysis::Matrix<TField> numerical_analysis::Matrix<TField>::operator*(const Matrix<TField> & _rhs) {
    // Check multiplication condition
    if (rows != _rhs.cols) throw std::logic_error("You must provide matrices mxn and nxp.");

    // Multiply
    Matrix<TField> prod {rows, _rhs.cols, 0};
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            for (int k = 0; k < cols; ++k) {
                prod[i][j] += this->data[i][k] * _rhs.get(k, j); 
            }   
        }
    }
    return prod;
}

template<typename TField>
numerical_analysis::Matrix<TField> numerical_analysis::operator*(const TField & _scalar, numerical_analysis::Matrix<TField> & _rhs) {
    Matrix<TField> prod {rows, cols, 0};
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            prod.set(i, j, data[i][j] * _scalar);
        } 
    }
    return prod;
}

template<typename TField>
numerical_analysis::Matrix<TField> numerical_analysis::operator*(numerical_analysis::Matrix<TField> & _rhs, const TField & _scalar) {
    return (_scalar)*(_rhs);
}

template<typename TField>
numerical_analysis::Matrix<TField> & numerical_analysis::operator*=(const TField & _scalar) {
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            data[i][j] =  data[i][j] * _scalar;
        } 
    }
    return *this;   
}
