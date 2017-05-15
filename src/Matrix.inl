#include "Matrix.h"
#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <cmath>

template<typename TField>
numerical_analysis::Matrix<TField>::Matrix(const int & _m, const int & _n, 
										   const TField & _diag, const TField & _others) : rows {_m}, cols {_n} {

    // Test nullity of number of cols and rows
    if (!_m or !_n) 
        throw std::logic_error("A matrix cannot have zero rows or columns.");

    // Create data and populate it
    data = new TField* [_m];
    for (int i = 0; i < _m; ++i) {
        *(data + i) = new TField[_n];
        for (int j = 0; j < _n; ++j) {
			if (i == j)
				data[i][j] = _diag;
			else
				data[i][j] = _others;
        }
    }
}

template<typename TField>
numerical_analysis::Matrix<TField>::Matrix(const int & _m, const int & _n,
		const TField & _initial) : Matrix(_m, _n, _initial, _initial) {/* empty */}

template<typename TField>
numerical_analysis::Matrix<TField>::Matrix(const Matrix<TField> & from) :  rows {from.rows}, cols {from.cols} {
    this->data = new TField * [from.rows]; 
    for (int i = 0; i < from.rows; ++i)
        this->data[i] = new TField[from.cols];
    for (int i = 0; i < from.rows; ++i) {
        for (int j = 0; j < from.cols; ++j) {
            this->data[i][j] = from[i][j];
        }
    }
}

template<typename TField>
numerical_analysis::Matrix<TField>::Matrix(const std::initializer_list<std::initializer_list<TField>> & l) : data{nullptr} {
    this->rows = l.size();
    
    if (this->rows == 0)
        throw std::invalid_argument("Zero lines not allowed!");
    else {
        this->cols = l.begin()->size();
        if (this->cols == 0)
            throw std::invalid_argument("Zero columns not allowed!");
        else {
            this->data = new TField * [this->rows]; 
            for (int i = 0; i < this->rows; ++i)
                this->data[i] = new TField[this->cols];
        }
    }

    int i = 0; 
    for (auto r : l) {
       int j = 0;
       if ((int) r.size() != this->cols)
           throw std::logic_error("Some element of the matrix is missing!");
       for (auto e : r) {
           data[i][j] = e;
           ++j;
       }
       ++i;
    } 
}

template<typename TField>
numerical_analysis::Matrix<TField>::~Matrix() {
    if (data != nullptr) {
        for (int i = 0; i < this->rows; ++i)
            delete [] this->data[i];
        delete [] this->data;
    }
}

template<typename TField>
void numerical_analysis::Matrix<TField>::set(const int & i, const int & j, const TField & value) {
    data[i][j] = value;
}

template<typename TField>
const TField & numerical_analysis::Matrix<TField>::at(const int & i, const int & j) const {
    return data[i][j];
}

template<typename TField>
void numerical_analysis::Matrix<TField>::swap_lines(const int & i, const int & j) {
    for (int c = 0; c < this->cols; ++c) {
        std::swap(data[i][c], data[j][c]);
    }
}

template<typename TField>
void numerical_analysis::Matrix<TField>::swap_lines(const int & i, const int & j, const int c1, const int c2) {
    for (int c = c1; c <= c2; ++c) {
        std::swap(data[i][c], data[j][c]);
    }
}

template<typename TField>
TField * & numerical_analysis::Matrix<TField>::operator[](const int & i) {
    if (i < 0 || i >= rows)
            throw std::logic_error("Accessing position out of bounds.");
    return data[i];
}

template<typename TField>
TField * numerical_analysis::Matrix<TField>::operator[](const int & i) const {
    if (i < 0 || i >= rows)
            throw std::logic_error("Accessing position out of bounds.");
    return data[i];
}


template<typename TField>
numerical_analysis::Matrix<TField> numerical_analysis::Matrix<TField>::operator+(const Matrix & _rhs) {
    if (rows != _rhs.rows || cols != _rhs.cols)
        throw std::logic_error("Matrices must have the same size!");
    Matrix<TField> sum {rows, cols, 0};
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            sum[i][j] = this->data[i][j] + _rhs[i][j];  
        }
    }
    return sum;
}

template<typename TField>
numerical_analysis::Matrix<TField> & numerical_analysis::Matrix<TField>::operator+=(const Matrix<TField> & _rhs) {
    if (rows != _rhs.rows || cols != _rhs.cols)
        throw std::logic_error("Matrices must have the same size!");
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            this->data[i][j] += _rhs[i][j];  
        }
    }   
    return *this;
}

template<typename TField>
numerical_analysis::Matrix<TField> numerical_analysis::Matrix<TField>::operator-(const Matrix<TField> & _rhs) {
    if (rows != _rhs.rows || cols != _rhs.cols)
        throw std::logic_error("Matrices must have the same size!");
    Matrix<TField> diff {rows, cols, 0};
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            diff[i][j] = this->data[i][j] - _rhs[i][j];  
        }
    }
    return diff;
}

template<typename TField>
numerical_analysis::Matrix<TField> numerical_analysis::Matrix<TField>::operator-() {
	Matrix<TField> opp {*this};
	for (int i = 0; i < opp.rows; ++i)
		for (int j = 0; j < opp.cols; ++j)
			opp[i][j] *= -1;
	return opp;
}

template<typename TField>
numerical_analysis::Matrix<TField> & numerical_analysis::Matrix<TField>::operator-=(const Matrix<TField> & _rhs) {
    if (rows != _rhs.rows || cols != _rhs.cols)
        throw std::logic_error("Matrices must have the same size!");
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            this->data[i][j] -= _rhs[i][j]; 
        }
    }   
    return *this;
}

template<typename TField>
numerical_analysis::Matrix<TField> numerical_analysis::Matrix<TField>::operator*(const Matrix<TField> & _rhs) {
    // Check multiplication condition
    if (cols != _rhs.rows) 
        throw std::logic_error("You must provide matrices mxn and nxp.");
    // Multiply
    Matrix<TField> prod {rows, _rhs.cols, 0};
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < _rhs.cols; ++j) {
            for (int k = 0; k < cols; ++k) {
                prod[i][j] = prod[i][j] + this->data[i][k] * _rhs[k][j]; 
            }   
        }
    }
    return prod;
}

template<typename TField>
numerical_analysis::Matrix<TField> & numerical_analysis::Matrix<TField>::operator=(Matrix<TField> m) {
    if (m.cols != this->cols || m.rows != this->rows) {
		if (data != nullptr) {
			for (int i = 0; i < this->rows; ++i)
				delete [] this->data[i];
			delete [] this->data;
		}
        this->data = new TField * [m.rows]; 
        for (int i = 0; i < m.rows; ++i)
            this->data[i] = new TField[m.cols];
    }
    for (int i = 0; i < m.rows; ++i) {
        for (int j = 0; j < m.cols; ++j) {
            this->data[i][j] = m[i][j];
        }
    }
    this->cols = m.cols;
    this->rows = m.rows;
    return *(this);
}

template<typename TField>
numerical_analysis::Matrix<TField> operator*(const TField & _scalar, numerical_analysis::Matrix<TField> & _rhs) {
    numerical_analysis::Matrix<TField> prod {_rhs.rows, _rhs.cols, 0};
    for (int i = 0; i < _rhs.rows; ++i) {
        for (int j = 0; j < _rhs.cols; ++j) {
            prod[i][j] = _rhs[i][j] * _scalar;
        } 
    }
    return prod;
}

template<typename TField>
numerical_analysis::Matrix<TField> operator*(numerical_analysis::Matrix<TField> & _rhs, const TField & _scalar) {
    return (_scalar)*(_rhs);
}


template<typename TField>
std::ostream& numerical_analysis::operator<<(std::ostream& os, const numerical_analysis::Matrix<TField>& matrix) {
    for (int i = 0; i < matrix.rows; ++i) {
        for (int j = 0; j < matrix.cols; ++j) {
            os << std::setprecision(6) << matrix[i][j] << " ";
        }
        os << std::endl;
    }
    os << std::endl;
    return os;
}

template<typename TField>
numerical_analysis::Matrix<TField> numerical_analysis::eval(
		const Matrix<std::function<TField(const Matrix<TField> &)>> & M, const Matrix<TField> & v) {

	if (v.cols > 1) 
		throw std::runtime_error("eval: You should pass a vector (just one column)."); 	

	Matrix<TField> result {M.rows, M.cols, 0};
	for (int i = 0; i < M.rows; ++i) 
		for (int j = 0; j < M.cols; ++j) 
			result[i][j] = M[i][j](v);
	return result;
}

template<typename TField>
numerical_analysis::Matrix<TField> numerical_analysis::Matrix<TField>::transpose() const{

    numerical_analysis::Matrix<TField> transposed {this->cols, this->rows, 0};
    for(int i = 0; i < rows; ++i){
        for(int j = 0; j < cols; ++j){
            transposed[j][i] = this->data[i][j];
        }
    }

    return transposed;
}

template<typename TField>
bool numerical_analysis::Matrix<TField>::isSymmetric() const{
    if(this->rows != this->cols) {
        return false;
    }else{
        for(int i = 0; i < this->rows; ++i){
            for(int j = i+1; j < this->cols; ++j){
                if(this->data[i][j] != this->data[j][i]) 
                    return false;
            }
        }
    }

    return true;
}

template<typename TField>
numerical_analysis::Matrix<TField> numerical_analysis::Matrix<TField>::diagonal() const{

    numerical_analysis::Matrix<TField> d {this->cols, this->rows, 0};
    for(int i = 0; i < rows; ++i){
        for(int j = 0; j < cols; ++j){
            if (i == j){
                d[i][j] = this->data[i][j];
            }
        }
    }

    return d;
}

template<typename TField>
numerical_analysis::Matrix<TField> numerical_analysis::Matrix<TField>::symmetric() const{

    numerical_analysis::Matrix<TField> s {this->cols, this->rows, 0};
    for(int i = 0; i < rows; ++i){
        for(int j = 0; j < cols; ++j){
                if (this->data[i][j] == 0){
                    s[i][j] = 0; // for some reason -1*0 = -0
                } else {
                    s[i][j] = -1*this->data[i][j];
                }
        }
    }
    return s;
}

template<typename TField>
numerical_analysis::Matrix<TField> numerical_analysis::Matrix<TField>::pow(int k) const{

    numerical_analysis::Matrix<TField> p {this->cols, this->rows, 0};
    for(int i = 0; i < rows; ++i){
        for(int j = 0; j < cols; ++j){
                p[i][j] = std::pow(this->data[i][j], k);
        }
    }
    return p;
}

template<typename TField>
double numerical_analysis::Matrix<TField>::norm_infinity(){
    double aux = 0;
    double max = 0;
    for (int i = 0; i < this->rows; i++){
        for (int j = 0; j < this->cols; j++){
            aux += std::abs(this->data[i][j]);
        }
        if (max < aux) max = aux;
        aux = 0;
    }
    return max;
}

template<typename TField>
double numerical_analysis::Matrix<TField>::norm_one(){
    double aux = 0;
    double max = 0;
    for (int i = 0; i < this->cols; i++){
        for (int j = 0; j < this->rows; j++){
            aux += std::abs(this->data[i][j]);
        }
        if (max < aux) max = aux;
        aux = 0;
    }
    return max;
}

template<typename TField>
bool numerical_analysis::Matrix<TField>::rowCriteria(){
    double sum = 0;
    for (int i = 0; i < this->rows; i++){
        for (int j = 0; j < this->cols; j++){
            if (i != j){
                sum += std::abs(this->data[i][j]);
            }
        }
        if (std::abs(this->data[i][i]) < sum){
            return false;
        }
        sum = 0;
    }
    return true;
}
