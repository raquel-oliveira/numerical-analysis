#include "LinearSystemsMethods.h"

template<typename TField>
void numerical_analysis::LinearSystemsMethods<TField>::getLinvU(const Matrix<TField> & source, 
        Matrix<TField> & b,
        Matrix<TField> & Linv,
        Matrix<TField> & U) {

    numerical_analysis::Matrix<double> UTemp {source};
    numerical_analysis::Matrix<double> LinvTemp {source.rows}; 

    // Calculate G0...Gi...Gn-2
    for (int i = 0; i < source.cols - 1; ++i) {
        // Partial pivoting
        int maxp = i;
        
        for (int j = i + 1; j < source.rows; ++j) {
            if (UTemp[j][i] > UTemp[maxp][i])
                maxp = j;
        }
        UTemp.swap_lines(i, maxp);
        b.swap_lines(i, maxp);

        numerical_analysis::Matrix<double> Gi {source.rows}; 
        for (int j = i + 1; j < source.rows; ++j) {
            Gi[j][i] = -(UTemp[j][i]/UTemp[i][i]);
        }
        LinvTemp = Gi * LinvTemp;
        UTemp = Gi * UTemp;
        b = Gi * b;
    }

    U = UTemp;
    Linv = LinvTemp;
}

template<typename TField>
void numerical_analysis::LinearSystemsMethods<TField>::backSubstitution(const Matrix<TField> A,
                      const Matrix<TField> b,
                      Matrix<TField> & x) {
    numerical_analysis::Matrix<double> xTemp {1, b.rows, 0};
    for (int i = A.rows - 1; i >= 0; --i) {
        xTemp[0][i] = b[i][0];
        for (int j = i + 1; j < A.cols; ++j) {
            xTemp[0][i] = xTemp[0][i] - xTemp[0][j]*A[i][j];
        } 
        xTemp[0][i] = xTemp[0][i]/A[i][i];
    }
    x = xTemp;
}

template<typename TField>
void numerical_analysis::LinearSystemsMethods<TField>::solveByLU(const Matrix<TField> A,
                      Matrix<TField> b,
                      Matrix<TField> & x) {

    numerical_analysis::Matrix<double> Linv {A.rows, 1};
    numerical_analysis::Matrix<double> U {A.rows, 1};
    numerical_analysis::LinearSystemsMethods<double>::getLinvU(A, b, Linv, U);
    numerical_analysis::LinearSystemsMethods<double>::backSubstitution(U, b, x);

}

template<typename TField>
void numerical_analysis::LinearSystemsMethods<TField>::solveByCholesky(const Matrix<TField> A,
                      Matrix<TField> b,
                      Matrix<TField> & x) {

    numerical_analysis::Matrix<TField> Linv {A.rows, 1};
    numerical_analysis::Matrix<TField> U {A.rows, 1};
    numerical_analysis::LinearSystemsMethods<TField>::getLinvU(A, b, Linv, U);

    numerical_analysis::Matrix<TField> Dinv {U.rows, U.cols, 0};

    for(int i = 0; i < U.cols; ++i) {
        if(U[i][i] < 0) 
            throw std::logic_error("A matrix cannot have negative values on diagonal."); 
        else 
            Dinv[i][i] =  1 / U[i][i];
    }

    x = Linv.transpose() * Dinv * b;
}
