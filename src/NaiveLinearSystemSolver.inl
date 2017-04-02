#include "NaiveLinearSystemSolver.h"
#include <math.h>

template<typename TField>
void numerical_analysis::NaiveLinearSystemSolver<TField>::get_linvu(const Matrix<TField> & source, 
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
void numerical_analysis::NaiveLinearSystemSolver<TField>::solve_by_lu(const Matrix<TField> A,
					  Matrix<TField> b,
					  Matrix<TField> & x) {

	numerical_analysis::Matrix<double> Linv {A.rows, 1};
	numerical_analysis::Matrix<double> U {A.rows, 1};
	numerical_analysis::NaiveLinearSystemSolver<double>::get_linvu(A, b, Linv, U);
	numerical_analysis::NaiveLinearSystemSolver<double>::back_substitution(U, b, x);

}

template<typename TField>
void numerical_analysis::NaiveLinearSystemSolver<TField>::solve_by_cholesky(const Matrix<TField> A,
						  Matrix<TField> b,
						  Matrix<TField> & x) {

	numerical_analysis::Matrix<TField> Linv {A.rows, 1};
	numerical_analysis::Matrix<TField> U {A.rows, 1};
    numerical_analysis::NaiveLinearSystemSolver<TField>::get_linvu(A, b, Linv, U);

    numerical_analysis::Matrix<TField> Dinv {U.rows, U.cols, 0};

    for(int i = 0; i < U.cols; ++i) {
        if(U[i][i] < 0) 
            throw std::logic_error("A matrix cannot have negative values on diagonal."); 
        else 
            Dinv[i][i] =  1 / U[i][i];
    }

    x = Linv.transpose() * Dinv * b;
}


template<typename TField>
double numerical_analysis::NaiveLinearSystemSolver<TField>::getNorm(Matrix<TField> c){
	//TODO:: verify if c has only 1 col or change to general norm
	double sum = 0;
	for (int i = 0; i < c.rows; ++i) {
        sum += std::pow(c[i][0], 2);
   	}
   	return sqrt(sum);
}

template<typename TField>
void numerical_analysis::NaiveLinearSystemSolver<TField>::solve_by_jacobi(const Matrix<TField> A,
		Matrix<TField> b,
		double c,
		Matrix<TField> & x) {
	if (c < 0 )
		throw std::logic_error("Varepsilon can not be < 0");
	else {
		numerical_analysis::Matrix<TField> aux {x.rows, x.cols};
		numerical_analysis::Matrix<TField> s {A.rows, 1};
		numerical_analysis::Matrix<TField> de {A.rows, 1};
		numerical_analysis::Matrix<TField> e = x;
		numerical_analysis::Matrix<TField> n {e.rows, e.cols}; //check correctness
		while (numerical_analysis::NaiveLinearSystemSolver<double>::getNorm(n) > c) {
			aux = e;
			de = (A.pow(-1)).diagonal();

			//s = A-A.diagonal();
			s = A;
			for (int i = 0; i < A.rows; ++i) {
				for (int j = 0; j < A.cols; ++j) {
					if (i == j) s[i][j] = 0;
				}
			}

			s = (de.symmetric()*s*e) + (de*b);
			e = s;

			n = e-aux;
		}
		x = s;
	}
}