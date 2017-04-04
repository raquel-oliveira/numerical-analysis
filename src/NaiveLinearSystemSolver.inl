template<typename TField>
void numerical_analysis::NaiveLinearSystemSolver<TField>::get_linvu(const Matrix<TField> & source, 
			Matrix<TField> & b,
			Matrix<TField> & Linv,
			Matrix<TField> & U) {

	numerical_analysis::Matrix<TField> UTemp {source};
	numerical_analysis::Matrix<TField> LinvTemp {source.rows}; 

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

		numerical_analysis::Matrix<TField> Gi {source.rows}; 
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

	numerical_analysis::Matrix<TField> Linv {A.rows, 1};
	numerical_analysis::Matrix<TField> U {A.rows, 1};
	numerical_analysis::NaiveLinearSystemSolver<TField>::get_linvu(A, b, Linv, U);
	numerical_analysis::NaiveLinearSystemSolver<TField>::back_substitution(U, b, x);

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
