template<typename TField>
void numerical_analysis::LinearSystemSolver<TField>::back_substitution (const numerical_analysis::Matrix<TField> A,
		const numerical_analysis::Matrix<TField> b,
		numerical_analysis::Matrix<TField> & x) {
	numerical_analysis::Matrix<TField> xTemp {b.rows, 1, 0};
	for (int i = A.rows - 1; i >= 0; --i) {
		xTemp[i][0] = b[i][0];
		for (int j = i + 1; j < A.cols; ++j) {
			xTemp[i][0] = xTemp[i][0] - xTemp[j][0]*A[i][j];
		} 
		xTemp[i][0] = xTemp[i][0]/A[i][i];
	}
	x = xTemp;
}

template<typename TField>
void numerical_analysis::LinearSystemSolver<TField>::forward_substitution (const numerical_analysis::Matrix<TField> A,
		const numerical_analysis::Matrix<TField> b,
		numerical_analysis::Matrix<TField> & x) {

	//numerical_analysis::Matrix<TField> xTemp {b.rows, 1, 0};
	for (int i = 0; i < A.rows; ++i) {
		x[i][0] = b[i][0];
		for (int j = 0; j < i; ++j) {
			x[i][0] = x[i][0] - x[j][0]*A[i][j];
		} 
		x[i][0] = x[i][0]/A[i][i];
	}
	//x = xTemp;
}


template<typename TField>
void numerical_analysis::LinearSystemSolver<TField>::solve_by_lu(const Matrix<TField> A,
					  Matrix<TField> b,
					  Matrix<TField> & x,
					  bool partial_piv) {
	numerical_analysis::Matrix<TField> L {A.rows};
	numerical_analysis::Matrix<TField> U {A};
	numerical_analysis::Matrix<TField> P {A.rows};
	numerical_analysis::Matrix<TField> y {A.rows, 1, 0.0};
	numerical_analysis::MatrixDecomposer<TField>::lu(A, L, U, P, partial_piv);
	std::cout << L << std::endl;
	std::cout << U << std::endl;
	std::cout << P << std::endl;
	std::cout << L*U << std::endl;
	numerical_analysis::LinearSystemSolver<TField>::forward_substitution(L, P*b, y);
	numerical_analysis::LinearSystemSolver<TField>::back_substitution(U, y, x);
}


template<typename TField>
void numerical_analysis::LinearSystemSolver<TField>::solve_by_cholesky(const Matrix<TField> A,
					  Matrix<TField> b,
					  Matrix<TField> &x) {
	
	numerical_analysis::Matrix<TField> L {A.rows};	
	numerical_analysis::MatrixDecomposer<TField>::cholesky(A, L);
	numerical_analysis::Matrix<TField> y {A.rows, 1, 0};
	numerical_analysis::LinearSystemSolver<TField>::forward_substitution(L, b, y);
	numerical_analysis::LinearSystemSolver<TField>::back_substitution(L.transpose(), y, x);
}
