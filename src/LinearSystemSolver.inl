template<typename TField>
void numerical_analysis::LinearSystemSolver<TField>::back_substitution (
		const numerical_analysis::Matrix<TField> A,
		const numerical_analysis::Matrix<TField> b,
		numerical_analysis::Matrix<TField> & x) {
	numerical_analysis::Matrix<TField> xTemp {1, b.rows, 0};
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
void numerical_analysis::LinearSystemSolver<TField>::forward_substitution (
		const numerical_analysis::Matrix<TField> A,
		const numerical_analysis::Matrix<TField> b,
		numerical_analysis::Matrix<TField> & x) {

	numerical_analysis::Matrix<TField> xTemp {1, b.rows, 0};
	for (int i = 0; i < A.rows; --i) {
		xTemp[0][i] = b[i][0];
		for (int j = 0; j < i; ++j) {
			xTemp[0][i] = xTemp[0][i] - xTemp[0][j]*A[i][j];
		} 
		xTemp[0][i] = xTemp[0][i]/A[i][i];
	}
	x = xTemp;
}


template<typename TField>
void numerical_analysis::LinearSystemSolver<TField>::solve_by_lu(const Matrix<TField> A,
					  Matrix<TField> b,
					  Matrix<TField> & x) {

}


template<typename TField>
void numerical_analysis::LinearSystemSolver<TField>::solve_by_cholesky(const Matrix<TField> A,
					  Matrix<TField> b,
					  Matrix<TField> &x) {
	
	numerical_analysis::Matrix<double> L {A.rows};	
	numerical_analysis::MatrixDecomposer<double>::cholesky(A, L);
	numerical_analysis::LinearSystemSolver<double>::forward_substitution(L, b, x);

}
