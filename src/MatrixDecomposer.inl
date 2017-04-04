template<typename TField>
void numerical_analysis::MatrixDecomposer<TField>::lu (const Matrix<TField> & source, 
		Matrix<TField> & L,
		Matrix<TField> & U) {

}

template<typename TField>
void numerical_analysis::MatrixDecomposer<TField>::ldlt (const Matrix<TField> & source,
		Matrix<TField> & L,
		Matrix<TField> & D) {

}

template<typename TField>
void numerical_analysis::MatrixDecomposer<TField>::cholesky (const Matrix<TField> & source,
		Matrix<TField> & L) {

	L[0][0] = std::sqrt(source[0][0]);
	
	for (int j = 1; j < source.rows; ++j) 
		L[j][0] = source[j][0]/L[0][0];

	for (int i = 1; i < source.rows - 1; ++i) {
		L[i][i] = source[i][i];
		for (int k = 0; k < i - 1; ++i)
			L[i][i] -= (L[i][k] * L[i][k]);
		L[i][i] = std::sqrt(L[i][i]);
		for (int j = i + 1; j < source.rows; ++j) {
			double sum = 0.0;
			for (int k = 0; k < i - 1; ++k) {
				sum += L[j][k]*L[i][k];
			}
			L[j][i] = (source[j][i] - sum) / L[i][i];
		}
	}

	double sum = 0.0;
	for (int k = 0; k < source.rows - 1; ++k) {
		sum += L[source.rows - 1][k]*L[source.rows - 1][k];
	L[source.rows-1][source.rows-1] = std::sqrt(source[source.rows-1][source.rows-1] - sum);
}
