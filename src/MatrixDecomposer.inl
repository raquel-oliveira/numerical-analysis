template<typename TField>
void numerical_analysis::MatrixDecomposer<TField>::lu (const Matrix<TField> & source, 
		Matrix<TField> & L,
		Matrix<TField> & U,
		Matrix<TField> & P,
		bool partial_piv) {
	int n = source.rows;
	for (int i = 0; i < n - 1; ++i) {
		// Perform partial pivoting
		if (partial_piv) {
			// Select max
			int imax = i;
			for (int j = i + 1; j < n; ++j) {
				if (std::abs(source[j][i]) > std::abs(source[imax][i]))
					imax = j;
			}
			// Exchange
			P.swap_lines(imax, i);		
			U.swap_lines(imax, i, i, source.rows - 1);		
			L.swap_lines(imax, i, 0, i - 1);		
		}
		// Gaussian
		for (int j = i + 1; j < source.rows; ++j) {
			L[j][i] = U[j][i]/U[i][i];
			for (int m = i; m < source.rows; ++m) {
				U[j][m] = U[j][m] - L[j][i]*U[i][m];
			}	
		}
		/*std::cout << "L("<<i<<"): \n" << L << std::endl;
		std::cout << "U("<<i<<"): \n" << U << std::endl;
		std::cout << "P("<<i<<"): \n" << P << std::endl;*/
	}
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

	/*std::cout << "L(0): " << std::endl;
	std::cout << L << std::endl;*/

	for (int i = 1; i < source.rows - 1; ++i) {
		TField sum = 0.0;
		for (int k = 0; k < i; ++k)
			sum += (L[i][k] * L[i][k]);
		L[i][i] = std::sqrt(source[i][i] - sum);
		for (int j = i + 1; j < source.rows; ++j) {
			TField sum = 0.0;
			for (int k = 0; k < i; ++k) {
				sum += L[j][k]*L[i][k];
			}
			L[j][i] = (source[j][i] - sum) / L[i][i];
		}

		/*std::cout << "L("<<i<<"): " << std::endl;
		std::cout << L << std::endl;*/
	}

	TField sum = 0.0;
	for (int k = 0; k < source.rows - 1; ++k) {
		sum += (L[source.rows - 1][k]*L[source.rows - 1][k]);
	}
	L[source.rows-1][source.rows-1] = std::sqrt(source[source.rows-1][source.rows-1] - sum);

	//std::cout << "L("<<source.rows-1<<"): " << std::endl;
	//std::cout << L << std::endl;
}
