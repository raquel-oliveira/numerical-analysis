template<typename TField>
void numerical_analysis::NaiveLinearSystemSolver<TField>::get_linvu(const Matrix<TField> & source, 
			Matrix<TField> & b,
			Matrix<TField> & Linv,
			Matrix<TField> & U) {

	numerical_analysis::Matrix<TField> UTemp {source};
	numerical_analysis::Matrix<TField> LinvTemp {source.rows, source.rows, 1, 0}; 

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

		numerical_analysis::Matrix<TField> Gi {source.rows, source.rows, 1, 0}; 
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

	numerical_analysis::Matrix<TField> Linv {A.rows, A.rows, 1, 0};
	numerical_analysis::Matrix<TField> U {A.rows, A.rows, 1, 0};
	numerical_analysis::NaiveLinearSystemSolver<TField>::get_linvu(A, b, Linv, U);
	numerical_analysis::NaiveLinearSystemSolver<TField>::back_substitution(U, b, x);

}

template<typename TField>
void numerical_analysis::NaiveLinearSystemSolver<TField>::solve_by_cholesky(const Matrix<TField> A,
						  Matrix<TField> b,
						  Matrix<TField> & x) {

	numerical_analysis::Matrix<TField> Linv {A.rows, A.rows, 1, 0};
	numerical_analysis::Matrix<TField> U {A.rows, A.rows, 1, 0};
    numerical_analysis::NaiveLinearSystemSolver<TField>::get_linvu(A, b, Linv, U);

    numerical_analysis::Matrix<TField> Dinv {U.rows, U.cols, 0, 0};

    for(int i = 0; i < U.cols; ++i) {
        if(U[i][i] < 0) 
            throw std::logic_error("A matrix cannot have negative values on diagonal."); 
        else 
            Dinv[i][i] =  1 / U[i][i];
    }

    x = Linv.transpose() * Dinv * b;
}

template<typename TField>
double numerical_analysis::NaiveLinearSystemSolver<TField>::norm_euclidean(Matrix<TField> m){
    double sum = 0;
    for (int i = 0; i < m.rows; ++i) {
        sum += std::pow(m[i][0], 2);
    }
    return sqrt(sum);
}


template<typename TField>
long numerical_analysis::NaiveLinearSystemSolver<TField>::solve_by_jacobi(Matrix<TField> A,
		Matrix<TField> b,
		double p,
		Matrix<TField> & x) {

	if (p < 0 )
		throw std::logic_error("Varepsilon can not be < 0");
	
	numerical_analysis::Matrix<TField> s_aux {A.rows, A.rows, 1};
	numerical_analysis::Matrix<TField> de {A.rows, A.rows, 1};
	numerical_analysis::Matrix<TField> s {A.rows, A.rows, 1};

	de = (A.pow(-1)).diagonal();
	s_aux = A-A.diagonal();
	
	s = s_aux*de;
	if (!(s.norm_infinity() < 1)) 
		throw std::logic_error("Not indicated to use Jacobi method\n");
	
	numerical_analysis::Matrix<TField> aux {x.rows, x.cols, 0};
	numerical_analysis::Matrix<TField> n {x.rows, x.cols, 0}; //to check correctness

	s_aux = de.symmetric()*s_aux;

	long count = 0;
	std::cout << "Iteration("<<count<<"): \n";
	std::cout << x << std::endl;

	while (n.norm_infinity() > p) {
		aux = x;
		s = (s_aux*x) + (de*b);
		x = s;

		count++;
		
		std::cout << "Iteration("<<count<<"): \n";
		std::cout << x << std::endl;
		n = x-aux;
	}
	return count;
}

template<typename TField>
long numerical_analysis::NaiveLinearSystemSolver<TField>::solve_by_seidel(Matrix<TField> A,
		Matrix<TField> b,
		double p,
		Matrix<TField> & xe) {

	if (p < 0 ) throw std::logic_error("Varepsilon can not be < 0\n");

	numerical_analysis::Matrix<double> s {xe.rows, 1, 0};
	numerical_analysis::Matrix<double> aux {xe.rows, 1, 0};
	double sum, n;
	long count = 0;
	bool checkConvergence = true;

	std::cout << "Iteration("<<count<<"): \n";
	std::cout << xe << std::endl;
	if (A.rowCriteria()){
		std::cout << "The system converges by the row criteria" << std::endl;
	} else{
		throw std::logic_error("The system won't converge by the row criteria\n");
	}
	while (checkConvergence) {
		aux = xe;
		for (int i = 0 ; i < xe.rows ; i++){
			sum = 0;
			for (int j = 0; j < i; j++){
				sum = sum + A[i][j] * s[j][0];
			}
			for (int j = i + 1; j < xe.rows; j++) {
				sum = sum + A[i][j] * xe[j][0];
			}
			sum = sum - b[i][0];
			s[i][0] =  (-1*sum)/ A[i][i];
		}
		xe = s;
		n = (xe-aux).norm_infinity();
		count++;

		std::cout << "Iteration("<<count<<"): \n";
		std::cout << xe << std::endl;

		if (n < p) { checkConvergence = false; }
	}
	return count;
}
