#include <iostream>
#include "NaiveLinearSystemSolver.h"
#include "LinearSystemSolver.h"
#include "Matrix.h"
#include <cstdlib>

using namespace std;

int main(int argn, char ** argc){
	
	int n = 3;
	int randomico;
	if(argn > 1)
		n = stoi(argc[1]);

	numerical_analysis::Matrix<double> M{n};
	
	for(int i = 0; i < n; ++i){
		for(int j = i; j < n; ++j){
			randomico = rand() % 10 + 1;
			if(i == j) {
				M[i][j] = randomico;
			}else{
				M[i][j] = M[j][i] = randomico;
			}
		}
	}

	numerical_analysis::Matrix<double> A{n};

	A = M*M.transpose();

	cout << A << endl;
	
	return 0;
}