#include <iostream>
#include <string>
#include <fstream>
#include <random>
#include "NaiveLinearSystemSolver.h"
#include "LinearSystemSolver.h"
#include "Matrix.h"
#include <cstdlib>

using namespace std;

int main(int argn, char ** argc){
	
	int n = 3;
	double randomico;
	if(argn > 1)
		n = stoi(argc[1]);

	numerical_analysis::Matrix<double> M{n};
	
	std::random_device rd;
	std::default_random_engine generator(rd()); // rd() provides a random seed
	std::uniform_real_distribution<double> distribution(0.1,1);

	for(int i = 0; i < n; ++i){
		for(int j = 0; j < n; ++j){
			randomico = distribution(generator);
			i != j ? 0 : M[i][j] = M[j][i] = randomico;
		}
	}

	numerical_analysis::Matrix<double> A{n};
	numerical_analysis::Matrix<double> b{n,1,0};

	A = M*M.transpose();

	cout << A << endl;

	for(int i = 0; i < n; ++i){
		for(int j = 0; j < n; ++j){
			b[i][0] += A[i][j];
		}
	}

	cout << b << endl;

	fstream input;
	input.open("inputs/indirect_input"+to_string(n)+"x"+to_string(n)+".txt", fstream::out);

	input << to_string(n) << endl;
	for(int i = 0; i < n; ++i){
		input << b[i][0];
		for(int j = 0; j < n; ++j){
			input << " " << A[i][j];
		}
		input << endl;
	}

	input.close();
	

	return 0;
}