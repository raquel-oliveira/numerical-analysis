#include "NaiveLinearSystemSolver.h"
#include "LinearSystemSolver.h"
#include "Matrix.h"
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>

using namespace std;

int main(int argn, char ** argc) {

    string path,input;    
    int n;
    double value;
    stringstream ss;

    if(argn > 1){
        path = argc[1];
    }

    ifstream matrix_txt;
    matrix_txt.open(path, std::ifstream::in);

    getline(matrix_txt, input);
    n = stoi(input);

    numerical_analysis::Matrix<double> matrix{n};
    numerical_analysis::Matrix<double> b{n,1,0};
    numerical_analysis::Matrix<double> x_lu {1, n, 0};
    numerical_analysis::Matrix<double> x_cholesky {1, n, 0};
    numerical_analysis::Matrix<double> x_jacobi {1, n, 0};
    numerical_analysis::Matrix<double> x_gs {1, n, 0};

    int line = 0;
    int col = 0;
    while(getline(matrix_txt, input)){
        ss << input;
        ss >> value;
        b[line][0] = value;
        
        while(ss >> value){
            matrix[line][col++] = value;
        }

        col = 0;
        ++line;
 
        ss.str(string());
        ss.clear();
    }

    numerical_analysis::NaiveLinearSystemSolver<double>::solve_by_lu(matrix, b, x_lu);
    numerical_analysis::NaiveLinearSystemSolver<double>::solve_by_cholesky(matrix, b, x_cholesky);
    //numerical_analysis::NaiveLinearSystemSolver<double>::solve_by_jacobi(matrix, b, x_jacobi);
    //numerical_analysis::NaiveLinearSystemSolver<double>::solve_by_gs(matrix, b, x_gs);

    cout <<"A: " << endl;
    cout << matrix;

    cout << "b: " << endl;
    cout << b; 

    cout << "Solve by lu: ";
    cout << x_lu;

    cout << "Solve by cholesky: ";
    cout << x_cholesky;

    return 0;
}
