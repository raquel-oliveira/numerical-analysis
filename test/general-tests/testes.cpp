#include "NaiveLinearSystemSolver.h"
#include "LinearSystemSolver.h"
#include "Matrix.h"
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <exception>

using namespace std;

int main(int argn, char ** argc) {

    long iterations;
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
    numerical_analysis::Matrix<double> x_jacobi {n, 1, 0};
    numerical_analysis::Matrix<double> x_gs {n, 1, 0};

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

    cout << "Matrix: " << endl;
    cout << matrix << endl;

    cout << b << endl;

    cout <<"A: " << endl;
    cout << matrix;

    cout << "b: " << endl;
    cout << b; 
    try {
        numerical_analysis::LinearSystemSolver<double>::solve_by_lu(matrix, b, x_lu);
        cout << "Solve by lu: ";
        cout << x_lu;
    } catch (exception& e){
        cout << e.what();
    }
    
    try {
        numerical_analysis::LinearSystemSolver<double>::solve_by_cholesky(matrix, b, x_cholesky);
        cout << "Solve by Cholesky: ";
        cout << x_cholesky;
    } catch (exception& e){
        cout << e.what();
    }


    cout << "Solve by Jacobi: ";
    try {
        iterations = numerical_analysis::NaiveLinearSystemSolver<double>::solve_by_jacobi(matrix, b, 0.00001, x_jacobi);
        cout << x_jacobi;
        cout << "Number of iterations " << iterations << endl;
    } catch (exception& e){
        cout << e.what();
    }

    cout << "Solve by Gauss: ";
    try {
        iterations = numerical_analysis::NaiveLinearSystemSolver<double>::solve_by_seidel(matrix, b, 0.000001, x_gs);
        cout << x_gs;
        cout << "Number of iterations " << iterations << endl;
    } catch (exception& e){
        cout << e.what();
    }
    

    return 0;
}
