#include "NaiveLinearSystemSolver.h"
#include "LinearSystemSolver.h"
#include "Matrix.h"
#include <chrono>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <exception>

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

    //cout << "Solve by lu: \n";
    try {
        auto start = std::chrono::steady_clock::now();
        numerical_analysis::LinearSystemSolver<double>::solve_by_lu(matrix, b, x_lu);
        auto end = std::chrono::steady_clock::now();
        std::chrono::duration<double, std::milli> elapsedLU = (end - start);
        std::cout << "& " << elapsedLU.count() << '\n';
        //cout << x_lu;
    } catch (exception& e){
        cout << e.what();
    }
    
    //cout << "Solve by Cholesky: \n";        
    try {
        auto start = std::chrono::steady_clock::now();
        numerical_analysis::LinearSystemSolver<double>::solve_by_cholesky(matrix, b, x_cholesky);
        auto end = std::chrono::steady_clock::now();
        std::chrono::duration<double, std::milli> elapsedCholesky = (end - start);
        std::cout << "& " << elapsedCholesky.count() << '\n';
        //cout << x_cholesky;
    } catch (exception& e){
        cout << e.what();
    }


    /*cout << "Solve by Jacobi: ";
    try {
        numerical_analysis::NaiveLinearSystemSolver<double>::solve_by_jacobi(matrix, b, 0.001, x_jacobi);
        cout << x_jacobi;
    } catch (exception& e){
        cout << e.what();
    }

    cout << "Solve by Gauss: ";
    try {
        numerical_analysis::NaiveLinearSystemSolver<double>::solve_by_seidel(matrix, b, 0.001, x_gs);
        cout << x_gs;
    } catch (exception& e){
        cout << e.what();
    }*/
    

    return 0;
}
