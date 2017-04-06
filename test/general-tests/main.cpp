#include "NaiveLinearSystemSolver.h"
#include "LinearSystemSolver.h"
#include "Matrix.h"
#include <chrono>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <exception>
#include <vector>

using namespace std;

int main(int argn, char ** argc) {

	string path,input;    
    int n;
    double value;
    stringstream ss;

    int v = 0;
    if(argn > 1){
        path = argc[1];
    }

    if(argn > 2){
    	value = stoi(argc[2]);
    }

    ifstream matrix_txt;
    matrix_txt.open(path, std::ifstream::in);

    getline(matrix_txt, input);
    n = stoi(input);

    numerical_analysis::Matrix<double> matrix{n};
    numerical_analysis::Matrix<double> b{n,1,0};
    numerical_analysis::Matrix<double> x_lu {1, n, 0};
    numerical_analysis::Matrix<double> x_cholesky {1, n, 0};
    numerical_analysis::Matrix<double> x_jacobi {n, 1, value};
    numerical_analysis::Matrix<double> x_gs {n, 1, value};

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

    cout << "Solve by lu: ";
    try {
        auto start = std::chrono::steady_clock::now();
        numerical_analysis::NaiveLinearSystemSolver<double>::solve_by_lu(matrix, b, x_lu);
        auto end = std::chrono::steady_clock::now();
        std::chrono::duration<double, std::milli> elapsedLU = (end - start);
        cout << x_lu.transpose();
        std::cout << "Time elapsed: " << elapsedLU.count() << std::endl;
        std::cout << "-------------------------------------------------\n";
    } catch (exception& e){
        cout << e.what();
    }
    
    cout << "Solve by Cholesky: ";        
    try {
        auto start = std::chrono::steady_clock::now();
        numerical_analysis::NaiveLinearSystemSolver<double>::solve_by_cholesky(matrix, b, x_cholesky);
        auto end = std::chrono::steady_clock::now();
        std::chrono::duration<double, std::milli> elapsedCholesky = (end - start);
        cout << x_cholesky.transpose();
        std::cout << "Time elapsed: " << elapsedCholesky.count() << std::endl;
        std::cout << "-------------------------------------------------\n";
    } catch (exception& e){
        cout << e.what();
    }


    cout << "Solve by Jacobi: ";
    try {
        auto start = std::chrono::steady_clock::now();
        int ite = numerical_analysis::NaiveLinearSystemSolver<double>::solve_by_jacobi(matrix, b, 0.01, x_jacobi);
        auto end = std::chrono::steady_clock::now();
        std::chrono::duration<double, std::milli> elapsedJacobi = (end - start);
        std::cout << x_jacobi.transpose();
        std::cout << "Time elapsed: " << elapsedJacobi.count() << std::endl;
        std::cout << "Iterations: "<< ite << std::endl;        
        std::cout << "-------------------------------------------------\n";
    } catch (exception& e){
        cout << e.what();
    }

    cout << "Solve by Gauss: ";
    try {
        auto start = std::chrono::steady_clock::now();
        int ite = numerical_analysis::NaiveLinearSystemSolver<double>::solve_by_seidel(matrix, b, 0.1, x_gs);
        auto end = std::chrono::steady_clock::now();
        std::chrono::duration<double, std::milli> elapsedSeidel = (end - start);
        std::cout << x_gs.transpose();
        std::cout << "Time elapsed: " << elapsedSeidel.count() << std::endl;
        std::cout << "Iterations: " << ite << std::endl;
        std::cout << "-------------------------------------------------\n";
    } catch (exception& e){
        cout << e.what();
    }
}