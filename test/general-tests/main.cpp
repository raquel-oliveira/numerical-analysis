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
    double value = 0, e = 0.1;
    string alg = "l";
    stringstream ss;

    if(argn > 1){
        path = argc[1];
    }

    if(argn > 2){
        alg = argc[2];
    }

    if(argn > 3){
    	value = stod(argc[3]);
    }

    if(argn > 4){
    	e = stod(argc[4]);
    }


    ifstream matrix_txt;
    matrix_txt.open(path, std::ifstream::in);

    getline(matrix_txt, input);
    n = stoi(input);

    numerical_analysis::Matrix<double> matrix{n, n, 1, 0};
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

    cout << "Matriz A: \n" << matrix;

    cout << "Matriz b: \n" << b;

    if(alg[0] == 'l'){
        cout << "\n---------- Solve by lu ----------\n\n";
        try {
            auto start = std::chrono::steady_clock::now();
            numerical_analysis::LinearSystemSolver<double>::solve_by_lu(matrix, b, x_lu);
            auto end = std::chrono::steady_clock::now();
            std::chrono::duration<double, std::milli> elapsedLU = (end - start);
            cout << "Solution: " << x_lu.transpose();
            std::cout << "Time elapsed: " << elapsedLU.count() << std::endl;
            std::cout << "-------------------------------------------------\n";
        } catch (exception& e){
            cout << e.what();
        }
    }else if(alg[0] == 'c'){
        cout << "\n---------- Solve by Cholesky ----------\n\n";        
        try {
            auto start = std::chrono::steady_clock::now();
            numerical_analysis::LinearSystemSolver<double>::solve_by_cholesky(matrix, b, x_cholesky);
            auto end = std::chrono::steady_clock::now();
            std::chrono::duration<double, std::milli> elapsedCholesky = (end - start);
            cout << "Solution: " << x_cholesky.transpose();
            std::cout << "Time elapsed: " << elapsedCholesky.count() << std::endl;
            std::cout << "-------------------------------------------------\n";
        } catch (exception& e){
            cout << e.what();
        }
    }else if(alg[0] == 'j'){
        cout << "Vetor de entrada: \n" << x_jacobi;
        cout << "Epsilon: " << e << endl;
        cout << "\n---------- Solve by Jacobi ----------\n\n";
        try {
            auto start = std::chrono::steady_clock::now();
            int ite = numerical_analysis::NaiveLinearSystemSolver<double>::solve_by_jacobi(matrix, b, e, x_jacobi);
            auto end = std::chrono::steady_clock::now();
            std::chrono::duration<double, std::milli> elapsedJacobi = (end - start);
            std::cout << "Solution: " << x_jacobi.transpose();
            std::cout << "Time elapsed: " << elapsedJacobi.count() << std::endl;
            std::cout << "Iterations: "<< ite << std::endl;        
            std::cout << "-------------------------------------------------\n";
        } catch (exception& e){
            cout << e.what();
        }
    }else if(alg[0] == 'g'){
        cout << "Vetor de entrada: \n" << x_jacobi;
        cout << "Epsilon: " << e << endl;
        cout << "\n---------- Solve by Gauss ----------\n\n";
        try {
            auto start = std::chrono::steady_clock::now();
            int ite = numerical_analysis::NaiveLinearSystemSolver<double>::solve_by_seidel(matrix, b, e, x_gs);
            auto end = std::chrono::steady_clock::now();
            std::chrono::duration<double, std::milli> elapsedSeidel = (end - start);
            std::cout << "Solution: " << x_gs.transpose();
            std::cout << "Time elapsed: " << elapsedSeidel.count() << std::endl;
            std::cout << "Iterations: " << ite << std::endl;
            std::cout << "-------------------------------------------------\n";
        } catch (exception& e){
            cout << e.what();
        }  
    }else{
        cout << "Choose a available option to test:" << endl;
        cout << "l : Solve by Lu" << endl;
        cout << "c : Solve by Cholesky" << endl;
        cout << "j : Solve by Jacobi" << endl;
        cout << "g : Solve by Gauss Seidel" << endl;
    }
    
}
