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

    //if(argn > 1){
    //    path = argc[1];
    //}

    /*vector<string> paths{ "inputs/input5x5.txt",
                          "inputs/input10x10.txt",
                          "inputs/input20x20.txt",
                          "inputs/input30x30.txt",
                          "inputs/input40x40.txt",
                          "inputs/input50x50.txt",
                          "inputs/input60x60.txt",
                          "inputs/input70x70.txt",
                          "inputs/input80x80.txt",
                          "inputs/input90x90.txt",
                          "inputs/input100x100.txt",
                          "inputs/input150x150.txt",
                          "inputs/input200x200.txt",
                          "inputs/input250x250.txt",
                          "inputs/input300x300.txt",
                          "inputs/input350x350.txt",
                          "inputs/input400x400.txt",
                          "inputs/input450x450.txt" };*/

    vector<string> paths{ "inputs/example1.txt",
                          "inputs/example2.txt",
                          "inputs/example3.txt",
                          "inputs/example4.txt",
                          "inputs/example5.txt" };

    for(int it = 0; it < paths.size(); ++it){
        path = paths[it];

        ifstream matrix_txt;
        matrix_txt.open(path, std::ifstream::in);

        getline(matrix_txt, input);
        n = stoi(input);

        numerical_analysis::Matrix<double> matrix{n};
        numerical_analysis::Matrix<double> b{n,1,0};
        numerical_analysis::Matrix<double> x_lu {1, n, 0};
        numerical_analysis::Matrix<double> x_cholesky {1, n, 0};
        numerical_analysis::Matrix<double> x_jacobi {n, 1, 1000000};
        numerical_analysis::Matrix<double> x_gs {n, 1, 1000000};

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

        /*//cout << "Solve by lu: \n";
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
        }*/


        //cout << "Solve by Jacobi: ";
        try {
            auto start = std::chrono::steady_clock::now();
            std::cout << "Jacobi: " << numerical_analysis::NaiveLinearSystemSolver<double>::solve_by_jacobi(matrix, b, 0.1, x_jacobi) << ", ";
            auto end = std::chrono::steady_clock::now();
            std::chrono::duration<double, std::milli> elapsedJacobi = (end - start);
            std::cout << elapsedJacobi.count() << " & ";
            //cout << x_jacobi;
        } catch (exception& e){
            cout << e.what();
        }

        //cout << "Solve by Gauss: ";
        try {
            auto start = std::chrono::steady_clock::now();
            std::cout << "GS: " << numerical_analysis::NaiveLinearSystemSolver<double>::solve_by_seidel(matrix, b, 0.1, x_gs) << ", ";
            auto end = std::chrono::steady_clock::now();
            std::chrono::duration<double, std::milli> elapsedSeidel = (end - start);
            std::cout << elapsedSeidel.count() << " \\\\ \n";
            //cout << x_gs;
        } catch (exception& e){
            cout << e.what();
        }
    }

    
    

    return 0;
}
