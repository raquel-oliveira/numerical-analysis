#include <complex>
#include "FunctionZeroesFinder.h"
#include <functional>
#include <set>
#include <unordered_set>
#include <string>
#include <iostream>
#include "Matrix.h"
#include <cstdlib>
#include <cmath>
#include <exception>
#include <iomanip>
#include "opencv2/core/core.hpp"
#include "opencv2/highgui/highgui.hpp"

using namespace numerical_analysis;
using namespace std::complex_literals;


double precision = 0.001;

struct CLess {
	bool operator()(std::complex<double> z, double e){
		return std::norm<double>(z) < e;
	};
};

struct CLessSet {
	bool operator()(std::complex<double> z1, std::complex<double> z2){
		return (std::abs(std::norm<double>(z1- z2)) < precision);
	};
};


int main(int argn, char * args[]) {

	int n = atoi(args[1]);
	precision = std::pow(10,-atoi(args[2]));

	std::function<std::complex<double> (const std::complex<double> &)> f = 
		[&](const std::complex<double> & z){return std::pow(z,n) - (1.);};

	std::function<std::complex<double> (const std::complex<double> &)> df = 
		[&](const std::complex<double> & z){return std::complex<double>(n + 0.0, 0.0) * std::pow(z,n-1);};
	
	std::vector<std::pair<std::complex<double>, int>> roots_of_unity;

	std::vector<cv::Scalar> colors = {
		{0, 183, 15},
		{154, 79, 0},
		{0,0,196},
		{0,255,20},
		{130,40,196},
		{239,187,54},
		{54,187,239},
		{255,255,0},
		{0,234, 30}
	};

	std::complex<double> root;
	int N0 = 50;

	// Image
	cv::Mat M {1000, 1000, CV_8UC3, cv::Scalar(0,0,0)};

	for (double i = 0; i < 1000; ++i) {
		for (double j = 0; j < 1000; ++j) {
			
			double ii = (2*i)/1000 - 1;
			double jj = (2*j)/1000 - 1;

			std::complex<double> initial {ii, jj};
			
			try{
				// Newton Method
				numerical_analysis::FunctionZeroesFinder<std::complex<double>, CLess>::newton(
					f,
					df, 
					initial, 
					numerical_analysis::FunctionZeroesFinder<std::complex<double>>::StopCriteria::IMAGE |
				   	numerical_analysis::FunctionZeroesFinder<std::complex<double>>::StopCriteria::DELTA_IMAGE,
					precision, 
					root, 
					N0 
				);
				
				std::pair<std::set<std::pair<std::complex<double>, int>>::iterator, bool> root_set;
				// Add to roots
				CLessSet comp;
				if (std::norm<double>(f(root)) < precision) {
					bool found = false;
					for (auto r : roots_of_unity) {
						if (comp(r.first, root)) {
							found = true;
							cv::Vec3b & color = M.at<cv::Vec3b>(i, j);
							color[0] = colors[r.second][0];
							color[1] = colors[r.second][1];
							color[2] = colors[r.second][2];
							break;
						}
					}
					if (!found) {
						roots_of_unity.push_back(std::make_pair(root, roots_of_unity.size()));
						cv::Vec3b & color = M.at<cv::Vec3b>(i, j);
						color[0] = colors[roots_of_unity.size() - 1][0];
						color[1] = colors[roots_of_unity.size() - 1][1];
						color[2] = colors[roots_of_unity.size() - 1][2];
					}
				}


			} catch (std::exception& e){
				cv::Vec3b & color = M.at<cv::Vec3b>(i, j);
				color[0] = 255;	
				color[1] = 255;	
				color[2] = 255;	
			}
		}
	}

	cv::imshow("test", M);
	cv::waitKey(0);

	cv::imwrite("/home/vitorgreati/nbasin.jpg", M);

	return 0;
}
