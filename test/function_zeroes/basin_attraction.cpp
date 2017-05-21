#include <complex>
#include "FunctionZeroesFinder.h"
#include <functional>
#include <iostream>
#include "Matrix.h"
#include <cmath>
#include <exception>
#include <iomanip>
#include "opencv2/core/core.hpp"
#include "opencv2/highgui/highgui.hpp"

using namespace numerical_analysis;
using namespace std::complex_literals;

struct CLess {
	bool operator()(std::complex<double> z, double e){
		return std::norm<double>(z) < e;
	};
};

int main(void) {

	std::function<std::complex<double> (const std::complex<double> &)> f = 
		[](const std::complex<double> & z){return std::pow(z,3) - (1.);};

	std::function<std::complex<double> (const std::complex<double> &)> df = 
		[](const std::complex<double> & z){return 3. * std::pow(z,2);};
	
	std::vector<std::complex<double>> roots = {
		std::complex<double>(1.0,0.0),
		std::complex<double>(-1.0/2, std::sqrt(3)/2),
		std::complex<double>(-1.0/2, -std::sqrt(3)/2)
	};

	std::vector<cv::Scalar> colors = {
		{0, 183, 15},
		{154, 79, 0},
		{0,0,196}
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
					0.001, 
					root, 
					N0 
				);
				
				// Check what root
				
				for (int k = 0; k < roots.size(); ++k) {
					if (std::abs(std::norm<double>(root - roots[k])) < 0.001) {
						cv::Vec3b & color = M.at<cv::Vec3b>(i, j);
						color[0] = colors[k][0];
						color[1] = colors[k][1];
						color[2] = colors[k][2];
						break;
					}
				}
			}catch (std::exception& e){
				cv::Vec3b & color = M.at<cv::Vec3b>(i, j);
				color[0] = 255;	
				color[1] = 255;	
				color[2] = 255;	
			}
		}
	}

	cv::imshow("test", M);
	cv::waitKey(0);

	cv::imwrite("~/n3basin.jpg", M);

	return 0;
}
