#pragma once
#include"opencv2/opencv.hpp"
using namespace cv;
void getGaus(Mat &gaus, double size_x, double size_y, double center_x, 
	         double center_y, double ro, double sigma_x, double sigma_y);
double* getGausParameter(Mat &gaus, double size_x, double size_y, double center_x, double center_y,
	double ro, double sigma_x, double sigma_y);
void getGaus2(Mat &gaus, double amp, double size_x, double size_y, double center_x, double center_y,
	double ro, double sigma_x, double sigma_y);