#pragma once
#include"opencv2/opencv.hpp"
using namespace std;
using namespace cv;
double* findArrayMaxPoint(Mat &gaus);
double* findConvMaxPoint(Mat sinc, Mat gaus);
double findMaxOfNPoint(vector<double> &a);