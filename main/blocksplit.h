#pragma once
#include"opencv2/opencv.hpp"
using namespace cv;
Mat getBlockSplit(Mat m, Vec<double, 5> v, int q);