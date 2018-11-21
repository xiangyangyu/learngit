#include"findMaxPoint.h"
#include"conv.h"
#include"iostream"
using namespace std;
double* findArrayMaxPoint(Mat &gaus)                         /*找一个Mat类型中最大的元素，并返回最大元素值、
                                                             最大值对应的横坐标、最大值对应的纵坐标*/
{
	double *p = new double[3];
	double minVal; double maxVal; Point minLoc; Point maxLoc;
	minMaxLoc(gaus, &minVal, &maxVal, &minLoc, &maxLoc);
	if (abs(maxVal) > abs(minVal))
	{
		p[0] = maxVal; p[1] = maxLoc.x; p[2] = maxLoc.y;
	}
	else
	{
		p[0] = minVal; p[1] = minLoc.x; p[2] = minLoc.y;
	}
	return p;
}

double* findConvMaxPoint(Mat sinc, Mat gaus)
{
	Mat sinc1 = sinc.clone();
	Mat gaus1 = gaus.clone();
	Mat conv1(sinc.rows + gaus.rows - 1, sinc.rows + gaus.rows - 1, CV_64FC1);
	convGet(sinc1, gaus1, conv1);
	return findArrayMaxPoint(conv1);
}

double findMaxOfNPoint(vector<double> &a)
{
	double p;
	vector<double>::iterator max = max_element(begin(a), end(a));
	p = *max;
	return p;
}