#include"getGaus.h"
void getGaus(Mat &gaus, double size_x, double size_y, double center_x, double center_y, 
	             double ro, double sigma_x, double sigma_y)
{
	const double PI = 4.0*atan(1.0); //圆周率π赋值
	double sector1 = sqrt(1 - ro * ro);
	double sector2 = 1 / (2 * sector1);
	for (int i = 0; i < size_x; i++)
	{
		for (int j = 0; j < size_y; j++)
		{
			gaus.at<double>(i, j) = (1 / (2 * PI*sigma_x*sigma_y*sector1)   *exp(-sector2 * ((i - center_x)*(i - center_x)
				/(sigma_x*sigma_x) - 2 * ro*(i - center_x)*(j - center_y) / (sigma_x*sigma_y)
				+ (j - center_y)*(j - center_y) / (sigma_y*sigma_y))));
		}
	}
	cv::normalize(gaus, gaus);                //将gaus矩阵归一化
}

double* getGausParameter(Mat &gaus, double size_x, double size_y, double center_x, double center_y, 
	              double ro, double sigma_x, double sigma_y)
{
	double *par = new double[5];
	const double PI = 4.0*atan(1.0); //圆周率π赋值
	double sector1 = sqrt(1 - ro * ro);
	double sector2 = 1 / (2 * sqrt(1 - ro * ro));
	for (int i = 0; i < size_x; i++)
	{
		for (int j = 0; j < size_y; j++)
		{
			gaus.at<double>(i, j) = (1 / (2 * PI*sigma_x*sigma_y*sector1)   *exp(-sector2 * ((i - center_x)*(i - center_x) /
				(sigma_x*sigma_x) - 2 * ro*(i - center_x)*(j - center_y) / (sigma_x*sigma_y) + (j - center_y)*(j - center_y) / (sigma_y*sigma_y))));
		}
	}
	cv::normalize(gaus, gaus);                //将gaus矩阵归一化
	par[0] = sigma_x; par[1] = sigma_y; par[2] = ro; par[3] = center_x; par[4] = center_y;
	return par;
}

void getGaus2(Mat &gaus, double amp, double size_x, double size_y, double center_x, double center_y,
	            double ro, double sigma_x, double sigma_y)
{
	const double PI = 4.0*atan(1.0); //圆周率π赋值
	double sector1 = sqrt(1 - ro * ro);
	double sector2 = 1 / (2 * sqrt(1 - ro * ro));
	for (int i = 0; i < size_x; i++)
	{
		for (int j = 0; j < size_y; j++)
		{
			gaus.at<double>(i, j) = amp *exp(-sector2 * ((i - center_x)*(i - center_x) /
				(sigma_x*sigma_x) - 2 * ro*(i - center_x)*(j - center_y) / (sigma_x*sigma_y) + (j - center_y)*(j - center_y) / (sigma_y*sigma_y)));
		}
	}
	//cv::normalize(gaus, gaus);                //将gaus矩阵归一化
}