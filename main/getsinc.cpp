#include"getsinc.h"
void getSinc(Mat &sinc, double size)
{
	const double PI = 4.0*atan(1.0); //‘≤÷‹¬ ¶–∏≥÷µ
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			if (i == size / 2 && j == size / 2)
				sinc.at<double>(i, j) = 1;
			else if (i == size / 2 && j != size / 2)
				sinc.at<double>(i, j) = (sin((j - size / 2)) / ((j - size / 2)));
			else if (i != size / 2 && j == size / 2)
				sinc.at<double>(i, j) = (sin((i - size / 2)) / ((i - size / 2)));
			else
				sinc.at<double>(i, j) = (sin((i - size / 2)) / ((i - size / 2)))*(sin((j - size / 2)) / ((j - size / 2)));

		}
	}
	cv::normalize(sinc, sinc);
}