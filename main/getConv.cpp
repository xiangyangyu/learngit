#include"conv.h"
void convGet(Mat gaus, Mat sinc, Mat &conv)
{
	Mat gaus1 = gaus.clone();
	Mat sinc1 = sinc.clone();
	flip(sinc1, sinc1, -1);//对要卷积的矩阵进行180度的旋转，成为卷积核
	copyMakeBorder(gaus1, gaus1, sinc1.rows - 1, 0, sinc1.rows - 1, 0, BORDER_CONSTANT, Scalar(0));//对要卷积的矩阵在它的上边和左边补零
	filter2D(gaus1, conv, -1, sinc1, cv::Point(0, 0), 0, BORDER_CONSTANT); //卷积
}