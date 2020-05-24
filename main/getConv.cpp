#include"conv.h"
void convGet(Mat gaus, Mat sinc, Mat &conv)
{
	Mat gaus1 = gaus.clone();
	Mat sinc1 = sinc.clone();
	flip(sinc1, sinc1, -1);//��Ҫ����ľ������180�ȵ���ת����Ϊ�����
	copyMakeBorder(gaus1, gaus1, sinc1.rows - 1, 0, sinc1.rows - 1, 0, BORDER_CONSTANT, Scalar(0));//��Ҫ����ľ����������ϱߺ���߲���
	filter2D(gaus1, conv, -1, sinc1, cv::Point(0, 0), 0, BORDER_CONSTANT); //���
}