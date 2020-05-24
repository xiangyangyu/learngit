#include"getInitialParameter.h"
#include"findMaxPoint.h"
#include"getGaus.h"
#include"iostream"
#include<time.h>
using namespace std;
double* getInitialParameter(Mat &sinc, double select)
{

	Mat m1(sinc.rows, sinc.cols, CV_64FC1), m2(sinc.rows, sinc.cols, CV_64FC1),
		m3(sinc.rows, sinc.cols, CV_64FC1), m4(sinc.rows, sinc.cols, CV_64FC1);
	double *ele = new double[7];
	double *p[8], *a1, *a2, *a3, *a4, center1, center2, ro, sigmax, sigmay, val;
	double *pp;
	pp = findArrayMaxPoint(sinc);

	if (select == 1)
	{
		//sigmax,sigmay设置太大了很容易发散
		p[0] = getGausParameter(m1, sinc.rows, sinc.cols, pp[1], pp[2], 0.0, 6, 6);
		p[1] = getGausParameter(m2, sinc.rows, sinc.cols, pp[1], pp[2], 0.0, 5, 5);
		p[2] = getGausParameter(m3, sinc.rows, sinc.cols, pp[1], pp[2], 0.0, 4, 4);
		p[3] = getGausParameter(m4, sinc.rows, sinc.cols, pp[1], pp[2], 0.0, 3, 3);
		a1 = findConvMaxPoint(m1, sinc);
		a2 = findConvMaxPoint(m2, sinc);
		a3 = findConvMaxPoint(m3, sinc);
		a4 = findConvMaxPoint(m4, sinc);
	}
	else if (select == 0)
	{
		double d[4];  int maxVal = 9, minVal = 1;
		srand((unsigned)time(NULL));
		for (int i = 0; i < 4; i++)
		{
			srand((int)time(0));
			d[i] = ((double)(rand() % 100) / 100)*(maxVal - minVal) + 1;
		}
		p[0] = getGausParameter(m1, sinc.rows, sinc.cols, pp[1], pp[2], 0.0, d[0], d[0]);
		p[1] = getGausParameter(m2, sinc.rows, sinc.cols, pp[1], pp[2], 0.0, d[1], d[1]);
		p[2] = getGausParameter(m3, sinc.rows, sinc.cols, pp[1], pp[2], 0.0, d[2], d[2]);
		p[3] = getGausParameter(m4, sinc.rows, sinc.cols, pp[1], pp[2], 0.0, d[3], d[3]);
		a1 = findConvMaxPoint(sinc, m1);
		a2 = findConvMaxPoint(sinc, m2);
		a3 = findConvMaxPoint(sinc, m3);
		a4 = findConvMaxPoint(sinc, m4);

	}
	double maxV;
	vector<double> v = { a1[0], a2[0], a3[0], a4[0] };
	maxV = findMaxOfNPoint(v);
	if (maxV == a1[0])
	{
		sigmax = p[0][0];
		sigmay = p[0][1];
		ro = p[0][2];
		center1 = p[0][3];
		center2 = p[0][4];
		val = a1[0];
	}
	else if (maxV == a2[0])
	{
		sigmax = p[1][0];
		sigmay = p[1][1];
		ro = p[1][2];
		center1 = p[1][3];
		center2 = p[1][4];
		val = a2[0];
	}
	else if (maxV == a3[0])
	{
		sigmax = p[2][0];
		sigmay = p[2][1];
		ro = p[2][2];
		center1 = p[2][3];
		center2 = p[2][4];
		val = a3[0];
	}
	else if (maxV == a4[0])
	{
		sigmax = p[3][0];
		sigmay = p[3][1];
		ro = p[3][2];
		center1 = p[3][3];
		center2 = p[3][4];
		val = a4[0];
	}
	ele[0] = sigmax; ele[1] = sigmay; ele[2] = ro; ele[3] = center1; ele[4] = center2; ele[5] = val;

	delete[]a1; delete[]a2; delete[]a3; delete[]a4; delete[]p[0]; delete[]p[1]; delete[]p[2]; delete[]p[3]; delete[]pp;
	//cout << q[0] << " " << q[1] << " " << q[2] << " " << q[3] << " " << q[4] << " " << q[5] << endl;
	return ele;
}
