#include"blocksplit.h"
#include"gausfitting.h"
#include"getGaus.h"
#include"iostream"
using namespace std;
Mat getBlockSplit(Mat m, Vec<double, 5> v, int q)
{
	Mat m1 = m.clone();
	Mat *m2 = new Mat[v[4]];
	Mat m3 = Mat::zeros(m1.rows, m1.cols, CV_64FC1);
	for (int i = 0; i < 4; i++)
	{
		double *fit;
		fit = gausFitting(m1, q);
		m2[i].create(m1.rows, m1.cols, CV_64FC1);
		getGaus2(m2[i], fit[5], m1.rows, m1.cols, fit[3], fit[4], fit[2], fit[0], fit[1]);
		m1 -= m2[i];
		m3 += m2[i];
		delete[]fit; 
	}
	delete[]m2;
	//cout << m3 << endl;
	return m3;
}