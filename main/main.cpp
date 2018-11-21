#include"blocksplit.h"
#include"iostream"
#include"findMaxPoint.h"
#include"getGaus.h"
#include"conv.h"
#include"getsinc.h"
using namespace std;

#define SIZE1 120
#define SIZE2 40
void main()
{
	Mat dirtyImage = imread("C:\\Users\\jlfu\\Desktop\\dirtyimage.png", 0);
	Mat PSF = imread("C:\\Users\\jlfu\\Desktop\\PSF.png", 0);
	dirtyImage.convertTo(dirtyImage, CV_64FC1, 1 / 255.0);
	PSF.convertTo(PSF, CV_64FC1, 1 / 255.0);
	Mat dirtyImage1 = dirtyImage(Range(201, 497), Range(95, 391));
	Mat PSF1 = PSF(Range(201, 497), Range(95, 391));

	Mat gaus1(SIZE1, SIZE1, CV_64FC1);
	Mat gaus2(SIZE1, SIZE1, CV_64FC1);
	Mat sinc(SIZE2, SIZE2, CV_64FC1);
	Mat conv(SIZE1 + SIZE2 - 1, SIZE1 + SIZE2 - 1, CV_64FC1);
	getGaus(gaus1, SIZE1, SIZE1, 60, 60, 0.2, 2, 3);
	getGaus(gaus2, SIZE1, SIZE1, 80, 80, 0.3, 1, 4);
	//gaus1 += gaus2;
	getSinc(sinc, SIZE2);
	convGet(gaus1, sinc, conv);

	Vec<double, 5> v = { 2,2,2,2,4 };
	Mat m = getBlockSplit(gaus1, v, 0);
	namedWindow("dirtyimage", WINDOW_NORMAL);
	cvNamedWindow("m", WINDOW_NORMAL);
	imshow("dirtyimage", gaus1*255);
	imshow("m", m*255);
	waitKey();

	//double p; vector<double> v = { 2,4.5,4,3,7 };
	//p = findMaxOfNPoint(v);
	//cout << p << endl;
	//system("pause");
}