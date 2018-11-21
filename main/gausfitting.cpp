#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>
#include"gausfitting.h"
#include"findMaxPoint.h"
#include"getGaus.h"
#include"getInitialParameter.h"
using namespace cv;

struct data {
	size_t n;
	int size_x;
	int size_y;
	double * y;
};

int expb_f(const gsl_vector * x, void *data,
	gsl_vector * f)
{
	size_t n = ((struct data *)data)->n;
	int size_x = ((struct data *)data)->size_x;
	int size_y = ((struct data *)data)->size_y;
	double *y = ((struct data *)data)->y;
	double sigma_x = gsl_vector_get(x, 0);
	double sigma_y = gsl_vector_get(x, 1);
	double ro = gsl_vector_get(x, 2);
	double center_x = gsl_vector_get(x, 3);
	double center_y = gsl_vector_get(x, 4);
	double max = gsl_vector_get(x, 5);

	/*高斯函数*/
	const double PI = 4.0*atan(1.0);
	double sector1 = 1 / (2 * sqrt(1 - ro * ro));
	for (size_t i = 0; i < size_x; i++)
		for (size_t j = 0; j < size_y; j++)
		{

			double t = i;
			double Yi = max * exp(-sector1 * ((i - center_x)*(i - center_x) /
				(sigma_x*sigma_x) - 2 * ro*(i - center_x)*(j - center_y) / (sigma_x*sigma_y) + (j - center_y)*(j - center_y) / (sigma_x*sigma_x)));
			gsl_vector_set(f, i*size_x + j, (Yi - y[i*size_x + j]));
		}
	return GSL_SUCCESS;
}

int expb_df(const gsl_vector * x, void *data,
	gsl_matrix * J)
{
	size_t n = ((struct data *)data)->n;
	int size_x = ((struct data *)data)->size_x;
	int size_y = ((struct data *)data)->size_y;

	double sigma_x = gsl_vector_get(x, 0);
	double sigma_y = gsl_vector_get(x, 1);
	double ro = gsl_vector_get(x, 2);
	double center_x = gsl_vector_get(x, 3);
	double center_y = gsl_vector_get(x, 4);
	double max = gsl_vector_get(x, 5);
	for (size_t i = 0; i < size_x; i++)
		for (size_t j = 0; j<size_y; j++)
		{
			const double PI = 4.0*atan(1.0); //圆周率π赋值
			double sector1 = sqrt(1 - ro * ro);
			double sector2 = 1 / (2 * sqrt(1 - ro * ro));
			double e = 1 / (2 * PI*sigma_x*sigma_y*sqrt(1 - ro * ro));
			double f = exp(-sector2 * ((i - center_x)*(i - center_x) /
				(sigma_x*sigma_x) - 2 * ro*(i - center_x)*(j - center_y) / (sigma_x*sigma_y) + (j - center_y)*(j - center_y) / (sigma_y*sigma_y)));
			double g = (i - center_x)*(i - center_x) /
				(sigma_x*sigma_x) - 2 * ro*(i - center_x)*(j - center_y) / (sigma_x*sigma_y) + (j - center_y)*(j - center_y) / (sigma_y*sigma_y);
			gsl_matrix_set(J, i*size_x + j, 0, max*f*((i - center_x)*(i - center_x) / (pow(sigma_x, 3)*sqrt(1 - ro * ro)) - ro * (i - center_x)*(j - center_y) / (sigma_x*sigma_x*sigma_y*sqrt(1 - ro * ro))));
			gsl_matrix_set(J, i*size_x + j, 1, max*f*((j - center_y)*(j - center_y) / (pow(sigma_y, 3)*sqrt(1 - ro * ro)) - ro * (i - center_x)*(j - center_y) / (sigma_y*sigma_y*sigma_x*sqrt(1 - ro * ro))));
			gsl_matrix_set(J, i*size_x + j, 2, max*f*(-1 / 2 * ro*pow(1 - ro * ro, -3 / 2)*g + (i - center_x)*(j - center_y) / (sigma_x*sigma_y*sqrt(1 - ro * ro))));
			gsl_matrix_set(J, i*size_x + j, 3, max*f*((i - center_x) / (sigma_x*sigma_x*sqrt(1 - ro * ro)) - ro * (j - center_y) / (sigma_x*sigma_y*sqrt(1 - ro * ro))));
			gsl_matrix_set(J, i*size_x + j, 4, max*f*((j - center_y) / (sigma_y*sigma_y*sqrt(1 - ro * ro)) - ro * (i - center_x) / (sigma_y*sigma_x*sqrt(1 - ro * ro))));
			gsl_matrix_set(J, i*size_x + j, 5, f);
		}
	return GSL_SUCCESS;
}

int expb_fdf(const gsl_vector * x, void *data,
	gsl_vector * f, gsl_matrix * J)
{
	expb_f(x, data, f);
	expb_df(x, data, J);
	return GSL_SUCCESS;
}

void print_state(size_t iter, gsl_multifit_fdfsolver * s);

double* gausFitting(Mat &gaus, int hh)
{
	Mat gaus1 = gaus.clone();
	const gsl_multifit_fdfsolver_type *T;
	gsl_multifit_fdfsolver *s;
	int status;
	size_t i, iter = 0;
	const int rows = gaus1.rows;
	const int cols = gaus1.cols;
	const size_t n1 = rows * cols;
	const size_t p = 6;

	double *y = new double[1000000]; double *sigma = new double[1000000];

	struct data d = { n1, rows, cols, y };
	gsl_multifit_function_fdf f;
	double *pp;
	pp = getInitialParameter(gaus1, hh);
	double x_init[6] = { pp[0],pp[1],pp[2],pp[3],pp[4],pp[5] };
	gsl_vector_view x = gsl_vector_view_array(x_init, p);

	f.f = &expb_f;
	f.df = &expb_df;
	f.fdf = &expb_fdf;
	f.n = n1;
	f.p = p;
	f.params = &d;
	/* This is the data to be fitted */
	const double PI = 4.0*atan(1.0);
	for (i = 0; i < rows; i++)
	{
		for (size_t j = 0; j < cols; j++)
		{
			y[i*rows + j] = gaus1.at<double>(i, j);
		}
	}
	T = gsl_multifit_fdfsolver_lmsder;
	s = gsl_multifit_fdfsolver_alloc(T, n1, p);
	gsl_multifit_fdfsolver_set(s, &f, &x.vector);
	print_state(iter, s);
	do
	{
		iter++;
		status = gsl_multifit_fdfsolver_iterate(s);
		printf("status = %s\n", gsl_strerror(status));
		print_state(iter, s);
		if (status)
			break;
		status = gsl_multifit_test_delta(s->dx, s->x,
			1e-4, 1e-4);
	} while (status == GSL_CONTINUE && iter < 500);

	double *ppp = new double[7];
#define FIT(i) gsl_vector_get(s->x, i)
#define ERR(i) sqrt(gsl_matrix_get(covar,i,i))
	{
		ppp[0] = FIT(0); ppp[1] = FIT(1); ppp[2] = FIT(2); ppp[3] = FIT(3); ppp[4] = FIT(4);
		ppp[5] = FIT(5);
		//double chi = gsl_blas_dnrm2(s->f);
		//double dof = n * n - p;
		//double c = GSL_MAX_DBL(1, chi / sqrt(dof));
		//printf("chisq/dof = %g\n", pow(chi, 2.0) / dof);
		printf("sigmax = %.5f\n", FIT(0));
		printf("sigmay = %.5f\n", FIT(1));
		printf("ro = %.5f\n", FIT(2));
		printf("center1 = %.5f\n", FIT(3));
		printf("center2 = %.5f\n", FIT(4));
		printf("max = %.5f\n", FIT(5));

	}
	printf("status = %s\n", gsl_strerror(status));
	gsl_multifit_fdfsolver_free(s);

	delete[]y; delete[]sigma; delete[]pp;
	return ppp;
}

void print_state(size_t iter, gsl_multifit_fdfsolver * s)
{
	printf("iter: %3u x = % 15.8f % 15.8f % 15.8f % 15.8f %15.8f %15.8f "
		"|f(x)| = %g\n",
		iter,
		gsl_vector_get(s->x, 0),
		gsl_vector_get(s->x, 1),
		gsl_vector_get(s->x, 2),
		gsl_vector_get(s->x, 3),
		gsl_vector_get(s->x, 4),
		gsl_vector_get(s->x, 5),
		gsl_blas_dnrm2(s->f));
}