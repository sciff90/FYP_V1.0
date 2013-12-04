#include<filter.h>

__device__
void filter(double *theta,double *y,const double *u,int num_samples,int order, long int N)
{
	int ii,jj;	
	double *a = &theta[0];
	double *b = &theta[(order+1)];

	for(ii=0;ii<num_samples;ii++)y[ii] = 0.0;	


	for (ii = 0; ii < (num_samples); ii++)
	{
		for (jj = 1; jj <= order; jj++)
		{
			if(ii-jj>=0)
				y[ii] += - a[jj]*y[ii-jj];
			else
				y[ii] += 0;
		}
		for (jj = 0; jj <= order; jj++)
		{
			if(ii-jj>=0)
				y[ii] += b[jj]*u[ii-jj];
			else
				y[ii] +=0;
		}

		y[ii] = y[ii]/a[0];		
	}
}
