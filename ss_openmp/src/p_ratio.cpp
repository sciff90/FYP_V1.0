#include<p_ratio.h>
double p_ratio(int num_samples,double elim,double *y,double *y_test)
{
	//cost function
	float max_diff = 0;
	int max_loc = 0;
	float diff;						

	for(int jj=0;jj<num_samples;jj++)
	{
		diff = std::abs(y[jj]-y_test[jj]);
		if(diff > max_diff) 
		{
			max_diff = diff;
			max_loc = jj;
		}
	}

	if(max_diff>elim)
		return -1e100;
	
	else return -num_samples*std::log(elim);
}
