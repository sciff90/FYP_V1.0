#include<mcmc.h>
void mcmc(double *u,double *y,double *theta,long int N,int order, int num_samples,double *theta_0,double elim)
{

	int u_size = num_samples*sizeof(double);
	long int theta_size = N*2*(order+1)*sizeof(double);
	int theta_0_size = 2*(order+1)*sizeof(double);

	printf("N = %ld\n",N);
	printf("order = %d\n",order);
	printf("num_samples = %d\n",num_samples);
	printf("theta_size = %ld\n",theta_size);

	std::cout << "Parallel MCMC" << std::endl;

	//Setup OMP part nthreads etc
	int nthreads,tid;	
	#pragma omp parallel private(tid)
	{
		tid = omp_get_thread_num();

		if(tid == 0)
		{
			nthreads = omp_get_num_threads();
		}
	}

	std::cout << "nthreads = " << nthreads << std::endl;	

	#pragma omp parallel private(tid) 
	{
		tid = omp_get_thread_num();	
		slicesample(u,y,theta,N,order,num_samples,theta_0,elim,nthreads,tid);
	}


	return ;
}
