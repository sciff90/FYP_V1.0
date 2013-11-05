#include<cmath>
#include<random>
#include<omp.h>
#include <iostream>
#include <fstream>

void filter_out(double *theta,double *y,const double *u,int num_samples,int order)
{
	int ii,jj;
	long int N = 80000000;
	double *a = &theta[0];
	double *b = &theta[N*(order+1)];

	for(ii=0;ii<num_samples;ii++)y[ii] = 0.0;	
	
	
	for (ii = 0; ii < (num_samples-order); ii++)
	{
		for (jj = 1; jj <= order; jj++)
		{
			if(ii-jj>=0)
				y[ii] += - a[jj*N]*y[ii-jj];
			else
				y[ii] += 0;
		}
		for (jj = 0; jj <= order; jj++)
		{
			if(ii-jj>=0)
				y[ii] += b[jj*N]*u[ii-jj];
			else
				y[ii] +=0;
		}

		y[ii] = y[ii]/a[0];		
	}
}

unsigned int hash(unsigned int a)
{
	a = (a+0x7ed55d16) + (a<<12);
	a = (a^0xc761c23c) ^ (a>>19);
	a = (a+0x165667b1) + (a<<5);
	a = (a+0xd3a2646c) ^ (a<<9);
	a = (a+0xfd7046c5) + (a<<3);
	a = (a^0xb55a4f09) ^ (a>>16);
	return a;
}


void mcmc_kernel(double *u,double *y,double *theta,long int N,int order,int num_samples,double *theta_0,double elim,int nthreads,int tid)
{
	const int tt = nthreads;	
	double y_test[50];
	if(tid<tt)
	{
		unsigned int seed_normal = hash(tid);
		unsigned int seed_uniform = hash(tid);
		//Random number generator
		std::mt19937 rng_normal(seed_normal);
		std::mt19937 rng_uniform(seed_uniform);

		//Random number distributions
		std::normal_distribution<double> dist_norm;
		std::uniform_real_distribution<double> dist_uniform(0.0,1.0);
		
		for(int ii=0;ii<2*(order+1);ii++)		
			theta[N*ii+tid] = theta_0[ii];		

		double sigma = 1.0;

		int kk=0;
		int accepted=0;
		int flg=0;
		
		//printf("Random number from thread %d = %f\n",tid,sigma*dist_norm(rng_normal));
		for(int ii=tid+tt;ii<N;ii +=tt)
		{
			//Generate Proposal
			
			for(int jj=0;jj<2*(order+1);jj++)
			{
				theta[N*jj+ii] = theta[N*jj+ii-tt] + sigma*dist_norm(rng_normal);
			}
			theta[ii] = 1.0;			

			filter_out(&theta[ii],y_test,u,num_samples,order);
			float max_diff = 0;
			int max_loc = 0;
			float diff;						
			
			for(int jj=0;jj<num_samples-order;jj++)
			{
				diff = abs(y[jj]-y_test[jj]);
				if(diff > max_diff) 
				{
					max_diff = diff;
					max_loc = jj;
				}
			}
			if(max_diff>elim)
			{
				//printf("a1=  %f\n",theta[N*1+ii]);
				for(int jj=0;jj<2*(order+1);jj++)
					theta[N*jj+ii] = theta[N*jj+ii-tt];
			}
			else
			{
				if(1/(2*elim)*(1-max_diff)>=dist_uniform(rng_uniform))
				{
					accepted++;				
				}
				else
				{
					for(int jj=0;jj<2*(order+1);jj++)
						theta[N*jj+ii] = theta[N*jj+ii-tt];
				}				
			}
			kk++;
			if(kk%1000==0 && kk!=0 && flg==0)
			{
				if((double)accepted/1000>0.3)
				{
					sigma=sigma*1.2;
					printf("accepted = %d sigma = %f at ii = %d arate = %f\n",accepted,sigma,ii,(double)accepted/1000);					
				}
				else if((double)accepted/1000<0.2)
				{
					sigma = sigma/1.2;
					printf("accepted = %d sigma = %f at ii = %d arate = %f\n",accepted,sigma,ii,(double)accepted/1000);
				}
				else
				{
					flg=1;				
					printf("sigma = %f at ii = %d arate = %f\n",sigma,ii,(double)accepted/1000);
					printf("tt = %d\n",tt);
				}
				//printf("a_rate = %f\n",(double)accepted/1000);
				kk=0;
				accepted=0;
				ii = tid+tt;
				
			}
			

		//if(tid==1)				
			//printf("accepted from thread %d = %d\n",tid,accepted);
		}
		

	}

}

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
		mcmc_kernel(u,y,theta,N,order,num_samples,theta_0,elim,nthreads,tid);
	}


	return ;
}
