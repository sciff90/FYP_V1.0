#include <thrust/random.h>
#include <thrust/device_vector.h>
#include <thrust/transform.h>
#include <thrust/iterator/counting_iterator.h>
#include <iostream>
#include <fstream>

#define Nthreads 256
#define PI 3.141592

	__device__
void filter_out(double *theta,double *y,const double *u,int num_samples,int order)
{
	int ii,jj;
	int N = 20000000;
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
	__device__
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

	__global__ 
void mcmc_kernel(double *u,double *y,double *theta,int N,int order,int num_samples,double *theta_0,double elim)
{
	const int tid = blockDim.x*blockIdx.x+threadIdx.x;
	const int tt = blockDim.x*gridDim.x;	
	double y_test[50];
	if(tid<tt)
	{
		unsigned int seed_normal = hash(tid);
		thrust::default_random_engine rng_normal(seed_normal);
		//thrust::default_random_engine rng_uniform(seed_uniform);

		thrust::random::experimental::normal_distribution<double> dist_norm(0, 1);
		//thrust::random::experimental::uniform_distribution<double> dist_uniform(0, 1);
		
		for(int ii=0;ii<2*(order+1);ii++)		
			theta[N*ii+tid] = theta_0[ii];		

		double sigma = 0.01;

		int kk=0;
		int accepted=0;
		int flg=0;
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
				//printf("flagged max_diff = %f at postion %d\n",max_diff,max_loc);
				//printf("y[max_loc] = %f y_test[max_loc] = %f\n",y[max_loc],y_test[max_loc]);
				//printf("a1=  %f\n",theta[N*1+ii]);
				for(int jj=0;jj<2*(order+1);jj++)
					theta[N*jj+ii] = theta[N*jj+ii-tt];
			}
			else
			{
				accepted++;				
			}
			kk++;
			if(kk%100==0 && kk!=0 && flg==0)
			{
				if((double)accepted/100>0.4)
				{
					sigma=sigma*1.2;					
				}
				else if((double)accepted/100<0.3)
				{
					sigma = sigma/1.2;
				}
				else
				{
					flg=1;				
					printf("sigma = %f\n",sigma);
				}
				//printf("a_rate = %f\n",(double)accepted/1000);
				kk=0;
				accepted=0;
				ii = tid+tt;
				
			}
			

						
		}

	}

}

void mcmc(double *u,double *y,double *theta,int N,int order, int num_samples,double *theta_0,double elim)
{

	double *d_u,*d_y,*d_y_test,*d_theta,*d_theta_0;

	int u_size = num_samples*sizeof(double);
	int theta_size = N*2*(order+1)*sizeof(double);
	int theta_0_size = 2*(order+1)*sizeof(double);

	printf("N = %d\n",N);
	printf("order = %d\n",order);
	printf("num_samples = %d\n",num_samples);
	printf("theta_size = %d\n",theta_size);

	cudaMalloc((void**)&d_u, u_size ); 
	cudaMalloc((void**)&d_y, u_size ); 
	cudaMalloc((void**)&d_theta, theta_size );
	cudaMalloc((void**)&d_y_test, u_size );
	cudaMalloc((void**)&d_theta_0,theta_0_size);


	cudaMemcpy(d_u, u, u_size, cudaMemcpyHostToDevice );
	cudaMemcpy(d_y, y, u_size, cudaMemcpyHostToDevice );
	cudaMemcpy(d_theta_0,theta_0,theta_0_size,cudaMemcpyHostToDevice);

	mcmc_kernel<<<8,128>>>(d_u,d_y,d_theta,N,order,num_samples,d_theta_0,elim);

	cudaMemcpy( theta, d_theta, theta_size, cudaMemcpyDeviceToHost ); 


	cudaFree(d_u);
	cudaFree(d_y);
	cudaFree(d_theta);
	cudaFree(d_theta_0);
	cudaFree(d_y_test);

	return ;
}
