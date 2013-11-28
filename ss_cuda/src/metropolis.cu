#include<metropolis.h>

	__global__ 
void metropolis(double *u,double *y,double *theta,long int N,int order,int num_samples,double *theta_0,double elim)
{
	const int tid = blockDim.x*blockIdx.x+threadIdx.x;
	const int tt = blockDim.x*gridDim.x;	
	double y_test[50];
	if(tid<tt)
	{
		unsigned int seed_normal = hash(tid);
		unsigned int seed_uniform = hash(tid);
		thrust::default_random_engine rng_normal(seed_normal);
		thrust::default_random_engine rng_uniform(seed_uniform);

		thrust::random::experimental::normal_distribution<double> dist_norm(0, 1);
		thrust::random::uniform_real_distribution<double> dist_uniform(0, 1);
		
		for(int ii=0;ii<2*(order+1);ii++)		
			theta[N*ii+tid] = theta_0[ii];		

		double sigma = 1.0;

		int kk=0;
		int win_width = 100;
		int win_accept = 0;
		int accepted=0;
		int flg=0;
		int ii=tid+tt;
				
		while(ii<N)
		{
			//Generate Proposal
		
			for(int jj=0;jj<2*(order+1);jj++)
			{
				theta[N*jj+ii] = theta[N*jj+ii-tt] + sigma*dist_norm(rng_normal);
			}
			theta[ii] = 1.0;			

			filter(&theta[ii],y_test,u,num_samples,order,N);
			//Calculate P_ratio
			//the ratio of the probability of the new theta guess to the old theta guess
			//P(theta_new)/P(theta_old)
			
			double Pr = p_ratio(num_samples,elim,y,y_test);


			if(Pr <= dist_uniform(rng_uniform))
			{
				//Don't accept the new theta value
				for(int jj=0;jj<2*(order+1);jj++)
					theta[N*jj+ii] = theta[N*jj+ii-tt];
			}
			else
			{	
				//Accept new theta value
				accepted++;				
				win_accept++;
			}
			kk++;
			//Check burnin
			if(kk%win_width ==0&&flg==0)
			{
				if((double)win_accept/win_width>0.25)
				{
					sigma=sigma*1.2;
				}
				else if((double)win_accept/win_width<0.2)
				{
					sigma = sigma/1.2;
				}
				if(kk%10000==0)
				{
					flg=1;
					printf("acceptance rate from thread %d = %f\n",tid,(double)win_accept/win_width);
					win_width = win_width*2;

				}
				for(int jj=0;jj<2*(order+1);jj++)
					theta[N*jj+tid] = theta[N*jj+ii-tt];
				win_accept=0;
				ii=tid;
				accepted=0;
				
			}
			
			ii +=tt;		
		}
		printf("acceptance rate from thread %d = %f\n",tid,(float)accepted/(N/(tt)));
		printf("sigma from thread %d = %f\n",tid,sigma);			


	}

}

