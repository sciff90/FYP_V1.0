#include<slicesample.h>

	__global__ 
void slicesample(double *u,double *y,double *theta,long int N,int order,int num_samples,double *theta_0,double elim)
{
	const int tid = blockDim.x*blockIdx.x+threadIdx.x;
	const int tt = blockDim.x*gridDim.x;	
		
	if(tid<tt)
	{
		unsigned int seed_normal = hash(tid);
		unsigned int seed_uniform = hash(tid);
		thrust::default_random_engine rng_normal(seed_normal);
		thrust::default_random_engine rng_uniform(seed_uniform);

		thrust::random::experimental::normal_distribution<double> dist_norm(0.0, 1.0);
		thrust::random::uniform_real_distribution<double> dist_uniform(0.0, 1.0);
		
		//Populate first theta value with initial theta_0 vector
		for(int jj=0;jj<2*(order+1);jj++)
		{		
			theta[N*jj+tid] = theta_0[jj];	
			theta[N*jj+(tid+tt)] = theta_0[jj];
			//printf("theta[%ld] = %f\n",N*jj+tid,theta[N*jj+tid]);		


		}		

		//Start proposal distribution with standard deviation sigma
		double width = 0.01;
		double y_test[100];
		double theta_new[4],theta_new_l[4],theta_new_r[4],theta_prime[4];
		
		int ii = tid+tt;
		for(int jj=0;jj<(2*(order+1));jj++)
		{
			theta_new[jj] = theta[N*jj+ii];
			//printf("theta_new[%d] = %f\n",jj,theta_new[jj]);
		}

		filter(theta_new,y_test,u,num_samples,order,N);
		double pstar = p_ratio(num_samples,elim,y,y_test);
				
		while(ii<N)
		{
			
			double Puprime = pstar+std::log(dist_uniform(rng_uniform));
		
			int kk=1;
			while(kk < 2*(order+1))
			{
				for(int jj=0;jj<2*(order+1);jj++)
				{
					theta_new[jj] = theta[N*jj+ii];
					theta_new_l[jj] = theta[N*jj+ii];
					theta_new_r[jj] = theta[N*jj+ii];
					theta_prime[jj] = theta[N*jj+ii];
				}
				
				double bit = dist_uniform(rng_uniform);
				theta_new_l[kk] = theta_new[kk]-bit*width;
				theta_new_r[kk] = theta_new[kk]+bit*width;
				
				//Step out to span target density
				filter(theta_new_l,y_test,u,num_samples,order,N);
				while(p_ratio(num_samples,elim,y,y_test)>Puprime)
				{
					theta_new_l[kk] = theta_new_l[kk]-width;
					filter(theta_new_l,y_test,u,num_samples,order,N);
				}

				filter(theta_new_r,y_test,u,num_samples,order,N);
				while(p_ratio(num_samples,elim,y,y_test)>Puprime)
				{
					theta_new_r[kk] = theta_new_r[kk]+width;
					filter(theta_new_l,y_test,u,num_samples,order,N);
				}

				int stepcount = 0;

				while(1)
				{
					stepcount++;
					
					thrust::random::uniform_real_distribution<double> dist_uniform_2(theta_new_l[kk],theta_new_r[kk]);
					theta_prime[kk] = dist_uniform_2(rng_uniform);
					filter(theta_prime,y_test,u,num_samples,order,N);
					pstar = p_ratio(num_samples,elim,y,y_test);

					if(pstar>Puprime)
					{
						//printf("BREAK\n");
						break;
					}
						
					else
					{
						if(theta_prime[kk]>theta_new[kk])
							theta_new_r[kk] = theta_prime[kk];
						else if(theta_prime[kk]<theta_new[kk])
							theta_new_l[kk] = theta_prime[kk];
						else
							printf("ERROR\n");
					}
				}
				theta[N*kk+ii] = theta_prime[kk];
				kk=kk+1;
				

			}

			
			ii = ii+tt;
			if(ii>=N)
				break;
			for(int jj=0;jj<2*(order+1);jj++)
				theta[N*jj+ii] = theta[N*jj+(ii-tt)];
		}
	}

			

}

