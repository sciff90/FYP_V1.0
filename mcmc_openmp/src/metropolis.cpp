#include<metropolis.h>

void metropolis(double *u,double *y,double *theta,long int N,int order,int num_samples,double *theta_0,double elim,int nthreads,int tid)
{
	const int tt = nthreads;		
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

		//Populate first theta value with initial theta_0 vector
		for(int ii=0;ii<2*(order+1);ii++)		
			theta[N*ii+tid] = theta_0[ii];		
		
		//Start proposal distribution with standard deviation sigma
		double sigma = 1.0;
		
		//kk is the window counter for checking acceptance rate
		//flg is the burnin flag 0 not burned in, 1 burned in
		//Different chains are strided by the nthreads tt

		int kk=0;
		int win_width = 1000;
		int win_accept = 0;
		int accepted=0;
		int flg=0;
		int ii = tid+tt;
		double y_test[num_samples];

		//Looped while the counter ii is less than the length of the combination of nthread chains N

		while(ii<N)
		{
			//Generate Proposal

			for(int jj=0;jj<2*(order+1);jj++)
			{
				theta[N*jj+ii] = theta[N*jj+ii-tt] + sigma*dist_norm(rng_normal);
			}
			//coefficient a_0 is set to be 1.0
			theta[ii] = 1.0;

			//Filter output for test theta coefficients
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
			if(kk%win_width ==0 && flg==0)
			{
				if((double)win_accept/win_width>0.25)
				{
					sigma=sigma*1.2;
				}
				else if((double)win_accept/win_width<0.2)
				{
					sigma = sigma/1.2;
				}
				if(kk%100000==0)
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
		printf("acceptance rate from thread %d = %f\n",tid,(float)accepted/(kk));
		printf("sigma from thread %d = %f\n",tid,sigma);



	}

}
