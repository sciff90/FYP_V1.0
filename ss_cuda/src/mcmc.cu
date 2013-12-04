#include<mcmc.h>

void mcmc(double *u,double *y,double *theta,long int N,int order, int num_samples,double *theta_0,double elim)
{

	double *d_u,*d_y,*d_y_test,*d_theta,*d_theta_0;

	int u_size = num_samples*sizeof(double);
	long int theta_size = N*2*(order+1)*sizeof(double);
	int theta_0_size = 2*(order+1)*sizeof(double);

	printf("N = %ld\n",N);
	printf("order = %d\n",order);
	printf("num_samples = %d\n",num_samples);
	printf("theta_size = %ld\n",theta_size);

	cudaMalloc((void**)&d_u, u_size ); 
	cudaMalloc((void**)&d_y, u_size ); 
	cudaMalloc((void**)&d_theta, theta_size );
	cudaMalloc((void**)&d_y_test, u_size );
	cudaMalloc((void**)&d_theta_0,theta_0_size);


	cudaMemcpy(d_u, u, u_size, cudaMemcpyHostToDevice );
	cudaMemcpy(d_y, y, u_size, cudaMemcpyHostToDevice );
	cudaMemcpy(d_theta_0,theta_0,theta_0_size,cudaMemcpyHostToDevice);

	slicesample<<<8,128>>>(d_u,d_y,d_theta,N,order,num_samples,d_theta_0,elim);

	cudaMemcpy( theta, d_theta, theta_size, cudaMemcpyDeviceToHost ); 


	cudaFree(d_u);
	cudaFree(d_y);
	cudaFree(d_theta);
	cudaFree(d_theta_0);
	cudaFree(d_y_test);

	cudaThreadExit();

	return ;
}
