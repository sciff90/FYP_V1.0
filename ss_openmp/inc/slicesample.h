#include<random>
#include<iostream>
#include<filter.h>
#include<hash.h>
#include<p_ratio.h>

void slicesample(double *u,double *y,double *theta,long int N,int order,int num_samples,double *theta_0,double elim,int nthreads,int tid);
