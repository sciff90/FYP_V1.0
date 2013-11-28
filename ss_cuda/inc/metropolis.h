#include <thrust/random.h>
#include<iostream>
#include<stdio.h>
#include<filter.h>
#include<hash.h>
#include<p_ratio.h>

__global__
void metropolis(double *u,double *y,double *theta,long int N,int order,int num_samples,double *theta_0,double elim);
