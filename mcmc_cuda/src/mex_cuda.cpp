/*********************************************************************
 * mex_cuda.cpp
 *
 * This file shows the basics of setting up a mex file to work with
 * Matlab. Linking matlab arrays into cpp for transfer onto cuda card.
 * 
 * Keep in mind:
 * <> Use 0-based indexing as always in C or C++
 * <> Indexing is column-based as in Matlab (not row-based as in C)
 * <> Use linear indexing.  [x*dimy+y] instead of [x][y]
 *
 *
 ********************************************************************/
#include <mex.h> 
#include <mcmc.cuh>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

//declare variables
    mxArray *u_m,*y_m,*theta_m,*theta_0_m;
    const mwSize *dims;
    int N,order,num_samples;
    double elim;
    double *u,*y,*theta,*theta_0;
       
//associate inputs
	
    u_m = mxDuplicateArray(prhs[0]);
    y_m = mxDuplicateArray(prhs[1]);  
    theta_0_m = mxDuplicateArray(prhs[4]);

    N = (int)mxGetScalar(prhs[2]);
    order = (int)mxGetScalar(prhs[3]);
    elim = (double)mxGetScalar(prhs[5]);

    dims = mxGetDimensions(prhs[0]);
    num_samples = (int)dims[1];  

//associate outputs
    theta_m = plhs[0] = mxCreateDoubleMatrix((mwSize)N,2*(order+1),mxREAL);    

//associate pointers
    u = mxGetPr(u_m);
    y = mxGetPr(y_m);
    theta = mxGetPr(theta_m);
    theta_0 = mxGetPr(theta_0_m);

    mcmc(u,y,theta,N,order,num_samples,theta_0,elim);  
//    mxFree(u_m);
//    mxFree(y_m);
//    mxFree(theta_m);
//    mxFree(theta_0_m);


    return;
}
