#include <cstdio>
#include <iostream>

#include <cstring>
#include <math.h>
#include <functional>
#include "pav.h"
#include "margin_pav.h"
#include "func_iterator.h"
#include "mex.h"
#include "helper.h"

/* solves min_{y,eps} 0.5*\|x-y\|_2^2 +0.5*\|eps-eps0\|_2^2
 *        s.t. Dy_j <=\eps_j, eps_j>=0, \sum_j eps_j=epsbar*/
extern "C" {
    void pav_double(double *in_array, size_t size, double *out_array) {
        pav::pav<std::less_equal, std::allocator>(in_array, in_array + size, out_array);
    }
    
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
        
    double *y, *eps, *z, *Jcol;
    
    const mwSize *dims;
    int n, d2;
    int c=0;
    const int verbose=0;
    
    
  
    //figure out dimensions vectors should be 1xn vextors
    dims = mxGetDimensions(prhs[0]);
    n = (int)dims[1];
    dims = mxGetDimensions(prhs[1]);
    if ((int)dims[1]==n)
        c=1;
    dims = mxGetDimensions(prhs[2]);
    d2 = (int)dims[1]-1;
    //mexPrintf("%d,%d,%d\n", n,d2,c);
    //associate outputs
    plhs[0] = mxCreateDoubleMatrix((int)dims[0],n,mxREAL);
    if (verbose>=1){
        mexPrintf("In\n");
        mexPrintf("%d %d\n", n, d2);
    }
    y = mxGetPr(prhs[0]);
    eps = mxGetPr(prhs[1]);
    Jcol = mxGetPr(prhs[2]);
    z = mxGetPr(plhs[0]);
    
    if (verbose==2){
        mexPrintf("y: ");
        print_const_size_array(y, n);
        mexPrintf("eps: n");
        print_const_size_array(eps, d2);
    }
    double v[n]; 
    double w[n]; 
    
    for(int j=0; j<d2; j++){
        int js=(int) Jcol[j];
        int nj=(int)Jcol[j+1]-js;           
        double norm_d_sq=nj*(nj*nj-1)/12.0;  
        
        for(int i=0; i<nj; i++)
            if (c)
                v[js+i] = eps[js+i];
            else
                v[js+i] = eps[j]*(-(nj-1)/2.0+i);
        
        double temp[nj];        
        vec_sub(&y[js], &v[js], nj, temp);
        pav_double(temp, nj, &w[js]);
    }
  
    if (verbose==2){
        mexPrintf("w: ");
        print_const_size_array(w, n);
    }
    
    if (verbose==2){
        mexPrintf("eps * v: ");
        print_const_size_array(v, n);
    }
    double temp[n];
    vec_add(w, v, n, z);
    if (verbose==2){
        mexPrintf("w + eps * v: ");
        print_const_size_array(z, n);
    }
    //std::memcpy(z, temp, sizeof w);
    //mexPrintf("chck");
    return;        
}