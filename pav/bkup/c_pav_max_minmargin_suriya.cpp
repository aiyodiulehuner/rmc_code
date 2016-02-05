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

extern "C" {
    void pav_double(double *in_array, size_t size, double *out_array) {
        pav::pav<std::less_equal, std::allocator>(in_array, in_array + size, out_array);
    }
    
}
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    double *y, *c, *z, *epsinit, *maxiters;
    
    const mwSize *dims;
    int dimx, dimy, numdims;
    const double tolerance = 0.00001;
    const int verbose=0;
    
    dims = mxGetDimensions(prhs[0]);
    numdims = mxGetNumberOfDimensions(prhs[0]);
    dimx = (int)dims[0];
    dimy = (int)dims[1];
    
    plhs[0]= mxCreateDoubleMatrix(dimx,dimy,mxREAL);   
    if (verbose==2){
        mexPrintf("In\n");
        mexPrintf("%d %d\n", dimx, dimy);
    }
    y = mxGetPr(prhs[0]);
    c = mxGetPr(prhs[1]);
    //zinit = mxGetPr(prhs[2]);
    epsinit = mxGetPr(prhs[2]);
    maxiters = mxGetPr(prhs[3]);
    z = mxGetPr(plhs[0]);
        
    double diff = 10.0;
    
    double d[dimy];
    double w[dimy];
    double w_old[dimy];
    double temp[dimy];
  
    int iter=0;
    double eps_0=1e-5;
    double eps = epsinit[0];
   
    // construct vector d
    for(int i=0; i<dimy; i++){
        d[i] = (-dimy/2)+i;
        w_old[i] = 0;
        //w[i]=zinit-epsinit*i;
    }
    
    double norm_d_sq = dotproduct(d, d, dimy);
    
    if (verbose==2){
        mexPrintf("maxiters: %f\n", maxiters[0]);
        mexPrintf("epsinit: %f\n", epsinit[0]);
        mexPrintf("y: ");
        print_const_size_array(y, dimy);
        mexPrintf("d: ");
        print_const_size_array(d, dimy);
        mexPrintf("c: %f\n", c[0]);
    }    

    while(diff > tolerance || iter <= 10){       
        vec_mulscalar(d, eps, dimy, temp);
        vec_sub(y, temp, dimy, temp);
        pav_double(temp , dimy, w);
        
        vec_sub(w, y, dimy, temp);
        eps = c[0] - dotproduct(d,temp , dimy);
        eps = eps/norm_d_sq;
        if(eps <eps_0){eps =eps_0;}
        
        vec_sub(w, w_old, dimy, temp);
        diff = norm(temp, dimy)/norm(w_old,dimy);
        
        if (verbose==2)
            mexPrintf("iter=%d diff=%f eps=%f\n", iter, diff, eps);
        iter = iter + 1;
        std::memcpy(w_old, w, sizeof w);
        if(iter >=maxiters[0]){
            //fprintf(cFile,"Max iters %f  reached! diff=%f eps=%f\n", maxiters[0], diff, eps);
            if (verbose>=1)
                mexPrintf("Max iters %f  reached! diff=%f eps=%f\n", maxiters[0], diff, eps);
            break;
        }
    }

    if ((iter<maxiters[0]) && (verbose>=1))
        //fprintf(cFile,"Exited in %d iterations. diff=%f eps=%f\n", iter,diff, eps);
        mexPrintf("Exited in %d iterations. diff=%f eps=%f\n", iter,diff, eps);
    
    if (verbose==2){
        mexPrintf("w: ");
        print_const_size_array(w, dimy);
    }
    
    vec_mulscalar(d, eps, dimy, temp);
    if (verbose==2){
        mexPrintf("eps * d: ");
        print_const_size_array( temp , dimy);
    }
    
    vec_add(w, temp, dimy, temp);
    if (verbose==2){
        mexPrintf("w + eps * d: ");
        print_const_size_array(temp, dimy);
    }
    
    std::memcpy(z, temp, sizeof w);
    //mexPrintf("chck");
    return;
//fclose(cFile);
}