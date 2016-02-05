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
    //FILE * cFile;
    //cFile = fopen ("PAVfile.txt","a");
    
    mxArray *y_mat, *c_mat, *winit_mat, *z_mat, *max_iters_mat;
    double *y, *c, *z, *winit, *maxiters;
    //float *y, *costs;
    //mxArray *z;
    
    const mwSize *dims;
    int dimx, dimy, numdims;
    const double tolerance = 0.0001;
    const int verbose=0;
    
    
    //associate inputs
    y_mat =     mxDuplicateArray(prhs[0]);   // input vector
    c_mat = mxDuplicateArray(prhs[1]);      // constant c (for tradeoff)
    winit_mat =  mxDuplicateArray(prhs[2]);  // initializing vector
    max_iters_mat =  mxDuplicateArray(prhs[3]);
    
// y = (float* ) mxGetData( mxDuplicateArray(prhs[0])) ;
// costs = (float* ) mxGetData(mxDuplicateArray(prhs[1]));
    
//figure out dimensions
    dims = mxGetDimensions(prhs[0]);
    numdims = mxGetNumberOfDimensions(prhs[0]);
    dimx = (int)dims[0];
    dimy = (int)dims[1];
    
//associate outputs
//z = plhs[0] = mxCreateNumericArray(2, dims, mxSINGLE_CLASS,  mxREAL);
    z_mat = plhs[0] = mxCreateDoubleMatrix(dimx,dimy,mxREAL);
    
//fprintf(cFile,"In\n");
//fprintf(cFile,"%d %d\n", dimx, dimy);
    if (verbose==2){
        mexPrintf("In\n");
        mexPrintf("%d %d\n", dimx, dimy);
    }
    y = mxGetPr(y_mat);
    c = mxGetPr(c_mat);
    winit = mxGetPr(winit_mat);
    maxiters = mxGetPr(max_iters_mat);
    z = mxGetPr(z_mat);
    
    
    double diff = 10.0;
    
    double d[dimy];
    double w[dimy];
    double w_old[dimy];
    
// construct vector d
    for(int i=0; i<dimy; i++){
        d[i] = i;
        w_old[i]=0;
        w[i]=winit[i];
    }
    
    double norm_d_sq = dotproduct(d, d, dimy);
    
    if (verbose==2){
        mexPrintf("maxiters: %f\n", maxiters[0]);
        mexPrintf("winit: ");
        print_const_size_array(winit, dimy);
        mexPrintf("y: ");
        print_const_size_array(y, dimy);
        mexPrintf("d: ");
        print_const_size_array(d, dimy);
        mexPrintf("c: %f\n", c[0]);
    }
    
    std::memcpy(w, winit, sizeof w);
    
    int iter=0;
    double eps_0=1e-5;
    double eps = eps_0;
    while(diff > tolerance || iter <= 1){
        eps = c[0] - dotproduct(d, vec_sub(w, y, dimy), dimy);
        eps = eps/norm_d_sq;
        if(eps <eps_0){eps =eps_0;}
        
        
        pav_double( vec_sub(y, vec_mulscalar(d, eps, dimy), dimy), dimy, w);
        
        diff = norm(vec_sub(w, w_old, dimy), dimy)/norm(w_old,dimy);
        if (verbose)
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
        mexPrintf("eps * d: ");
        print_const_size_array( vec_mulscalar(d, eps, dimy) , dimy);
        mexPrintf("w + eps * d: ");
        print_const_size_array(  vec_add(w,  vec_mulscalar(d, eps, dimy), dimy) , dimy);
    }
    std::memcpy(z, vec_add(w, vec_mulscalar(d, eps, dimy), dimy), sizeof w);
    
    
    mxDestroyArray(y_mat);
    mxDestroyArray(c_mat);
    mxDestroyArray(winit_mat);
    mxDestroyArray(max_iters_mat);
    
    return;
    

//fclose(cFile);
}