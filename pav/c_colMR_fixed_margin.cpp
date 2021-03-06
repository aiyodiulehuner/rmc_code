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

/* solves min_{y} 0.5*\|x-y\|_2^2 +0.5*\|eps-eps0\|_2^2
 *        s.t. Dy_j <=\eps_j, eps_j>=0, \sum_j eps_j=epsbar*/
extern "C" {
    void pav_double(double *in_array, size_t size, double *out_array) {
        pav::pav<std::less_equal, std::allocator>(in_array, in_array + size, out_array);
    }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    double *y, *eps, *z, *Jcol;//, *ii_in, *ii_out;
    
    const mwSize *dims;
    int n, d2;
    int c=0;
    const int verbose=1;
    //figure out dimensions vectors should be 1xn vextors
    dims = mxGetDimensions(prhs[0]);
    n = (int) dims[1];
    dims = mxGetDimensions(prhs[1]);
    if ((int) dims[1]==n)
        c=1;
    dims = mxGetDimensions(prhs[2]);
    d2 = (int) dims[1]-1;
    //mexPrintf("%d,%d,%d\n", n,d2,c);
    //associate outputs
    plhs[0] = mxCreateDoubleMatrix(1,n,mxREAL);
    if (verbose>=1){       
        mexPrintf("\t\t PAV: n=%d d2=%d 1=%d c=%d. ", n, d2, dims[0],c);
    }

    y = mxGetPr(prhs[0]);
    eps = mxGetPr(prhs[1]);
    Jcol = mxGetPr(prhs[2]);
    //ii_in = mxGetPr(prhs[3]);
    z = mxGetPr(plhs[0]);
    //ii_out = mxGetPr(plhs[1]);
    if (verbose==2){
        mexPrintf("y: \n");
        print_const_size_array(y, n);
        mexPrintf("eps: \n");
        print_const_size_array(eps, d2);
    }
    double change=0; 
    for(int j=0; j<d2; j++){
        int js=(int) Jcol[j];
        int nj=(int) Jcol[j+1]-js;                    
        double w[nj]; 
        double temp[nj];         

        if (!c){
            double v[nj];
            for(int i=0; i<nj; i++)
                v[i] = eps[j]*(-(nj-1)/2.0+i);               
                vec_sub(&y[js], v, nj, temp);
                pav_double(temp, nj, w);
                vec_add(w, v, nj, &z[js]);}
        else{
            vec_sub(&y[js], &eps[js], nj, temp);          
            pav_double(temp, nj, w);            
            vec_add(w, &eps[js], nj, &z[js]);
            vec_sub(&y[js],&z[js],nj,temp);
            change+=dotproduct(temp,temp,nj);
        }
//        mxDestroyArray(w);
//        mxDestroyArray(temp);
    }
//    double temp[n];
//    vec_sub(y,z,n,temp);
//    mexPrintf("ch=%f\n",dotproduct(temp,temp,n));
    mexPrintf(" ch=%f ",change);
    return;        
}
