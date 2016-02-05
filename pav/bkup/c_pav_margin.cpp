#include <cstdio>
#include <iostream>

#include <cstring>
#include <math.h>
#include <functional>
#include "pav.h"
#include "margin_pav.h"
#include "func_iterator.h"
#include "mex.h"

extern "C" {
void pav_float(float *in_array, size_t size, float *out_array) {
  pav::pav<std::less_equal, std::allocator>(in_array, in_array + size, out_array);
  }

void pav_double(double *in_array, size_t size, double *out_array) {
  pav::pav<std::less_equal, std::allocator>(in_array, in_array + size, out_array);
  }

void lbound_pav_float(float* in_strt, size_t size, float* margins, float* out ) {
  pav::lbound_margin_pav<std::allocator>(in_strt, in_strt + size, out, margins); 
  }

void ubound_pav_float(float* in_strt, size_t size, float* margins, float* out ) {
  pav::ubound_margin_pav<std::allocator>(in_strt, in_strt + size, out, margins); 
  }

  void lbound_mmpav_float(float* in_strt, size_t size, float* out, const float* c ) {
    pav::lbound_maxmargin_pav<std::allocator>(in_strt, in_strt + size, out, c);
  }
  
 void lbound_mmpav_double(double* in_strt, size_t size, double* out, const double* c ) {
    pav::lbound_maxmargin_pav<std::allocator>(in_strt, in_strt + size, out, c);
  }
  
  void lbound_mmfpav_float(float* in_strt, size_t size, float* out) {
    pav::lbound_maxmargin_pav<std::allocator>(in_strt, in_strt + size, out,
       pav::FuncIterator<pav::Arith<float>,
                         pav::greater_<0> >(size, 0.0) );
  }
}


void print_const_size_array(double *a_out, size_t size) {
  for (int i = 0; i < size; ++i){
    mexPrintf("%f ", a_out[i]);    
  }
  mexPrintf("\n");
}


double* vec_add(double* v1, double* v2, int d){
double* outp = new double[d];

for(int i=0; i<d; i++){
outp[i] = v1[i] + v2[i];

}
return (outp);
}


double* vec_sub(double* v1, double* v2, int d){
double* outp = new double[d];

for(int i=0; i<d; i++){
outp[i] = v1[i] - v2[i];

}
return (outp);
}

double norm(double*v, int d){
 double nn = 0.0;
 for(int i=0; i<d; i++){
   nn = nn + v[i]*v[i];
 }

 nn = nn/d;
 nn = sqrt(nn);
 return(nn);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) { 
 mxArray *y_mat, *costs_mat, *winit_mat, *z_mat, *max_iters_mat;   
 double *y, *costs, *z, *winit, *maxiters;
 //float *y, *costs;
 //mxArray *z;   
 
 const mwSize *dims;   
 int dimx, dimy, numdims;      
 const double tolerance = 0.0001;

 //associate inputs   
 y_mat =     mxDuplicateArray(prhs[0]);   
 costs_mat = mxDuplicateArray(prhs[1]);      
 winit_mat =  mxDuplicateArray(prhs[2]);
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
mexPrintf("In\n");
mexPrintf("%d %d\n", dimx, dimy);

y = mxGetPr(y_mat);
costs = mxGetPr(costs_mat);
winit = mxGetPr(winit_mat);
maxiters = mxGetPr(max_iters_mat);
z = mxGetPr(z_mat);


double diff_d = 10.0;
double diff_w = 10.0;

double d_old[dimy];
double w_old[dimy];
double d[dimy];
double w[dimy];
double *y_minus_costs;


y_minus_costs = vec_sub(y, costs, dimy);
mexPrintf("maxiters: %f\n", maxiters[0]);
mexPrintf("winit: ");
print_const_size_array(winit, dimy);
mexPrintf("y: ");
print_const_size_array(y, dimy);
mexPrintf("costs: ");
print_const_size_array(costs, dimy);
mexPrintf("y_minus_costs: ");
print_const_size_array(y_minus_costs, dimy);


std::memcpy(w, winit, sizeof w);
std::memcpy(w, winit, sizeof w);

mexPrintf("w after copy: ");
print_const_size_array(w, dimy);

int iter=0;
while(diff_d > tolerance || diff_w > tolerance || iter <= 1){
 pav_double( vec_sub(y_minus_costs, w, dimy), dimy, d); 
 pav_double( vec_sub(y, d, dimy), dimy, w);
 
 diff_d = norm(vec_sub(d, d_old, dimy), dimy);
 diff_w = norm(vec_sub(w, w_old, dimy), dimy);

 mexPrintf("iter=%d diff_d=%f diff_w=%f\n", iter, diff_d, diff_w);
 iter = iter + 1; 
 std::memcpy(d_old, d, sizeof d);
 std::memcpy(w_old, w, sizeof w);
 if (iter>maxiters[0]){
     mexPrintf("Reached maximum iterations: %d\n",maxiters[0]);
     break;
 }

}
mexPrintf("d: ");
print_const_size_array(d, dimy);
mexPrintf("w ");
print_const_size_array(w, dimy);

std::memcpy(z, vec_add(d, w, dimy), sizeof d);
}


