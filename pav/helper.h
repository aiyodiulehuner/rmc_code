void print_const_size_array(double *a_out, size_t size) {
  for (int i = 0; i < size; ++i)
    mexPrintf("%f ", a_out[i]);
  mexPrintf("\n");
}


double dotproduct(double *v1, double* v2, int d){
	double dp = 0.0;
	
	for(int i=0; i< d; i++){
		dp = dp + v1[i]*v2[i];
		
	}
	return(dp);
}

void vec_mulscalar(double* v1, double c, int d, double* outp){
//double* outp = new double[d];

for(int i=0; i<d; i++){
outp[i] = v1[i] * c;

}
//return (outp);
}



void vec_add(double* v1, double* v2, int d, double* outp){
//double* outp = new double[d];

for(int i=0; i<d; i++){
outp[i] = v1[i] + v2[i];

}
//return (outp);
}


void vec_sub(double* v1, double* v2, int d, double* outp){
//double* outp = new double[d];

for(int i=0; i<d; i++){
outp[i] = v1[i] - v2[i];

}
//return (outp);
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

