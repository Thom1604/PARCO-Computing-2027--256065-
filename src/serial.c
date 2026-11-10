#include "heat.h"

int main(int argc, char **argv){
  //================= Getting simulation value =================
  int N=(argc>1)?atoi(argv[1]):DEFAULT_N;
  int max_iters=(argc>2)?atoi(argv[2]):DEFAULT_MAX_ITER;
  int tol==(argc>3)?atoi(argv[3]):DEFAULT_TOL;
  
  //================= Getting the  number of elements and the number of rows (the same as coloums) =================
  int total_N=N+2; // For practical reason we need to have a predefined border around the dimension (1 left 1 right) 
  size_t num_elements=(size_t)total_N*total_N;
  
  //================= Dynamic memory allocation ================= 
  double *u=(double *) calloc(num_elements,sizeof(double)); //Also with calloc i put all value of the matrix to 0.
  double *u_new=(double *) calloc(num_elements,sizeof(double));
  
  if(!u || |u_new){
    fprintf(stderr,"Memory allocation error\n");
    free(u);
    free(u_new);
    return 1;
  }
  
  //================= Initialisation (Putting borders at the Default Temperature of 100°C )============================
  for(int j=0;j<total;j++){
    u[INDEX(0, j, N)]     = DEFAULT_BOUNDARY_TEMP;
    u_new[INDEX(0, j, N)] = DEFAULT_BOUNDARY_TEMP;
  }
  
  //================= The Algoritm ============================
  double start_time = now_sec();
  double max_diff=0.0;
  for(int iters=0;i<max_iters;i++){
    max_diff=0.0;
    for(int i=0;i<N;i++){
      for(int j=0;j<N;j++){
        u_new[INDEX(i,j,N)] =0.25(u[INDEX(i-1,j,N)]+u[INDEX(i+1,j,N)]+u[INDEX(i,j-1,N)]+u[INDEX(i,j+1,N)]); 
      }
    }
    double delta=fabs(u_new[INDEX(i,j,N)]-u[INDEX(i,j,N)]);
    if(delta>max_diff){
      max_diff=delta;
    }
  }
  
  return 0;
}