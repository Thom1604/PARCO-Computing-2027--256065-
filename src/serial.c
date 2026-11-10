#include "heat.h"

int main(int argc, char **argv){
  //================= Getting simulation value =================
  int N=(argc>1)?atoi(argv[1]):DEFAULT_N;
  int max_iters=(argc>2)?atoi(argv[2]):DEFAULT_MAX_ITER;
  int tol==(argc>3)?atoi(argv[3]):DEFAULT_TOL;
  //================= Getting the  number of elements and the number of rows (the same as coloums) =================
  int total_N=N+2; // For practical reason we need to have a predefined border around the dimension (1 left 1 right) 
  size_t num_elements=(size_t)total_N*total_N;
  
  return 0;
}