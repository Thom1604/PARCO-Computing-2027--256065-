#include "heat.h"

int main(int argc, char **argv){
  //================= Getting simulation value =================
  int N=(argc>1)?atoi(argv[1]):DEFAULT_N;
  int max_iters=(argc>2)?atoi(argv[2]):DEFAULT_MAX_ITER;
  double tol=(argc>3)?atof(argv[3]):DEFAULT_TOL; 
  int warm_up=(argc>4)?atoi(argv[4]):15;
  int reps=(argc>5)?atoi(argv[5]):30;
  
  //================= Getting the  number of elements and the number of rows (the same as coloums) =================
  int total_N=N+2; // For practical reason we need to have a predefined border around the dimension (1 left 1 right) 
  size_t num_elements=(size_t)total_N*total_N;
  double mem_size_mb=(double)(num_elements*sizeof(double)*2)/(1024.0*1024.0); // Calulation Memory Needed

  //================= Print Simulation Start Banner =================
  fprintf(stdout, "\n======================================================================\n");
  fprintf(stdout, "          2D HEAT EQUATION SIMULATION (JACOBI STENCIL)                \n");
  fprintf(stdout, "======================================================================\n");
  fprintf(stdout, " [*] Status: Program started successfully.\n");
  fprintf(stdout, " [*] Simulation Parameters:\n");
  fprintf(stdout, "     - Grid Size (N x N)           : %d x %d (Total with borders: %d x %d)\n", N, N, total_N, total_N);
  fprintf(stdout, "     - Max Iterations              : %d\n", max_iters);
  fprintf(stdout, "     - Convergence Tolerance       : %.1e\n", tol);
  fprintf(stdout, "     - Warm-up / Repetitions       : %d / %d\n", warm_up, reps);
  fprintf(stdout, "     - Allocated Memory (u + u_new): %.2f MB\n", mem_size_mb);
  fprintf(stdout, "----------------------------------------------------------------------\n");
  
  //================= Dynamic memory allocation ================= 
  double *u=(double *) malloc(num_elements*sizeof(double)); //Also with calloc i put all value of the matrix to 0.
  double *u_new=(double *) malloc(num_elements*sizeof(double));
  
  if(!u || !u_new){
    fprintf(stderr,"Memory allocation error\n");
    free(u);
    free(u_new);
    return 1;
  }
  
  fprintf(stdout, "Numeric calculation in progress...\n");
  fflush(stdout);
  
  //================= The Algoritm ============================
  for(int run=0;run<warm_up+reps;run++){
  
    //================= Initialisation (Putting borders at 100°C and interior at 0°C) ============================
    for(int i=0;i<total_N;i++){
      for(int j=0;j<total_N;j++){
        if(i==0||i==total_N-1||j==0||j==total_N-1){
          u[INDEX(i,j,N)]=DEFAULT_BOUNDARY_TEMP;
          u_new[INDEX(i,j,N)]=DEFAULT_BOUNDARY_TEMP;
        }else{
          u[INDEX(i,j,N)]=0.0;
          u_new[INDEX(i,j,N)]=0.0;
        }
      }
    }

    double start_time = now_sec();
    double max_diff=0.0;
    int iters=0;
    int converged=0; // Variabile to track convergence status
  
    for(iters=0;iters<max_iters;iters++){
      max_diff=0.0;
      for(int i=1;i<=N;i++){
        for(int j=1;j<=N;j++){
          u_new[INDEX(i,j,N)] =0.25*(u[INDEX(i-1,j,N)]+u[INDEX(i+1,j,N)]+u[INDEX(i,j-1,N)]+u[INDEX(i,j+1,N)]); 
          // I determine under the maximum error to be confronted outside the cycle with the tolerance
          double delta=fabs(u_new[INDEX(i,j,N)]-u[INDEX(i,j,N)]);
          if(delta>max_diff){
            max_diff=delta;
          }
        }
      }
      // Pointer changing 
      double *temp=u;
      u=u_new;
      u_new=temp;
      
      // Verify the convergence (if the value is under a tolerance it means that the system converges) 
      if(max_diff<tol){
        iters++; //I increment to get the last value update
        converged=1;
        break;
      }
    }

    // Calculation time of execution of the algoritm
    double elapsed_time=now_sec()-start_time;
    
    if(run>=warm_up){
      // Performance values
      fprintf(stdout, "Status: %s | Max diff: %.2e | Iters: %d | Time: %.6f s\n",converged?"CONVERGED":"NOT_CONVERGED",max_diff,iters,elapsed_time);
    }
  } 
  
  fprintf(stdout, "Numeric calculation finished...\n");
  fflush(stdout);
  
  //Deallocation
  free(u);
  free(u_new);
  
  return 0;
}