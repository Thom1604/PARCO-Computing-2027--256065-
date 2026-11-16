#include "heat.h"

int main(int argc, char **argv){
  //================= Getting simulation value =================
  int N=(argc>1)?atoi(argv[1]):DEFAULT_N;
  int max_iters=(argc>2)?atoi(argv[2]):DEFAULT_MAX_ITER;
  double tol=(argc>3)?atof(argv[3]):DEFAULT_TOL; 
  
  //================= Getting the  number of elements and the number of rows (the same as coloums) =================
  int total_N=N+2; // For practical reason we need to have a predefined border around the dimension (1 left 1 right) 
  size_t num_elements=(size_t)total_N*total_N;
  double mem_size_mb=(double)(num_elements*sizeof(double)*2)/(1024.0*1024.0); // Calulation Memory Needed

  //================= Print Simulation Start Banner =================
  fprintf(stdout, "\n======================================================================\n");
  fprintf(stdout, "          SIMULAZIONE EQUAZIONE DEL CALORE 2D (JACOBI STENCIL)        \n");
  fprintf(stdout, "======================================================================\n");
  fprintf(stdout, " [*] Stato: Programma avviato correttamente.\n");
  fprintf(stdout, " [*] Parametri di Simulazione:\n");
  fprintf(stdout, "     - Dimensione Griglia (N x N)    : %d x %d (Totale con bordi: %d x %d)\n", N, N, total_N, total_N);
  fprintf(stdout, "     - Iterazioni Massime            : %d\n", max_iters);
  fprintf(stdout, "     - Tolleranza Convergenza        : %.1e\n", tol);
  fprintf(stdout, "     - Memoria Allocata (u + u_new)  : %.2f MB\n", mem_size_mb);
  fprintf(stdout, "----------------------------------------------------------------------\n");
  
  //================= Dynamic memory allocation ================= 
  double *u=(double *) calloc(num_elements,sizeof(double)); //Also with calloc i put all value of the matrix to 0.
  double *u_new=(double *) calloc(num_elements,sizeof(double));
  
  if(!u || !u_new){
    fprintf(stderr," [!] Memory allocation error\n");
    free(u);
    free(u_new);
    return 1;
  }
  
  //================= Initialisation (Putting borders at the Default Temperature of 100°C )============================
  for(int j=0;j<total_N;j++){
    u[INDEX(0, j, N)]               = DEFAULT_BOUNDARY_TEMP; // Top
    u[INDEX(total_N - 1, j, N)]     = DEFAULT_BOUNDARY_TEMP; // Bottom
    u[INDEX(j, 0, N)]               = DEFAULT_BOUNDARY_TEMP; // Left
    u[INDEX(j, total_N - 1, N)]     = DEFAULT_BOUNDARY_TEMP; // Right

    u_new[INDEX(0, j, N)]           = DEFAULT_BOUNDARY_TEMP;
    u_new[INDEX(total_N - 1, j, N)] = DEFAULT_BOUNDARY_TEMP;
    u_new[INDEX(j, 0, N)]           = DEFAULT_BOUNDARY_TEMP;
    u_new[INDEX(j, total_N - 1, N)] = DEFAULT_BOUNDARY_TEMP;
  }
  
  fprintf(stdout, " [*] Numeric calculation in progress...");
  fflush(stdout);
  
  //================= The Algoritm ============================
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
  double elapsed_time = now_sec() - start_time;

  fprintf(stdout, " Completed!\n");

  // Performance calculation
  double total_gflops = ((double)iters * N * N * 5.0) / 1e9;
  double performance  = total_gflops / elapsed_time;

  fprintf(stdout, "======================================================================\n");
  fprintf(stdout, "                         PERFORMANCE REPORT                           \n");
  fprintf(stdout, "======================================================================\n");
  fprintf(stdout, " %-30s : %s (Max diff: %.2e)\n", 
          "Convergence Status", converged ? "CONVERGED" : "NOT CONVERGED", max_diff);
  fprintf(stdout, " %-30s : %d\n", "Iterations Executed", iters);
  fprintf(stdout, " %-30s : %.4f s\n", "Execution Time", elapsed_time);
  fprintf(stdout, " %-30s : %.4f GFLOP\n", "Total Work", total_gflops);
  fprintf(stdout, " %-30s : %.3f GFLOP/s\n", "Performance", performance);
  fprintf(stdout, "======================================================================\n\n");
  
  //Deallocation
  free(u);
  free(u_new);
  
  return 0;
}