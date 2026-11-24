#define _MPI
#include "heat.h"

int main(int argc, char **argv){
  //================= MPI Inizialisation =================
  MPI_Init(&argc, &argv);
  int rank=0,size=0;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);

  //================= Getting simulation value =================
  int N=(argc>1)?atoi(argv[1]):DEFAULT_N;
  int max_iters=(argc>2)?atoi(argv[2]):DEFAULT_MAX_ITER;
  double tol=(argc>3)?atof(argv[3]):DEFAULT_TOL; 
  int warm_up=(argc>4)?atoi(argv[4]):15;
  int reps=(argc>5)?atoi(argv[5]):30;
  
  
  //================= Getting the  number of elements and the number of rows (the same as coloums) =================
  int total_N=N+2; // For practical reason we need to have a predefined border around the dimension (1 left 1 right)
  int local_N=(N/size)+(rank<(N%size)?1:0); 
  int local_rows=local_N+2; 
  size_t num_elements=(size_t)local_rows*total_N;
  double mem_size_mb=(double)(num_elements*sizeof(double)*2*size)/(1024.0*1024.0); // Calulation Memory Needed

  //================= Print Simulation Start Banner =================
  if(rank==0){
    fprintf(stdout, "\n======================================================================\n");
    fprintf(stdout, "          2D HEAT EQUATION SIMULATION (JACOBI STENCIL)                \n");
    fprintf(stdout, "======================================================================\n");
    fprintf(stdout, " [*] Status: Program started successfully (%d MPI processes).\n", size);
    fprintf(stdout, " [*] Simulation Parameters:\n");
    fprintf(stdout, "     - Grid Size (N x N)           : %d x %d (Total with borders: %d x %d)\n", N, N, total_N, total_N);
    fprintf(stdout, "     - Max Iterations              : %d\n", max_iters);
    fprintf(stdout, "     - Convergence Tolerance       : %.1e\n", tol);
    fprintf(stdout, "     - Allocated Memory (u + u_new): %.2f MB\n", mem_size_mb);
    fprintf(stdout, "----------------------------------------------------------------------\n");
  }
    
  //================= Dynamic memory allocation ================= 
  double *u=(double *) calloc(num_elements,sizeof(double)); //Also with calloc i put all value of the matrix to 0.
  double *u_new=(double *) calloc(num_elements,sizeof(double));
  
  if(!u || !u_new){
    fprintf(stderr," [!] Memory allocation error\n");
    free(u);
    free(u_new);
    MPI_Abort(MPI_COMM_WORLD,1);
    return 1;
  }
  
  // Neighbour control , I assign if there is a neighbour or if I'm on the edges (Halo exchange)
  int top_neighbor=(rank>0)?(rank-1):MPI_PROC_NULL;
  int bottom_neighbor=(rank<size-1)?(rank+1):MPI_PROC_NULL;
  if(rank==0){
    fprintf(stdout, " [*] Numeric calculation in progress...");
    fflush(stdout);
  }
    
  // Algoritm
  for(int run=0;run<warm_up+reps;run++){
    //================= Initialisation (Putting borders at the Default Temperature of 100°C )============================
    for(int i=0;i<local_rows;i++){
      for(int j=0;j<total_N;j++){
        if((rank==0 && i==0) || (rank==size-1 && i==local_rows-1) || (j==0) || (j==total_N-1)){
            u[INDEX(i,j,N)]=DEFAULT_BOUNDARY_TEMP;
            u_new[INDEX(i,j,N)]=DEFAULT_BOUNDARY_TEMP;
        }
        else{
          u[INDEX(i,j,N)]=0.0;
          u_new[INDEX(i,j,N)]=0.0;
        }
      }
    }
  
    MPI_Barrier(MPI_COMM_WORLD); // To synchronize all the process before the execution
    double start_time=now_sec();
    double max_diff=0.0;
    int iters=0;
    int converged=0;
  
    for(iters=0;iters<max_iters;iters++){
      // Halo exchange (Changing of the shared value on != rows)
      MPI_Sendrecv(&u[INDEX(1,0,N)],N+2,MPI_DOUBLE,top_neighbor,0,
        &u[INDEX(local_N+1,0,N)],N+2,MPI_DOUBLE,bottom_neighbor,0,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  
      MPI_Sendrecv(&u[INDEX(local_N, 0, N)],N+2,MPI_DOUBLE,bottom_neighbor,1,
        &u[INDEX(0,0,N)],N+2,MPI_DOUBLE,top_neighbor,1,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  
      // Calculation of the average of the 4 near neighbours
      double local_max_diff = 0.0;
      for (int i = 1; i <= local_N; i++) {
        for (int j = 1; j <= N; j++) {
          double val = 0.25 * (u[INDEX(i - 1, j, N)] + u[INDEX(i + 1, j, N)] +
                               u[INDEX(i, j - 1, N)] + u[INDEX(i, j + 1, N)]);
          u_new[INDEX(i, j, N)] = val;
  
          // I determine under the maximum error to be confronted outside the cycle with the tolerance
          double delta=fabs(u_new[INDEX(i,j,N)]-u[INDEX(i,j,N)]);
          if(delta>local_max_diff){
            local_max_diff=delta;
          }
        }
      }
  
      // Global maximum error reduction across all processes
      MPI_Allreduce(&local_max_diff, &max_diff, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  
      // Pointer changing 
      double *temp=u;
      u=u_new;
      u_new=temp;
      
      // Verify the convergence (if the value is under a tolerance it means that the system converges) Ma 
      if(max_diff<tol){
        iters++; //I increment to get the last value update
        converged=1;
        break;
      }
    }
    
    // Calculation time of execution of the algoritm
    MPI_Barrier(MPI_COMM_WORLD);
    double elapsed_time = now_sec() - start_time;
    if(run>=warm_up){  
      if (rank == 0) {
        fprintf(stdout, "Status: %s | Max diff: %.2e | Iters: %d | Time: %.6f s\n",converged?"CONVERGED":"NOT_CONVERGED",max_diff,iters,elapsed_time);
      }
    }
  }
  if(rank==0){
    fprintf(stdout, "Numeric calculation finished...\n");
    fflush(stdout);
  }
  
  //Deallocation
  free(u);
  free(u_new);
  
  //MPI Finalisation
  MPI_Finalize();
  
  return 0;
}