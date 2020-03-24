#include "ParallelPoissonSolver.h"

ParallelPoissonSolver::ParallelPoissonSolver(double *v_intial, double dx, double dy, int Nx, int Ny, int rank, MPI_Comm grid_comm){
  //Constructing A Matrix, row-major
  // b = v = [v11, v12, v13 ... v1(Ny-1), v21, v22....]'

  A_width = (Nx - 2)*(Ny-2);
  A_total = A_width*A_width;
  A = new double[A_total];

  double dx_sq = dx*dx;
  double dy_sq = dy*dy;

  double diag = 2.0/dx_sq + 2.0/dy_sq;
  double neg_recip_dx_sq = -1.0/dx_sq;
  double neg_recip_dy_sq = -1.0/dy_sq;

  for (int i = 0;i<(Nx-2);++i){
    for (int j = 0; j<(Ny-2); ++j){
      // i,j represent the coordinate of the stream function for which the row (field address) represents
      int field_address = i*(Ny-2)+j; //Will essentially just count up 0,1,2,3,4...
      int A_diag_address = field_address*(A_width + 1); //Will give addresses for diagonal
      // Add the diagonal values

      A[A_diag_address] = diag;

      // Add the horizontal components
      if (i > 0){
        //not on left column
        A[A_diag_address-(Ny-2)] = neg_recip_dx_sq;
      }

      if (i < (Nx-3)){
        //not on the right
        A[A_diag_address+(Ny-2)] = neg_recip_dx_sq;
      }

      if (j > 0){
        //not on the botton
        A[A_diag_address-1] = neg_recip_dy_sq;
      }

      if (j < (Ny-3)){
        //not on the top
        A[A_diag_address+1] = neg_recip_dy_sq;
      }
    }
  }
  printMatrixRM(A,A_width,A_width);
}


ParallelPoissonSolver::~ParallelPoissonSolver(){
}

double ParallelPoissonSolver::Solve(){
  double t = 3;
  return t;
}
