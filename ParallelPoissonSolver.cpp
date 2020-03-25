#include "ParallelPoissonSolver.h"

ParallelPoissonSolver::ParallelPoissonSolver(double *v_initial, double dx, double dy, int Nx, int Ny, int rank, MPI_Comm grid_comm){
  //Constructing A Matrix, row-major
  // b = v = [v11, v12, v13 ... v1(Ny-1), v21, v22....]'

  A_width = (Nx - 2)*(Ny-2);
  A_total = A_width*A_width;
  A = new double[A_total];
  b = new double[A_width];
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

      //Generate b matrix
      b[i*(Ny-2)+j] = v_initial[(i+1)*(Ny)+(j+1)];
    }
  }
  cout << "Printing A Matrix:" << endl;
  printMatrixRM(A,A_width,A_width);

  // Generate b matrix:
  //b = v_initial;
  cout << "Printing B Vector:" << endl;
  printMatrixRM(b,A_width,1);

  x = new double[A_width];
}


ParallelPoissonSolver::~ParallelPoissonSolver(){
}

<<<<<<< HEAD
double* ParallelPoissonSolver::GetX(){
=======
double ParallelPoissonSolver::GetX(){
>>>>>>> 3cd712f99b7ccdeecc668b09844ef0aaa22e0ee9
  return x;

}
void ParallelPoissonSolver::Solve(){

  double* r = new double[A_width];
  double* p = new double[A_width];
  double* temp = new double[A_width];



  int k;
  double alpha;
  double beta;
  double eps;
  double tol = 0.00001;


  for (int i = 0;i<A_width;++i){x[i] = 0.1;}
  cout << "x vector: " << endl;
  printMatrixRM(x,A_width,1);

  cblas_dcopy(A_width, b, 1, r, 1);        // r_0 = b (i.e. b)
  cblas_dsymv(CblasRowMajor, CblasUpper, A_width, -1.0, A, A_width, x, 1, 1.0, r, 1);  // r_0 = b - A x_0
  cout << "r vector: " << endl;
  printMatrixRM(r,A_width,1);
  cblas_dcopy(A_width, r, 1, p, 1);        // p_0 = r_0

  cout << "r vector: " << endl;
  printMatrixRM(r,A_width,1);

  cout << "p vector: " << endl;
  printMatrixRM(p,A_width,1);
  k = 0;
  do {
      cblas_dsymv(CblasRowMajor, CblasUpper, A_width, 1.0, A, A_width, p, 1, 0.0, temp, 1); // temp= A p_k

      //cout << "t vector: " << endl;
      //printMatrixRM(temp,A_width,1);

      alpha = cblas_ddot(A_width, temp, 1, p, 1);  // alpha = p_k^T A p_k
      //cout << "alphav1: " << alpha << endl;
      alpha = cblas_ddot(A_width, r, 1, r, 1) / alpha; // compute alpha_k
      //cout << "alphav2: " << alpha << endl;
      beta  = cblas_ddot(A_width, r, 1, r, 1);  // r_k^T r_k

      cblas_daxpy(A_width, alpha, p, 1, x, 1);  // x_{k+1} = x_k + alpha_k p_k
      cblas_daxpy(A_width, -alpha, temp, 1, r, 1); // r_{k+1} = r_k - alpha_k A p_k

      eps = cblas_dnrm2(A_width, r, 1);
      cout << "Iteration " << k << ": eps=" << eps << endl;
      if (eps < tol*tol) {
          break;
      }
      beta = cblas_ddot(A_width, r, 1, r, 1) / beta;

      cblas_dcopy(A_width, r, 1, temp, 1);
      cblas_daxpy(A_width, beta, p, 1, temp, 1);
      cblas_dcopy(A_width, temp, 1, p, 1);

      k++;
    } while (k < 5000); // Set a maximum number of iterations

  //printMatrixCM(x,A_width,1);
}
