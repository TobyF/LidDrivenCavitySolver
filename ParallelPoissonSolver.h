#pragma once

#include <string>
#include <iostream>
#include <mpi.h>
#include <cmath>
#include <iomanip>
#include "matrix_print.h"
#include <cblas.h>

using namespace std;

class ParallelPoissonSolver
{
public:
  ParallelPoissonSolver(double *v_intial, double dx, double dy, int Nx, int Ny, int rank, MPI_Comm grid_comm);
  ~ParallelPoissonSolver();
  void Solve();
  double GetX();

private:
  double* A = nullptr;
  double* b = nullptr;
  double* x = nullptr;
  double diag;
  double dx;
  double dy;
  int Nx;
  int Ny;
  int rank;
  MPI_Comm grid_comm;
  int A_width;
  int A_total;



};
