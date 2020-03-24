#pragma once

#include <string>
#include <iostream>
#include <mpi.h>
#include <cmath>
#include <iomanip>
#include "matrix_print.h"

using namespace std;

class ParallelPoissonSolver
{
public:
  ParallelPoissonSolver();
  ~ParallelPoissonSolver();
  double Solve();
private:


};
