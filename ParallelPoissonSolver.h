#pragma once

#include <string>
#include <iostream>
#include <mpi.h>
#include <cmath>
#include <iomanip>

using namespace std;

class ParallelPoissonSolver
{
public:
  ParallelPoissonSolver();
  ~ParallelPoissonSolver();
  double Solve();
private:


};
