#pragma once

#include <string>
#include <iostream>
#include <mpi.h>
#include <cmath>
#include "ParallelPoissonSolver.h"
#include "matrix_print.h"
using namespace std;

class LidDrivenCavity
{
public:
    // Useage: LidDrivenCavity(domain_grid, grid_rank, neighbours, grid_points, dx, dy, dt, T, Re);
    LidDrivenCavity(MPI_Comm grid_comm, int rank, int neighbours[4], int grid_size[2], double dx, double dy, double dt, double T, double Re); //Constructor
    ~LidDrivenCavity(); //Destructor

    void SetDomainSize(double xlen, double ylen);
    void SetGridSize(int nx, int ny);
    void SetTimeStep(double deltat);
    void SetFinalTime(double finalt);
    void SetReynoldsNumber(double Re);

    void Initialise();
    void Integrate();
    void Test();
    void Gather();
    // Add any other public functions

private:
    void UpdateBoundaryConditions();
    void CalculateInteriorVorticity();
    void CalculateFutureInteriorVorticity();
    void UpdateSharedInterfaces(double *field);

    double dt;
    double dy;
    double dx;
    double T;
    int    Nx;
    int    Ny;
    double Lx;
    double Ly;
    double Re;
    int neighbours[4];
    int rank;
    MPI_Comm grid_comm;
    int grid_size[2];
    double* v = nullptr;
    double* s = nullptr;
    double* s_new = nullptr;
    double* v_new = nullptr;
    double* u_vel = nullptr;
    double* v_vel = nullptr;
    double U;

};
