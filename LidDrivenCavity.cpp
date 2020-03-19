#include "LidDrivenCavity.h"

LidDrivenCavity::LidDrivenCavity(MPI_Comm grid_comm, int rank, int neighbours[4], int grid_size[2], double dx, double dy, double dt, double T, double Re)
{
  this -> grid_comm = grid_comm;

  this -> rank = rank;

  this -> neighbours[0] = neighbours[0];
  this -> neighbours[1] = neighbours[1];
  this -> neighbours[2] = neighbours[2];
  this -> neighbours[3] = neighbours[3];

  this -> grid_size[0] = grid_size[0];
  this -> grid_size[1] = grid_size[1];

  this -> Nx = grid_size[0];
  this -> Ny = grid_size[1];

  this -> dx = dx;

  this -> dy = dy;

  this -> dt = dt;

  this -> T = T;

  this -> Re = Re;

  this -> U = 1;

}

LidDrivenCavity::~LidDrivenCavity()
{
}

void LidDrivenCavity::Test()
{
  cout << "I am a process, rank: " << rank << " and my neighbours are:" << neighbours[0] << neighbours[1] << neighbours[2] << neighbours[3] << endl;
}

void LidDrivenCavity::SetDomainSize(double xlen, double ylen)
{
}

void LidDrivenCavity::SetGridSize(int nx, int ny)
{
}

void LidDrivenCavity::SetTimeStep(double deltat)
{
}

void LidDrivenCavity::SetFinalTime(double finalt)
{
}

void LidDrivenCavity::SetReynoldsNumber(double re)
{
}

void LidDrivenCavity::Initialise()
{
  //Create flow grid of stream functions, vorticitys and velocitys. COULD BE MOVED TO constructor
  this -> s = new double[grid_size[0]*grid_size[1]];
  this -> s_new = new double[grid_size[0]*grid_size[1]];
  this -> v = new double[grid_size[0]*grid_size[1]];
  this -> v_new = new double[grid_size[0]*grid_size[1]];
  this -> u_vel = new double[grid_size[0]*grid_size[1]];
  this -> v_vel = new double[grid_size[0]*grid_size[1]];


  // Set boundary conditions:
  // If this set of points has a boundary on an edge
  // set the boundary points as per the first equation step.
  int i;
  if (neighbours[0] == -2){
    for (i=0;i<Ny;i++){
      v[i*Ny] = (s[i*Ny]-s[i*Ny+1])*2/(pow(dy,2));
    }
  }

}

void LidDrivenCavity::Integrate()
{
}
