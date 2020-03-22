#include "LidDrivenCavity.h"

LidDrivenCavity::LidDrivenCavity(MPI_Comm grid_comm, int rank, int neighbours[4], int grid_size[2], double dx, double dy, double dt, double T, double Re)
{
  this -> grid_comm = grid_comm;

  this -> rank = rank;

  this -> neighbours[0] = neighbours[0];
  this -> neighbours[1] = neighbours[1];
  this -> neighbours[2] = neighbours[2];
  this -> neighbours[3] = neighbours[3];

  this -> Nx = grid_size[0];
  this -> Ny = grid_size[1];


  //this -> Nx = grid_size[0];
  //this -> Ny = grid_size[1];

  // Alter grid size to account for a shared row/column with every non-boundary edge.
  if (neighbours[0] == -2 || neighbours[1] == -2) {Nx++;}
  else {Nx+=2;}

  if (neighbours[2] == -2 || neighbours[3] == -2) {Ny++;}
  else {Ny+=2;}


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
  //cout << "(Rank" << rank << ") I have " << Nx << "x points and " << Ny << "y points" << " my neighbours are:" << neighbours[0] << neighbours[1] << neighbours[2] << neighbours[3] << endl;

  //cout << "I am a process, rank: " << rank << " and my neighbours are:" << neighbours[0] << neighbours[1] << neighbours[2] << neighbours[3] << endl;
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
  this -> s = new double[Nx*Ny];
  this -> s_new = new double[Nx*Ny];
  this -> v = new double[Nx*Ny];
  this -> v_new = new double[Nx*Ny];
  this -> u_vel = new double[Nx*Ny];
  this -> v_vel = new double[Nx*Ny];
}

void LidDrivenCavity::UpdateBoundaryConditions()
{
  // Step 1 - Set boundary conditions:
  // If this set of points has a boundary on an edge
  // set the boundary points as per the first equation step.

  // Top
  int i;
  if (neighbours[0] == -2){ //Top boundary
    for (i=1;i<=Nx;i++){
      v[i*Ny-1] = (s[i*Ny-1]-s[i*Ny-2])*2/pow(dy,2) - 2*U/dy;
    }
  }
  // Bottom
  if (neighbours[1] == -2){
    for (i=0;i<Nx;i++){
      v[Ny*i] = (s[Ny*i] - s[Ny*i+1])*2/pow(dy,2);
    }
  }
  // Left
  if (neighbours[2] == -2){
    for (i=0;i<Ny;i++){
      v[i] = (s[i] - s[Ny+i])*2/pow(dx,2);
    }
  }
  // Right
  if (neighbours[2] == -2){
    for (i=0;i<Nx;i++){
      v[i] = (s[i] - s[Ny+i])*2/pow(dx,2);
    }
  }
}

void LidDrivenCavity::CalculateInteriorVorticity(){
  int i;
  int j;
  // i =1 --> Nx -1 is to only calculate interior points
  for (i=1;i<Nx-1;i++){
    for (j=1;j<Ny-1;j++){
      v[i*Ny + j] = (s[(i-1)*Ny + j] - 2*s[i*Ny + j] + s[(i+1)*Ny + j])/pow(dx,2) + (s[i*Ny + j - 1] - 2*s[i*Ny + j] + s[i*Ny + j + 1]);
    }
  }
}

void LidDrivenCavity::CalculateFutureInteriorVorticity(){
  //This is a long equation which has been split into 3 partitions
  //Details follow coursework handout

  double a;
  double b;
  double rhs;
  int i;
  int j;
  //For interior points
  for (i=1;i<Nx-1;i++){
    for (j=1;j<Ny-1;j++){
      a = (v[(i+1)*Ny +j]-v[(i-1)*Ny + j])*(s[i*Ny+j+1]-s[i*Ny+j-1])/(4*dx*dy);
      b = (s[(i+1)*Ny +j]-s[(i-1)*Ny + j])*(v[i*Ny+j+1]-v[i*Ny+j-1])/(4*dx*dy);
      rhs = (1/Re)*((v[(i+1)*Ny+j] - 2*v[i*Ny+j] + v[(i-1)*Ny+j])/pow(dx,2) + (v[i*Ny+j+1] - 2*v[i*Ny+j] + v[i*Ny+j-1])/pow(dy,2));

      v_new[i*Ny+j] = (rhs+b-a)*dt + v[i*Ny+j];

    }
  }
}

void LidDrivenCavity::UpdateSharedInterfaces(){

    int i;

    //Horizontal
    double send_left[Ny];
    double send_right[Ny];
    double recv_left[Ny];
    double recv_right[Ny];


    //Horizontal
    if ((neighbours[0] != -2) && (neighbours[1] != -2)){ //no boundary on left or right
      for (i=0;i<Ny;i++){
        send_left[i] = v[Ny + i];
        send_right[i] = v[(Nx-2)*Ny + i];
      }
      cout << rank << " send/recv data flowing left <---" << endl;
      MPI_Sendrecv(&send_left, Ny, MPI_DOUBLE, neighbours[0], 0,
                &recv_right, Ny , MPI_DOUBLE, neighbours[1], MPI_ANY_TAG,
                grid_comm,  MPI_STATUS_IGNORE);

      cout << rank << " send/recv data flowing left --->" << endl;
      MPI_Sendrecv(&send_right, Ny, MPI_DOUBLE, neighbours[1], 1,
                &recv_left, Ny , MPI_DOUBLE, neighbours[0], MPI_ANY_TAG,
                grid_comm,  MPI_STATUS_IGNORE);
    }
    else if ((neighbours[0] != -2)){ //attached to right wall and left dree
      cout << rank << " send data flowing left <---" << endl;
      MPI_Send(&send_left, Ny, MPI_DOUBLE, neighbours[0],0,grid_comm);
      cout << rank << " recv data flowing right --->" << endl;
      MPI_Recv(&recv_left, Ny, MPI_DOUBLE, neighbours[0],MPI_ANY_TAG,grid_comm, MPI_STATUS_IGNORE);

    }
    else{ //Attached to left wall
      cout << rank << " recv data flowing left <---" << endl;
      MPI_Recv(&recv_right, Ny, MPI_DOUBLE, neighbours[1],MPI_ANY_TAG,grid_comm, MPI_STATUS_IGNORE);
      cout << rank << " send data flowing right --->" << endl;
      MPI_Send(&send_right, Ny, MPI_DOUBLE, neighbours[1],1,grid_comm);
    }


    //Unpack any recvievd files, updating the stored arrays
    for (i=0;i<Ny;i++){
        if (neighbours[0] != -2) { v[i] = recv_left[i];}
        if (neighbours[1] != -2) {v[(Nx-1)*Ny + i] = recv_right[i];}
    }
    cout << rank << " Unpacked horizontal data" << endl;

    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0){cout << "----Vertical----" << endl;}
    //Vertical
    double send_up[Ny];
    double send_down[Ny];
    double recv_up[Ny];
    double recv_down[Ny];

    if ((neighbours[2] != -2) && (neighbours[3] != -2)){ //no boundary on left or right
      for (i=0;i<Nx;i++){
        send_down[i] = v[Ny*i+1]; //Second from bottom row
        send_up[i] = v[Ny*(i+1)-2]; //Second to top row
      }
      cout << rank << "Yep im here - 1" << endl;
      MPI_Sendrecv(&send_down, Nx, MPI_DOUBLE, neighbours[2], 2,
                &recv_up, Nx , MPI_DOUBLE, neighbours[3], MPI_ANY_TAG,
                grid_comm,  MPI_STATUS_IGNORE);
      cout << rank << "Now im here - 2" << endl;

      MPI_Sendrecv(&send_up, Nx, MPI_DOUBLE, neighbours[3], 3,
                &recv_down, Nx , MPI_DOUBLE, neighbours[2], MPI_ANY_TAG,
                grid_comm,  MPI_STATUS_IGNORE);

      cout << rank << "Now im finished" << endl;
    }
    else if ((neighbours[2] != -2)){ //attached to top lid
      MPI_Send(&send_down, Nx, MPI_DOUBLE, neighbours[2],2,grid_comm);
      MPI_Recv(&recv_down, Nx, MPI_DOUBLE, neighbours[2],MPI_ANY_TAG,grid_comm, MPI_STATUS_IGNORE);

    }
    else{ //Attached to bottom wall
      MPI_Recv(&recv_up, Nx, MPI_DOUBLE, neighbours[3],MPI_ANY_TAG,grid_comm, MPI_STATUS_IGNORE);
      MPI_Send(&send_up, Nx, MPI_DOUBLE, neighbours[3], 3,grid_comm);
    }

    //Unpack any recvievd files, updating the stored arrays
    for (i=0;i<Nx;i++){
        if (neighbours[2] != -2) { v[Ny*i] = recv_down[i];}
        if (neighbours[3] != -2) {v[Ny*(i+1)-1] = recv_up[i];}
    }
    cout << rank << " Unpacked Vertical Data" << endl;
    //MPI_Barrier(grid_comm);

}

void LidDrivenCavity::PoissonSolver(){

}

void LidDrivenCavity::Integrate()
{
  if (rank == 0){cout << "Integrating..." << endl;}
  //Sets the boundary values of vorticity based on current stream function (if attached to wall
  UpdateBoundaryConditions();
  MPI_Barrier(grid_comm);

  // Sets the interior vorticity values.
  if (rank == 0){cout << "Step 2: Calculating Interior Vorticity" << endl;}
  CalculateInteriorVorticity();
  MPI_Barrier(grid_comm);

  // GET VORTICITY AT DOMAIN INTERFACES
  if (rank == 0){cout << "Update the shared boundary points" << endl;}
  UpdateSharedInterfaces();
  MPI_Barrier(grid_comm);

  if (rank == 0){cout << "Step 3: Calculate Future Interior Vorticity" << endl;}
  CalculateFutureInteriorVorticity();
  MPI_Barrier(grid_comm);

  if (rank == 0){cout << "Step 4: Solve for stream functions" << endl;}
  PoissonSolver();


}
