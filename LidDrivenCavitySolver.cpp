#include <iostream>
#include <mpi.h>
using namespace std;

#include "LidDrivenCavity.h"

int main(int argc, char* argv[])
{
	MPI_Init(&argc, &argv); //Init MPI

	//Get rank and overall size
	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	//Display a message from each process
	cout << "Yo! I am a process and my rank is:" << rank <<" of the " << size << endl;

	MPI_Comm domain_grid; //Define a communicator for the grid
	const int dims = 2; //Working on a 2D flow
	int sizes[dims] = {4,4}; //Needs to be changed by the arguments

	if ((sizes[0]*sizes[1]) != size) {
		if (rank==0)cout << "No. of processes does not match domain allocation (n != Ny * Nx)" << endl;
		return 1;
	}

	int periods[dims] = {0,0};
	int reorder = 1;
	int coords[dims];
	int grid_rank;
	MPI_Cart_create(MPI_COMM_WORLD, dims, sizes, periods, reorder, &domain_grid); //Create a cartesian grid
	MPI_Comm_rank(domain_grid, &grid_rank);
	MPI_Cart_coords(domain_grid, grid_rank, dims, coords);

	cout << "I live at coords:" << coords << endl;

   	 // Create a new instance of the LidDrivenCavity class
    	LidDrivenCavity* solver = new LidDrivenCavity();

    	// Configure the solver here...
    	solver->Initialise();

   	 // Run the solver
    	solver->Integrate();


	MPI_Finalize();
	return 0;
}
