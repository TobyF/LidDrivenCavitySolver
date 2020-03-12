#include <iostream>
#include <mpi.h>
using namespace std;

#include "LidDrivenCavity.h"

int main(int argc, char* argv[])
{   
	MPI_Init(&argc, &argv); //Init MPI
	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	

	cout << "Yo! I am a process and my rank is:" << rank << endl;
   	 // Create a new instance of the LidDrivenCavity class
    	LidDrivenCavity* solver = new LidDrivenCavity();

    	// Configure the solver here...
    	solver->Initialise();

   	 // Run the solver
    	solver->Integrate();


	MPI_Finalize();
	return 0;
}
