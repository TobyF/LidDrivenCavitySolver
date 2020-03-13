#include <iostream>

//Multi-Process Interface
#include <mpi.h>

using namespace std;

//Allows program options to be read from cmd line
#include <boost/program_options.hpp>
namespace po = boost::program_options;

//User written class
#include "LidDrivenCavity.h"

int main(int argc, char* argv[])
{
	//Get command line arguments
	po::options_description desc("Allowed options");
  desc.add_options()
      ("help", "Produce help message")
			("Lx", po::value<int>(), "Length of the domain in the x-direction.")
			("Ly", po::value<int>(), "Length of the domain in the y-direction.")
      ("Nx", po::value<int>(), "Number of grid points in x-direction.")
			("Ny", po::value<int>(), "Number of grid points in y-direction.")
			("Px", po::value<int>(), "Number of partitions in the x-direction. (parallel)")
			("Py", po::value<int>(), "Number of partitions in the y-direction. (parallel)")
			("dt", po::value<int>(), "Time step size.")
			("T", po::value<int>(), "Final time.")
			("Re", po::value<int>(), "Reynolds number.")

  ;

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);

  if (vm.count("help")) {
      cout << desc << "\n";
      return 0;
  }

	//Initialise MPI
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
