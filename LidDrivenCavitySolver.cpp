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
			("Lx", po::value<double>()->default_value(1), "Length of the domain in the x-direction.")
			("Ly", po::value<double>()->default_value(1), "Length of the domain in the y-direction.")
      ("Nx", po::value<int>()->default_value(16), "Number of grid points in x-direction.")
			("Ny", po::value<int>()->default_value(16), "Number of grid points in y-direction.")
			("Px", po::value<int>()->default_value(4), "Number of partitions in the x-direction. (parallel)")
			("Py", po::value<int>()->default_value(4), "Number of partitions in the y-direction. (parallel)")
			("dt", po::value<double>()->default_value(0.1), "Time step process_count.")
			("T", po::value<double>()->default_value(10), "Final time.")
			("Re", po::value<double>()->default_value(1), "Reynolds number.")

  ;

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);

  if (vm.count("help")) {
      cout << desc << "\n";
      return 0;
  }
	// Initialise domain parameters
	const int dims = 2; //Working on a 2D flow
	double cavity_size[dims];
	int grid_points[dims], sub_domains[dims];

	// Assign domain parameters.
	cavity_size[0] = vm["Lx"].as<double>();
	cavity_size[1] = vm["Ly"].as<double>();
 	grid_points[0] = vm["Nx"].as<int>();
	grid_points[1] = vm["Ny"].as<int>();
	sub_domains[0] = vm["Px"].as<int>();
	sub_domains[1] = vm["Py"].as<int>();

	// Load sim parameters
	const double dt = vm["dt"].as<double>();
	const double T = vm["T"].as<double>();
	const double Re = vm["Re"].as<double>();

	//Initialise MPI
	MPI_Init(&argc, &argv); //Init MPI

	//Get rank and overall process_count
	int rank, process_count;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &process_count);

	//Display a message from each process
	cout << "Yo! I am a process and my rank is:" << rank <<" of the " << process_count << endl;

	MPI_Comm domain_grid; //Define a communicator for the grid

	if ((sub_domains[0]*sub_domains[1]) != process_count) {
		if (rank==0)cout << "No. of processes does not match domain allocation (n != Py * Px)" << endl;
		return 1;
	}

	int periods[dims] = {0,0};
	int reorder = 1;
	int coords[dims];
	int grid_rank;
	MPI_Cart_create(MPI_COMM_WORLD, dims, sub_domains, periods, reorder, &domain_grid); //Create a cartesian grid
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
