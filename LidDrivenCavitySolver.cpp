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
			("dt", po::value<double>()->default_value(0.01), "Time step process_count.")
			("T", po::value<double>()->default_value(10), "Final time.")
			("Re", po::value<double>()->default_value(1000), "Reynolds number.")

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
	int global_grid_points[dims], sub_domains[dims];

	// Assign domain parameters.
	cavity_size[0] = vm["Lx"].as<double>();
	cavity_size[1] = vm["Ly"].as<double>();
 	global_grid_points[0] = vm["Nx"].as<int>();
	global_grid_points[1] = vm["Ny"].as<int>();
	sub_domains[0] = vm["Px"].as<int>();
	sub_domains[1] = vm["Py"].as<int>();

	// Load sim parameters
	const double dt = vm["dt"].as<double>();
	const double T = vm["T"].as<double>();
	const double Re = vm["Re"].as<double>();
	const double dx = cavity_size[0]/(global_grid_points[0]-1);
	const double dy = cavity_size[1]/(global_grid_points[1]-1);

	// Validate dt
	if (dt > Re*dx*dy/4){
		cout << "dt too large for Re and number of points - change dt < Re*dx*dy/4" << endl;
		MPI_Finalize();
		return 1;
	}

	//Initialise MPI
	MPI_Init(&argc, &argv); //Init MPI

	//Get rank and overall process_count
	int global_rank, process_count;
	MPI_Comm_rank(MPI_COMM_WORLD, &global_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &process_count);

	//Display a message from each process
	//cout << "Yo! I am a process and my global rank is:" << global_rank <<" of the " << process_count << endl;

	// Sort out grid and parallelisms

	MPI_Comm domain_grid; //Define a communicator for the grid

	if ((sub_domains[0]*sub_domains[1]) != process_count) {
		if (global_rank==0)cout << "No. of processes does not match domain allocation (n != Py * Px)" << endl;
		MPI_Finalize();
		return 1;
	}

	int periods[dims] = {0,0}; //Dont allow any wrap-around
	int reorder = 1; //Allow re-ordering
	int coords[dims];
	int grid_rank;
	MPI_Cart_create(MPI_COMM_WORLD, dims, sub_domains, periods, reorder, &domain_grid); //Create a cartesian grid
	MPI_Comm_rank(domain_grid, &grid_rank); //Get rank in grid
	MPI_Cart_coords(domain_grid, grid_rank, dims, coords); //Get coordinated n grid

	//cout << "I am rank" << grid_rank <<" and I live at coords:" << coords[0] << coords[1] << endl;

	int neighbours[4] = {grid_rank, grid_rank, grid_rank, grid_rank};
  MPI_Cart_shift(domain_grid, 0, 1, &neighbours[0], &neighbours[1]);
  MPI_Cart_shift(domain_grid, 1, 1, &neighbours[2], &neighbours[3]);

	//for (int i = 0; i<4; i++) cout << neighbours[i] << endl;

	// Calculate gridding - tries to evenely distribute grid points.
	int grid_points[dims];

	// Get base points
	grid_points[0] = global_grid_points[0]/sub_domains[0];
	grid_points[1] = global_grid_points[1]/sub_domains[1];

	// Top up with any 'remainers'
	if (global_grid_points[0]%sub_domains[0] > grid_rank) grid_points[0]++;
	if (global_grid_points[1]%sub_domains[1] > grid_rank) grid_points[1]++;

	cout << "(Rank" << grid_rank << ") I have " << grid_points[0] << "x points and " << grid_points[1] << "y points" << ",my coords are:" << coords[0] << coords[1] << " my neighbours are:" << neighbours[0] << neighbours[1] << neighbours[2] << neighbours[3] << endl;

	// Create a new instance of the LidDrivenCavity class
	LidDrivenCavity* solver = new LidDrivenCavity(domain_grid, grid_rank, neighbours, grid_points, dx, dy, dt, T, Re);

	// Print test
	solver->Test();

	// Configure the solver here...
	solver->Initialise();

	 // Run the solver
	solver->Integrate();


	MPI_Finalize();
	return 0;
}
