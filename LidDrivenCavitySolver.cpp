#include <iostream>
using namespace std;

#include "LidDrivenCavity.h"

int main(int argc, char **argv)
{
    // Create a new instance of the LidDrivenCavity class
    LidDrivenCavity* solver = new LidDrivenCavity();

    // Configure the solver here...
    // ...
    solver->Initialise();

    // Run the solver
    solver->Integrate();

	return 0;
}
