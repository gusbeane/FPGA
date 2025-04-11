#include <iostream>
#include <vector>
#include <iomanip>  // For std::setprecision
#include "CellState.hpp"   // Make sure these header paths match your project
#include "HydroSolver.hpp"

int main() {
    // Define the number of cells and cell width.
    const int Ncell = 100;
    const double dx = 1.0 / Ncell;

    // Create vectors for the conserved variables.
    // Each vector is filled with the same value for the test.
    std::vector<double> rho(Ncell, 1.0);
    std::vector<double> momentum(Ncell, 0.1);
    std::vector<double> energy(Ncell, 0.1);

    rho[0] = 100.0;  // Set the density in the first cell to a higher value.

    // Instantiate the cell state. The last parameter is gamma (default here is 5/3).
    CellState cellstate(Ncell, dx, rho, momentum, energy, 5.0/3.0);

    // Instantiate the hydro solver, passing a pointer to the cell state.
    HydroSolver hydrosolver(&cellstate);

    // Define the simulation end time and timestep
    double t_end = 1.0;
    double dt = 0.001;

    // Run the solver until time t_end with timestep dt.
    hydrosolver.run(t_end, dt);

    // Optionally, print some results from the simulation.
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "Final simulation time: " << hydrosolver.t << std::endl;
    std::cout << "Final density in cell 0: " << cellstate.rho[0] << std::endl;
    std::cout << "Final momentum in cell 0: " << cellstate.momentum[0] << std::endl;
    std::cout << "Final energy in cell 0: " << cellstate.energy[0] << std::endl;

    return 0;
}