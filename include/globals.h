#pragma once
#include <vector>
#include <string>
#include "atom.h"
#include "force.h" 
using namespace std;

// simulation parameters
extern double dt, T_target, Tdamp, mass, epsilon, sigma, rc;
extern int nsteps, thermo, M;
extern string units;
extern double Lx, Ly, Lz;
extern string coords_file, output_file, plot_file;

// NHC chain
extern vector<double> xi, eta, Q;

// Boltzmann constant
extern double kB;

<<<<<<< HEAD
// Add random seed and neighbor list parameters
extern unsigned int rand_seed;  // Random seed
extern double skin;             // Neighbor list buffer
extern NeighborList nlist;
=======
// Add random seed and neighbor list parameters
extern unsigned int rand_seed;  // Random seed
extern double skin;             // Neighbor list buffer
extern NeighborList nlist;
// current simulation step (set in main loop) for diagnostics
extern int current_step;
>>>>>>> 62b2587 (Initial commit: English comments, cleaned workspace, ready for GitHub)
