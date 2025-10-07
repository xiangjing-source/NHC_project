#include "../include/globals.h"
#include "../include/atom.h"
#include "../include/force.h"
#include <vector>
#include <string>

using namespace std;

// ------------------- Global atom data -------------------
vector<Atom> atoms;
int N;

// ------------------- Simulation parameters -------------------
double dt = 1.0;
double T_target = 180.0;
double Tdamp = 5.0;
double mass = 39.948;
double epsilon = 0.234;
double sigma = 3.504;
double rc = 8.76;
double skin = 0.3;

int nsteps = 100000;
int thermo = 100;
int M = 3;

string units = "real";
double Lx=21.04, Ly=21.04, Lz=21.04;

string coords_file = "coords.data";
string output_file = "output.log";
string plot_file = "md_energy_temp.png";

// ------------------- NHC chain -------------------
vector<double> xi, eta, Q;

// ------------------- Physical constants -------------------
double kB = 0.0019872041;

// ------------------- Random seed -------------------
unsigned int rand_seed = 0;

// ------------------- Global neighbor list -------------------
NeighborList nlist;

// current simulation step (for diagnostics)
int current_step = 0;

