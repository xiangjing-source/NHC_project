#pragma once
#include <vector>
using std::vector;

struct NeighborList {
    double rc;
    double skin;
    vector<vector<int>> neigh;

    NeighborList() : rc(0.0), skin(0.0) {}            
    NeighborList(double rc_, double skin_) : rc(rc_), skin(skin_) {}
    void build();
};

extern NeighborList nlist;

double compute_kinetic();
double compute_temperature(double ekin);
double compute_forces();
