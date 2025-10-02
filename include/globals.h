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

// 添加随机种子和邻居列表参数
extern unsigned int rand_seed;  // 随机种子
extern double skin;             // 邻居列表缓冲
extern NeighborList nlist;