#include "../include/io.h"
#include "../include/atom.h"
#include "../include/globals.h"
<<<<<<< HEAD
=======
#include "../include/force.h"
>>>>>>> 62b2587 (Initial commit: English comments, cleaned workspace, ready for GitHub)
#include <fstream>
#include <sstream>
#include <random>
#include <chrono>
#include <iostream>
#include <cmath>

using namespace std;

<<<<<<< HEAD
// ------------------- 声明所有用到的全局变量 -------------------
=======
// ------------------- Declare all the global variables used -------------------
>>>>>>> 62b2587 (Initial commit: English comments, cleaned workspace, ready for GitHub)
extern int N, nsteps, thermo, M;
extern double dt, T_target, Tdamp, mass, Lx, Ly, Lz, epsilon, sigma, rc, kB;
extern string units, coords_file, output_file, plot_file;
extern vector<Atom> atoms;

<<<<<<< HEAD
// ------------------- 读取 input.in -------------------
=======
// ------------------- Read input.in -------------------
>>>>>>> 62b2587 (Initial commit: English comments, cleaned workspace, ready for GitHub)
void read_input(const string &filename) {
    ifstream fin(filename);
    if(!fin.is_open()) {
        cerr << "Error: cannot open input file " << filename << endl;
        exit(1);
    }

    string key;
    while(fin >> key) {
        if(key[0]=='#') { fin.ignore(1e6,'\n'); continue; }

        if(key=="units") fin >> units;
        else if(key=="dt") fin >> dt;
        else if(key=="steps") fin >> nsteps;
        else if(key=="thermo") fin >> thermo;
        else if(key=="T_target") fin >> T_target;
        else if(key=="nhc_chain") fin >> M;
        else if(key=="Tdamp") fin >> Tdamp;
        else if(key=="Lx") fin >> Lx;
        else if(key=="Ly") fin >> Ly;
        else if(key=="Lz") fin >> Lz;
        else if(key=="epsilon") fin >> epsilon;
        else if(key=="sigma") fin >> sigma;
        else if(key=="rc") fin >> rc;
        else if(key=="mass") fin >> mass;
        else if(key=="coords") fin >> coords_file;
        else if(key=="output") fin >> output_file;
        else if(key=="plot") fin >> plot_file;
        else if(key=="rand_seed") fin >> rand_seed;
        else if(key=="skin") fin >> skin;
        else { fin.ignore(1e6,'\n'); 
        }
    }

    fin.close();

<<<<<<< HEAD
    // 设置单位制
=======
    // Set unit system
>>>>>>> 62b2587 (Initial commit: English comments, cleaned workspace, ready for GitHub)
    if(units=="lj") kB = 1.0;
    else kB = 0.0019872041;
}

<<<<<<< HEAD
// ------------------- 从 coords.data 读取盒子尺寸 -------------------
=======
// ------------------- Read box size from coords.data -------------------
>>>>>>> 62b2587 (Initial commit: English comments, cleaned workspace, ready for GitHub)
void read_box(const string &filename) {
    ifstream fin(filename);
    if(!fin.is_open()) { 
        cerr << "Cannot open file " << filename << endl; 
        exit(1); 
    }

    string line;
<<<<<<< HEAD
    getline(fin, line); // 第一行: 原子数，跳过
    getline(fin, line); // 第二行: Lattice
    size_t pos1 = line.find("Lattice=\"");
    if(pos1 != string::npos){
        pos1 += 9; // 跳过 Lattice="
=======
    getline(fin, line); // First line: number of atoms, skip
    getline(fin, line); // Second line: Lattice
    size_t pos1 = line.find("Lattice=\"");
    if(pos1 != string::npos){
        pos1 += 9; // Skip Lattice="
>>>>>>> 62b2587 (Initial commit: English comments, cleaned workspace, ready for GitHub)
        size_t pos2 = line.find("\"", pos1);
        string lattice_str = line.substr(pos1, pos2 - pos1);
        double lx, m1, m2, m3, ly, m4, m5, m6, lz;
        istringstream iss(lattice_str);
        iss >> lx >> m1 >> m2 >> m3 >> ly >> m4 >> m5 >> m6 >> lz;
        Lx = lx;
        Ly = ly;
        Lz = lz;
    } else {
        cerr << "Cannot parse lattice information, using default box sizes\n";
    }

    fin.close();
}

<<<<<<< HEAD
// ------------------- 读取初始结构 -------------------
=======
// ------------------- Read initial structure -------------------
>>>>>>> 62b2587 (Initial commit: English comments, cleaned workspace, ready for GitHub)
void read_xyz_extended(const string &filename) {
    ifstream fin(filename);
    if(!fin.is_open()) { 
        cerr << "Cannot open coords file " << filename << endl; 
        exit(1); 
    }

    string line;
    getline(fin, line);
    N = stoi(line);
    if(N <= 0){
        cerr << "Error: number of atoms <=0 in coords file\n";
        exit(1);
    }

    getline(fin, line); // comment line

    atoms.resize(N);

    string sym;
    for(int i=0; i<N; i++){
        fin >> sym >> atoms[i].x[0] >> atoms[i].x[1] >> atoms[i].x[2];
        for(int d=0; d<3; d++){
            atoms[i].v[d] = 0.0;
            atoms[i].f[d] = 0.0;
        }
    }

    fin.close();
}

#include <random>
#include <chrono>
#include <iostream>
#include <cmath>       // sqrt
#include <algorithm>   // std::max

void init_velocities() {
<<<<<<< HEAD
=======
    // Use the seed from input file, if not set then use current time
>>>>>>> 62b2587 (Initial commit: English comments, cleaned workspace, ready for GitHub)
    unsigned seed = (rand_seed > 0) ? rand_seed : std::chrono::system_clock::now().time_since_epoch().count();
    std::mt19937 gen(seed);
    double stdv = sqrt(kB*T_target/mass);
    std::normal_distribution<double> dist(0.0, stdv);

    for(int i=0;i<N;i++)
        for(int d=0;d<3;d++)
            atoms[i].v[d] = dist(gen);

<<<<<<< HEAD
    // 去质心速度
=======
    // Remove center-of-mass velocity
>>>>>>> 62b2587 (Initial commit: English comments, cleaned workspace, ready for GitHub)
    double vcm[3]={0,0,0};
    for(int i=0;i<N;i++)
        for(int d=0;d<3;d++) vcm[d]+=atoms[i].v[d];
    for(int d=0;d<3;d++) vcm[d]/=N;
    for(int i=0;i<N;i++)
        for(int d=0;d<3;d++) atoms[i].v[d]-=vcm[d];

<<<<<<< HEAD
    // 打印一个原子速度模长最大值
=======
    // Print the maximum velocity magnitude of one atom
>>>>>>> 62b2587 (Initial commit: English comments, cleaned workspace, ready for GitHub)
    double vmax = 0;
    for(int i=0;i<N;i++){
        double vlen = sqrt(atoms[i].v[0]*atoms[i].v[0] +
                           atoms[i].v[1]*atoms[i].v[1] +
                           atoms[i].v[2]*atoms[i].v[2]);
        if(vlen>vmax) vmax = vlen;
    }
    std::cout << "Random seed = " << seed << ", max initial velocity = " << vmax << std::endl;
<<<<<<< HEAD
=======
    // --- Rescale velocities so that initial temperature matches T_target ---
    // compute current kinetic energy
    double ekin = compute_kinetic();
    double temp_now = compute_temperature(ekin);
    if(temp_now > 0){
        double scale = sqrt((T_target)/temp_now);
        for(int i=0;i<N;i++)
            for(int d=0;d<3;d++) atoms[i].v[d] *= scale;
        // recompute vmax for reporting
        double vmax2 = 0;
        for(int i=0;i<N;i++){
            double vlen = sqrt(atoms[i].v[0]*atoms[i].v[0] +
                               atoms[i].v[1]*atoms[i].v[1] +
                               atoms[i].v[2]*atoms[i].v[2]);
            if(vlen>vmax2) vmax2 = vlen;
        }
        std::cout << "After rescale, max initial velocity = " << vmax2 << ", initial temp = " << compute_temperature(compute_kinetic()) << std::endl;
    }
>>>>>>> 62b2587 (Initial commit: English comments, cleaned workspace, ready for GitHub)
}
