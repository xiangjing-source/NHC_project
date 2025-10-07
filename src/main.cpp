#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <chrono>      
#include "../include/atom.h"
#include "../include/io.h"
#include "../include/force.h"
#include "../include/thermostat_nhc.h"
#include "../include/globals.h"

using namespace std;
using namespace std::chrono;

extern NeighborList nlist;

int main(int argc, char* argv[]){
    if(argc < 2){
        cerr << "Usage: ./md input.in" << endl;
        return 1;
    }

    auto start_time = high_resolution_clock::now();  // Start timing

    // ------------------- Read input -------------------
    read_input(argv[1]);           // Read input.in, update rc and skin
    read_xyz_extended(coords_file);
    read_box(coords_file);
    init_velocities();

    // ------------------- Print input parameters -------------------
    cout << "================= Input Parameters =================" << endl;
    cout << "units      = " << units << endl;
    cout << "dt         = " << dt << endl;
    cout << "steps      = " << nsteps << endl;
    cout << "thermo     = " << thermo << endl;
    cout << "rand_seed  = " << rand_seed << endl;
    cout << "T_target   = " << T_target << endl;
    cout << "nhc_chain  = " << M << endl;
    cout << "Tdamp      = " << Tdamp << endl;
    cout << "Lx, Ly, Lz = " << Lx << " " << Ly << " " << Lz << endl;
    cout << "epsilon    = " << epsilon << endl;
    cout << "sigma      = " << sigma << endl;
    cout << "rc         = " << rc << endl;
    cout << "skin       = " << skin << endl;
    cout << "mass       = " << mass << endl;
    cout << "coords     = " << coords_file << endl;
    cout << "output     = " << output_file << endl;
    cout << "plot       = " << plot_file << endl;
    cout << "N (atoms)  = " << N << endl;
    cout << "====================================================" << endl;

    if(N <= 0){
        cerr << "Error: N <= 0 after reading coords file!" << endl;
        return 1;
    }

    if(N > 0){
        cout << "First atom coords    = "
             << atoms[0].x[0] << " "
             << atoms[0].x[1] << " "
             << atoms[0].x[2] << endl;
        cout << "First atom velocity  = "
             << atoms[0].v[0] << " "
             << atoms[0].v[1] << " "
             << atoms[0].v[2] << endl;
    }

    // ------------------- Set neighbor list -------------------
    nlist.rc = rc;
    nlist.skin = skin;
    nlist.build();

    // ------------------- Initialize NHC -------------------
    xi.assign(M,0.0);
    eta.assign(M,0.0);
    Q.assign(M,0.0);
    int Nf = 3*N-3;
    for(int i=0;i<M;i++) Q[i]=kB*T_target*Tdamp*Tdamp;
    Q[0]=Nf*kB*T_target*Tdamp*Tdamp;

    // ------------------- Main loop -------------------
    double epot = compute_forces();
    double ekin = compute_kinetic();

    ofstream fout(output_file);
    fout << "# Step   E_pot   E_kin   E_tot   Temp\n";

    for(int step=0; step<=nsteps; step++){
        current_step = step;
    // Rebuild neighbor list every step (written as modulo 1 to make intent explicit)
        if(step%1==0) nlist.build();

        nhc_halfstep(ekin);

        for(int i=0;i<N;i++)
            for(int d=0;d<3;d++)
                atoms[i].v[d] += 0.5*dt*atoms[i].f[d]/mass;

        for(int i=0;i<N;i++)
            for(int d=0;d<3;d++){
                atoms[i].x[d] += dt*atoms[i].v[d];
                if(atoms[i].x[d] < 0) atoms[i].x[d] += (d==0?Lx:(d==1?Ly:Lz));
                if(atoms[i].x[d] > (d==0?Lx:(d==1?Ly:Lz))) atoms[i].x[d] -= (d==0?Lx:(d==1?Ly:Lz));
            }

        epot = compute_forces();

        for(int i=0;i<N;i++)
            for(int d=0;d<3;d++)
                atoms[i].v[d] += 0.5*dt*atoms[i].f[d]/mass;

        ekin = compute_kinetic();
        nhc_halfstep(ekin);

        double temp = compute_temperature(ekin);
        double etot = epot + ekin;

        if(step % thermo == 0){
            fout << setw(8) << step
                 << " " << setw(12) << epot
                 << " " << setw(12) << ekin
                 << " " << setw(12) << etot
                 << " " << setw(12) << temp
                 << endl;
        }
    }

    fout.close();

    // ------------------- Call Python for plotting -------------------
    string cmd = "python3 data.py ";
    cmd += argv[1];   // Pass input filename
    int ret = system(cmd.c_str());
    (void)ret;  // Ignore return value, do not output warning

    // ------------------- End timing -------------------
    auto end_time = high_resolution_clock::now();
    auto duration = duration_cast<seconds>(end_time - start_time);
    cout << "Simulation finished. Results saved in "
         << output_file << " and " << plot_file << endl;
    cout << "Total simulation time: " << duration.count() << " seconds" << endl;

    return 0;
}
