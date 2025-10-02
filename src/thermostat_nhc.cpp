#include "../include/thermostat_nhc.h"
#include "../include/force.h"
#include "../include/atom.h"
#include "../include/globals.h"
#include <cmath>

extern int N, M;
extern double dt, kB, T_target;
extern std::vector<double> xi, eta, Q;

void nhc_halfstep(double &ekin){
    double dt2 = dt * 0.5;
    double dt4 = dt2 * 0.5;
    double dt8 = dt4 * 0.5;

    double G = xi[M-2]*xi[M-2]/Q[M-2] - kB*T_target;
    xi[M-1] += dt4 * G;

    for(int m = M-2; m >= 0; m--){
        double tmp = exp(-dt8 * xi[m+1]/Q[m+1]);
        if(m==0) G = 2.0*ekin - (3.0*N-3)*kB*T_target;
        else     G = xi[m-1]*xi[m-1]/Q[m-1] - kB*T_target;
        xi[m] = tmp * (tmp * xi[m] + dt4*G);
    }
    for(int m = M-1; m >=0; m--) eta[m] += dt2*xi[m]/Q[m];

    double scale = exp(-dt2*xi[0]/Q[0]);
    for(int i=0;i<N;i++) for(int d=0;d<3;d++) atoms[i].v[d] *= scale;

    for(int m =0; m<M-1; m++){
        double tmp = exp(-dt8*xi[m+1]/Q[m+1]);
        if(m==0) G = 2.0*ekin*scale*scale - (3.0*N-3)*kB*T_target;
        else     G = xi[m-1]*xi[m-1]/Q[m-1] - kB*T_target;
        xi[m] = tmp*(tmp*xi[m]+dt4*G);
    }

    G = xi[M-2]*xi[M-2]/Q[M-2] - kB*T_target;
    xi[M-1] += dt4*G;

    ekin = compute_kinetic();
}

