#include "../include/force.h"
#include "../include/globals.h"
#include "../include/atom.h"
#include <cmath>

// ------------------- Neighbor List -------------------
void NeighborList::build() {
    neigh.assign(N, vector<int>());
    for(int i=0;i<N;i++){
        for(int j=i+1;j<N;j++){
            double dx = pbc(atoms[i].x[0]-atoms[j].x[0], Lx);
            double dy = pbc(atoms[i].x[1]-atoms[j].x[1], Ly);
            double dz = pbc(atoms[i].x[2]-atoms[j].x[2], Lz);
            double r2 = dx*dx + dy*dy + dz*dz;
            if(r2 < (rc + skin)*(rc + skin)){
                neigh[i].push_back(j);
                neigh[j].push_back(i);
            }
        }
    }
}

// ------------------- 力计算 -------------------
double compute_forces() {
    for(int i=0;i<N;i++) for(int d=0;d<3;d++) atoms[i].f[d]=0;
    double epot=0;
    double shift = 4*epsilon*(pow(sigma/rc,12)-pow(sigma/rc,6));

    for(int i=0;i<N;i++){
        for(int j:nlist.neigh[i]){
            if(j>i){
                double dx = pbc(atoms[i].x[0]-atoms[j].x[0],Lx);
                double dy = pbc(atoms[i].x[1]-atoms[j].x[1],Ly);
                double dz = pbc(atoms[i].x[2]-atoms[j].x[2],Lz);
                double r2 = dx*dx + dy*dy + dz*dz;
                if(r2<rc*rc){
                    double r2i = 1.0/r2;
                    double r6i = pow(sigma*sigma*r2i,3);
                    double ff = 48*epsilon*r6i*(r6i-0.5)*r2i;
                    atoms[i].f[0]+=ff*dx; atoms[i].f[1]+=ff*dy; atoms[i].f[2]+=ff*dz;
                    atoms[j].f[0]-=ff*dx; atoms[j].f[1]-=ff*dy; atoms[j].f[2]-=ff*dz;
                    epot += 4*epsilon*r6i*(r6i-1) - shift;
                }
            }
        }
    }
    return epot;
}

// ------------------- 动能和温度 -------------------
double compute_kinetic() {
    double vcm[3]={0.0,0.0,0.0};
    for(int i=0;i<N;i++) for(int d=0;d<3;d++) vcm[d]+=atoms[i].v[d];
    for(int d=0;d<3;d++) vcm[d]/=N;

    double ekin=0.0;
    for(int i=0;i<N;i++){
        for(int d=0;d<3;d++){
            double vv = atoms[i].v[d]-vcm[d];
            ekin += 0.5*mass*vv*vv;
        }
    }
    return ekin;
}

double compute_temperature(double ekin){
    return (2.0*ekin)/((3.0*N-3)*kB);
}
