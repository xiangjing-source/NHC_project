#pragma once
#include <vector>
using namespace std;

struct Atom {
    double x[3];
    double v[3];
    double f[3];
};

extern vector<Atom> atoms;
extern int N;

inline double pbc(double x, double L){
    while (x>0.5*L) x-=L;
    while (x<-0.5*L) x+=L;
    return x;
}

