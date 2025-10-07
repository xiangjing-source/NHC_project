#include "../include/thermostat_nhc.h"
#include "../include/force.h"
#include "../include/atom.h"
#include "../include/globals.h"
#include <cmath>
#include <iostream>
#include <fstream>

extern int N, M;
extern double dt, kB, T_target;
extern std::vector<double> xi, eta, Q;

void nhc_halfstep(double &ekin){
    double dt2 = dt * 0.5;
    double dt4 = dt2 * 0.5;
    double dt8 = dt4 * 0.5;

    double G = xi[M-2]*xi[M-2]/Q[M-2] - kB*T_target;
    xi[M-1] += dt4 * G;

    // --- runtime diagnostics: check for NaN/inf or extremely large xi/eta/Q ---
    // soft-protection parameters
    const double XI_SOFT_MAX = 1e5;   // soft clamp for xi
    const double EKIN_SOFT_MAX = 1e7; // kinetic energy threshold for recovery

    for(int m=0;m<M;m++){
        if(!std::isfinite(xi[m]) || !std::isfinite(eta[m]) || !std::isfinite(Q[m])){
            std::cerr << "NHC diagnostic: non-finite xi/eta/Q at step " << current_step << "\n";
            std::ofstream snap("crash_snapshot.xyz");
            snap << N << "\n" << "Step " << current_step << "\n";
            for(int i=0;i<N;i++) snap << "Ar " << atoms[i].x[0] << " " << atoms[i].x[1] << " " << atoms[i].x[2] << "\n";
            snap.close();
            // try to recover: reset chain variables and continue
            for(int ii=0; ii<M; ++ii){ xi[ii]=0.0; eta[ii]=0.0; }
            std::ofstream rlog("nhc_recovery.log", std::ios::app);
            rlog << "Recovery: non-finite xi/eta/Q at step "<<current_step<<"; chain reset.\n";
            rlog.close();
            break;
        }
        if(std::abs(xi[m]) > XI_SOFT_MAX){
            std::cerr << "NHC recovery: large xi["<<m<<"]="<<xi[m]<<" at step "<<current_step<<" - applying soft clamp\n";
            std::string fname = "nhc_diagnostic_" + std::to_string(current_step) + ".log";
            std::ofstream dout(fname);
            dout << "step "<<current_step<<" xi_index="<<m<<" xi_value="<<xi[m]<<"\n";
            dout << "xi: "; for(int ii=0;ii<M;ii++) dout<<xi[ii]<<" "; dout<<"\n";
            dout << "eta: "; for(int ii=0;ii<M;ii++) dout<<eta[ii]<<" "; dout<<"\n";
            dout << "Q: "; for(int ii=0;ii<M;ii++) dout<<Q[ii]<<" "; dout<<"\n";
            dout << "ekin_before="<<ekin<<"\n";
            // compute epot and min distance
            double epot_tmp = 0.0;
            double min_r2 = 1e12;
            for(int i=0;i<N;i++){
                for(int j=i+1;j<N;j++){
                    double dx = pbc(atoms[i].x[0]-atoms[j].x[0], Lx);
                    double dy = pbc(atoms[i].x[1]-atoms[j].x[1], Ly);
                    double dz = pbc(atoms[i].x[2]-atoms[j].x[2], Lz);
                    double r2 = dx*dx+dy*dy+dz*dz;
                    if(r2<min_r2) min_r2 = r2;
                    if(r2<rc*rc){
                        double r2i = 1.0/std::max(r2,1e-12);
                        double r6i = pow(sigma*sigma*r2i,3);
                        epot_tmp += 4*epsilon*r6i*(r6i-1);
                    }
                }
            }
            dout << "epot_approx="<<epot_tmp<<" min_r="<<sqrt(min_r2)<<"\n";
            // neighbor counts
            for(int i=0;i<N;i++) dout << "neigh_count["<<i<<"]="<<nlist.neigh[i].size()<<"\n";
            dout.close();

            // soft-clamp xi and reset eta to avoid runaway
            for(int ii=0; ii<M; ++ii){ xi[ii] = std::copysign(std::min(std::abs(xi[ii]), XI_SOFT_MAX), xi[ii]); eta[ii]=0.0; }
            std::ofstream rlog("nhc_recovery.log", std::ios::app);
            rlog << "Recovery at step "<<current_step<<": clamped xi to "<<XI_SOFT_MAX<<" and zeroed eta. ekin_before="<<ekin<<"\n";
            rlog.close();
            break; // continue integration after soft recovery
        }
    }

    for(int m = M-2; m >= 0; m--){
        double tmp = exp(-dt8 * xi[m+1]/Q[m+1]);
        if(m==0) G = 2.0*ekin - (3.0*N-3)*kB*T_target;
        else     G = xi[m-1]*xi[m-1]/Q[m-1] - kB*T_target;
        xi[m] = tmp * (tmp * xi[m] + dt4*G);
    }
    for(int m = M-1; m >=0; m--) eta[m] += dt2*xi[m]/Q[m];

    double scale = exp(-dt2*xi[0]/Q[0]);
        // strict per-step scale clamp to avoid instant huge velocity multipliers
        const double SCALE_MIN = 0.5; // do not shrink velocities by more than 2x in one half-step
        const double SCALE_MAX = 2.0; // do not expand velocities by more than 2x in one half-step
        if(!std::isfinite(scale) || std::abs(scale) > 1e6){
            std::cerr << "NHC diagnostic: bad scale="<<scale<<" at step "<<current_step<<" - applying recovery\n";
        std::ofstream snap("crash_snapshot.xyz");
        snap << N << "\n" << "Step " << current_step << "\n";
        for(int i=0;i<N;i++) snap << "Ar " << atoms[i].x[0] << " " << atoms[i].x[1] << " " << atoms[i].x[2] << "\n";
        snap.close();
        // reset chain and rescale velocities to target temperature
        for(int ii=0; ii<M; ++ii){ xi[ii]=0.0; eta[ii]=0.0; }
        double target_ekin = 0.5*(3.0*N-3)*kB*T_target;
        double cur_ekin = compute_kinetic();
        if(std::isfinite(cur_ekin) && cur_ekin > 0){
            double s = sqrt(std::max(1e-12, target_ekin / cur_ekin));
            for(int i=0;i<N;i++) for(int d=0; d<3; d++) atoms[i].v[d] *= s;
        }
            // clamp scale to [SCALE_MIN, SCALE_MAX]
            double orig_scale = scale;
            if(scale < SCALE_MIN) scale = SCALE_MIN;
            if(scale > SCALE_MAX) scale = SCALE_MAX;
            if(scale != orig_scale){
                std::ofstream rlog("nhc_recovery.log", std::ios::app);
                rlog << "Scale clamped at step "<<current_step<<": orig="<<orig_scale<<" clamped="<<scale<<"\n";
                rlog.close();
            }
    }
    for(int i=0;i<N;i++) for(int d=0;d<3;d++) atoms[i].v[d] *= scale;
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
    if(!std::isfinite(ekin) || ekin > EKIN_SOFT_MAX){
        std::cerr << "NHC recovery: bad ekin="<<ekin<<" at step "<<current_step<<" - rescaling to target temperature\n";
        std::ofstream rlog("nhc_recovery.log", std::ios::app);
        rlog << "Recovery: ekin="<<ekin<<" at step "<<current_step<<"; rescaling velocities and resetting chain.\n";
        rlog.close();
        // write snapshot for post-mortem
        std::ofstream snap("crash_snapshot.xyz");
        snap << N << "\n" << "Step " << current_step << "\n";
        for(int i=0;i<N;i++) snap << "Ar " << atoms[i].x[0] << " " << atoms[i].x[1] << " " << atoms[i].x[2] << "\n";
        snap.close();
        // reset thermostat chain
        for(int ii=0; ii<M; ++ii){ xi[ii]=0.0; eta[ii]=0.0; }
        // rescale velocities to target temperature
        double target_ekin = 0.5*(3.0*N-3)*kB*T_target;
        if(std::isfinite(ekin) && ekin > 0){
            double s = sqrt(std::max(1e-12, target_ekin / ekin));
            for(int i=0;i<N;i++) for(int d=0; d<3; d++) atoms[i].v[d] *= s;
        }
        ekin = compute_kinetic();
    }
}

