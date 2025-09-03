#include <vector>
#include <iostream>
#include <cmath>
#include "force.h"
#include "particles.h"
#include "kdtree.h" 
#include "density.h" 
#include "config.h"
#include "utils.h"

using namespace std; 

void GetAcc(Particles& pts, Tree& tree){
    int n = pts.endid; 
    const double alpha = Config::getInstance().get("alpha"); 
    const double beta = Config::getInstance().get("beta");

    for (int ia=0; ia<n; ia++){
        double acc1=0, acc2=0, acc3=0, accu=0;
        vector<int> ibs = tree.search(pts.x1[ia], pts.x2[ia], pts.x3[ia], 2*pts.len[ia]); 
        for (int ib: ibs){
            if (ia == ib){
                continue;
            }
            double qa=0, qb=0, wha=0, whb=0; 

            double dis_ab = pts.dist(ia, ib); 
            double inv_dis_ab = 1.0 / (dis_ab+1e-10); 
            vector<double> vr = pts.vecr(ia, ib); 
            vector<double> vv = pts.vecv(ia, ib); 
            double vdotr = dot(vr, vv); 
            vdotr = vdotr*inv_dis_ab; 

            // viscosity terms 
            if (vdotr<0){
                qa = -0.5 * pts.den[ia] * (alpha*pts.cs(ia) - beta*vdotr) * vdotr; 
                qb = -0.5 * pts.den[ib] * (alpha*pts.cs(ib) - beta*vdotr) * vdotr; 
            }
            
            // dv/dt = - Sig{ mb [ (Pa + qa)/ρa^2*div(Wa_ab) 
            //                   + (Pb + qb)/ρb^2*div(Wb_ab) ] }
            wha = wpq_m4(dis_ab, pts.len[ia]); // wha*dx_i/dist(a,b) = div(W)
            whb = wpq_m4(dis_ab, pts.len[ib]); 
            double term_a = ((pts.pres(ia)+qa)/(pts.den[ia]*pts.den[ia])) * wha; 
            double term_b = ((pts.pres(ib)+qb)/(pts.den[ib]*pts.den[ib])) * whb; 
            double term = - pts.mas[ib] * (term_a + term_b) * inv_dis_ab; 
            acc1 += term * vr[0]; 
            acc2 += term * vr[1]; 
            acc3 += term * vr[2]; 
            // du/dt = Sig[ mb (Pa + qa)/ρa^2 v_ab * div(W_ab) ]
            accu += pts.mas[ib] * term_a * vdotr;
        }
        pts.acc1[ia] = acc1; 
        pts.acc2[ia] = acc2; 
        pts.acc3[ia] = acc3; 
        pts.accu[ia] = accu; 
    } 
    return; 
}