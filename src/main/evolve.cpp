#include <vector>
#include "evolve.h"
#include "particles.h"
#include "kdtree.h" 
#include "density.h"
#include "force.h"

using namespace std;

void Derivs(Particles& pts, Tree& tree){
    const int n_smooth=3; 
    pts.set_boundary(); 
    tree = Tree(pts); 
    for (int i=0; i<n_smooth; i++){
        GetDensity(pts, tree); 
        Smoothing(pts); 
    }
    pts.set_boundary();
    GetAcc(pts, tree); 
    return; 
}

void Integral(Particles& pts, Tree& tree, double dt){
    int n = pts.endid; 
    for (int a=0; a<n; a++){
        pts.x1[a] += pts.vel1[a]*dt + 0.5*dt*dt*pts.acc1[a]; 
        pts.x2[a] += pts.vel2[a]*dt + 0.5*dt*dt*pts.acc2[a]; 
        pts.x3[a] += pts.vel3[a]*dt + 0.5*dt*dt*pts.acc3[a]; 
        pts.vel1[a] += 0.5*dt*pts.acc1[a]; 
        pts.vel2[a] += 0.5*dt*pts.acc2[a]; 
        pts.vel3[a] += 0.5*dt*pts.acc3[a]; 
        pts.ene[a] += 0.5*dt*pts.accu[a]; 
    } 
    Derivs(pts, tree); 
    for (int a=0; a<n; a++){
        pts.vel1[a] += 0.5*dt*pts.acc1[a]; 
        pts.vel2[a] += 0.5*dt*pts.acc2[a]; 
        pts.vel3[a] += 0.5*dt*pts.acc3[a]; 
        pts.ene[a] += 0.5*dt*pts.accu[a]; 
    } 
    return; 
}