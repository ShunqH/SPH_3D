#include <vector> 
#include "density.h" 
#include "utils.h"
#include "particles.h"
#include "kdtree.h"
#include "config.h"

using namespace std; 

double w_m4(double q, double h){
    double result = 1./pi; 
    if (q>=0 && q<1){
        result *= 1-1.5*q*q + 0.75*q*q*q ;
    }else if (q>=1 && q<2){
        result *= 0.25*(2-q)*(2-q)*(2-q); 
    }else{
        result = 0; 
    } 
    return result/(h*h*h); 
}

double wpq_m4(double dist, double h){
    double result = 1./pi; 
    double q = dist/h; 
    if (q>=0 && q<1){
        result *= -3.*q + 3*0.75*q*q;
    }else if (q>=1 && q<2){
        result *= -0.75*(2-q)*(2-q); 
    }else{
        result = 0; 
    } 
    return result/(h*h*h*h); 
}

double w_m6(double q, double h){
    double result = 1./(120*pi); 
    if (q>=0 && q<1){
        result *= power(3-q, 5) - 6*power(2-q, 5) + 15*power(1-q, 5); 
    }else if (q>=1 && q<2){
        result *= power(3-q, 5) - 6*power(2-q, 5); 
    }else if (q>=2 && q<3){
        result *= power(3-q, 5); 
    }else{
        result = 0; 
    } 
    return result/(h*h*h); 
}

void GetDensity(Particles& pts, Tree& tree){
    int n = pts.endid; 
    #pragma omp parallel for schedule(dynamic)
    for (int ia=0; ia<n; ia++){
        double dens = 0; 
        double h = pts.len[ia]; 
        vector<int> ibs = tree.search(pts.x1[ia], pts.x2[ia], pts.x3[ia], 2.0*h); 
        for (int ib : ibs){
            double q = pts.dist(ia, ib)/h; 
            dens += pts.mas[ib]*w_m4(q, h); 
        }
        pts.den[ia] = dens; 
    }
}

void Smoothing(Particles& pts){
    double hfact = Config::getInstance().get("hfact"); 
    for (int i=0; i<pts.endid; i++){
        double new_len = hfact*cbrt(pts.mas[i]/pts.den[i]);
        pts.len[i] = max(new_len, 1e-5);  
    }
} 
