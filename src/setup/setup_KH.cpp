#include <iostream>  
#include <vector>    
#include <stdexcept> 
#include <cmath>
#include "setup.h"
#include "config.h"
#include "random.h"
#include "utils.h"
#include "eos.h"

using namespace std;

double Coef_R(double y, double delta=0.25){
    double coef_f = 1./(1 + exp(2*(y+0.25)/delta)); 
    double coef_g = 1./(1 + exp(2*(0.25-y)/delta)); 
    double R = (1-coef_f)*(1-coef_g); 
    return R; 
}

double inverse_cdf(double u,
                   const std::vector<double>& xs,
                   const std::vector<double>& ys){
    auto it = std::lower_bound(ys.begin(), ys.end(), u);
    if (it == ys.begin()) return xs.front();
    if (it == ys.end())   return xs.back();
    int i = std::distance(ys.begin(), it);
    double x0 = xs[i-1], x1 = xs[i];
    double y0 = ys[i-1], y1 = ys[i];
    return x0 + (u-y0)*(x1-x0)/(y1-y0);
}

void Setup(Particles& pts){
    
    int nx = (int)Config::getInstance().get("nx"); 
    int ny = (int)Config::getInstance().get("ny"); 
    int nz = (int)Config::getInstance().get("nz"); 
    double xmin = Config::getInstance().get("xmin");  
    double xmax = Config::getInstance().get("xmax"); 
    double ymin = Config::getInstance().get("ymin");  
    double ymax = Config::getInstance().get("ymax"); 
    double zmin = Config::getInstance().get("zmin");  
    double zmax = Config::getInstance().get("zmax"); 
    double rho1 = Config::getInstance().get("rho1"); 
    double rho2 = Config::getInstance().get("rho2"); 
    double v1 = Config::getInstance().get("v1"); 
    double v2 = Config::getInstance().get("v2"); 
    double press = Config::getInstance().get("p"); 
    double hfact = Config::getInstance().get("hfact"); 
    int eos_type = Config::getInstance().get("EoS"); 
    int bdx1l_type = Config::getInstance().get("X1BDLeft"); 
    int bdx1r_type = Config::getInstance().get("X1BDRight"); 
    int bdx2l_type = Config::getInstance().get("X2BDLeft"); 
    int bdx2r_type = Config::getInstance().get("X2BDRight"); 
    int bdx3l_type = Config::getInstance().get("X3BDLeft"); 
    int bdx3r_type = Config::getInstance().get("X3BDRight"); 

    double gamma = 1.4; 
    Particle pt; 

    double eos_parm = 0; 
    switch (eos_type)
    {
        case 0:
            eos_parm = Config::getInstance().get("gamma"); 
            gamma = eos_parm;
            break;
        case 1:
            eos_parm = sqrt(press/rho1); 
            break;
        default:
            throw invalid_argument("unknow equation of state."); 
            break;
    }
    pts.eos = set_eos(eos_type, eos_parm); 

    double dx = 1.0/nx;
    // double dy = 1.0/ny;
    double dz = 1.0/nz;

    int Nygrid = 2000;
    vector<double> ylist(Nygrid), pdf(Nygrid), cdf(Nygrid);
    for (int i=0; i<Nygrid; i++){
        ylist[i] = ymin + (ymax-ymin) * i / (Nygrid-1);
        pdf[i]   = rho1 + Coef_R(ylist[i])*(rho2 - rho1);; 
    }
    double sum = 0.0;
    for (int i=0; i<Nygrid; i++){ sum += pdf[i]; cdf[i] = sum; }
    for (int i=0; i<Nygrid; i++){ cdf[i] /= sum; }

    double total_mass = 0.0;
    for (int i=0; i<Nygrid-1; i++){
        double dy = ylist[i+1]-ylist[i];
        total_mass += 0.5*(pdf[i]+pdf[i+1])*dy;
    }
    total_mass *= (xmax-xmin)*(zmax-zmin); 
    int Ntot = nx*ny*nz;
    double m0 = total_mass / Ntot;

    // RandomGenerator ran(-0.05*dx, 0.05*dx); 
    // int idx = 0; 
    for (int i = 0; i < nx; i++) {
        double xnow = xmin + (i + 0.5) * dx;
        for (int j = 0; j < ny; j++) {
            // double ynow = ymin + (j + 0.5) * dy;
            for (int k = 0; k < nz; k++) {
                double u = (j+0.5)/double(ny);
                double ynow = inverse_cdf(u, ylist, cdf); 
                double vel1now = v1 + Coef_R(ynow)*(v2-v1);
                double znow = zmin + (k + 0.5) * dz;
                pt.x1 = xnow; 
                pt.x2 = ynow; 
                pt.x3 = znow; 
                pt.den = rho1 + Coef_R(ynow)*(rho2-rho1); 
                pt.mas = m0; 
                double pervel = 0.1*sin(2*pi*pt.x1/(xmax-xmin)); 
                pt.vel1 = vel1now + pervel; 
                pt.vel2 = 0; 
                pt.vel3 = 0; 
                pt.len = hfact*dx ;     // 1.2*dxleft
                pt.ene = 1e-5 ;
                pts.add_particle(pt); 
                // idx++; 
            }
        }
    }

    // int n = nx*ny*nz; 
    // double dV = (xmax-xmin)*(ymax-ymin)*(zmax-zmin)/n ;
    // double dr = cbrt(dV); 
    // double m0 = rho*dV; 
    // RandomGenerator ranx(xmin, xmax); 
    // RandomGenerator rany(ymin, ymax); 
    // RandomGenerator ranz(zmin, zmax); 
    // // RandomGenerator ranvel(-0.1*v1, 0.1*v1); 
    // for (int i=0; i<n; i++){
    //     pt.x1 = ranx(); 
    //     pt.x2 = rany(); 
    //     pt.x3 = ranz(); 
        
    //     double pervel = 0.1*sin(2*pi*pt.x1/(xmax-xmin)); 
    //     if (fabs(pt.x2)<=0.25) {
    //         pt.den = rho; 
    //         pt.mas = m0; 
    //         pt.vel1 = v1+pervel; 
    //         pt.status = 1; 
    //     }else{
    //         pt.den = 2*rho; 
    //         pt.mas = 2*m0; 
    //         pt.vel1 = v2+pervel; 
    //         pt.status = 2; 
    //     }
    //     pt.vel2 = 0; 
    //     pt.vel3 = 0; 
    //     pt.len = hfact*dr ;     // 1.2*dxleft
    //     pt.ene = press/((gamma-1.)*rho) ;
    //     pts.add_particle(pt);
    // }

    cout << "Particles initialized: " << pts.size() << endl;

    // Binding the boundary conditions
    if (bdx1l_type != 0){
        double dx1l = Config::getInstance().get("dx1l"); 
        pts.BDX1L = Set_BDX1L(bdx1l_type, xmin, xmax, dx1l); 
    }
    if (bdx1r_type != 0){
        double dx1r = Config::getInstance().get("dx1r"); 
        pts.BDX1R = Set_BDX1R(bdx1r_type, xmin, xmax, dx1r); 
    }
    if (bdx2l_type != 0){
        double dx2l = Config::getInstance().get("dx2l"); 
        pts.BDX2L = Set_BDX2L(bdx2l_type, ymin, ymax, dx2l); 
    }
    if (bdx2r_type != 0){
        double dx2r = Config::getInstance().get("dx2r"); 
        pts.BDX2R = Set_BDX2R(bdx2r_type, ymin, ymax, dx2r); 
    }
    if (bdx3l_type != 0){
        double dx3l = Config::getInstance().get("dx3l"); 
        pts.BDX3L = Set_BDX3L(bdx3l_type, zmin, zmax, dx3l); 
    }
    if (bdx3r_type != 0){
        double dx3r = Config::getInstance().get("dx3r"); 
        pts.BDX3R = Set_BDX3R(bdx3r_type, zmin, zmax, dx3r); 
    }
    
    return; 
}


