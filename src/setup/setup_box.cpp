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

void Setup(Particles& pts){
    
    int nxunify = (int)Config::getInstance().get("nxunify"); 
    int nyunify = (int)Config::getInstance().get("nyunify"); 
    int nzunify = (int)Config::getInstance().get("nzunify"); 
    double xmin = Config::getInstance().get("xmin");  
    double xmax = Config::getInstance().get("xmax"); 
    double ymin = Config::getInstance().get("ymin");  
    double ymax = Config::getInstance().get("ymax"); 
    double zmin = Config::getInstance().get("zmin");  
    double zmax = Config::getInstance().get("zmax"); 
    double rho = Config::getInstance().get("rho"); 
    double press = Config::getInstance().get("p"); 
    double hfact = Config::getInstance().get("hfact"); 
    int eos_type = Config::getInstance().get("EoS"); 
    int bdx1l_type = Config::getInstance().get("X1BDLeft"); 
    int bdx1r_type = Config::getInstance().get("X1BDRight"); 
    int bdx2l_type = Config::getInstance().get("X2BDLeft"); 
    int bdx2r_type = Config::getInstance().get("X2BDRight"); 
    int bdx3l_type = Config::getInstance().get("X3BDLeft"); 
    int bdx3r_type = Config::getInstance().get("X3BDRight"); 

    // Binding equation of state
    double gamma = 1.4; 
    double eos_parm = 0; 
    switch (eos_type)
    {
        case 0:
            eos_parm = Config::getInstance().get("gamma"); 
            gamma = eos_parm;
            break;
        case 1:
            eos_parm = sqrt(press/rho); 
            break;
        default:
            throw invalid_argument("unknow equation of state."); 
            break;
    }
    pts.eos = set_eos(eos_type, eos_parm); 
    
    double dx = 1.0/nxunify;
    double dy = 1.0/nyunify;
    double dz = 1.0/nzunify;
    int nx = (xmax-xmin)/dx;
    int ny = (ymax-ymin)/dy;
    int nz = (zmax-zmin)/dz;

    Particle pt; 
    
    /************ unify grid particle setup ************/ 
    // double m0 = rho*dx*dy*dz; 
    // RandomGenerator ran(-0.05*dx, 0.05*dx); 
    // for (int i = 0; i < nx; i++) {
    //     double xnow = xmin + (i + 0.5) * dx;
    //     for (int j = 0; j < ny; j++) {
    //         double ynow = ymin + (j + 0.5) * dy;
    //         for (int k = 0; k < nz; k++) {
    //             double znow = zmin + (k + 0.5) * dz;
    //             pt.x1 = xnow+ran(); 
    //             pt.x2 = ynow+ran(); 
    //             pt.x3 = znow+ran(); 
    //             pt.den = rho; 
    //             pt.mas = m0; 
    //             pt.vel1 = 1.0*cos(pt.x1*4*pi/(xmax-xmin)); 
    //             pt.vel2 = 0; 
    //             pt.vel3 = 0; 
    //             pt.len = hfact*dx ;     // 1.2*dxleft
    //             pt.ene = press/((gamma-1.)*rho) ;
    //             pts.add_particle(pt); 
    //         }
    //     }
    // }

    /************ random particle setup ************/ 
    int n = nx*ny*nz; 
    double dV = (xmax-xmin)*(ymax-ymin)*(zmax-zmin)/n ;
    double dr = cbrt(dV); 
    double m0 = rho*dV; 
    RandomGenerator ranx(xmin, xmax); 
    RandomGenerator rany(ymin, ymax); 
    RandomGenerator ranz(zmin, zmax); 
    for (int i=0; i<n; i++){
        pt.x1 = ranx(); 
        pt.x2 = rany(); 
        pt.x3 = ranz(); 
        pt.den = rho; 
        pt.mas = m0; 
        pt.vel1 = 1.0 + 0.5*cos(pt.x1*4*pi/(xmax-xmin));
        pt.vel2 = 0; 
        pt.vel3 = 0; 
        pt.len = hfact*dr ;     // 1.2*dxleft
        pt.ene = press/((gamma-1.)*rho) ;
        pts.add_particle(pt);
    }

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
