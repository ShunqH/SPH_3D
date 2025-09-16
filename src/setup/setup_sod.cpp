#include <iostream>  
#include <vector>    
#include <stdexcept> 
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
    double x0 = Config::getInstance().get("x0");  
    double xmin = Config::getInstance().get("xmin");  
    double xmax = Config::getInstance().get("xmax"); 
    double ymin = Config::getInstance().get("ymin");  
    double ymax = Config::getInstance().get("ymax"); 
    double zmin = Config::getInstance().get("zmin");  
    double zmax = Config::getInstance().get("zmax"); 
    double rholeft = Config::getInstance().get("rholeft"); 
    double rhoright = Config::getInstance().get("rhoright"); 
    double pleft = Config::getInstance().get("pleft"); 
    double pright = Config::getInstance().get("pright"); 
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
            eos_parm = sqrt(pleft/rholeft); 
            break;
        default:
            throw invalid_argument("unknow equation of state."); 
            break;
    }
    pts.eos = set_eos(eos_type, eos_parm); 
    
    double dx = 1.0/nxunify;
    double dy = 1.0/nyunify;
    double dz = 1.0/nzunify;
    int nx = (x0-xmin)/dx;
    int ny = (ymax-ymin)/dy;
    int nz = (zmax-zmin)/dz;
    
    double m0 = rholeft*dx*dy*dz; 
    double dr = cbrt(dx*dy*dz); 
    // double dr = dy; 
    RandomGenerator ran(-0.05*dx, 0.05*dx); 
    for (int i = 0; i < nx; i++) {
        double xnow = xmin + (i + 0.5) * dx;
        for (int j = 0; j < ny; j++) {
            double ynow = ymin + (j + 0.5) * dy;
            for (int k = 0; k < nz; k++) {
                double znow = zmin + (k + 0.5) * dz;
                pt.x1 = xnow+ran(); 
                pt.x2 = ynow+ran(); 
                pt.x3 = znow+ran(); 
                pt.den = rholeft; 
                pt.mas = m0; 
                pt.vel1 = 0; 
                pt.vel2 = 0; 
                pt.vel3 = 0; 
                pt.len = hfact*dr ;     // 1.2*dxleft
                pt.ene = pleft/((gamma-1.)*rholeft) ;
                pts.add_particle(pt); 
            }
        }
    }

    dx = cbrt(rholeft/rhoright)*dx;
    dy = cbrt(rholeft/rhoright)*dy;
    dz = cbrt(rholeft/rhoright)*dz;
    dr = cbrt(rholeft/rhoright)*dr;
    nx = (xmax-x0)/dx; 
    ny = (ymax-ymin)/dy;
    nz = (zmax-zmin)/dz;
    // cout << nright << endl;
    for (int i = 0; i < nx; i++) {
        double xnow = 0 + (i + 0.5) * dx;
        for (int j = 0; j < ny; j++) {
            double ynow = ymin + (j + 0.5) * dy;
            for (int k = 0; k < nz; k++) {
                double znow = zmin + (k + 0.5) * dz;
                pt.x1 = xnow+ran(); 
                pt.x2 = ynow+ran(); 
                pt.x3 = znow+ran(); 
                pt.den = rhoright; 
                pt.mas = m0; 
                pt.vel1 = 0; 
                pt.vel2 = 0; 
                pt.vel3 = 0; 
                pt.len = hfact*dr ;     // 1.2*dxright
                pt.ene = pright/((gamma-1.)*rhoright) ;
                pts.add_particle(pt); 
            }
        }
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
