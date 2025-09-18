#include <iostream>  
#include <vector>   
#include <array> 
#include <stdexcept> 
#include "setup.h"
#include "config.h"
#include "random.h"
#include "utils.h"
#include "eos.h"

using namespace std;

// vector<vector<double>> Setup(int nleft, double x0, double xmin, double xmax, double rholeft, double rhoright, double pleft, double pright, double hfact, int nghost){
void Setup(Particles& pts){
    
    int N = Config::getInstance().get("N"); 
    double xmin = Config::getInstance().get("xmin");  
    double xmax = Config::getInstance().get("xmax"); 
    double ymin = Config::getInstance().get("ymin");  
    double ymax = Config::getInstance().get("ymax"); 
    double zmin = Config::getInstance().get("zmin");  
    double zmax = Config::getInstance().get("zmax"); 
    double r1 = Config::getInstance().get("r1");  
    double rho = Config::getInstance().get("rho"); 
    double totale = Config::getInstance().get("totale"); 
    double elow = Config::getInstance().get("elow"); 
    double hfact = Config::getInstance().get("hfact"); 
    int bdx1l_type = Config::getInstance().get("X1BDLeft"); 
    int bdx1r_type = Config::getInstance().get("X1BDRight"); 
    int bdx2l_type = Config::getInstance().get("X2BDLeft"); 
    int bdx2r_type = Config::getInstance().get("X2BDRight"); 
    int bdx3l_type = Config::getInstance().get("X3BDLeft"); 
    int bdx3r_type = Config::getInstance().get("X3BDRight"); 

    Particle pt; 

    double eos_parm = Config::getInstance().get("gamma"); 
    pts.eos = set_eos(0, eos_parm);             // adiabatic

    double V = (xmax - xmin) * (ymax - ymin) * (zmax - zmin);
    double m0 = rho*V/N; 
    double a = cbrt(4.0 * V / N);
    // RandomGenerator ranx(-0.05*dx, 0.05*dx); 
    // RandomGenerator rany(-0.05*dy, 0.05*dy); 
    // RandomGenerator ranz(-0.05*dz, 0.05*dz); 

    array<array<double,3>,4> basis = {{
        {0.0, 0.0, 0.0},
        {0.5, 0.5, 0.0},
        {0.5, 0.0, 0.5},
        {0.0, 0.5, 0.5}
    }};
    int count = 0; 
    for (int i = 0; i * a <= (xmax - xmin); i++) {
        for (int j = 0; j * a <= (ymax - ymin); j++) {
            for (int k = 0; k * a <= (zmax - zmin); k++) {
                double shift_x = xmin + i * a;
                double shift_y = ymin + j * a;
                double shift_z = zmin + k * a;
                for (auto &b : basis) {
                    double xnow = shift_x + b[0] * a;
                    double ynow = shift_y + b[1] * a;
                    double znow = shift_z + b[2] * a;
                    pt.x1 = xnow; 
                    pt.x2 = ynow; 
                    pt.x3 = znow; 
                    pt.den = rho; 
                    pt.mas = m0; 
                    pt.vel1 = 0; 
                    pt.vel2 = 0; 
                    pt.vel3 = 0; 
                    pt.len = hfact*a;     // 1.2*dxleft
                    pt.ene = elow;
                    if (pt.x1*pt.x1+pt.x2*pt.x2+pt.x3*pt.x3<=r1*r1) {
                            count++; 
                        }
                    pts.add_particle(pt); 
                }
            }
        }
    }

    if (count == 0) {
        throw runtime_error("No particles found inside r1, check your setup!");
    }
    double ehigh = totale/count; 
    for (int i=0; i<pts.endid; i++){
        if (pts.x1[i]*pts.x1[i] + pts.x2[i]*pts.x2[i] + pts.x3[i]*pts.x3[i]<=r1*r1){
            pts.ene[i] = ehigh/pts.mas[i]; 
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