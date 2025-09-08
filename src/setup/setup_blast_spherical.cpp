#include <iostream>  
#include <vector>    
#include <stdexcept> 
#include "setup.h"
#include "config.h"
#include "random.h"
#include "utils.h"
#include "eos.h"

using namespace std;

// vector<vector<double>> Setup(int nleft, double x0, double xmin, double xmax, double rholeft, double rhoright, double pleft, double pright, double hfact, int nghost){
void Setup(Particles& pts){
    
    int nr = (int)Config::getInstance().get("Nr"); 
    int nthe = (int)Config::getInstance().get("Ntheta"); 
    int nphi = (int)Config::getInstance().get("Nphi"); 
    double r1 = Config::getInstance().get("r1");  
    double r2 = Config::getInstance().get("r2");  
    double rho = Config::getInstance().get("rho"); 
    double totale = Config::getInstance().get("totale"); 
    double elow = Config::getInstance().get("elow"); 
    double hfact = Config::getInstance().get("hfact"); 

    Particle pt; 

    double eos_parm = Config::getInstance().get("gamma"); 
    pts.eos = set_eos(0, eos_parm);             // adiabatic

    int ntot = nr*nthe*nphi; 
    double dV = 4./3.*pi*r2*r2*r2/ntot;
    double dx = cbrt(dV); 
    double m0 = rho*dV;
    RandomGenerator ranr(-0.05/nr, 0.05/nr); 
    RandomGenerator ranthe(-0.05/nthe, 0.05/nthe); 
    RandomGenerator ranphi(-0.05/nphi, 0.05/nphi); 
    int count = 0; 
    for (int i = 0; i < nr; i++) {
        for (int j = 0; j < nthe; j++) {
            for (int k = 0; k < nphi; k++) {
                double rnow = r2*cbrt((i+0.5)/nr + ranr());
                double thetanow = acos(2*(j+0.5)/nthe-1 + ranthe());
                double phinow = 2*pi*(k+0.5)/nphi + ranphi();
                pt.x1 = rnow*sin(thetanow)*cos(phinow); 
                pt.x2 = rnow*sin(thetanow)*sin(phinow); 
                pt.x3 = rnow*cos(thetanow); 
                pt.den = rho; 
                pt.mas = m0; 
                pt.vel1 = 0; 
                pt.vel2 = 0; 
                pt.vel3 = 0; 
                pt.len = hfact*dx ;     // 1.2*dxleft
                pt.ene = elow ;
                if (rnow<=r1) {
                    count++; 
                }
                pts.add_particle(pt); 
            }
        }
    }
    double ehight = totale/count; 
    for (int i=0; i<pts.endid; i++){
        if (pts.x1[i]*pts.x1[i] + pts.x2[i]*pts.x2[i] + pts.x3[i]*pts.x3[i]<=r1*r1){
            pts.ene[i] = ehight/pts.mas[i]; 
        }
    }

    cout << "Particles initialized: " << pts.size() << endl;
    return; 
}