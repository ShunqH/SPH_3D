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
    
    int Nuniform = (int)Config::getInstance().get("Nuniform"); 
    double r1 = Config::getInstance().get("r1");  
    double r2 = Config::getInstance().get("r2");  
    double rho = Config::getInstance().get("rho"); 
    double totale = Config::getInstance().get("totale"); 
    double elow = Config::getInstance().get("elow"); 
    double hfact = Config::getInstance().get("hfact"); 

    Particle pt; 

    double eos_parm = Config::getInstance().get("gamma"); 
    pts.eos = set_eos(0, eos_parm);             // adiabatic

    double dx = 1.0/Nuniform;
    int npt = 2*(r2/dx + 1); 
    double dV = dx*dx*dx; 
    double m0 = rho*dV;
    RandomGenerator ran(-0.05*dx, 0.05*dx); 
    int count = 0; 
    for (int i = 0; i < npt; i++) {
        double xnow = -r2 + (i + 0.5) * dx;
        for (int j = 0; j < npt; j++) {
            double ynow = -r2 + (j + 0.5) * dx;
            for (int k = 0; k < npt; k++) {
                double znow = -r2 + (k + 0.5) * dx;
                if (xnow*xnow+ynow*ynow+znow*znow<=r2*r2) {
                    pt.x1 = xnow+ran(); 
                    pt.x2 = ynow+ran(); 
                    pt.x3 = znow+ran(); 
                    pt.den = rho; 
                    pt.mas = m0; 
                    pt.vel1 = 0; 
                    pt.vel2 = 0; 
                    pt.vel3 = 0; 
                    pt.len = hfact*dx ;     // 1.2*dxleft
                    pt.ene = elow ;
                    if (xnow*xnow+ynow*ynow+znow*znow<=r1*r1) {
                        count++; 
                    }
                    pts.add_particle(pt); 
                }
                
            }
        }
    }
    double ehight = totale/count; 
    for (int i=0; i<pts.endid; i++){
        if (pts.x1[i]*pts.x1[i] + pts.x2[i]*pts.x2[i] + pts.x3[i]*pts.x3[i]<=r1*r1){
            pts.ene[i] = ehight; 
            cout<<ehight<<endl; 
        }
    }

    cout << "Particles initialized: " << pts.size() << endl;
    return; 
}