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
    double rhohigh = Config::getInstance().get("rhohigh"); 
    double rholow = Config::getInstance().get("rholow"); 
    double phigh = Config::getInstance().get("phigh"); 
    double plow = Config::getInstance().get("plow"); 
    double hfact = Config::getInstance().get("hfact"); 
    int eos_type = Config::getInstance().get("EoS"); 

    const double gamma = 1.4; 
    Particle pt; 

    double eos_parm = 0; 
    switch (eos_type)
    {
        case 0:
            eos_parm = Config::getInstance().get("gamma"); 
            break; 
        case 1: 
            eos_parm = sqrt(phigh/rhohigh); 
            break; 
        default: 
            throw invalid_argument("unknow equation of state."); 
            break; 
    } 
    pts.eos = set_eos(eos_type, eos_parm); 

    double dx = 1.0/Nuniform;
    int ninside = 2*(r1/dx + 1); 
    double dV = dx*dx*dx; 
    double m0 = rhohigh*dV;
    RandomGenerator ran_in(-0.2*dx, 0.2*dx); 
    for (int i = 0; i < ninside; i++) {
        double xnow = -r1 + (i + 0.5) * dx;
        for (int j = 0; j < ninside; j++) {
            double ynow = -r1 + (j + 0.5) * dx;
            for (int k = 0; k < ninside; k++) {
                double znow = -r1 + (k + 0.5) * dx;
                if (xnow*xnow+ynow*ynow+znow*znow<=r1*r1) {
                    double theta = k*pi/4; 
                    pt.x1 = xnow*cos(theta)-ynow*sin(theta) + ran_in(); 
                    pt.x2 = xnow*sin(theta)+ynow*cos(theta) + ran_in(); 
                    pt.x3 = znow+ran_in(); 
                    pt.den = rhohigh; 
                    pt.mas = m0; 
                    pt.vel1 = 0; 
                    pt.vel2 = 0; 
                    pt.vel3 = 0; 
                    pt.len = hfact*dx ;     // 1.2*dxleft
                    pt.ene = phigh/((gamma-1.)*rhohigh) ;
                    pts.add_particle(pt); 
                }
            }
        }
    }

    dV = m0/rholow; 
    dx = cbrt(dV);
    int noutside = 2*(r2/dx + 1); 
    RandomGenerator ran_out(-0.2*dx, 0.2*dx); 
    for (int i = 0; i < noutside; i++) {
        double xnow = -r2 + (i + 0.5) * dx;
        for (int j = 0; j < noutside; j++) {
            double ynow = -r2 + (j + 0.5) * dx;
            for (int k = 0; k < noutside; k++) {
                double znow = -r2 + (k + 0.5) * dx;
                if (xnow*xnow+ynow*ynow+znow*znow<=r2*r2 &&
                    xnow*xnow+ynow*ynow+znow*znow>r1*r1) {
                    double theta = k*pi/4; 
                    pt.x1 = xnow*cos(theta)-ynow*sin(theta) + ran_out(); 
                    pt.x2 = xnow*sin(theta)+ynow*cos(theta) + ran_out(); 
                    pt.x3 = znow+ran_out(); 
                    pt.den = rholow; 
                    pt.mas = m0; 
                    pt.vel1 = 0; 
                    pt.vel2 = 0; 
                    pt.vel3 = 0; 
                    pt.len = hfact*dx ;     // 1.2*dxleft
                    pt.ene = plow/((gamma-1.)*rholow) ;
                    pts.add_particle(pt); 
                }
            }
        }
    }
    cout << "Particles initialized: " << pts.size() << endl;
    return; 
}