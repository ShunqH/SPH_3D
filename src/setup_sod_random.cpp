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
    
    int nleft = (int)Config::getInstance().get("nleft"); 
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
    // int nghost = (int)Config::getInstance().get("nghost"); 

    const double gamma = 1.4; 
    Particle pt; 

    double eos_parm = 0; 
    switch (eos_type)
    {
        case 0:
            eos_parm = Config::getInstance().get("gamma"); 
            break;
        case 1:
            eos_parm = sqrt(pleft/rholeft); 
            break;
        default:
            throw invalid_argument("unknow equation of state."); 
            break;
    }
    pts.eos = set_eos(eos_type, eos_parm); 

    double dx = 1.0/nleft;
    double dV = dx*dx*dx; 
    double m0 = rholeft*dV; 
    RandomGenerator ran(-0.5*dx, 0.5*dx); 
    for (int i = 0; i < nleft; i++) {
        double xnow = xmin + (i + 0.5) * dx;
        for (int j = 0; j < nleft; j++) {
            double ynow = ymin + (j + 0.5) * dx;
            for (int k = 0; k < nleft; k++) {
                double znow = zmin + (k + 0.5) * dx;
                pt.x1 = xnow+ran(); 
                pt.x2 = ynow+ran(); 
                pt.x3 = znow+ran(); 
                pt.den = rholeft; 
                pt.mas = m0; 
                pt.vel1 = 0; 
                pt.vel2 = 0; 
                pt.vel3 = 0; 
                pt.len = hfact*dx ;     // 1.2*dxleft
                pt.ene = pleft/((gamma-1.)*rholeft) ;
                pts.add_particle(pt); 
            }
        }
    }

    dV = m0/rhoright; 
    dx = cbrt(dV);
    int nright = (xmax-x0)/dx; 
    // cout << nright << endl;
    for (int i = 0; i < nright; i++) {
        double xnow = 0 + (i + 0.5) * dx;
        for (int j = 0; j < nright; j++) {
            double ynow = ymin + (j + 0.5) * dx;
            for (int k = 0; k < nright; k++) {
                double znow = zmin + (k + 0.5) * dx;
                pt.x1 = xnow+ran(); 
                pt.x2 = ynow+ran(); 
                pt.x3 = znow+ran(); 
                pt.den = rhoright; 
                pt.mas = m0; 
                pt.vel1 = 0; 
                pt.vel2 = 0; 
                pt.vel3 = 0; 
                pt.len = hfact*dx ;     // 1.2*dxright
                pt.ene = pright/((gamma-1.)*rhoright) ;
                pts.add_particle(pt); 
            }
        }
    }
    cout << "Particles initialized: " << pts.size() << endl;
    return; 
}