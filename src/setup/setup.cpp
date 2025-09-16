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

    double gamma = 1.4; 
    Particle pt; 
    RandomGenerator ranx_left(xmin, x0);
    RandomGenerator ranx_right(x0, xmax);
    RandomGenerator rany(ymin, ymax); 
    RandomGenerator ranz(zmin, zmax); 

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

    double dV = (x0-xmin)*(ymax-ymin)*(zmax-zmin)/nleft; 
    double dr = cbrt(dV);
    double m0 = rholeft*dV; 
    for (int i=0; i<nleft; i++){
        pt.x1 = ranx_left(); 
        pt.x2 = rany(); 
        pt.x3 = ranz(); 
        pt.den = rholeft; 
        pt.mas = m0; 
        pt.vel1 = 0; 
        pt.vel2 = 0; 
        pt.vel3 = 0; 
        pt.len = hfact*dr ;     // 1.2*dxleft
        pt.ene = pleft/((gamma-1.)*rholeft) ;
        pts.add_particle(pt); 
    }

    dV = m0/rhoright; 
    dr = cbrt(dV);
    int nright = (xmax-x0)*(ymax-ymin)*(zmax-zmin)/dV; 
    for (int i=0; i<nright; i++){
        pt.x1 = ranx_right(); 
        pt.x2 = rany(); 
        pt.x3 = ranz(); 
        pt.den = rhoright; 
        pt.mas = m0; 
        pt.vel1 = 0; 
        pt.vel2 = 0; 
        pt.vel3 = 0; 
        pt.len = hfact*dr ;     // 1.2*dxleft
        pt.ene = pright/((gamma-1.)*rhoright) ;
        pts.add_particle(pt);
    }
    cout << "Particles initialized: " << pts.size() << endl;
    return; 
}