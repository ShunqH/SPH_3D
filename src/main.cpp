#include <iostream>  
#include <vector> 
#include "particles.h"
#include "setup.h"
#include "config.h"
#include "utils.h"
#include "kdtree.h"
#include "density.h"
#include "force.h"
#include "evolve.h"

using namespace std;

int main(int argc, char* argv[]){
    double dt, t, tmax, dtoutput; 
    int step = 0; 
    // int nghost; 

    if (argc < 3 || std::string(argv[1]) != "-i") {
        cerr << "Usage: " << argv[0] << " -i input.in" << endl;
        return 1;
    }
    Config::getInstance().loadFromFile(argv[2]);

    const int NMAX = Config::getInstance().get("NMAX"); 
    Particles pts(NMAX); 

    Setup(pts); 
    Tree tree(pts); 
    
    dt = Config::getInstance().get("dt"); 
    tmax = Config::getInstance().get("tmax"); 
    dtoutput = Config::getInstance().get("dtoutput"); 

    t = 0; 
    int interval = dtoutput/dt; 
    WriteParticles(pts, step); 
    cout<<"Start running"<<endl; 
    while (t<tmax){
        Integral(pts, tree, dt); 
        t = t+dt; 
        step++; 
        if (step%5 == 0){
            Tree tree(pts); 
        }
        if (step%interval == 0){
            cout << "\rProgress: t = " << t << "       " << flush;
            WriteParticles(pts, step/interval); 
        }            
    }
    cout << endl; 
    cout<< "finish!" <<endl; 
    
    return 0; 
}

