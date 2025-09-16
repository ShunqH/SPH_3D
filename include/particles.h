#ifndef PARTICLES_H
#define PARTICLES_H
#include <iostream>  
#include <vector>   
#include <cmath>
#include <memory> 
#include "eos.h" 
#include "boundary.h"

// a single particle; 
class Particle{
public:
    double x1, x2, x3; 
    double vel1, vel2, vel3;
    double acc1, acc2, acc3, accu;
    double den; 
    double mas; 
    double len; 
    double ene; 
    int status; 

    Particle(); 
};

// particles' quantity arrays/structrue
class Particles{
public: 
    int NMAX;                                   // maximum length of the data array
    std::vector<double> x1, x2, x3;             // positions x, y, z
    std::vector<double> vel1, vel2, vel3;       // velocity vx, vy, xz 
    std::vector<double> acc1, acc2, acc3, accu; // accleration for x, y, z, and energy (viscoty heating) 
    std::vector<double> den;                    // density 
    std::vector<double> mas;                    // mass of the particles
    std::vector<double> len;                    // smooth length (size) of the particle
    std::vector<double> ene;                    // energy 
    std::vector<int> status;                    // status: 0 alive, 1 dead
    int endid;                                  // the number of activated particles (dead included)
    int ghostid; 
    std::shared_ptr<EoS> eos = nullptr;
    std::shared_ptr<Boundary> BDX1L = nullptr;
    std::shared_ptr<Boundary> BDX1R = nullptr;
    std::shared_ptr<Boundary> BDX2L = nullptr;
    std::shared_ptr<Boundary> BDX2R = nullptr;
    std::shared_ptr<Boundary> BDX3L = nullptr;
    std::shared_ptr<Boundary> BDX3R = nullptr;

    Particles(int nmax); 

    void add_particle(const Particle& pt);

    Particle extract(int index);

    int size(); 

    void check_valid(int i) const;

    double dist(int ia, int ib); 
    std::vector<double> vecr(int ia, int ib); 
    std::vector<double> vecv(int ia, int ib); 

    double pres(int i) const;                 
    double cs(int i) const;

    void clean_ghost(); 
    void set_boundary(); 
    

private: 
    
}; 


#endif