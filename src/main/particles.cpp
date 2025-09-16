#include <iostream>  
#include <vector>   
#include <cmath>
#include <memory> 
#include <stdexcept>
#include "eos.h"   
#include "particles.h"
#include "boundary.h"

using namespace std; 

Particle::Particle(): 
    x1(0), x2(0), x3(0),
    vel1(0), vel2(0), vel3(0),
    acc1(0), acc2(0), acc3(0), accu(0),
    den(0), mas(0), len(0), ene(0), status(1)
{}

// allocate arrarys by nmax input
Particles::Particles(int nmax): 
NMAX(nmax), endid(0), ghostid(0), 
eos(nullptr), 
BDX1L(std::make_unique<BoundaryX1L_Free>()),
BDX1R(std::make_unique<BoundaryX1R_Free>()),
BDX2L(std::make_unique<BoundaryX2L_Free>()),
BDX2R(std::make_unique<BoundaryX2R_Free>()),
BDX3L(std::make_unique<BoundaryX3L_Free>()),
BDX3R(std::make_unique<BoundaryX3R_Free>())
{
    x1.resize(nmax); 
    x2.resize(nmax); 
    x3.resize(nmax); 
    vel1.resize(nmax); 
    vel2.resize(nmax); 
    vel3.resize(nmax); 
    acc1.resize(nmax); 
    acc2.resize(nmax); 
    acc3.resize(nmax); 
    accu.resize(nmax); 
    den.resize(nmax); 
    mas.resize(nmax);
    len.resize(nmax); 
    ene.resize(nmax); 
    status.resize(nmax); 
}

// Particles method: add particle, input: class Particle
// note that you have to do boundary again if you add an activated particle
void Particles::add_particle(const Particle& pt){
    int id = 0; 
    if (pt.status == 0 || pt.status == 1 || pt.status == 2){
        id = endid; 
        if (endid >= NMAX) {
            cerr << "Error: Too many particles! activated particles = "
                 << endid 
                 << ", maximum =" 
                 << NMAX << endl;
            return;
        }
        endid ++; 
        ghostid = endid; 
    }else if (pt.status == -1){
        id = ghostid; 
        if (ghostid >= NMAX) {
            cerr << "Error: Too many particles! activated particles = "
                 << endid 
                 << ", ghost particles = "
                 << ghostid - endid
                 << ", maximum =" 
                 << NMAX << endl;
            return;
        }
        ghostid ++; 
    }else{
        cerr << "unknown particle types, -1 for ghost, 1 for gas, "
             << endl;
        return;
    }
    x1[id] = pt.x1; 
    x2[id] = pt.x2; 
    x3[id] = pt.x3; 
    vel1[id] = pt.vel1;
    vel2[id] = pt.vel2; 
    vel3[id] = pt.vel3; 
    acc1[id] = pt.acc1;
    acc2[id] = pt.acc2; 
    acc3[id] = pt.acc3; 
    accu[id] = pt.accu; 
    den[id] = pt.den; 
    mas[id] = pt.mas;
    len[id] = pt.len; 
    ene[id] = pt.ene; 
    status[id] = pt.status; 
    
    return; 
}

// Particles method: extract a single Particle form Particles 
Particle Particles::extract(int index){
    // error if index out of number of particles (endid)
    if (index < 0 || index >= ghostid) {
        cerr << "Error: extract index out of range (index = " 
                  << index << ", ghostid = " << ghostid << ").\n";
        return Particle();  
    }
    Particle pt; 
    // extract particle at index
    pt.x1 = x1[index]; 
    pt.x2 = x2[index]; 
    pt.x3 = x3[index]; 
    pt.vel1 = vel1[index]; 
    pt.vel2 = vel2[index]; 
    pt.vel3 = vel3[index]; 
    pt.acc1 = acc1[index]; 
    pt.acc2 = acc2[index]; 
    pt.acc3 = acc3[index]; 
    pt.accu = accu[index]; 
    pt.den = den[index]; 
    pt.mas = mas[index]; 
    pt.len = len[index]; 
    pt.ene = ene[index]; 
    pt.status = status[index]; 
    return pt; 
}

int Particles::size(){
    return endid; 
}

double Particles::dist(int ia, int ib){
    return sqrt((x1[ia]-x1[ib])*(x1[ia]-x1[ib]) 
               +(x2[ia]-x2[ib])*(x2[ia]-x2[ib]) 
               +(x3[ia]-x3[ib])*(x3[ia]-x3[ib])); 
}

vector<double> Particles::vecr(int ia, int ib){
    vector<double> vec(3,0); 
    if (ia ==ib){
        return vec; 
    }
    vec[0] = x1[ia]-x1[ib]; 
    vec[1] = x2[ia]-x2[ib]; 
    vec[2] = x3[ia]-x3[ib]; 
    return vec; 
}

vector<double> Particles::vecv(int ia, int ib){
    vector<double> vec(3,0); 
    if (ia ==ib){
        return vec; 
    }
    vec[0] = vel1[ia]-vel1[ib]; 
    vec[1] = vel2[ia]-vel2[ib]; 
    vec[2] = vel3[ia]-vel3[ib]; 
    return vec; 
}

void Particles::check_valid(int i) const {
    if (i < 0 || i >= ghostid) 
      throw std::out_of_range("Particle index out of range");
  }

double Particles::pres(int i) const{
    check_valid(i); 
    if (!eos) {
        throw runtime_error("eos was not initialized");
    }
    return eos->pressure(den[i], ene[i]); 
}

double Particles::cs(int i) const{
    check_valid(i); 
    if (!eos) {
        throw runtime_error("eos was not initialized");
    }
    return eos->cs(den[i], ene[i]); 
}

void Particles::clean_ghost(){
    ghostid = endid; 
    return; 
}

void Particles::set_boundary(){
    clean_ghost(); 
    BDX1L->SetBoundary(*this); 
    BDX1R->SetBoundary(*this); 
    BDX2L->SetBoundary(*this); 
    BDX2R->SetBoundary(*this); 
    BDX3L->SetBoundary(*this); 
    BDX3R->SetBoundary(*this); 
    return; 
}
