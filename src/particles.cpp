#include <iostream>  
#include <vector>   
#include <cmath>
#include <memory> 
#include <stdexcept>
#include "eos.h"   
#include "particles.h"

using namespace std; 

Particle::Particle()
    : x1(0), x2(0), x3(0),
      vel1(0), vel2(0), vel3(0),
      acc1(0), acc2(0), acc3(0), accu(0),
      den(0), mas(0), len(0), ene(0)
{}

// allocate arrarys by nmax input
Particles::Particles(int nmax): NMAX(nmax), endid(0), eos(nullptr){
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
void Particles::add_particle(Particle pt){
    x1[endid] = pt.x1; 
    x2[endid] = pt.x2; 
    x3[endid] = pt.x3; 
    vel1[endid] = pt.vel1;
    vel2[endid] = pt.vel2; 
    vel3[endid] = pt.vel3; 
    acc1[endid] = pt.acc1;
    acc2[endid] = pt.acc2; 
    acc3[endid] = pt.acc3; 
    accu[endid] = pt.accu; 
    den[endid] = pt.den; 
    mas[endid] = pt.mas;
    len[endid] = pt.len; 
    ene[endid] = pt.ene; 
    status[endid] = 0; 
    endid ++; 
    if (endid >= NMAX) {
        cerr << "Error: Too many particles! Max = " << NMAX << endl;
        return;
    }
    return; 
}

// Particles method: extract a single Particle form Particles 
Particle Particles::extract(int index){
    // error if index out of number of particles (endid)
    if (index < 0 || index >= endid) {
        cerr << "Error: extract index out of range (index = " 
                  << index << ", endid = " << endid << ").\n";
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
    vec[0] = x1[ib]-x1[ia]; 
    vec[1] = x2[ib]-x2[ia]; 
    vec[2] = x3[ib]-x3[ia]; 
    return vec; 
}

vector<double> Particles::vecv(int ia, int ib){
    vector<double> vec(3,0); 
    if (ia ==ib){
        return vec; 
    }
    vec[0] = vel1[ib]-vel1[ia]; 
    vec[1] = vel2[ib]-vel2[ia]; 
    vec[2] = vel3[ib]-vel3[ia]; 
    return vec; 
}

void Particles::check_valid(int i) const {
    if (i < 0 || i >= endid) 
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

