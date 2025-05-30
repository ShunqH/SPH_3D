#ifndef DENSITY_H
#define DENSITY_H
#include <vector>
#include "density.h" 
#include "utils.h"
#include "particles.h"
#include "kdtree.h"

double w_m4(double q, double h); 

double wpq_m4(double dist, double h); 

double w_m6(double q, double h); 

void GetDensity(Particles& pts, Tree& tree); 

void Smoothing(Particles& pts); 

#endif 