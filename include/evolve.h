#ifndef EVOLVE_H
#define EVOLVE_H
#include "particles.h"
#include "kdtree.h"

void Derivs(Particles& pts, Tree& tree); 

void Integral(Particles& pts, Tree& tree, double dt); 

#endif