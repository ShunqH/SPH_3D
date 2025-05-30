#ifndef UTILS_H
#define UTILS_H
#include <vector>
#include "particles.h"

const double pi = 3.141592653589793238462643383; 

void WriteParticles(Particles& pts, int step); 

double dot(const std::vector<double>& va, const std::vector<double>& vb); 

double power(double base, int order); 

#endif 