#include "random.h"

using namespace std; 

RandomGenerator::RandomGenerator(double low, double high):
    gen(random_device{}()), dis(low, high) {}

double RandomGenerator::operator()() {
    return dis(gen);
}