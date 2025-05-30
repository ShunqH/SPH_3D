#ifndef RANDOM_H
#define RANDOM_H
#include <ctime>
#include <random>


class RandomGenerator {
    std::mt19937 gen;
    std::uniform_real_distribution<> dis;

public:
    RandomGenerator(double low, double high);

    double operator()(); 
};

#endif

