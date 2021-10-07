#ifndef GLOBAL_HPP
#define GLOBAL_HPP
#include "Polynomial.hpp"
#include <random>

class CGLOBAL{
    public:
    static unsigned int seed;
    static std::default_random_engine RAND_generator;
    static CPolynomial GXP;
    //static std::default_random_engine RAND_generatorA;
    //static std::default_random_engine RAND_generatorB;
    //static std::default_random_engine RAND_generatorC;
    //static std::default_random_engine RAND_generatorD;
    //static std::default_random_engine RAND_generatorE;
};
#endif