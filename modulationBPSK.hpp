#ifndef MODULATIONBPSK_HPP
#define MODULATIONBPSK_HPP
#include "modulation.hpp"

class CModulationBPSK:public CModulation{
public:
    static CModulation* new_gen_noise(std::normal_distribution<float> &norm,std::default_random_engine &generator);
    CModulationBPSK(Code);
    CModulationBPSK();
    virtual void set_code(Code in);
    //virtual Code get_code();
    //virtual Code get_code2();
    virtual CModulation* new_add(CModulation *b);
    virtual modulation modulation_type();
    virtual double reliability();
    virtual double get_N0();
    virtual double distance_pow2(CModulation *b);
    virtual double distance(CModulation *b);
    virtual void computing_probability(double noise);
    virtual void show();
    double value[_m_];
};
#endif