#ifndef MODULATION16QAM_HPP
#define MODULATION16QAM_HPP
#include "modulation.hpp"

class CModulation16QAM:public CModulation{
public:
    const Code g_v[4]={0,1,3,2};
    static CModulation* new_gen_noise(std::normal_distribution<float> &norm,std::default_random_engine &generator);
    CModulation16QAM(Code);
    CModulation16QAM();
    virtual void set_code(Code in);
    virtual double get_N0();
    virtual modulation modulation_type();
};
#endif