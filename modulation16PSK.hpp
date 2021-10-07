#ifndef MODULATION16PSK_HPP
#define MODULATION16PSK_HPP
#include "modulation.hpp"

class CModulation16PSK:public CModulation{
public:
    //const Code garry[16]={0x0,0x8,0x9,0xb,0xa,0xe,0xf,0xd,0xc,0x4,0x5,0x7,0x6,0x2,0x3,0x1};
    const Code bin[16]={0,15,13,14,9,10,12,11,1,2,4,3,8,7,5,6};
    static CModulation* new_gen_noise(std::normal_distribution<float> &norm,std::default_random_engine &generator);
    CModulation16PSK(Code);
    CModulation16PSK();
    virtual void set_code(Code in);
    virtual double get_N0();
    virtual modulation modulation_type();
};
#endif