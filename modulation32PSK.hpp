#ifndef MODULATION32PSK_HPP
#define MODULATION32PSK_HPP
#include "modulation.hpp"

class CModulation32PSK:public CModulation{
public:
    //const Code garry[16]={0x0,0x8,0x9,0xb,0xa,0xe,0xf,0xd,0xc,0x4,0x5,0x7,0x6,0x2,0x3,0x1};
    const Code bin[32]={0,1,3,2,7,6,4,5,15,14,12,13,8,9,11,10,31,30,28,29,24,25,27,26,16,17,19,18,23,22,20,21};
    static CModulation* new_gen_noise(std::normal_distribution<float> &norm,std::default_random_engine &generator);
    CModulation32PSK(Code);
    CModulation32PSK();
    virtual void set_code(Code in);
    virtual double get_N0();
    virtual modulation modulation_type();
};
#endif