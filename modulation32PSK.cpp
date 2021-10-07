#include "modulation32PSK.hpp"
#include <iostream>
#include <cmath>
CModulation* CModulation32PSK::new_gen_noise(std::normal_distribution<float> &norm,std::default_random_engine &generator){
    CModulation32PSK* symbol=new CModulation32PSK();
    symbol->real=norm(generator);
    symbol->inag=norm(generator);
    return symbol;
}

CModulation32PSK::CModulation32PSK(){}

CModulation32PSK::CModulation32PSK(Code in){
    set_code(in);
}

void CModulation32PSK::set_code(Code in){
    total_probability=0;
    double sita=PI*bin[in]/16;
    inag=1.0*sin(sita);
    real=1.0*cos(sita);
}

double CModulation32PSK::get_N0(){return 1;}

modulation CModulation32PSK::modulation_type(){return psk32;}