#include "modulation16PSK.hpp"
#include <iostream>
#include <cmath>
CModulation* CModulation16PSK::new_gen_noise(std::normal_distribution<float> &norm,std::default_random_engine &generator){
    CModulation16PSK* symbol=new CModulation16PSK();
    symbol->real=norm(generator);
    symbol->inag=norm(generator);
    return symbol;
}

CModulation16PSK::CModulation16PSK(){}

CModulation16PSK::CModulation16PSK(Code in){
    set_code(in);
}

void CModulation16PSK::set_code(Code in){
    total_probability=0;
    double sita=PI*bin[in]/8;
    inag=1.0*sin(sita);
    real=1.0*cos(sita);
}

double CModulation16PSK::get_N0(){return 1;}

modulation CModulation16PSK::modulation_type(){return psk16;}