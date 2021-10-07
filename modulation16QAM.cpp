#include "modulation16QAM.hpp"
#include <iostream>
#include <cmath>
CModulation* CModulation16QAM::new_gen_noise(std::normal_distribution<float> &norm,std::default_random_engine &generator){
    CModulation16QAM* symbol=new CModulation16QAM();
    symbol->real=norm(generator);
    symbol->inag=norm(generator);
    return symbol;
}

CModulation16QAM::CModulation16QAM(){}

CModulation16QAM::CModulation16QAM(Code in){
    set_code(in);
}

void CModulation16QAM::set_code(Code in){
    total_probability=0;
    Code msb=(in>>2)&3;
    Code lsb=in&3;
    inag = -3.0 +g_v[lsb]*2;
    real = -3.0 +g_v[msb]*2;
}

double CModulation16QAM::get_N0(){
    return 10;//sqrt(10);
}

modulation CModulation16QAM::modulation_type(){return qam16;}