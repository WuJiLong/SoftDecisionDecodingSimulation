#include "modulation32QAM.hpp"
#include <iostream>
#include <cmath>
CModulation* CModulation32QAM::new_gen_noise(std::normal_distribution<float> &norm,std::default_random_engine &generator){
    CModulation32QAM* symbol=new CModulation32QAM();
    symbol->real=norm(generator);
    symbol->inag=norm(generator);
    return symbol;
}

CModulation32QAM::CModulation32QAM(){}

CModulation32QAM::CModulation32QAM(Code in){
    set_code(in);
}

void CModulation32QAM::set_code(Code in){
    total_probability=0;
    Code MSB=in>>2;
    Code LSB=in&0x3;
    if(MSB==0b000){
        if((LSB&0x2)==0) inag=5.0;
        else inag=-5.0;
        if((LSB&0x1)==0) real=-3.0;
        else real=-1.0;
    }else if(MSB==0b100){
        if((LSB&0x2)==0) inag=5.0;
        else inag=-5.0;
        if((LSB&0x1)==0) real=3.0;
        else real=1.0;
    }else{
        if(LSB==0) inag=3.0;
        else if(LSB==1) inag=1.0;
        else if(LSB==3) inag=-1.0;
        else inag=-3.0;
        if(MSB==1) real=-5.0;
        else if(MSB==3) real=-3.0;
        else if(MSB==2) real=-1.0;
        else if(MSB==6) real=1.0;
        else if(MSB==7) real=3.0;
        else real=5.0;//MSB==5
    }
    /*Code msb=(in>>2)&3;
    Code lsb=in&3;
    inag = -3.0 +g_v[lsb]*2;
    real = -3.0 +g_v[msb]*2;*/
}

double CModulation32QAM::get_N0(){
    return 20;//sqrt(20);
}

modulation CModulation32QAM::modulation_type(){return qam32;}