#include "modulationBPSK.hpp"
#include "function.hpp"
#include <iostream>
#include <cmath>
CModulation* CModulationBPSK::new_gen_noise(std::normal_distribution<float> &norm,std::default_random_engine &generator){
    CModulationBPSK* symbol=new CModulationBPSK();
    for(int i=0;i<_m_;i++)
        symbol->value[i]=norm(generator);
    return symbol;
}

CModulationBPSK::CModulationBPSK(){}

CModulationBPSK::CModulationBPSK(Code in){
    set_code(in);
}

void CModulationBPSK::set_code(Code in){
    total_probability=0;
    for(int i=0;i<_m_;i++){
        if(in&1){
            value[i]=-1.0;
        }else{
            value[i]=1.0;
        }
        in=in>>1;
    }
}
/*
Code CModulationBPSK::get_code(){
    Code ans=0;
    for(int i=_m_-1;i>=0;i--){
        ans=ans<<1;
        if(value[i]<0){
            ans=ans|0x01ull;
        }
    }
    return ans;
}

Code CModulationBPSK::get_code2(){//not use
    std::cerr<<"The type is BPSK! function\"get_code2()\" can't to use!!"<<std::endl;
    return 0;
}*/

CModulation* CModulationBPSK::new_add(CModulation *b){
    //std::cout<<"sadd"<<std::endl;
    CModulationBPSK* ans=new CModulationBPSK();
    CModulationBPSK* bbb=(CModulationBPSK*)b;
    for(int i=0;i<_m_;i++)
        ans->value[i]=this->value[i]+bbb->value[i];
    return ans;
}
double CModulationBPSK::get_N0(){return 1;}

double CModulationBPSK::reliability(){
    double r=ABS(value[0]);
    for(int i=1;i<_m_;i++){
        if(ABS(value[i])<r) r=ABS(value[i]);
    }
    return r;
}
double CModulationBPSK::distance_pow2(CModulation *b){// is not original distance;
    CModulationBPSK* bbb=(CModulationBPSK*)b;
    double sum=0;
    for(int i=0;i<_m_;i++){
        //if the signed is not equal,then add the bit's reliability;
        if(value[i]*bbb->value[i] < 0) sum+=ABS(value[i]);
    }
    return sum;
}
double CModulationBPSK::distance(CModulation *b){// this is distance_pow2;
    CModulationBPSK* bbb=(CModulationBPSK*)b;
    double sum=0;
    for(int i=0;i<_m_;i++){  
        sum+=pow(value[i]-bbb->value[i],2);
    }
    return sum;
}
void CModulationBPSK::computing_probability(double noise){
    total_probability=0;
    double bitstochastic[_m_];
    double sigama =  get_sigma2(this->modulation_type(),noise)/2;
    for(int j=0;j<_m_;j++){
        bitstochastic[j]=1.0/(1.0+exp(value[j]*2/sigama));
    }
    double max1=-1;double max2=-1;
    for(int i=0;i<=MASK;i++){
        code[i]=i;
        probability[i]=1.0;
        Code in=i;
        for(int j=0;j<_m_;j++){
            if(in&1){
                probability[i]*=bitstochastic[j];
            }else{
                probability[i]*=(1-bitstochastic[j]);
            }
            in=in>>1;
        }
        total_probability+=probability[i];
        /*if(probability[i]>max1){
            max2=max1;
            max1=probability[i];
            code2=code1;
            code1=i;
        }else if(probability[i]>max2){
            max2=probability[i];
            code2=i;
        }*/
    }
    for(int i=0;i<=MASK;i++){
        for(int j=i+1;j<=MASK;j++){
            if(probability[code[i]]<probability[code[j]]){
                Code tmp=code[i];
                code[i]=code[j];
                code[j]=tmp;
            }
        }
    }

}
modulation CModulationBPSK::modulation_type(){return bpsk;}
void CModulationBPSK::show(){
    for(int i=0;i<_m_;i++){
        std::cout<<"\t"<<i<<":"<<value[i]<<std::endl;
    }
}