#include "modulation.hpp"
#include "modulationBPSK.hpp"
#include "modulation16PSK.hpp"
#include "modulation16QAM.hpp"
#include "modulation32PSK.hpp"
#include "modulation32QAM.hpp"
#include "function.hpp"
#include <iostream>
#include <cmath>

//
// EDIT!!!!!!!!!
//
CSignal* CSignal::new_gen_noise(modulation type,double r,std::default_random_engine &generator){
    CSignal* Noise=new CSignal();
    double sigma = 0;
    //Es/N0
    /*if(type==bpsk){
        sigma = sqrt((  pow(10,(-(r/10)))  )/2);
    }else if(type==psk16){
        sigma = sqrt((  pow(10,(-(r/10)))  )/2);
    }else if(type==qam16){
        sigma = sqrt((  pow(10,(-(r/10)))  ) * 5);
    }else if(type==psk32){
        sigma = sqrt((  pow(10,(-(r/10)))  )/2);
    }else if(type==qam32){
        sigma = sqrt((  pow(10,(-(r/10)))  ) * 10);
    }
    */
    //Eb/N0
    /*if(type==bpsk){
        sigma = sqrt((  pow(10,(-(r/10)))  )/2);
    }else if(type==psk16){
        sigma = sqrt((  pow(10,(-(r/10)))  )/8);
    }else if(type==qam16){
        sigma = sqrt((  pow(10,(-(r/10)))  ) * 1.25);
    }else if(type==psk32){
        sigma = sqrt((  pow(10,(-(r/10)))  )/10);
    }else if(type==qam32){
        sigma = sqrt((  pow(10,(-(r/10)))  ) * 2);
    }*/
    sigma=sqrt(get_sigma2(type,r)/2);
    //std::cout<<sigma<<std::endl;
    std::normal_distribution<float> norm(0.0, sigma);
    for(int i=0;i<MASK;i++){
        if(type==bpsk){
            Noise->signal[i]=CModulationBPSK::new_gen_noise(norm,generator);
        }else if(type==psk16){
            Noise->signal[i]=CModulation16PSK::new_gen_noise(norm,generator);
        }else if(type==qam16){
            Noise->signal[i]=CModulation16QAM::new_gen_noise(norm,generator);
        }else if(type==psk32){
            Noise->signal[i]=CModulation32PSK::new_gen_noise(norm,generator);
        }else if(type==qam32){
            Noise->signal[i]=CModulation32QAM::new_gen_noise(norm,generator);
        }
    }
    return Noise;
}
CModulation* new_gen_modulation(modulation type,Code value){
    if(type==bpsk){
        return new CModulationBPSK(value);
    }else if(type==psk16){
        return new CModulation16PSK(value);
    }else if(type==qam16){
        return new CModulation16QAM(value);
    }else if(type==psk32){
        return new CModulation32PSK(value);
    }else if(type==qam32){
        return new CModulation32QAM(value);
    }
    return NULL;
}

void CSignal::computing_probability(double noise){
    for(int i=0;i<MASK;i++){
        signal[i]->computing_probability(noise);
    }
}
CSignal::CSignal(CPolynomial R,modulation type){  
    for(int i=0;i<MASK;i++){
        signal[i] =new_gen_modulation(type,R[i].getNum());
    }
}
double CSignal::euclidean_distance(CPolynomial &C){
    modulation type=modulation_type();
    double ed=0;
    CModulation* aa=new_gen_modulation(type,0);
    for(int i=0;i<MASK;i++){
        aa->set_code(C[i].getNum());
        if(type==bpsk){
            ed+=signal[i]->distance(aa);//this is the BPSK distance_pow2
        }else{
            ed+=signal[i]->distance_pow2(aa);
        }
    }
    delete aa;
    return ed;
}
CSignal::CSignal(){
    for(int i=0;i<MASK;i++)
        signal[i]=NULL;
}
CSignal::~CSignal(){
    for(int i=0;i<MASK;i++)
        if(signal[i]) delete signal[i];
}
CPolynomial CSignal::get_data(){
    CPolynomial R;
    for(int i=0;i<MASK;i++){
        R[i]=signal[i]->get_code();
    }
    return R;
}
CModulation* CSignal::operator[](unsigned int index){
    return signal[index];
}
CSignal* CSignal::new_add(CSignal* b){
    CSignal* message=new CSignal();
    for(int i=0;i<MASK;i++){
        //std::cout<<signal[i]<<std::endl;
        message->signal[i] = signal[i]->new_add((*b)[i]);
    }
    return message;
}
double CSignal::euclidean_distance(CSignal* b){
    double ed=0;
    for(int i=0;i<MASK;i++){
        ed+=signal[i]->distance_pow2(b->signal[i]);
    }
    return ed;
}

modulation CSignal::modulation_type(){
    if(signal[0]) return signal[0]->modulation_type();
    return nulltype;
}
CPolynomial CSignal::get_random_data(std::default_random_engine &generator){
    CPolynomial R;
    for(int i=0;i<MASK;i++){
        R[i]=signal[i]->random_code(generator);
    }
    return R;
}
//virtual static
//CModulation* CModulation::gen_noise(double r,std::default_random_engine &generator){return NULL;}

CModulation::CModulation(){
    real=0.0;inag=0.0;
    total_probability=0;
    //code1=0;code2=0;
    for(int i=0;i<=MASK;i++) code[i]=0;
}

CModulation::CModulation(double r,double i){
    real=r;inag=i;
    total_probability=0;
    //code1=0;code2=0;
    for(int i=0;i<=MASK;i++) code[i]=0;
}
modulation CModulation::modulation_type(){return nulltype;}//virtual
void CModulation::set_code(Code in){std::cerr<<"[CModulation:set_code]the class is not valid!"<<std::endl;exit(1);}//virtual
double CModulation::get_N0(){
    std::cerr<<"[CModulation:get_N0]the class is not valid!"<<std::endl;exit(1);
    return 0;
}
Code CModulation::get_code(){
    if(total_probability<=0){std::cerr<<"[get_code()]you don't to computing the probability of the symbol!!!"<<std::endl;exit(1);}
    //return code1;
    return code[0];
}//virtual
Code CModulation::get_code2(){
    if(total_probability<=0){std::cerr<<"[get_code2()]you don't to computing the probability of the symbol!!!"<<std::endl;exit(1);}
    //return code2;
    return code[1];
}//virtual

double CModulation::reliability(){
    modulation type= this->modulation_type();
    CModulation* a=new_gen_modulation(type,this->get_code());
    CModulation* b=new_gen_modulation(type,this->get_code2());
    double r=this->distance_pow2(b) - this->distance_pow2(a);
    delete a;delete b;
    if(r<0){ std::cout<<r<<std::endl;exit(87);}
    //std::cout<<"?";
    return r;
}//virtual
CModulation* CModulation::new_add(CModulation *b){
    CModulation* ans=new_gen_modulation(this->modulation_type(),0);
    ans->real=this->real+b->real;
    ans->inag=this->inag+b->inag;
    return ans;
}//virtual
void CModulation::computing_probability(double noise){// "this" 適用多型
    CModulation *test=new_gen_modulation(this->modulation_type(),0);
    total_probability=0;
    //double max1=-1;double max2=-1;
    double sigama =  get_sigma2(this->modulation_type(),noise);
    for(int i=0;i<=MASK;i++){
        code[i]=i;
        test->set_code(i);
        //std::cout<<test->modulation_type()<<std::endl;
        double d=this->distance_pow2(test);
        //probability[i]=exp(-d/(2));//---------------------------
        probability[i]=exp(-d/(2*sigama));
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
    //std::cout<<"total:"<<total_probability<<std::endl;
    //for(int i=0;i<=MASK;i++){
    //   std::cout<<code[i]<<" "<<probability[code[i]]<<std::endl;
    //}
    //std::cout<<std::endl;
    delete test;
}
Code CModulation::random_code(std::default_random_engine &generator){
    if(total_probability<=0){
        std::cerr<<"[random_code()]you don't to computing the probability of the symbol!!!"<<std::endl;
        exit(1);
    }
    std::uniform_real_distribution<double> unif(0.0,total_probability);
    double rand=unif(generator);
    for(int i=0;i<=MASK;i++){
        if(rand<probability[i]) return i;
        rand-=probability[i];
    }
    if(rand>0) exit(78);
    return MASK;
}
double CModulation::distance(CModulation *b){
    double dx = real-b->real;
    double dy = inag-b->inag;
    return sqrt(dx*dx+dy*dy);
}

double CModulation::distance_pow2(CModulation *b){
    double dx = real-b->real;
    double dy = inag-b->inag;
    return dx*dx+dy*dy;
}

void CModulation::show(){
    std::cout<<"\treal:"<<real<<std::endl;
    std::cout<<"\tinag:"<<inag<<std::endl;
}
void CSignal::showSignal(){
    for(int i=0;i<MASK;i++){
        std::cout<<i<<std::endl;
        signal[i]->show();
    }
}