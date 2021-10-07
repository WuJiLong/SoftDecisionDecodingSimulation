#ifndef MODULATION_HPP
#define MODULATION_HPP
#include "define.hpp"
#include "Galois.hpp"
#include "Polynomial.hpp"
#include <random>
//#include <random>
//#define PI 3.14159265358979323846

enum modulation{nulltype=0,bpsk=1,psk16=2,qam16=3,psk32=4,qam32=5};
class CModulation{
public:
    //virtual CModulation* gen_noise(double r,std::default_random_engine &generator);
    CModulation();
    CModulation(double,double);
    virtual void set_code(Code in);
    virtual Code get_code();
    virtual Code get_code2();
    virtual CModulation* new_add(CModulation *b);
    virtual modulation modulation_type();
    virtual double reliability();
    virtual double distance(CModulation *b);
    virtual double distance_pow2(CModulation *b);
    virtual void computing_probability(double noise);
    virtual double get_N0();
    Code random_code(std::default_random_engine &generator);
    virtual void show();
    //CModulation operator+(CModulation &b);
//protected:
    double real;
    double inag;
    double probability[MASK+1];
    double total_probability;
    //Code code1;
    //Code code2;
    Code code[MASK+1];
};
class CSignal{
public:
    static CSignal* new_gen_noise(modulation type,double r,std::default_random_engine &generator);
    CSignal();
    CSignal(CPolynomial,modulation type);
    ~CSignal();
    CPolynomial get_data();
    CPolynomial get_random_data(std::default_random_engine &generator);
    double euclidean_distance(CSignal*);
    double euclidean_distance(CPolynomial&);
    CSignal* new_add(CSignal*);
    //CSignal operator+(CSignal&);
    CModulation* operator[](unsigned int);
    modulation modulation_type();
    void computing_probability(double noise);
    void showSignal();
//private:
    CModulation* signal[MASK];
};
CModulation* new_gen_modulation(modulation,Code);


#endif
