#ifndef ALGORITHM_HPP
#define ALGORITHM_HPP
#include "modulation.hpp"

CPolynomial algorithm_EC2(CPolynomial &R,CSignal &y,unsigned long long &decode_count,bool *pass=NULL);
CPolynomial algorithm_CHU(CPolynomial &R,CSignal &y,unsigned long long &decode_count,bool *pass=NULL);
CPolynomial algorithm_CHU_and_EC2(CPolynomial &R,CSignal &y,unsigned long long  &decode_count,bool *pass=NULL);
CPolynomial algorithm_CHU2(CPolynomial &R,CSignal &y,unsigned long long &decode_count,bool *pass=NULL);
//CPolynomial algorithm_CHU22(CPolynomial &R,CSignal &y,unsigned long long &decode_count,bool *pass=NULL);
CPolynomial algorithm_CHU2_and_EC2(CPolynomial &R,CSignal &y,unsigned long long &decode_count,bool *pass=NULL);


CPolynomial algorithm_CHU_SORT(CPolynomial &R,CSignal &y,unsigned long long &decode_count,bool *pass=NULL);
CPolynomial algorithm_CHU_and_EC2_SORT(CPolynomial &R,CSignal &y,unsigned long long  &decode_count,bool *pass=NULL);
CPolynomial algorithm_CHU2_SORT(CPolynomial &R,CSignal &y,unsigned long long &decode_count,bool *pass=NULL);
CPolynomial algorithm_CHU2_and_EC2_SORT(CPolynomial &R,CSignal &y,unsigned long long &decode_count,bool *pass=NULL);


CPolynomial algorithm_SCA_EC2(Code pattern,CPolynomial &R,CSignal &y,unsigned long long &decode_count,bool *pass=NULL);
CPolynomial algorithm_SCA_CHU(Code pattern,CPolynomial &R,CSignal &y,unsigned long long &decode_count,bool *pass=NULL);
CPolynomial algorithm_SCA_CHU_and_EC2(Code pattern,CPolynomial &R,CSignal &y,unsigned long long  &decode_count,bool *pass=NULL);
CPolynomial algorithm_SCA_CHU2(Code pattern,CPolynomial &R,CSignal &y,unsigned long long &decode_count,bool *pass=NULL);
CPolynomial algorithm_SCA_CHU2_and_EC2(Code pattern,CPolynomial &R,CSignal &y,unsigned long long &decode_count,bool *pass=NULL);


CPolynomial algorithm_SSCA(Code pattern,CPolynomial &R,CSignal &y,unsigned long long &decode_count,bool *pass=NULL);

CPolynomial algorithm_SSCA_EC2(Code pattern,CPolynomial &R,CSignal &y,unsigned long long &decode_count,bool *pass=NULL);
CPolynomial algorithm_SSCA_CHU(Code pattern,CPolynomial &R,CSignal &y,unsigned long long &decode_count,bool *pass=NULL);
CPolynomial algorithm_SSCA_CHU_and_EC2(Code pattern,CPolynomial &R,CSignal &y,unsigned long long  &decode_count,bool *pass=NULL);
CPolynomial algorithm_SSCA_CHU2(Code pattern,CPolynomial &R,CSignal &y,unsigned long long &decode_count,bool *pass=NULL);
CPolynomial algorithm_SSCA_CHU2_and_EC2(Code pattern,CPolynomial &R,CSignal &y,unsigned long long &decode_count,bool *pass=NULL);



CPolynomial algorithm_CHU3(CPolynomial &R,CSignal &y,unsigned long long &decode_count,bool *pass=NULL);
CPolynomial algorithm_CHU3_and_EC2(CPolynomial &R,CSignal &y,unsigned long long &decode_count,bool *pass=NULL);
CPolynomial algorithm_CHU3_SORT(CPolynomial &R,CSignal &y,unsigned long long &decode_count,bool *pass=NULL);
CPolynomial algorithm_CHU3_and_EC2_SORT(CPolynomial &R,CSignal &y,unsigned long long &decode_count,bool *pass=NULL);
CPolynomial algorithm_SCA_CHU3(Code pattern,CPolynomial &R,CSignal &y,unsigned long long &decode_count,bool *pass=NULL);
CPolynomial algorithm_SCA_CHU3_and_EC2(Code pattern,CPolynomial &R,CSignal &y,unsigned long long &decode_count,bool *pass=NULL);
CPolynomial algorithm_SSCA_CHU3(Code pattern,CPolynomial &R,CSignal &y,unsigned long long &decode_count,bool *pass=NULL);
CPolynomial algorithm_SSCA_CHU3_and_EC2(Code pattern,CPolynomial &R,CSignal &y,unsigned long long &decode_count,bool *pass=NULL);

CPolynomial algorithm_ENCODE_CHU3(Code pattern,CPolynomial &R,CSignal &y,unsigned long long &decode_count,bool *pass=NULL);
#endif