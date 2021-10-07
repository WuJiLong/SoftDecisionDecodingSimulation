#ifndef DECODING_HPP
#define DECODING_HPP
#include "define.hpp"
#include "Galois.hpp"
#include "Polynomial.hpp"

CPolynomial generateGX(int t,int b=1,CGalois a=2);
CPolynomial getMessage(int n);
CPolynomial RSencode(CPolynomial M,CPolynomial &GX,int t=T);
CPolynomial getNose(int n);
CPolynomial Euclidean_Algorithm(CPolynomial &R,bool *pass=NULL);
CPolynomial Berlekamp_Massey(CPolynomial &R,bool* pass=NULL);

#endif