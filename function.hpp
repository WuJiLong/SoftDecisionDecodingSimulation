#ifndef FUNCTION_HPP
#define FUNCTION_HPP
#include "define.hpp"
#include "Galois.hpp"
#include "Polynomial.hpp"
#include "modulation.hpp"
int get_weight(Code c);
unsigned int computing_BER(CPolynomial &C,CPolynomial &DC);
unsigned int computing_SER(CPolynomial &C,CPolynomial &DC);
unsigned int computing_CER(CPolynomial &C,CPolynomial &DC);
double get_sigma2(modulation type,double r);

#endif