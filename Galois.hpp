#ifndef GALOIS_HPP
#define GALOIS_HPP
#include <iostream>
#include "define.hpp"
typedef unsigned long long Code;
class CPolynomial;
class CGalois{
    friend std::ostream &operator<<(std::ostream&,CGalois);
    friend CPolynomial;
    public:
    CGalois(Code init=0,Code gx=GFX,Code m=MASK);
    CGalois operator+(CGalois);
    //CGalois operator+(Code);
    CGalois operator*(CGalois);
    //CGalois operator*(Code);
    CGalois operator*=(CGalois);
    CGalois operator/=(CGalois);
    //CGalois operator*=(Code);
    CGalois operator/(CGalois);
    CGalois operator^(int);
    //CGalois operator/(Code);
    bool operator!=(CGalois);
    bool operator==(CGalois);
    Code getNum();
    private:
    Code mul(Code,Code);
    Code pow(Code,int);
    Code number;
    Code GFx;
    Code mask;
};
std::ostream &operator<<(std::ostream&,CGalois);
#endif