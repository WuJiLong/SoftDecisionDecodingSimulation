#ifndef POLYNOMIAL_HPP
#define POLYNOMIAL_HPP
#include <iostream>
#include "define.hpp"
#include "Galois.hpp"
class CPolynomial{
    friend std::ostream &operator<<(std::ostream&,CPolynomial);
    public:
        CPolynomial();
        CPolynomial(Code*,int);
        //CPolynomial(Code[],int);
        void initial(Code*,int);
        int get_max();
        int get_max(int);
        CPolynomial operator+(CPolynomial);
        CPolynomial operator*(CGalois);
        CPolynomial operator/(CGalois);
        CPolynomial operator*(CPolynomial);
        CPolynomial operator/(CPolynomial);
        CPolynomial operator%(CPolynomial);
        CPolynomial operator<<(int);
        CPolynomial operator>>(int);
        CGalois &operator[](Code);
        CPolynomial differential();
        bool operator==(CPolynomial);
        bool operator!=(CPolynomial);
        CGalois inputx(CGalois);
        unsigned int distance(CPolynomial);
        bool empty();
        unsigned int weight();
    private:
        CGalois Coefficient[PolynomialSize];
};
std::ostream &operator<<(std::ostream&,CPolynomial);
#endif
