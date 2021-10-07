#include "Galois.hpp"


CGalois::CGalois(Code init,Code gx,Code m){
    
    this->GFx=gx;
    this->mask=m;
    int x=0;
    int a=1;
    if(init>m){
        //std::cout<<"debug";
        while(init!=0){
            if(init%2==1)
                x=x^a;
            a=a<<1;
            if(a>mask)
                a=a^gx;
            init=init>>1;
        }
        //std::cout<<x<<std::endl;
    }else{
        x=init;
    }
    this->number=x;
}
Code CGalois::getNum(){
    return number;
}
/*
CGalois CGalois::operator+(Code b){
    CGalois R(this->number,this->GFx,this->mask);
    R.number^=b;
    R.number&=mask;
    return R;
}*/
CGalois CGalois::operator+(CGalois b){
    CGalois R(this->number,this->GFx,this->mask);
    R.number^=b.number;
    R.number&=mask;
    return R;
}
CGalois CGalois::operator*(CGalois b){
    return CGalois(mul(this->number,b.number),this->GFx,this->mask);
}
/*
CGalois CGalois::operator*(Code b){
    return CGalois(mul(this->number, b),this->GFx,this->mask);
}*/
CGalois CGalois::operator*=(CGalois b){
    this->number=mul(this->number, b.number);
    return *this;
}
CGalois CGalois::operator/=(CGalois b){
    this->number=mul(this->number, pow(b.number,-1));
    return *this;
}
/*
CGalois CGalois::operator*=(Code b){
    this->number=mul(this->number, b);
    return *this;
}*/
CGalois CGalois::operator/(CGalois  b){
    return CGalois(mul(this->number, pow(b.number,-1)),this->GFx,this->mask);
}
CGalois CGalois::operator^(int  b){
    return CGalois(pow(this->number, b),this->GFx,this->mask);
}
/*
CGalois CGalois::operator/(Code b){
    return CGalois(mul(this->number, pow(b,-1)),this->GFx,this->mask);
}*/
Code CGalois::mul(Code a,Code b){
    Code ans=0;
    if(a&1) ans^=b;
    for(int i=1;a!=0;i++){
        a>>=1;b<<=1;
        if(b&(~mask)) b=b^GFx;
        if(a&1) ans^=b;
    }
    return ans;
}
Code CGalois::pow(Code a,int i) {
    if(a==0) return 0;
    while(i<0)
        i=mask+i;
        //i=(1<<GFbit)-1+i;
    Code ans=1;
    while(i!=0){
        if(i&1) ans=mul(ans,a);
        a=mul(a,a);
        i>>=1;
    }
    return ans;
}
bool CGalois::operator!=(CGalois b){
    return number!=b.number;
}
bool CGalois::operator==(CGalois b){
    return number==b.number;
}
std::ostream &operator<<(std::ostream &out,CGalois b){
    out<<b.number;
    return out;
}