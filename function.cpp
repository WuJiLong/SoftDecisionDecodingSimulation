#include "function.hpp"





int get_weight(Code c){
    int cnt=0;
    while(c){
        c=c&(c-1);
        cnt++;
    }
    return cnt;
}
unsigned int computing_BER(CPolynomial &C,CPolynomial &DC){
    unsigned int count=0;
    for(int i=0;i<MASK;i++){
        count+=get_weight((C[i]+DC[i]).getNum());
    }
    return count;
}
unsigned int computing_SER(CPolynomial &C,CPolynomial &DC){
    unsigned int count=0;
    for(int i=0;i<MASK;i++){
        if(C[i]!=DC[i]) count+=1;
    }
    return count;
}
unsigned int computing_CER(CPolynomial &C,CPolynomial &DC){
    if(C!=DC) return 1;
    return 0;
}

double get_sigma2(modulation type,double r){
    //Eb/N0
    if(type==bpsk){
        return pow(10,(-(r/10)));
    }else if(type==psk16){
        return (pow(10,(-(r/10))))/4;
    }else if(type==qam16){
        return (pow(10,(-(r/10)))) * 2.5;
    }else if(type==psk32){
        return (pow(10,(-(r/10))))/5;
    }else if(type==qam32){
        return (pow(10,(-(r/10)))) * 4;
    }
    return 0;
}