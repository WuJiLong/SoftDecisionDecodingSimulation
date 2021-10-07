#include "Polynomial.hpp"

CPolynomial::CPolynomial(){
    for(int i=0;i<PolynomialSize;i++)
        this->Coefficient[i].number=0;
}
CPolynomial::CPolynomial(Code *b,int size){
    initial(b,size);
}
/*CPolynomial::CPolynomial(Code b[],int size){
    initial(b,size);
}*/
void CPolynomial::initial(Code *b,int size){
    int i=0;
    for(;i<size;i++){
        Coefficient[i]=b[i];
    }
    for(;i<PolynomialSize;i++){
        Coefficient[i]=0;
    }
}
int CPolynomial::get_max(){//deg
    for(int i=PolynomialSize-1;i>=0;i--)
        if(Coefficient[i]!=0) return i;
    return -1;
}
int CPolynomial::get_max(int m){//deg
    for(int i=m;i>=0;i--)
        if(Coefficient[i]!=0) return i;
    return -1;
}
CPolynomial CPolynomial::operator+(CPolynomial b){
    CPolynomial R=*this;
    for(int i=0;i<PolynomialSize;i++)
        R[i]=R[i]+b[i];
    return R;
}
CPolynomial CPolynomial::operator<<(int b){
    CPolynomial N=*this;
    int i=PolynomialSize-1;
    for(;i>=b;i--)
        N[i]=N[i-b];
    for(;i>=0;i--)
        N[i]=0;
    return N;
}
CPolynomial CPolynomial::operator>>(int b){
    CPolynomial N=*this;
    int i=0;
    for(;i+b<PolynomialSize;i++)
        N[i]=N[i+b];
    for(;i<PolynomialSize;i++)
        N[i]=0;
    return N;
}
CGalois &CPolynomial::operator[](Code b){
    //if(b>PolynomialSize) return Coefficient[0];
    return Coefficient[b];
}
CPolynomial CPolynomial::operator*(CGalois b){
    CPolynomial N=*this;
    for(int i=0;i<PolynomialSize;i++){
        N[i]*=b;
    }
    return N;
}
CPolynomial CPolynomial::operator/(CGalois b){
    CPolynomial N=*this;
    for(int i=0;i<PolynomialSize;i++){
        N[i]/=b;
    }
    return N;
}
CPolynomial CPolynomial::operator*(CPolynomial b){
    CPolynomial N;
    for(int i=0;i<PolynomialSize;i++){
        if(Coefficient[i]!=0){
            N=N+((b<<i)*Coefficient[i]);
        }
    }
    return N;
}
CPolynomial CPolynomial::operator%(CPolynomial b){
    CPolynomial N=*this;
    int nmax=N.get_max();
    int bmax=b.get_max();
    int de=nmax-bmax;
    while(de>=0){
        CGalois n=N[nmax]/b[bmax];
        N=N+(b<<de)*n;
        nmax=N.get_max(nmax);
        //bmax=b.get_max(bmax);
        de=nmax-bmax;
    }
    return N;
}
CPolynomial CPolynomial::operator/(CPolynomial b){
    CPolynomial N=*this;
    CPolynomial ANS;
    int nmax=N.get_max();
    int bmax=b.get_max();
    int de=nmax-bmax;
    while(de>=0){
        CGalois n=N[nmax]/b[bmax];
        ANS[de]=n;
        N=N+(b<<de)*n;
        nmax=N.get_max(nmax);
        //bmax=b.get_max(bmax);
        de=nmax-bmax;
    }
    return ANS;
}
bool CPolynomial::operator==(CPolynomial b){
    for(int i=0;i<PolynomialSize;i++)
        if(this->Coefficient[i]!=b[i]) return false;
    return true;
}
bool CPolynomial::operator!=(CPolynomial b){
    for(int i=0;i<PolynomialSize;i++)
        if(this->Coefficient[i]!=b[i]) return true;
    return false;
}
CGalois CPolynomial::inputx(CGalois b){
    CGalois X=0;
    CGalois a=1;
    for(int i=0;i<PolynomialSize;i++){
        X=X+Coefficient[i]*a;
        a=a*b;
    }
    return X;
}
CPolynomial CPolynomial::differential(){
    CPolynomial D=*this;
    for(int i=0;i<PolynomialSize;i+=2)
        D[i]=0;
    for(int i=0;i<PolynomialSize-1;i++){
        D[i]=D[i+1];
    }
    D[PolynomialSize-1]=0;
    return D;
}
unsigned int CPolynomial::distance(CPolynomial in){
    unsigned int count=0;
    for(int i=0;i<PolynomialSize;i++){
        if(in[i]!=Coefficient[i]) count++;
    }
    return count;
}
unsigned int CPolynomial::weight(){
    unsigned int count=0;
    for(int i=0;i<PolynomialSize;i++){
        if(Coefficient[i]!=0) count++;
    }
    return count;
}
bool CPolynomial::empty(){
    for(int i=0;i<PolynomialSize;i++){
        if(Coefficient[i]!=0) return false;
    }
    return true;
}
std::ostream &operator<<(std::ostream &out,CPolynomial b){
    bool nzeor=false;
    for(int i=PolynomialSize-1;i>0;i--){
        if(nzeor){
            if(b[i]==0)
                out<<"\033[46m"<<b[i]<<"\033[49m"<<"X^"<<i<<"+";
            else
                out<<"\033[44m"<<b[i]<<"\033[49m"<<"X^"<<i<<"+";
        }else if(b[i]!=0){
            nzeor=true;
            out<<"\033[44m"<<b[i]<<"\033[49m"<<"X^"<<i<<"+";
        }
    }
    if(b[0]==0)
        out<<"\033[46m"<<b[0]<<"\033[49m";
    else
        out<<"\033[44m"<<b[0]<<"\033[49m";
    return out;
}

