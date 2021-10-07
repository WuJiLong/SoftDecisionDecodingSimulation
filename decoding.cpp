#include "decoding.hpp"

CPolynomial generateGX(int t,int b,CGalois a){
    CPolynomial GXPT;GXPT[0]=1;
    CGalois m=a^b;
    for(int i=b;i<(t*2)+b;i++){
        CPolynomial MT;MT[0] = m;MT[1] = 1;
        GXPT=GXPT*MT;
        m=m*a;
    }
    return GXPT;
}
CPolynomial getMessage(int n){
    CPolynomial M;
    for(int i=0;i<n;i++){
        M[i]=rand()%(MASK+1);
    }
    return M;
}
CPolynomial RSencode(CPolynomial M,CPolynomial &GX,int t){
    CPolynomial C=M<<(t*2);
    return C+C%GX;
}
CPolynomial getNose(int n){
    CPolynomial E;
    for(int i=0;i<n;i++){
        int randindex=rand()%MASK;
        if(E[randindex]==0)
            E[randindex]=rand()%(MASK+1);
        else
            i--;
    }
    return E;
}
CPolynomial Euclidean_Algorithm(CPolynomial &R,bool *pass){
    CGalois Si[T*2+1];
    for(int i=0;i<T*2+1;i++){
        Si[i]=R.inputx(1ull<<i);
    }
    bool ckerr=true;
    for(int i=1;i<=T*2;i++){
        if(Si[i]!=0){
            ckerr=false;
            break;
        }
    }
    if(pass) *pass=true;
    if(ckerr)
        return R;
    CPolynomial S;
    for(int i=0;i<2*T;i++)
        S[i]=Si[i+1];
    CPolynomial r_n1;r_n1[2*T]=1;
    CPolynomial r_0=S;
    CPolynomial t_n1;
    CPolynomial t_0;t_0[0]=1;
    while(r_0.get_max()>=T){
        CPolynomial q=r_n1/r_0;
        CPolynomial r_i=r_n1+(r_0*q);
        CPolynomial t_i=t_n1+(t_0*q);
        r_n1=r_0; r_0=r_i;
        t_n1=t_0; t_0=t_i;
    }
    CPolynomial Lambda_=t_0.differential();
    CPolynomial E;
    CGalois a=1;
    for(int i=0;i<MASK;i++,a=a*2){
        if(t_0.inputx(a)==0){
            int errorloc=(MASK-i)%MASK;
            CGalois errorvalue=r_0.inputx(a)/Lambda_.inputx(a);
             
            E[errorloc]=errorvalue;
        }
    }
    #ifdef DEBUG
    cout<<"Euclidean Algorithm."<<endl;
    for(int i=0;i<T*2+1;i++){
        cout<<"S"<<i<<"="<<Si[i]<<" ";
    }cout<<endl;
    cout<<"S(x)="<<S<<endl;
    cout<<"A(x)="<<t_0<<endl;
    cout<<"O(x)="<<r_0<<endl;
    cout<<"A'(x)="<<Lambda_<<endl;
    cout<<"DE(x)="<<E<<endl;
    #endif
    CPolynomial DC=R+E;
    CGalois Sii[T*2+1];
    for(int i=1;i<T*2+1;i++){
        Sii[i]=DC.inputx(1ull<<i);
    }
    for(int i=1;i<=T*2;i++){
        if(Sii[i]!=0){
            if(pass) *pass=false;
            return R;
        }
    }
    if(pass) *pass=true;
    return DC;
}
CPolynomial Berlekamp_Massey(CPolynomial &R,bool* pass){
    CGalois Si[T*2+1];
    for(int i=0;i<T*2+1;i++){
        Si[i]=R.inputx(1ull<<i);
    }
    bool ckerr=true;
    for(int i=1;i<=T*2;i++){
        if(Si[i]!=0){
            ckerr=false;
            break;
        }
    }

    if(pass) *pass=true;
    if(ckerr)
        return R;
    CPolynomial Lambda;Lambda[0]=1;//A(x)=1
    Code Lu=0;
    CPolynomial B;B[0]=1;
    for(int u=0;u<=2*T-1;u++){
         CGalois du=Si[u+1];
         for(int i=u,j=1;i>=u-Lu+1;i--,j++){
             du=du+Lambda[j]*Si[i];
         }
         if(du==0){
             B=B<<1;
         }else if(2*Lu>u){//du!=0
            Lambda=Lambda+(B<<1)*du;
            B=B<<1;
         }else{//du!=0 and 2*Lu<=u
            CPolynomial Lam=Lambda;
            Lambda=Lambda+(B<<1)*du;
            Lu=u+1-Lu;
            B=Lam*(du^-1);
         }
    }


    CPolynomial S;
    S[0]=1;
    for(int i=0;i<2*T;i++)
        S[i]=Si[i+1];
    CPolynomial keyequation;
    keyequation[0]=1;
    keyequation=Lambda*S%(keyequation<<(2*T));
    CPolynomial Lambda_=Lambda.differential();

    //Chien's search
    CPolynomial E;
    CGalois a=1;
    int ecount=0;
    for(int i=0;i<MASK;i++,a=a*2){
        std::cout << a <<" "<<Lambda.inputx(a)<< std::endl;
        if(Lambda.inputx(a)==0){
            int errorloc=(MASK-i)%MASK;
            CGalois errorvalue=keyequation.inputx(a)/Lambda_.inputx(a);
            //std::cout << Lambda.inputx(a) <<" "<< errorloc << std::endl;
 
            //cout<<"error localtion is "<<errorloc<<",error value is "<<errorvalue<<endl;
            E[errorloc]=errorvalue;
            if(errorvalue!=0)
                ecount++;
        }
    }
    #ifdef DEBUG
    cout<<"Berlekamp Massey"<<endl;
    for(int i=0;i<T*2+1;i++){
        cout<<"S"<<i<<"="<<Si[i]<<" ";
    }cout<<endl;
    cout<<"S(x)="<<S<<endl;
    cout<<"A(x)="<<Lambda<<endl;
    cout<<"O(x)="<<keyequation<<endl;
    cout<<"A'(x)="<<Lambda_<<endl;
    cout<<"DE(x)="<<E<<endl;
    #endif
    CPolynomial DC=R+E;
    CGalois Sii[T*2+1];
    for(int i=1;i<T*2+1;i++){
        Sii[i]=DC.inputx(1ull<<i);
    }
    for(int i=1;i<=T*2;i++){
        if(Sii[i]!=0){
            if(pass) *pass=false;
            return R;
        }
    }
    if(pass) *pass=true;
    return DC;
}