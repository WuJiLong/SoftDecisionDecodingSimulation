#include "algorithm.hpp"
#include "decoding.hpp"
#include "global.hpp"
#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <algorithm>

double MLE(CPolynomial &R,CPolynomial &C,CSignal &y){
    double sum=0;
    CModulation* p1=new_gen_modulation(y.modulation_type(),0);
    CModulation* p2=new_gen_modulation(y.modulation_type(),0);
    for(int i=0;i<MASK;i++){
        if(R[i]!=C[i]){
            //std::cout<<i<<"  ";
            p1->set_code(C[i].getNum());
            p2->set_code(R[i].getNum());
            sum += (y[i]->distance_pow2(p1) - y[i]->distance_pow2(p2));
            //std::cout<<(y[i]->distance_pow2(p1))<<std::endl;
            //if(y[i]->distance_pow2(p2)!=0) exit(87);
            //if(y[i]->distance_pow2(p1)<0) exit(78);
            //std::cout<<i<<":"<<(y[i]->distance_pow2(p1) - y[i]->distance_pow2(p2))<<std::endl;
        }
    }
    delete p1;
    delete p2;
    return sum;
}
double ErwinChu_B(CPolynomial &R,CPolynomial &C,CSignal &y,int *index){
    double sum=0;
    int d = R.distance(C);
    int c = 0;
    int i = 0;
    while(d+c<T*2+1){
        if(R[index[i]]==C[index[i]]){
            sum+=y[index[i]]->reliability();
            c++;
        }
        i++;
    }
    return sum;
}
double ErwinChu2_B(CPolynomial &R,CPolynomial V1,CPolynomial V2,CSignal &y,int *index){
    double sum=0;
    int d1 = (T*2+1) - R.distance(V1);
    int d2 = (T*2+1) - R.distance(V2);
    if(d1<d2){
        CPolynomial tmp = V1;V1=V2;V2=tmp;
        int t=d1;d1=d2;d2=t;
    }
    //std::cout<<d1<<" "<<d2<<std::endl;
    int x = d1-d2;
    int cnt=0;
    int i=0;
    while(cnt<d1 && i<MASK){
        int ind=index[i];
        if(R[ind]==V1[ind] && R[ind]==V2[ind]){
            sum+=y[ind]->reliability();
            cnt++;
        }else if(R[ind]==V1[ind] && R[ind]!=V2[ind]&& x>0){
            sum+=y[ind]->reliability();
            cnt++;
            x--;
        }
        i++;
    }
    //std::cout<<V1<<std::endl;
    //std::cout<<V2<<std::endl;
    return sum;
}

double ErwinChu3_B(CPolynomial &R,CPolynomial V1,CPolynomial V2,CSignal &y,int *index){
    double sum=0;
    int d1 = (T*2+1) - R.distance(V1);
    int d2 = (T*2+1) - R.distance(V2);
    if(d1<d2){
        CPolynomial tmp = V1;V1=V2;V2=tmp;
        int t=d1;d1=d2;d2=t;
    }
    //std::cout<<d1<<" "<<d2<<std::endl;
    int x = d1-d2;
    int x2 = (d1-d2)/2;
    int cnt=0;
    int i=0;
    while(cnt<d1 && i<MASK){
        int ind=index[i];
        if(R[ind]==V1[ind] && R[ind]==V2[ind]){
            sum+=y[ind]->reliability();
            cnt++;
        }else if(R[ind]==V1[ind] && R[ind]!=V2[ind]&& x>0){
            if(y[ind]->get_code2()==V2[ind].getNum() && x2>0){
                sum+=y[ind]->reliability();
                x2--;
                cnt++;
                x--;
            }else if(y[ind]->get_code2()!=V2[ind].getNum() ){
                sum+=y[ind]->reliability();
                cnt++;
                x--;
            }
        }
        i++;
    }
    //std::cout<<V1<<std::endl;
    //std::cout<<V2<<std::endl;
    return sum;
}

CPolynomial algorithm_EC2(CPolynomial &R,CSignal &y,unsigned long long &decode_count,bool *pass){
    std::vector<CPolynomial> all_dc;
    bool hddpass=false;
    decode_count=decode_count+1;
    CPolynomial Algebraic = Berlekamp_Massey(R,&hddpass);
    if(hddpass) all_dc.push_back(Algebraic);
    int index[MASK];
    for(int i=0;i<MASK;i++) index[i]=i;
    for(int i=0;i<MASK;i++){
        for(int j=i+1;j<MASK;j++)
            if(y[index[i]]->reliability()>y[index[j]]->reliability()){
                int tmp=index[i];
                index[i]=index[j];
                index[j]=tmp;
            }
    } 
    CPolynomial N;
    for(Code i=1;i<pow(MASK+1,T);i++){
        for(int j=0;j<T;j++){
           N[index[j]] = (i>>(j*_m_))&MASK;
        }
        CPolynomial RpN=R+N;
        bool ispass=false;
        for(std::vector<CPolynomial>::iterator it=all_dc.begin();it!=all_dc.end();it++){
            if(it->distance(RpN)<=T){
                ispass=true;
                break;
            }
        }
        if(ispass)
            continue;
        decode_count=decode_count+1;
        CPolynomial ddc=Berlekamp_Massey(RpN,&ispass);
        if(ispass) all_dc.push_back(ddc);
    }
    if(pass) *pass=true;
    if(all_dc.size()==0){
        if(pass) *pass=false;
        return R;
    }
    std::vector<CPolynomial>::iterator minDC;
    double min=-1;
    for(std::vector<CPolynomial>::iterator it=all_dc.begin();it!=all_dc.end();it++){
        double Ed=MLE(R,*it,y);//y.euclidean_distance(*it);
        if(min==-1 || min > Ed){
            min=Ed;
            minDC=it;
        }
    }
    return *minDC;
}
CPolynomial algorithm_CHU(CPolynomial &R,CSignal &y,unsigned long long &decode_count,bool *pass){
    std::vector<CPolynomial> all_dc;
    bool hddpass=false;
    decode_count=decode_count+1;
    CPolynomial Algebraic = Berlekamp_Massey(R,&hddpass);
    if(hddpass) all_dc.push_back(Algebraic);

    int index[MASK];
    for(int i=0;i<MASK;i++) index[i]=i;
    for(int i=0;i<MASK;i++){
        for(int j=i+1;j<MASK;j++)
            if(y[index[i]]->reliability()>y[index[j]]->reliability()){
                int tmp=index[i];
                index[i]=index[j];
                index[j]=tmp;
            }
    }
    if(hddpass){
        if(MLE(R,Algebraic,y)<=ErwinChu_B(R,Algebraic,y,index)){
            if(pass) *pass=true;
            return Algebraic;
        }
    }
    CPolynomial N;
    for(Code i=1;i<pow(MASK+1,T);i++){
        for(int j=0;j<T;j++){
           N[index[j]] = (i>>(j*_m_))&MASK;
        }
        CPolynomial RpN=R+N;
        bool ispass=false;
        decode_count=decode_count+1;
        CPolynomial ddc=Berlekamp_Massey(RpN,&ispass);
        if(ispass){
            if(MLE(R,ddc,y)<=ErwinChu_B(R,ddc,y,index)){
                if(pass) *pass=true;
                return ddc;
            }
            all_dc.push_back(ddc);
        }
    }
    if(pass) *pass=true;
    if(all_dc.size()==0){
        if(pass) *pass=false;
        return R;
    }
    std::vector<CPolynomial>::iterator minDC;
    double min=-1;
    for(std::vector<CPolynomial>::iterator it=all_dc.begin();it!=all_dc.end();it++){
        double Ed=MLE(R,*it,y);//y.euclidean_distance(*it);
        if(min==-1 || min > Ed){
            min=Ed;
            minDC=it;
        }
    }
    return *minDC;
}
CPolynomial algorithm_CHU_and_EC2(CPolynomial &R,CSignal &y,unsigned long long  &decode_count,bool *pass){
    std::vector<CPolynomial> all_dc;
    bool hddpass=false;
    decode_count=decode_count+1;
    CPolynomial Algebraic = Berlekamp_Massey(R,&hddpass);
    if(hddpass) all_dc.push_back(Algebraic);

    int index[MASK];
    for(int i=0;i<MASK;i++) index[i]=i;
    for(int i=0;i<MASK;i++){
        for(int j=i+1;j<MASK;j++)
            if(y[index[i]]->reliability()>y[index[j]]->reliability()){
                int tmp=index[i];
                index[i]=index[j];
                index[j]=tmp;
            }
    }

    if(hddpass){
        if(MLE(R,Algebraic,y)<=ErwinChu_B(R,Algebraic,y,index)){
            if(pass) *pass=true;
            return Algebraic;
        }
    }
    CPolynomial N;
    for(Code i=1;i<pow(MASK+1,T);i++){
        for(int j=0;j<T;j++){
           N[index[j]] = (i>>(j*_m_))&MASK;
        }
        CPolynomial RpN=R+N;
        bool ispass=false;
        for(std::vector<CPolynomial>::iterator it=all_dc.begin();it!=all_dc.end();it++){
            if(it->distance(RpN)<=T){
                ispass=true;
                break;
            }
        }
        if(ispass)
            continue;
        decode_count=decode_count+1;
        ispass=true;
        CPolynomial ddc=Berlekamp_Massey(RpN,&ispass);
        if(ispass){
            if(MLE(R,ddc,y)<=ErwinChu_B(R,ddc,y,index)){
                if(pass) *pass=true;
                return ddc;
            }
            all_dc.push_back(ddc);
        }
    }
    if(pass) *pass=true;
    if(all_dc.size()==0){
        if(pass) *pass=false;
        return R;
    }
    std::vector<CPolynomial>::iterator minDC;
    double min=-1;
    for(std::vector<CPolynomial>::iterator it=all_dc.begin();it!=all_dc.end();it++){
        //PSK16_Signal D(*it);
        double Ed=MLE(R,*it,y);//y.euclidean_distance(*it);
        if(min==-1 || min > Ed){
            min=Ed;
            minDC=it;
        }
    }
    return *minDC;
}
CPolynomial algorithm_CHU2(CPolynomial &R,CSignal &y,unsigned long long &decode_count,bool *pass){
    int index[MASK];
    for(int i=0;i<MASK;i++) index[i]=i;
    for(int i=0;i<MASK;i++){
        for(int j=i+1;j<MASK;j++)
            if(y[index[i]]->reliability()>y[index[j]]->reliability()){
                int tmp=index[i];
                index[i]=index[j];
                index[j]=tmp;
            }
    }
    bool V1_valid=false;
    decode_count++;
    CPolynomial V1=Berlekamp_Massey(R,&V1_valid);
    double V1_lambda=MLE(R,V1,y);
    CPolynomial N;
    if(V1_valid){
        if(V1_lambda<=ErwinChu_B(R,V1,y,index)){
            if(pass) *pass=true;
            return V1;
        }
    }
    for(Code i=1;i<pow(MASK+1,T);i++){
        for(int j=0;j<T;j++){
           N[index[j]] = (i>>(j*_m_))&MASK;
        }
        CPolynomial RpN=R+N;
        bool ispass=true;
        decode_count++;
        CPolynomial V2=Berlekamp_Massey(RpN,&ispass);
        if(ispass){
            if(V1_valid){
                double V2_lambda=MLE(R,V2,y);
                if(V2_lambda<V1_lambda){//swap
                    V1_lambda=V2_lambda;
                    CPolynomial tmp=V1; V1=V2; V2=tmp;
                }
                if(  V1_lambda <= ErwinChu2_B(R,V1,V2,y,index) ){
                    return V1;
                }
            }else{
                V1=V2;
                V1_lambda=MLE(R,V1,y);
                V1_valid=true;
                if(V1_lambda<=ErwinChu_B(R,V1,y,index)){
                    if(pass) *pass=true;
                    return V1;
                }
            }
        }
    }
    if(V1_valid){
        return V1;
    }else{
        return R;
    }
}
// CPolynomial algorithm_CHU22(CPolynomial &R,CSignal &y,unsigned long long &decode_count,bool *pass){
//     int index[MASK];
//     for(int i=0;i<MASK;i++) index[i]=i;
//     for(int i=0;i<MASK;i++){
//         for(int j=i+1;j<MASK;j++)
//             if(y[index[i]]->reliability()>y[index[j]]->reliability()){
//                 int tmp=index[i];
//                 index[i]=index[j];
//                 index[j]=tmp;
//             }
//     }
//     bool V1_valid=false;
//     decode_count++;
//     CPolynomial V1=Berlekamp_Massey(R,&V1_valid);
//     double V1_lambda=MLE(R,V1,y);
//     CPolynomial N;
//     if(V1_valid){
//         if(V1_lambda<=ErwinChu_B(R,V1,y,index)){
//             if(pass) *pass=true;
//             return V1;
//         }
//     }
//     for(Code i=1;i<pow(MASK+1,T);i++){
//         for(int j=0;j<T;j++){
//            N[index[j]] = (i>>(j*_m_))&MASK;
//         }
//         CPolynomial RpN=R+N;
//         bool ispass=true;
//         decode_count++;
//         CPolynomial V2=Berlekamp_Massey(RpN,&ispass);
//         if(ispass){
//             if(V1_valid){
//                 double V2_lambda=MLE(R,V2,y);
//                 if(V2_lambda<V1_lambda){//swap
//                     V1_lambda=V2_lambda;
//                     CPolynomial tmp=V1; V1=V2; V2=tmp;
//                 }
//                 if(  V1_lambda <= ErwinChu2_B2(R,V1,V2,y,index) ){
//                     return V1;
//                 }
//             }else{
//                 V1=V2;
//                 V1_lambda=MLE(R,V1,y);
//                 V1_valid=true;
//                 if(V1_lambda<=ErwinChu_B(R,V1,y,index)){
//                     if(pass) *pass=true;
//                     return V1;
//                 }
//             }
//         }
//     }
//     if(V1_valid){
//         return V1;
//     }else{
//         return R;
//     }
// }
CPolynomial algorithm_CHU2_and_EC2(CPolynomial &R,CSignal &y,unsigned long long &decode_count,bool *pass){
    std::vector<CPolynomial> all_dc;
    int index[MASK];
    for(int i=0;i<MASK;i++) index[i]=i;
    for(int i=0;i<MASK;i++){
        for(int j=i+1;j<MASK;j++)
            if(y[index[i]]->reliability()>y[index[j]]->reliability()){
                int tmp=index[i];
                index[i]=index[j];
                index[j]=tmp;
            }
    }

    bool V1_valid=false;
    decode_count++;
    CPolynomial V1=Berlekamp_Massey(R,&V1_valid);
    double V1_lambda=MLE(R,V1,y);

    if(V1_valid){
        if(V1_lambda<=ErwinChu_B(R,V1,y,index)){
            if(pass) *pass=true;
            return V1;
        }
        all_dc.push_back(V1);
    }
    CPolynomial N;
    for(Code i=1;i<pow(MASK+1,T);i++){
        for(int j=0;j<T;j++){
           N[index[j]] = (i>>(j*_m_))&MASK;
        }
        CPolynomial RpN=R+N;
        bool ispass=false;
        CPolynomial V2;
        for(std::vector<CPolynomial>::iterator it=all_dc.begin();it!=all_dc.end();it++){
            if(it->distance(RpN)<=T){
                ispass=true;
                V2=*it;
                break;
            }
        }
        if(!ispass){
            ispass=true;
            decode_count++;
            V2=Berlekamp_Massey(RpN,&ispass);
        }
        if(ispass){
            all_dc.push_back(V2);
            if(V1_valid){
                double V2_lambda=MLE(R,V2,y);
                if(V2_lambda<V1_lambda){//swap
                    V1_lambda=V2_lambda;
                    CPolynomial tmp=V1; V1=V2; V2=tmp;
                }
                if(  V1_lambda <= ErwinChu2_B(R,V1,V2,y,index) ){
                    return V1;
                }
            }else{
                V1=V2;
                V1_lambda=MLE(R,V1,y);
                V1_valid=true;
                if(V1_lambda<=ErwinChu_B(R,V1,y,index)){
                    if(pass) *pass=true;
                    return V1;
                }
            }
        }
    }
    if(V1_valid){
        return V1;
    }else{
        return R;
    }
}

CPolynomial algorithm_CHU_SORT(CPolynomial &R,CSignal &y,unsigned long long &decode_count,bool *pass){
    std::vector<CPolynomial> all_dc;
    bool hddpass=false;
    decode_count=decode_count+1;
    CPolynomial Algebraic = Berlekamp_Massey(R,&hddpass);
    if(hddpass) all_dc.push_back(Algebraic);

    int index[MASK];
    for(int i=0;i<MASK;i++) index[i]=i;
    for(int i=0;i<MASK;i++){
        for(int j=i+1;j<MASK;j++)
            if(y[index[i]]->reliability()>y[index[j]]->reliability()){
                int tmp=index[i];
                index[i]=index[j];
                index[j]=tmp;
            }
    }
    if(hddpass){
        if(MLE(R,Algebraic,y)<=ErwinChu_B(R,Algebraic,y,index)){
            if(pass) *pass=true;
            //std::cout<<MLE(R,Algebraic,y)<<std::endl;
            return Algebraic;
        }
    }
    //CPolynomial N;
    for(Code i=1;i<pow(MASK+1,T);i++){
        CPolynomial RpN=R;
        for(int j=0;j<T;j++){
           RpN[index[j]] = y[index[j]]->code[(i>>(j*_m_))&MASK];
        }
        bool ispass=false;
        decode_count=decode_count+1;
        CPolynomial ddc=Berlekamp_Massey(RpN,&ispass);
        if(ispass){
            if(MLE(R,ddc,y)<=ErwinChu_B(R,ddc,y,index)){
                if(pass) *pass=true;
                return ddc;
                //std::cout<<MLE(R,ddc,y)<<std::endl;
            }
            all_dc.push_back(ddc);
        }
    }
    if(pass) *pass=true;
    if(all_dc.size()==0){
        if(pass) *pass=false;
        return R;
    }
    std::vector<CPolynomial>::iterator minDC;
    double min=-1;
    for(std::vector<CPolynomial>::iterator it=all_dc.begin();it!=all_dc.end();it++){
        double Ed=MLE(R,*it,y);//y.euclidean_distance(*it);
        if(min==-1 || min > Ed){
            min=Ed;
            minDC=it;
        }
    }
    //std::cout<<MLE(R,*minDC,y)<<std::endl;
    return *minDC;
}
CPolynomial algorithm_CHU_and_EC2_SORT(CPolynomial &R,CSignal &y,unsigned long long  &decode_count,bool *pass){
    std::vector<CPolynomial> all_dc;
    bool hddpass=false;
    decode_count=decode_count+1;
    CPolynomial Algebraic = Berlekamp_Massey(R,&hddpass);
    if(hddpass) all_dc.push_back(Algebraic);

    int index[MASK];
    for(int i=0;i<MASK;i++) index[i]=i;
    for(int i=0;i<MASK;i++){
        for(int j=i+1;j<MASK;j++)
            if(y[index[i]]->reliability()>y[index[j]]->reliability()){
                int tmp=index[i];
                index[i]=index[j];
                index[j]=tmp;
            }
    }

    if(hddpass){
        if(MLE(R,Algebraic,y)<=ErwinChu_B(R,Algebraic,y,index)){
            if(pass) *pass=true;
            return Algebraic;
        }
    }
    //CPolynomial N;
    for(Code i=1;i<pow(MASK+1,T);i++){
        CPolynomial RpN=R;
        for(int j=0;j<T;j++){
           RpN[index[j]] = y[index[j]]->code[(i>>(j*_m_))&MASK];
        }
        bool ispass=false;
        for(std::vector<CPolynomial>::iterator it=all_dc.begin();it!=all_dc.end();it++){
            if(it->distance(RpN)<=T){
                ispass=true;
                break;
            }
        }
        if(ispass)
            continue;
        decode_count=decode_count+1;
        ispass=true;
        CPolynomial ddc=Berlekamp_Massey(RpN,&ispass);
        if(ispass){
            if(MLE(R,ddc,y)<=ErwinChu_B(R,ddc,y,index)){
                if(pass) *pass=true;
                return ddc;
            }
            all_dc.push_back(ddc);
        }
    }
    if(pass) *pass=true;
    if(all_dc.size()==0){
        if(pass) *pass=false;
        return R;
    }
    std::vector<CPolynomial>::iterator minDC;
    double min=-1;
    for(std::vector<CPolynomial>::iterator it=all_dc.begin();it!=all_dc.end();it++){
        //PSK16_Signal D(*it);
        double Ed=MLE(R,*it,y);//y.euclidean_distance(*it);
        if(min==-1 || min > Ed){
            min=Ed;
            minDC=it;
        }
    }
    return *minDC;
}
CPolynomial algorithm_CHU2_SORT(CPolynomial &R,CSignal &y,unsigned long long &decode_count,bool *pass){
    int index[MASK];
    for(int i=0;i<MASK;i++) index[i]=i;
    for(int i=0;i<MASK;i++){
        for(int j=i+1;j<MASK;j++)
            if(y[index[i]]->reliability()>y[index[j]]->reliability()){
                int tmp=index[i];
                index[i]=index[j];
                index[j]=tmp;
            }
    }
    bool V1_valid=false;
    decode_count++;
    CPolynomial V1=Berlekamp_Massey(R,&V1_valid);
    double V1_lambda=MLE(R,V1,y);
    //CPolynomial N;
    if(V1_valid){
        if(V1_lambda<=ErwinChu_B(R,V1,y,index)){
            if(pass) *pass=true;
            //std::cout<<MLE(R,V1,y)<<std::endl;
            return V1;
        }
    }
    //CPolynomial N;
    for(Code i=1;i<pow(MASK+1,T);i++){
        CPolynomial RpN=R;
        for(int j=0;j<T;j++){
           RpN[index[j]] = y[index[j]]->code[(i>>(j*_m_))&MASK];
        }
        bool ispass=true;
        decode_count++;
        CPolynomial V2=Berlekamp_Massey(RpN,&ispass);
        if(ispass){
            if(V1_valid){
                double V2_lambda=MLE(R,V2,y);
                if(V2_lambda<V1_lambda){//swap
                    V1_lambda=V2_lambda;
                    CPolynomial tmp=V1; V1=V2; V2=tmp;
                }
                if(  V1_lambda <= ErwinChu2_B(R,V1,V2,y,index) ){
                    //std::cout<<V2<<std::endl;
                    //std::cout<<"FUCK:"<<MLE(R,V1,y)<<std::endl;
                    //std::cout<<"FUCK:"<<MLE(R,V2,y)<<std::endl;
                    //std::cout<<"FUCK:"<< ErwinChu2_B(R,V1,V2,y,index)<<std::endl;
                    return V1;
                }
            }else{
                V1=V2;
                V1_lambda=MLE(R,V1,y);
                V1_valid=true;
                if(V1_lambda<=ErwinChu_B(R,V1,y,index)){
                    if(pass) *pass=true;
                    //std::cout<<MLE(R,V1,y)<<std::endl;
                    return V1;
                }
            }
        }
    }
    if(V1_valid){
        //std::cout<<MLE(R,V1,y)<<std::endl;
        return V1;
    }else{
        return R;
    }
}
CPolynomial algorithm_CHU2_and_EC2_SORT(CPolynomial &R,CSignal &y,unsigned long long &decode_count,bool *pass){
    std::vector<CPolynomial> all_dc;
    int index[MASK];
    for(int i=0;i<MASK;i++) index[i]=i;
    for(int i=0;i<MASK;i++){
        for(int j=i+1;j<MASK;j++)
            if(y[index[i]]->reliability()>y[index[j]]->reliability()){
                int tmp=index[i];
                index[i]=index[j];
                index[j]=tmp;
            }
    }
    decode_count++;
    bool V1_valid=false;
    CPolynomial V1=Berlekamp_Massey(R,&V1_valid);
    double V1_lambda=MLE(R,V1,y);

    if(V1_valid){
        if(V1_lambda<=ErwinChu_B(R,V1,y,index)){
            if(pass) *pass=true;
            return V1;
        }
        all_dc.push_back(V1);
    }
    //CPolynomial N;
    for(Code i=1;i<pow(MASK+1,T);i++){
        CPolynomial RpN=R;
        for(int j=0;j<T;j++){
           RpN[index[j]] = y[index[j]]->code[(i>>(j*_m_))&MASK];
        }
        bool ispass=false;
        CPolynomial V2;
        for(std::vector<CPolynomial>::iterator it=all_dc.begin();it!=all_dc.end();it++){
            if(it->distance(RpN)<=T){
                ispass=true;
                V2=*it;
                break;
            }
        }
        if(!ispass){
            ispass=true;
            decode_count++;
            V2=Berlekamp_Massey(RpN,&ispass);
        }
        if(ispass){
            all_dc.push_back(V2);
            if(V1_valid){
                double V2_lambda=MLE(R,V2,y);
                if(V2_lambda<V1_lambda){//swap
                    V1_lambda=V2_lambda;
                    CPolynomial tmp=V1; V1=V2; V2=tmp;
                }
                if(  V1_lambda <= ErwinChu2_B(R,V1,V2,y,index) ){
                    return V1;
                }
            }else{
                V1=V2;
                V1_lambda=MLE(R,V1,y);
                V1_valid=true;
                if(V1_lambda<=ErwinChu_B(R,V1,y,index)){
                    if(pass) *pass=true;
                    return V1;
                }
            }
        }
    }
    if(V1_valid){
        return V1;
    }else{
        return R;
    }
}





CPolynomial algorithm_SCA_EC2(Code pattern,CPolynomial &R,CSignal &y,unsigned long long &decode_count,bool *pass){
    std::vector<CPolynomial> all_dc;
    bool hddpass=false;
    decode_count=decode_count+1;
    CPolynomial Algebraic = Berlekamp_Massey(R,&hddpass);
    if(hddpass) all_dc.push_back(Algebraic);
    int index[MASK];
    for(int i=0;i<MASK;i++) index[i]=i;
    for(int i=0;i<MASK;i++){
        for(int j=i+1;j<MASK;j++)
            if(y[index[i]]->reliability()>y[index[j]]->reliability()){
                int tmp=index[i];
                index[i]=index[j];
                index[j]=tmp;
            }
    }
    CPolynomial N;
    for(Code i=1;i<pattern;i++){
        for(int j=0;j<T;j++){
           N[index[j]] = (i>>(j*_m_))&MASK;
        }
        CPolynomial RpN=R+N;
        bool ispass=false;
        for(std::vector<CPolynomial>::iterator it=all_dc.begin();it!=all_dc.end();it++){
            if(it->distance(RpN)<=T){
                ispass=true;
                break;
            }
        }
        if(ispass)
            continue;
        decode_count=decode_count+1;
        CPolynomial ddc=Berlekamp_Massey(RpN,&ispass);
        if(ispass) all_dc.push_back(ddc);
    }
    if(pass) *pass=true;
    if(all_dc.size()==0){
        if(pass) *pass=false;
        return R;
    }
    std::vector<CPolynomial>::iterator minDC;
    double min=-1;
    for(std::vector<CPolynomial>::iterator it=all_dc.begin();it!=all_dc.end();it++){
        double Ed=MLE(R,*it,y);//y.euclidean_distance(*it);
        if(min==-1 || min > Ed){
            min=Ed;
            minDC=it;
        }
    }
    return *minDC;
}
CPolynomial algorithm_SCA_CHU(Code pattern,CPolynomial &R,CSignal &y,unsigned long long &decode_count,bool *pass){
    std::vector<CPolynomial> all_dc;
    bool hddpass=false;
    decode_count=decode_count+1;
    CPolynomial Algebraic = Berlekamp_Massey(R,&hddpass);
    if(hddpass) all_dc.push_back(Algebraic);

    int index[MASK];
    for(int i=0;i<MASK;i++) index[i]=i;
    for(int i=0;i<MASK;i++){
        for(int j=i+1;j<MASK;j++)
            if(y[index[i]]->reliability()>y[index[j]]->reliability()){
                int tmp=index[i];
                index[i]=index[j];
                index[j]=tmp;
            }
    }
    if(hddpass){
        if(MLE(R,Algebraic,y)<=ErwinChu_B(R,Algebraic,y,index)){
            if(pass) *pass=true;
            return Algebraic;
        }
    }
    CPolynomial N;
    for(Code i=1;i<pattern;i++){
        for(int j=0;j<T;j++){
           N[index[j]] = (i>>(j*_m_))&MASK;
        }
        CPolynomial RpN=R+N;
        bool ispass=false;
        decode_count=decode_count+1;
        CPolynomial ddc=Berlekamp_Massey(RpN,&ispass);
        if(ispass){
            if(MLE(R,ddc,y)<=ErwinChu_B(R,ddc,y,index)){
                if(pass) *pass=true;
                return ddc;
            }
            all_dc.push_back(ddc);
        }
    }
    if(pass) *pass=true;
    if(all_dc.size()==0){
        if(pass) *pass=false;
        return R;
    }
    std::vector<CPolynomial>::iterator minDC;
    double min=-1;
    for(std::vector<CPolynomial>::iterator it=all_dc.begin();it!=all_dc.end();it++){
        double Ed=MLE(R,*it,y);//y.euclidean_distance(*it);
        if(min==-1 || min > Ed){
            min=Ed;
            minDC=it;
        }
    }
    return *minDC;
}
CPolynomial algorithm_SCA_CHU_and_EC2(Code pattern,CPolynomial &R,CSignal &y,unsigned long long  &decode_count,bool *pass){
    std::vector<CPolynomial> all_dc;
    bool hddpass=false;
    decode_count=decode_count+1;
    CPolynomial Algebraic = Berlekamp_Massey(R,&hddpass);
    if(hddpass) all_dc.push_back(Algebraic);

    int index[MASK];
    for(int i=0;i<MASK;i++) index[i]=i;
    for(int i=0;i<MASK;i++){
        for(int j=i+1;j<MASK;j++)
            if(y[index[i]]->reliability()>y[index[j]]->reliability()){
                int tmp=index[i];
                index[i]=index[j];
                index[j]=tmp;
            }
    }

    if(hddpass){
        if(MLE(R,Algebraic,y)<=ErwinChu_B(R,Algebraic,y,index)){
            if(pass) *pass=true;
            return Algebraic;
        }
    }
    CPolynomial N;
    for(Code i=1;i<pattern;i++){
        for(int j=0;j<T;j++){
           N[index[j]] = (i>>(j*_m_))&MASK;
        }
        CPolynomial RpN=R+N;
        bool ispass=false;
        for(std::vector<CPolynomial>::iterator it=all_dc.begin();it!=all_dc.end();it++){
            if(it->distance(RpN)<=T){
                ispass=true;
                break;
            }
        }
        if(ispass)
            continue;
        decode_count=decode_count+1;
        ispass=true;
        CPolynomial ddc=Berlekamp_Massey(RpN,&ispass);
        if(ispass){
            if(MLE(R,ddc,y)<=ErwinChu_B(R,ddc,y,index)){
                if(pass) *pass=true;
                return ddc;
            }
            all_dc.push_back(ddc);
        }
    }
    if(pass) *pass=true;
    if(all_dc.size()==0){
        if(pass) *pass=false;
        return R;
    }
    std::vector<CPolynomial>::iterator minDC;
    double min=-1;
    for(std::vector<CPolynomial>::iterator it=all_dc.begin();it!=all_dc.end();it++){
        //PSK16_Signal D(*it);
        double Ed=MLE(R,*it,y);//y.euclidean_distance(*it);
        if(min==-1 || min > Ed){
            min=Ed;
            minDC=it;
        }
    }
    return *minDC;
}
CPolynomial algorithm_SCA_CHU2(Code pattern,CPolynomial &R,CSignal &y,unsigned long long &decode_count,bool *pass){
    int index[MASK];
    for(int i=0;i<MASK;i++) index[i]=i;
    for(int i=0;i<MASK;i++){
        for(int j=i+1;j<MASK;j++)
            if(y[index[i]]->reliability()>y[index[j]]->reliability()){
                int tmp=index[i];
                index[i]=index[j];
                index[j]=tmp;
            }
    }
    bool V1_valid=false;
    decode_count++;
    CPolynomial V1=Berlekamp_Massey(R,&V1_valid);
    double V1_lambda=MLE(R,V1,y);
    CPolynomial N;
    if(V1_valid){
        if(V1_lambda<=ErwinChu_B(R,V1,y,index)){
            if(pass) *pass=true;
            return V1;
        }
    }
    for(Code i=1;i<pattern;i++){
        for(int j=0;j<T;j++){
           N[index[j]] = (i>>(j*_m_))&MASK;
        }
        CPolynomial RpN=R+N;
        bool ispass=true;
        decode_count++;
        CPolynomial V2=Berlekamp_Massey(RpN,&ispass);
        if(ispass){
            if(V1_valid){
                double V2_lambda=MLE(R,V2,y);
                if(V2_lambda<V1_lambda){//swap
                    V1_lambda=V2_lambda;
                    CPolynomial tmp=V1; V1=V2; V2=tmp;
                }
                if(  V1_lambda <= ErwinChu2_B(R,V1,V2,y,index) ){
                    return V1;
                }
            }else{
                V1=V2;
                V1_lambda=MLE(R,V1,y);
                V1_valid=true;
                if(V1_lambda<=ErwinChu_B(R,V1,y,index)){
                    if(pass) *pass=true;
                    return V1;
                }
            }
        }
    }
    if(V1_valid){
        return V1;
    }else{
        return R;
    }
}
CPolynomial algorithm_SCA_CHU2_and_EC2(Code pattern,CPolynomial &R,CSignal &y,unsigned long long &decode_count,bool *pass){
    std::vector<CPolynomial> all_dc;
    int index[MASK];
    for(int i=0;i<MASK;i++) index[i]=i;
    for(int i=0;i<MASK;i++){
        for(int j=i+1;j<MASK;j++)
            if(y[index[i]]->reliability()>y[index[j]]->reliability()){
                int tmp=index[i];
                index[i]=index[j];
                index[j]=tmp;
            }
    }

    bool V1_valid=false;
    decode_count++;
    CPolynomial V1=Berlekamp_Massey(R,&V1_valid);
    double V1_lambda=MLE(R,V1,y);

    if(V1_valid){
        if(V1_lambda<=ErwinChu_B(R,V1,y,index)){
            if(pass) *pass=true;
            return V1;
        }
        all_dc.push_back(V1);
    }
    CPolynomial N;
    for(Code i=1;i<pattern;i++){
        for(int j=0;j<T;j++){
           N[index[j]] = (i>>(j*_m_))&MASK;
        }
        CPolynomial RpN=R+N;
        bool ispass=false;
        CPolynomial V2;
        for(std::vector<CPolynomial>::iterator it=all_dc.begin();it!=all_dc.end();it++){
            if(it->distance(RpN)<=T){
                ispass=true;
                V2=*it;
                break;
            }
        }
        if(!ispass){
            ispass=true;
            decode_count++;
            V2=Berlekamp_Massey(RpN,&ispass);
        }
        if(ispass){
            all_dc.push_back(V2);
            if(V1_valid){
                double V2_lambda=MLE(R,V2,y);
                if(V2_lambda<V1_lambda){//swap
                    V1_lambda=V2_lambda;
                    CPolynomial tmp=V1; V1=V2; V2=tmp;
                }
                if(  V1_lambda <= ErwinChu2_B(R,V1,V2,y,index) ){
                    return V1;
                }
            }else{
                V1=V2;
                V1_lambda=MLE(R,V1,y);
                V1_valid=true;
                if(V1_lambda<=ErwinChu_B(R,V1,y,index)){
                    if(pass) *pass=true;
                    return V1;
                }
            }
        }
    }
    if(V1_valid){
        return V1;
    }else{
        return R;
    }
}

CPolynomial algorithm_SSCA(Code pattern,CPolynomial &R,CSignal &y,unsigned long long &decode_count,bool *pass){
    std::vector<CPolynomial> all_dc;
    std::vector<CPolynomial> SWAP;
    bool hddpass=false;
    SWAP.push_back(R);
    decode_count=decode_count+1;
    CPolynomial Algebraic = Berlekamp_Massey(R,&hddpass);
    if(hddpass) all_dc.push_back(Algebraic);
    /* int index[MASK];
    for(int i=0;i<MASK;i++) index[i]=i;
    for(int i=0;i<MASK;i++){
        for(int j=i+1;j<MASK;j++)
            if(y[index[i]]->reliability()>y[index[j]]->reliability()){
                int tmp=index[i];
                index[i]=index[j];
                index[j]=tmp;
            }
    }*/
    //std::default_random_engine RAND_generator = std::default_random_engine( 0 );//seed=0
    for(Code i=1;i<pattern;i++){
        CPolynomial RpN=y.get_random_data(CGLOBAL::RAND_generator);//R+N;
        if(find(SWAP.begin(),SWAP.end(),RpN) != SWAP.end()){
            continue;
        }
        SWAP.push_back(RpN);
        bool ispass=false;
     
        decode_count=decode_count+1;
        CPolynomial ddc=Berlekamp_Massey(RpN,&ispass);
        if(ispass) all_dc.push_back(ddc);
    }
    if(pass) *pass=true;
    if(all_dc.size()==0){
        if(pass) *pass=false;
        return R;
    }
    std::vector<CPolynomial>::iterator minDC;
    double min=-1;
    for(std::vector<CPolynomial>::iterator it=all_dc.begin();it!=all_dc.end();it++){
        double Ed=MLE(R,*it,y);//y.euclidean_distance(*it);
        if(min==-1 || min > Ed){
            min=Ed;
            minDC=it;
        }
    }
    return *minDC;
}

CPolynomial algorithm_SSCA_EC2(Code pattern,CPolynomial &R,CSignal &y,unsigned long long &decode_count,bool *pass){
    std::vector<CPolynomial> all_dc;
    std::vector<CPolynomial> SWAP;
    bool hddpass=false;
    decode_count=decode_count+1;
    SWAP.push_back(R);
    CPolynomial Algebraic = Berlekamp_Massey(R,&hddpass);
    if(hddpass) all_dc.push_back(Algebraic);
    int index[MASK];
    for(int i=0;i<MASK;i++) index[i]=i;
    for(int i=0;i<MASK;i++){
        for(int j=i+1;j<MASK;j++)
            if(y[index[i]]->reliability()>y[index[j]]->reliability()){
                int tmp=index[i];
                index[i]=index[j];
                index[j]=tmp;
            }
    }
    //std::default_random_engine RAND_generator = std::default_random_engine( 0 );//seed=0
    for(Code i=1;i<pattern;i++){
        CPolynomial RpN=y.get_random_data(CGLOBAL::RAND_generator);//R+N;
        if(find(SWAP.begin(),SWAP.end(),RpN) != SWAP.end()){
            continue;
        }
        SWAP.push_back(RpN);
        bool ispass=false;
        for(std::vector<CPolynomial>::iterator it=all_dc.begin();it!=all_dc.end();it++){
            if(it->distance(RpN)<=T){
                ispass=true;
                break;
            }
        }
        if(ispass)
            continue;
        decode_count=decode_count+1;
        CPolynomial ddc=Berlekamp_Massey(RpN,&ispass);
        if(ispass) all_dc.push_back(ddc);
    }
    if(pass) *pass=true;
    if(all_dc.size()==0){
        if(pass) *pass=false;
        return R;
    }
    std::vector<CPolynomial>::iterator minDC;
    double min=-1;
    for(std::vector<CPolynomial>::iterator it=all_dc.begin();it!=all_dc.end();it++){
        double Ed=MLE(R,*it,y);//y.euclidean_distance(*it);
        if(min==-1 || min > Ed){
            min=Ed;
            minDC=it;
        }
    }
    return *minDC;
}

CPolynomial algorithm_SSCA_CHU(Code pattern,CPolynomial &R,CSignal &y,unsigned long long &decode_count,bool *pass){
    std::vector<CPolynomial> all_dc;
    std::vector<CPolynomial> SWAP;
    bool hddpass=false;
    decode_count=decode_count+1;
    SWAP.push_back(R);
    CPolynomial Algebraic = Berlekamp_Massey(R,&hddpass);
    if(hddpass) all_dc.push_back(Algebraic);

    int index[MASK];
    for(int i=0;i<MASK;i++) index[i]=i;
    for(int i=0;i<MASK;i++){
        for(int j=i+1;j<MASK;j++)
            if(y[index[i]]->reliability()>y[index[j]]->reliability()){
                int tmp=index[i];
                index[i]=index[j];
                index[j]=tmp;
            }
    }
    if(hddpass){
        if(MLE(R,Algebraic,y)<=ErwinChu_B(R,Algebraic,y,index)){
            if(pass) *pass=true;
            return Algebraic;
        }
    }
    for(Code i=1;i<pattern;i++){
        CPolynomial RpN=y.get_random_data(CGLOBAL::RAND_generator);//R+N;
        if(find(SWAP.begin(),SWAP.end(),RpN) != SWAP.end()){
            continue;
        }
        SWAP.push_back(RpN);
        bool ispass=false;
        decode_count=decode_count+1;
        CPolynomial ddc=Berlekamp_Massey(RpN,&ispass);
        if(ispass){
            if(MLE(R,ddc,y)<=ErwinChu_B(R,ddc,y,index)){
                if(pass) *pass=true;
                //for(i++;i<pattern;i++) RpN=y.get_random_data(CGLOBAL::RAND_generatorB);
                return ddc;
            }
            all_dc.push_back(ddc);
        }
    }
    if(pass) *pass=true;
    if(all_dc.size()==0){
        if(pass) *pass=false;
        return R;
    }
    std::vector<CPolynomial>::iterator minDC;
    double min=-1;
    for(std::vector<CPolynomial>::iterator it=all_dc.begin();it!=all_dc.end();it++){
        double Ed=MLE(R,*it,y);//y.euclidean_distance(*it);
        if(min==-1 || min > Ed){
            min=Ed;
            minDC=it;
        }
    }
    return *minDC;
}
CPolynomial algorithm_SSCA_CHU_and_EC2(Code pattern,CPolynomial &R,CSignal &y,unsigned long long  &decode_count,bool *pass){
    std::vector<CPolynomial> all_dc;
    std::vector<CPolynomial> SWAP;
    bool hddpass=false;
    decode_count=decode_count+1;
    SWAP.push_back(R);
    CPolynomial Algebraic = Berlekamp_Massey(R,&hddpass);
    if(hddpass) all_dc.push_back(Algebraic);

    int index[MASK];
    for(int i=0;i<MASK;i++) index[i]=i;
    for(int i=0;i<MASK;i++){
        for(int j=i+1;j<MASK;j++)
            if(y[index[i]]->reliability()>y[index[j]]->reliability()){
                int tmp=index[i];
                index[i]=index[j];
                index[j]=tmp;
            }
    }

    if(hddpass){
        if(MLE(R,Algebraic,y)<=ErwinChu_B(R,Algebraic,y,index)){
            if(pass) *pass=true;
            return Algebraic;
        }
    }
    for(Code i=1;i<pattern;i++){
        CPolynomial RpN=y.get_random_data(CGLOBAL::RAND_generator);//R+N;
        if(find(SWAP.begin(),SWAP.end(),RpN) != SWAP.end()){
            continue;
        }
        SWAP.push_back(RpN);
        bool ispass=false;
        for(std::vector<CPolynomial>::iterator it=all_dc.begin();it!=all_dc.end();it++){
            if(it->distance(RpN)<=T){
                ispass=true;
                break;
            }
        }
        if(ispass)
            continue;
        decode_count=decode_count+1;
        ispass=true;
        CPolynomial ddc=Berlekamp_Massey(RpN,&ispass);
        if(ispass){
            if(MLE(R,ddc,y)<=ErwinChu_B(R,ddc,y,index)){
                if(pass) *pass=true;
                return ddc;
            }
            all_dc.push_back(ddc);
        }
    }
    if(pass) *pass=true;
    if(all_dc.size()==0){
        if(pass) *pass=false;
        return R;
    }
    std::vector<CPolynomial>::iterator minDC;
    double min=-1;
    for(std::vector<CPolynomial>::iterator it=all_dc.begin();it!=all_dc.end();it++){
        //PSK16_Signal D(*it);
        double Ed=MLE(R,*it,y);//y.euclidean_distance(*it);
        if(min==-1 || min > Ed){
            min=Ed;
            minDC=it;
        }
    }
    return *minDC;
}
CPolynomial algorithm_SSCA_CHU2(Code pattern,CPolynomial &R,CSignal &y,unsigned long long &decode_count,bool *pass){
    std::vector<CPolynomial> SWAP;
    SWAP.push_back(R);
    int index[MASK];
    for(int i=0;i<MASK;i++) index[i]=i;
    for(int i=0;i<MASK;i++){
        for(int j=i+1;j<MASK;j++)
            if(y[index[i]]->reliability()>y[index[j]]->reliability()){
                int tmp=index[i];
                index[i]=index[j];
                index[j]=tmp;
            }
    }
    bool V1_valid=false;
    decode_count++;
    CPolynomial V1=Berlekamp_Massey(R,&V1_valid);
    double V1_lambda=MLE(R,V1,y);
    if(V1_valid){
        if(V1_lambda<=ErwinChu_B(R,V1,y,index)){
            if(pass) *pass=true;
            return V1;
        }
    }
    for(Code i=1;i<pattern;i++){
        CPolynomial RpN=y.get_random_data(CGLOBAL::RAND_generator);//R+N;
        if(find(SWAP.begin(),SWAP.end(),RpN) != SWAP.end()){
            continue;
        }
        SWAP.push_back(RpN);
        bool ispass=false;
        decode_count++;
        CPolynomial V2=Berlekamp_Massey(RpN,&ispass);
        if(ispass){
            if(V1_valid){
                double V2_lambda=MLE(R,V2,y);
            
                if(V2_lambda<V1_lambda){//swap
                    V1_lambda=V2_lambda;
                    CPolynomial tmp=V1; V1=V2; V2=tmp;
                }
                if(  V1_lambda <= ErwinChu2_B(R,V1,V2,y,index) ){
                    return V1;
                }
            }else{
                V1=V2;
                V1_lambda=MLE(R,V1,y);
                V1_valid=true;
                if(V1_lambda<=ErwinChu_B(R,V1,y,index)){
                    if(pass) *pass=true;
                    return V1;
                }
            }
        }
    }
    if(V1_valid){
        //std::cout<<"fuck"<<std::endl;
        return V1;
    }else{
        return R;
    }
}
CPolynomial algorithm_SSCA_CHU2_and_EC2(Code pattern,CPolynomial &R,CSignal &y,unsigned long long &decode_count,bool *pass){
    std::vector<CPolynomial> all_dc;
    std::vector<CPolynomial> SWAP;
    SWAP.push_back(R);
    int index[MASK];
    for(int i=0;i<MASK;i++) index[i]=i;
    for(int i=0;i<MASK;i++){
        for(int j=i+1;j<MASK;j++)
            if(y[index[i]]->reliability()>y[index[j]]->reliability()){
                int tmp=index[i];
                index[i]=index[j];
                index[j]=tmp;
            }
    }

    bool V1_valid=false;
    decode_count++;
    CPolynomial V1=Berlekamp_Massey(R,&V1_valid);
    double V1_lambda=MLE(R,V1,y);

    if(V1_valid){
        if(V1_lambda<=ErwinChu_B(R,V1,y,index)){
            if(pass) *pass=true;
            return V1;
        }
        all_dc.push_back(V1);
    }
    for(Code i=1;i<pattern;i++){
        CPolynomial RpN=y.get_random_data(CGLOBAL::RAND_generator);//R+N;
        if(find(SWAP.begin(),SWAP.end(),RpN) != SWAP.end()){
            continue;
        }
        SWAP.push_back(RpN);
        bool ispass=false;
        CPolynomial V2;
        for(std::vector<CPolynomial>::iterator it=all_dc.begin();it!=all_dc.end();it++){
            if(it->distance(RpN)<=T){
                ispass=true;
                V2=*it;
                break;
            }
        }
        if(!ispass){
            ispass=true;
            decode_count++;
            V2=Berlekamp_Massey(RpN,&ispass);
        }
        if(ispass){
            all_dc.push_back(V2);
            if(V1_valid){
                double V2_lambda=MLE(R,V2,y);
                if(V2_lambda<V1_lambda){//swap
                    V1_lambda=V2_lambda;
                    CPolynomial tmp=V1; V1=V2; V2=tmp;
                }
                if(  V1_lambda <= ErwinChu2_B(R,V1,V2,y,index) ){
                    return V1;
                }
            }else{
                V1=V2;
                V1_lambda=MLE(R,V1,y);
                V1_valid=true;
                if(V1_lambda<=ErwinChu_B(R,V1,y,index)){
                    if(pass) *pass=true;
                    return V1;
                }
            }
        }
    }
    if(V1_valid){
        return V1;
    }else{
        return R;
    }
}

//-----------------------------------------------------------------
//-----------------------------------------------------------------
//-----------------------------------------------------------------


CPolynomial algorithm_CHU3(CPolynomial &R,CSignal &y,unsigned long long &decode_count,bool *pass){
    int index[MASK];
    for(int i=0;i<MASK;i++) index[i]=i;
    for(int i=0;i<MASK;i++){
        for(int j=i+1;j<MASK;j++)
            if(y[index[i]]->reliability()>y[index[j]]->reliability()){
                int tmp=index[i];
                index[i]=index[j];
                index[j]=tmp;
            }
    }
    bool V1_valid=false;
    decode_count++;
    CPolynomial V1=Berlekamp_Massey(R,&V1_valid);
    double V1_lambda=MLE(R,V1,y);
    CPolynomial N;
    if(V1_valid){
        if(V1_lambda<=ErwinChu_B(R,V1,y,index)){
            if(pass) *pass=true;
            return V1;
        }
    }
    for(Code i=1;i<pow(MASK+1,T);i++){
        for(int j=0;j<T;j++){
           N[index[j]] = (i>>(j*_m_))&MASK;
        }
        CPolynomial RpN=R+N;
        bool ispass=true;
        decode_count++;
        CPolynomial V2=Berlekamp_Massey(RpN,&ispass);
        if(ispass){
            if(V1_valid){
                double V2_lambda=MLE(R,V2,y);
                if(V2_lambda<V1_lambda){//swap
                    V1_lambda=V2_lambda;
                    CPolynomial tmp=V1; V1=V2; V2=tmp;
                }
                if(  V1_lambda <= ErwinChu3_B(R,V1,V2,y,index) ){
                    return V1;
                }
            }else{
                V1=V2;
                V1_lambda=MLE(R,V1,y);
                V1_valid=true;
                if(V1_lambda<=ErwinChu_B(R,V1,y,index)){
                    if(pass) *pass=true;
                    return V1;
                }
            }
        }
    }
    if(V1_valid){
        return V1;
    }else{
        return R;
    }
}

CPolynomial algorithm_CHU3_and_EC2(CPolynomial &R,CSignal &y,unsigned long long &decode_count,bool *pass){
    std::vector<CPolynomial> all_dc;
    int index[MASK];
    for(int i=0;i<MASK;i++) index[i]=i;
    for(int i=0;i<MASK;i++){
        for(int j=i+1;j<MASK;j++)
            if(y[index[i]]->reliability()>y[index[j]]->reliability()){
                int tmp=index[i];
                index[i]=index[j];
                index[j]=tmp;
            }
    }

    bool V1_valid=false;
    decode_count++;
    CPolynomial V1=Berlekamp_Massey(R,&V1_valid);
    double V1_lambda=MLE(R,V1,y);

    if(V1_valid){
        if(V1_lambda<=ErwinChu_B(R,V1,y,index)){
            if(pass) *pass=true;
            return V1;
        }
        all_dc.push_back(V1);
    }
    CPolynomial N;
    for(Code i=1;i<pow(MASK+1,T);i++){
        for(int j=0;j<T;j++){
           N[index[j]] = (i>>(j*_m_))&MASK;
        }
        CPolynomial RpN=R+N;
        bool ispass=false;
        CPolynomial V2;
        for(std::vector<CPolynomial>::iterator it=all_dc.begin();it!=all_dc.end();it++){
            if(it->distance(RpN)<=T){
                ispass=true;
                V2=*it;
                break;
            }
        }
        if(!ispass){
            ispass=true;
            decode_count++;
            V2=Berlekamp_Massey(RpN,&ispass);
        }
        if(ispass){
            all_dc.push_back(V2);
            if(V1_valid){
                double V2_lambda=MLE(R,V2,y);
                if(V2_lambda<V1_lambda){//swap
                    V1_lambda=V2_lambda;
                    CPolynomial tmp=V1; V1=V2; V2=tmp;
                }
                if(  V1_lambda <= ErwinChu3_B(R,V1,V2,y,index) ){
                    return V1;
                }
            }else{
                V1=V2;
                V1_lambda=MLE(R,V1,y);
                V1_valid=true;
                if(V1_lambda<=ErwinChu_B(R,V1,y,index)){
                    if(pass) *pass=true;
                    return V1;
                }
            }
        }
    }
    if(V1_valid){
        return V1;
    }else{
        return R;
    }
}

CPolynomial algorithm_CHU3_SORT(CPolynomial &R,CSignal &y,unsigned long long &decode_count,bool *pass){
    int index[MASK];
    for(int i=0;i<MASK;i++) index[i]=i;
    for(int i=0;i<MASK;i++){
        for(int j=i+1;j<MASK;j++)
            if(y[index[i]]->reliability()>y[index[j]]->reliability()){
                int tmp=index[i];
                index[i]=index[j];
                index[j]=tmp;
            }
    }
    bool V1_valid=false;
    decode_count++;
    CPolynomial V1=Berlekamp_Massey(R,&V1_valid);
    double V1_lambda=MLE(R,V1,y);
    //CPolynomial N;
    if(V1_valid){
        if(V1_lambda<=ErwinChu_B(R,V1,y,index)){
            if(pass) *pass=true;
            //std::cout<<MLE(R,V1,y)<<std::endl;
            return V1;
        }
    }
    //CPolynomial N;
    for(Code i=1;i<pow(MASK+1,T);i++){
        CPolynomial RpN=R;
        for(int j=0;j<T;j++){
           RpN[index[j]] = y[index[j]]->code[(i>>(j*_m_))&MASK];
        }
        bool ispass=true;
        decode_count++;
        CPolynomial V2=Berlekamp_Massey(RpN,&ispass);
        if(ispass){
            if(V1_valid){
                double V2_lambda=MLE(R,V2,y);
                if(V2_lambda<V1_lambda){//swap
                    V1_lambda=V2_lambda;
                    CPolynomial tmp=V1; V1=V2; V2=tmp;
                }
                if(  V1_lambda <= ErwinChu3_B(R,V1,V2,y,index) ){
                    //std::cout<<V2<<std::endl;
                    //std::cout<<"FUCK:"<<MLE(R,V1,y)<<std::endl;
                    //std::cout<<"FUCK:"<<MLE(R,V2,y)<<std::endl;
                    //std::cout<<"FUCK:"<< ErwinChu2_B(R,V1,V2,y,index)<<std::endl;
                    return V1;
                }
            }else{
                V1=V2;
                V1_lambda=MLE(R,V1,y);
                V1_valid=true;
                if(V1_lambda<=ErwinChu_B(R,V1,y,index)){
                    if(pass) *pass=true;
                    //std::cout<<MLE(R,V1,y)<<std::endl;
                    return V1;
                }
            }
        }
    }
    if(V1_valid){
        //std::cout<<MLE(R,V1,y)<<std::endl;
        return V1;
    }else{
        return R;
    }
}
CPolynomial algorithm_CHU3_and_EC2_SORT(CPolynomial &R,CSignal &y,unsigned long long &decode_count,bool *pass){
    std::vector<CPolynomial> all_dc;
    int index[MASK];
    for(int i=0;i<MASK;i++) index[i]=i;
    for(int i=0;i<MASK;i++){
        for(int j=i+1;j<MASK;j++)
            if(y[index[i]]->reliability()>y[index[j]]->reliability()){
                int tmp=index[i];
                index[i]=index[j];
                index[j]=tmp;
            }
    }
    decode_count++;
    bool V1_valid=false;
    CPolynomial V1=Berlekamp_Massey(R,&V1_valid);
    double V1_lambda=MLE(R,V1,y);

    if(V1_valid){
        if(V1_lambda<=ErwinChu_B(R,V1,y,index)){
            if(pass) *pass=true;
            return V1;
        }
        all_dc.push_back(V1);
    }
    //CPolynomial N;
    for(Code i=1;i<pow(MASK+1,T);i++){
        CPolynomial RpN=R;
        for(int j=0;j<T;j++){
           RpN[index[j]] = y[index[j]]->code[(i>>(j*_m_))&MASK];
        }
        bool ispass=false;
        CPolynomial V2;
        for(std::vector<CPolynomial>::iterator it=all_dc.begin();it!=all_dc.end();it++){
            if(it->distance(RpN)<=T){
                ispass=true;
                V2=*it;
                break;
            }
        }
        if(!ispass){
            ispass=true;
            decode_count++;
            V2=Berlekamp_Massey(RpN,&ispass);
        }
        if(ispass){
            all_dc.push_back(V2);
            if(V1_valid){
                double V2_lambda=MLE(R,V2,y);
                if(V2_lambda<V1_lambda){//swap
                    V1_lambda=V2_lambda;
                    CPolynomial tmp=V1; V1=V2; V2=tmp;
                }
                if(  V1_lambda <= ErwinChu3_B(R,V1,V2,y,index) ){
                    return V1;
                }
            }else{
                V1=V2;
                V1_lambda=MLE(R,V1,y);
                V1_valid=true;
                if(V1_lambda<=ErwinChu_B(R,V1,y,index)){
                    if(pass) *pass=true;
                    return V1;
                }
            }
        }
    }
    if(V1_valid){
        return V1;
    }else{
        return R;
    }
}

CPolynomial algorithm_SCA_CHU3(Code pattern,CPolynomial &R,CSignal &y,unsigned long long &decode_count,bool *pass){
    int index[MASK];
    for(int i=0;i<MASK;i++) index[i]=i;
    for(int i=0;i<MASK;i++){
        for(int j=i+1;j<MASK;j++)
            if(y[index[i]]->reliability()>y[index[j]]->reliability()){
                int tmp=index[i];
                index[i]=index[j];
                index[j]=tmp;
            }
    }
    bool V1_valid=false;
    decode_count++;
    CPolynomial V1=Berlekamp_Massey(R,&V1_valid);
    double V1_lambda=MLE(R,V1,y);
    CPolynomial N;
    if(V1_valid){
        if(V1_lambda<=ErwinChu_B(R,V1,y,index)){
            if(pass) *pass=true;
            return V1;
        }
    }
    for(Code i=1;i<pattern;i++){
        for(int j=0;j<T;j++){
           N[index[j]] = (i>>(j*_m_))&MASK;
        }
        CPolynomial RpN=R+N;
        bool ispass=true;
        decode_count++;
        CPolynomial V2=Berlekamp_Massey(RpN,&ispass);
        if(ispass){
            if(V1_valid){
                double V2_lambda=MLE(R,V2,y);
                if(V2_lambda<V1_lambda){//swap
                    V1_lambda=V2_lambda;
                    CPolynomial tmp=V1; V1=V2; V2=tmp;
                }
                if(  V1_lambda <= ErwinChu3_B(R,V1,V2,y,index) ){
                    return V1;
                }
            }else{
                V1=V2;
                V1_lambda=MLE(R,V1,y);
                V1_valid=true;
                if(V1_lambda<=ErwinChu_B(R,V1,y,index)){
                    if(pass) *pass=true;
                    return V1;
                }
            }
        }
    }
    if(V1_valid){
        return V1;
    }else{
        return R;
    }
}
CPolynomial algorithm_SCA_CHU3_and_EC2(Code pattern,CPolynomial &R,CSignal &y,unsigned long long &decode_count,bool *pass){
    std::vector<CPolynomial> all_dc;
    int index[MASK];
    for(int i=0;i<MASK;i++) index[i]=i;
    for(int i=0;i<MASK;i++){
        for(int j=i+1;j<MASK;j++)
            if(y[index[i]]->reliability()>y[index[j]]->reliability()){
                int tmp=index[i];
                index[i]=index[j];
                index[j]=tmp;
            }
    }

    bool V1_valid=false;
    decode_count++;
    CPolynomial V1=Berlekamp_Massey(R,&V1_valid);
    double V1_lambda=MLE(R,V1,y);

    if(V1_valid){
        if(V1_lambda<=ErwinChu_B(R,V1,y,index)){
            if(pass) *pass=true;
            return V1;
        }
        all_dc.push_back(V1);
    }
    CPolynomial N;
    for(Code i=1;i<pattern;i++){
        for(int j=0;j<T;j++){
           N[index[j]] = (i>>(j*_m_))&MASK;
        }
        CPolynomial RpN=R+N;
        bool ispass=false;
        CPolynomial V2;
        for(std::vector<CPolynomial>::iterator it=all_dc.begin();it!=all_dc.end();it++){
            if(it->distance(RpN)<=T){
                ispass=true;
                V2=*it;
                break;
            }
        }
        if(!ispass){
            ispass=true;
            decode_count++;
            V2=Berlekamp_Massey(RpN,&ispass);
        }
        if(ispass){
            all_dc.push_back(V2);
            if(V1_valid){
                double V2_lambda=MLE(R,V2,y);
                if(V2_lambda<V1_lambda){//swap
                    V1_lambda=V2_lambda;
                    CPolynomial tmp=V1; V1=V2; V2=tmp;
                }
                if(  V1_lambda <= ErwinChu3_B(R,V1,V2,y,index) ){
                    return V1;
                }
            }else{
                V1=V2;
                V1_lambda=MLE(R,V1,y);
                V1_valid=true;
                if(V1_lambda<=ErwinChu_B(R,V1,y,index)){
                    if(pass) *pass=true;
                    return V1;
                }
            }
        }
    }
    if(V1_valid){
        return V1;
    }else{
        return R;
    }
}
CPolynomial algorithm_SSCA_CHU3(Code pattern,CPolynomial &R,CSignal &y,unsigned long long &decode_count,bool *pass){
    std::vector<CPolynomial> SWAP;
    SWAP.push_back(R);
    int index[MASK];
    for(int i=0;i<MASK;i++) index[i]=i;
    for(int i=0;i<MASK;i++){
        for(int j=i+1;j<MASK;j++)
            if(y[index[i]]->reliability()>y[index[j]]->reliability()){
                int tmp=index[i];
                index[i]=index[j];
                index[j]=tmp;
            }
    }
    bool V1_valid=false;
    decode_count++;
    CPolynomial V1=Berlekamp_Massey(R,&V1_valid);
    double V1_lambda=MLE(R,V1,y);
    if(V1_valid){
        if(V1_lambda<=ErwinChu_B(R,V1,y,index)){
            if(pass) *pass=true;
            return V1;
        }
    }
    for(Code i=1;i<pattern;i++){
        CPolynomial RpN=y.get_random_data(CGLOBAL::RAND_generator);//R+N;
        if(find(SWAP.begin(),SWAP.end(),RpN) != SWAP.end()){
            continue;
        }
        SWAP.push_back(RpN);
        bool ispass=false;
        decode_count++;
        CPolynomial V2=Berlekamp_Massey(RpN,&ispass);
        if(ispass){
            if(V1_valid){
                double V2_lambda=MLE(R,V2,y);
            
                if(V2_lambda<V1_lambda){//swap
                    V1_lambda=V2_lambda;
                    CPolynomial tmp=V1; V1=V2; V2=tmp;
                }
                if(  V1_lambda <= ErwinChu3_B(R,V1,V2,y,index) ){
                    return V1;
                }
            }else{
                V1=V2;
                V1_lambda=MLE(R,V1,y);
                V1_valid=true;
                if(V1_lambda<=ErwinChu_B(R,V1,y,index)){
                    if(pass) *pass=true;
                    return V1;
                }
            }
        }
    }
    if(V1_valid){
        //std::cout<<"fuck"<<std::endl;
        return V1;
    }else{
        return R;
    }
}
CPolynomial algorithm_SSCA_CHU3_and_EC2(Code pattern,CPolynomial &R,CSignal &y,unsigned long long &decode_count,bool *pass){
    std::vector<CPolynomial> all_dc;
    std::vector<CPolynomial> SWAP;
    SWAP.push_back(R);
    int index[MASK];
    for(int i=0;i<MASK;i++) index[i]=i;
    for(int i=0;i<MASK;i++){
        for(int j=i+1;j<MASK;j++)
            if(y[index[i]]->reliability()>y[index[j]]->reliability()){
                int tmp=index[i];
                index[i]=index[j];
                index[j]=tmp;
            }
    }

    bool V1_valid=false;
    decode_count++;
    CPolynomial V1=Berlekamp_Massey(R,&V1_valid);
    double V1_lambda=MLE(R,V1,y);

    if(V1_valid){
        if(V1_lambda<=ErwinChu_B(R,V1,y,index)){
            if(pass) *pass=true;
            return V1;
        }
        all_dc.push_back(V1);
    }
    for(Code i=1;i<pattern;i++){
        CPolynomial RpN=y.get_random_data(CGLOBAL::RAND_generator);//R+N;
        if(find(SWAP.begin(),SWAP.end(),RpN) != SWAP.end()){
            continue;
        }
        SWAP.push_back(RpN);
        bool ispass=false;
        CPolynomial V2;
        for(std::vector<CPolynomial>::iterator it=all_dc.begin();it!=all_dc.end();it++){
            if(it->distance(RpN)<=T){
                ispass=true;
                V2=*it;
                break;
            }
        }
        if(!ispass){
            ispass=true;
            decode_count++;
            V2=Berlekamp_Massey(RpN,&ispass);
        }
        if(ispass){
            all_dc.push_back(V2);
            if(V1_valid){
                double V2_lambda=MLE(R,V2,y);
                if(V2_lambda<V1_lambda){//swap
                    V1_lambda=V2_lambda;
                    CPolynomial tmp=V1; V1=V2; V2=tmp;
                }
                if(  V1_lambda <= ErwinChu3_B(R,V1,V2,y,index) ){
                    return V1;
                }
            }else{
                V1=V2;
                V1_lambda=MLE(R,V1,y);
                V1_valid=true;
                if(V1_lambda<=ErwinChu_B(R,V1,y,index)){
                    if(pass) *pass=true;
                    return V1;
                }
            }
        }
    }
    if(V1_valid){
        return V1;
    }else{
        return R;
    }
}

//CPolynomial algorithm_ENCODE_CHU3(Code pattern,CPolynomial &R,CSignal &y,unsigned long long &decode_count,bool *pass=NULL);
CPolynomial algorithm_ENCODE_CHU3(Code pattern,CPolynomial &R,CSignal &y,unsigned long long &decode_count,bool *pass){
    if(pass) *pass=true;

    int index[MASK];
    int index2[MASK-T*2];
    for(int i=0;i<MASK;i++)index[i]=i;
    for(int i=0;i<MASK-T*2;i++)index2[i]=i;

    for(int i=0;i<MASK;i++){
        for(int j=i+1;j<MASK;j++)
            if(y[index[i]]->reliability()>y[index[j]]->reliability()){
                int tmp=index[i];
                index[i]=index[j];
                index[j]=tmp;
            }
    }

    for(int i=0;i<MASK-T*2;i++){
        for(int j=i+1;j<MASK-T*2;j++)
            if(y[index2[i]+T*2]->reliability()>y[index2[j]+T*2]->reliability()){
                int tmp=index2[i];
                index2[i]=index2[j];
                index2[j]=tmp;
            }
    }
    //CPolynomial RpN=y.get_random_data(CGLOBAL::RAND_generator);
    //CPolynomial M=R>>(T*2);
    decode_count++;
    CPolynomial V1 = Berlekamp_Massey(R);
    //M=M>>(T*2);
    double V1_lambda=MLE(R,V1,y);

    if(V1_lambda<=ErwinChu_B(R,V1,y,index)){
        return V1;
    }
    
    CPolynomial N;
    for(Code i=0;i<pattern;i++){
        CPolynomial MpN=y.get_random_data(CGLOBAL::RAND_generator)>>(T*2);
        decode_count++;
        CPolynomial V2 = RSencode( MpN,CGLOBAL::GXP);
 
        double V2_lambda=MLE(R,V2,y);

        if(V2_lambda<V1_lambda){//swap
            V1_lambda=V2_lambda;
            CPolynomial tmp=V1; V1=V2; V2=tmp;
        }

        if(  V1_lambda <= ErwinChu3_B(R,V1,V2,y,index) ){
            return V1;
        }
        
    }
    return V1;
    /*CPolynomial M=Berlekamp_Massey(R);
    decode_count++;
    CPolynomial V1 = M;//Berlekamp_Massey(R,&V1_valid);
    M=M>>(T*2);
    double V1_lambda=MLE(R,V1,y);

    if(V1_lambda<=ErwinChu_B(R,V1,y,index)){
        return V1;
    }
    
    CPolynomial N;
    for(Code i=0;i<pattern;i++){
        for(int j=0; j<T && T*2+j<MASK;j++){
           N[index2[j]] = (i>>(j*_m_))&MASK;
        }
        CPolynomial MpN=M+N;
        decode_count++;
        CPolynomial V2 = RSencode( MpN,CGLOBAL::GXP);
 
        double V2_lambda=MLE(R,V2,y);

        if(V2_lambda<V1_lambda){//swap
            V1_lambda=V2_lambda;
            CPolynomial tmp=V1; V1=V2; V2=tmp;
        }

        if(  V1_lambda <= ErwinChu3_B(R,V1,V2,y,index) ){
            return V1;
        }
        
    }
    return V1;*/
}