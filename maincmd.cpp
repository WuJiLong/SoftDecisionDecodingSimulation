#include <iostream>
#include <string>
#include <ctime>
#include <random>
#include <vector>
#include "define.hpp"
#include "Galois.hpp"
#include "Polynomial.hpp"
#include "function.hpp"
#include "decoding.hpp"
#include "modulation.hpp"
#include "algorithm.hpp"
#include <fstream>
#include <string>
#include <iomanip>
using namespace std;

//#define TESTNAME "16PSKSORT_EbNo"
#include "global.hpp"
//Global val inital
unsigned int CGLOBAL::seed = 0;
CPolynomial CGLOBAL::GXP = generateGX(T,1);
std::default_random_engine CGLOBAL::RAND_generator=std::default_random_engine( CGLOBAL::seed );

string a_name[30]={"HDD","EC2","STC1","STC2","STC3","STC1EC2","STC2EC2","STC3EC2","SORTSTC1","SORTSTC2","SORTSTC3","SORTSTC1EC2","SORTSTC2EC2","SORTSTC3EC2","SCAEC2","SCASTC1","SCASTC2","SCASTC3","SCASTC1EC2","SCASTC2EC2","SCASTC3EC2","SSCA","SSCAEC2","SSCASTC1","SSCASTC2","SSCASTC3","SSCASTC1EC2","SSCASTC2EC2","SSCASTC3EC2","ENcodeSTC3"};
/*
0:HDD
1~7: class1
8~13: class2 : SORT
14~20 class3 : SCA
20~28 class4 : SSCA
29 encode
*/

void showinfo(int i=0){
    cerr<<"error! :command A B C D E F"<<endl;
    if(i==0||i==1){
        cerr<<"A:調變技術"<<endl;
        cerr<<"  1:BPSK  2:16PSK  3:16QAM"<<endl;
    }
    if(i==0||i==2){
        cerr<<"B:雜訊比"<<endl;
        cerr<<"  浮點數(double)數值"<<endl;
    }
    if(i==0||i==3){
        cerr<<"C:模擬次數"<<endl;
        cerr<<"  長整數(unsigned long long)數值"<<endl;
    }
    if(i==0||i==4){
        cerr<<"D-1:模擬演算法";//<<endl;
        for(int i=0;i<30;i++){
            if(i%7==0)cerr<<endl<<"  ";
            cerr<<i<<":"<<a_name[i]<<"\t";
        }
        cerr<<endl<<"D-2:次數"<<endl;
        cerr<<"  長整數(unsigned long long)數值"<<endl;
        cerr<<"  當選擇SCA和SSCA系列時才需填寫"<<endl;
    }
    if(i==0||i==5){
        cerr<<"E:檔案註解(選填)"<<endl;
        cerr<<"  無空格字串"<<endl;
        
    }
    if(i==0||i==6){
        cerr<<"F:亂數種子(選填)"<<endl;
        cerr<<"  長整數(unsigned long long)數值"<<endl;
    }
}


int main(int n,char **argv){
    cout<< CGLOBAL::GXP <<endl;
    CPolynomial C;
    C[1]=7;C[2]=8;C[4]=4;
    CPolynomial R = RSencode(C,CGLOBAL::GXP);
    R[5]=3;R[6]=4;
    bool t=true;
    CPolynomial cc=Berlekamp_Massey(R,&t);
    cout << t <<endl;
    cout<< cc<<endl;
    CGalois a = CGalois(9);
    for(int i=0;i<15;i++)
    cout << (a^i) <<endl;
    return 0;
    /*srand( time(NULL) );
    CPolynomial GXP=generateGX(T,1);
    CPolynomial M;//=getMessage(MASK-T*2);//rand
    M[2]=4;M[1]=2;M[0]=3;
    CPolynomial C = RSencode(M,GXP);
    CPolynomial R=C;
    R[1]=7;
    R[4]=5;
    R[6]=6;
    cout<<"G(x)="<<GXP<<endl;
    cout<<"M(x)="<<M<<endl;
    cout<<"C(x)="<<C<<endl;
    cout<<"R(x)="<<R<<endl;
    bool t=true;
    CPolynomial cc=Berlekamp_Massey(R,&t);
    cout<<"C'(x)="<< cc<<endl;
    cout<<""<<t<<endl;
    return 0;*/
    //command noise count filename seed
    if(n<5){
        showinfo(0);
        return 1;
    }
    string TESTNAME="";
    string title="調變：";
    modulation type=bpsk;
    int switchtype=atoi(argv[1]);
    switch(switchtype){
        case 1:
            type=bpsk;
            title=title+"BPSK";
            TESTNAME=TESTNAME+"BPSK";
            break;
        case 2:
            type=psk16;
            title=title+"16PSK";TESTNAME=TESTNAME+"16PSK";
            break;
        case 3:
            type=qam16;
            title=title+"16QAM";TESTNAME=TESTNAME+"16QAM";
            break;
        default:
            showinfo(1);
            exit(1);
    }
    double Noise = strtod(argv[2],NULL);
    title=title+"  Eb/No:"+string(argv[2]);
    unsigned long long count = atol(argv[3]);
    int switchalg=atoi(argv[4]);
    int randnum=16;
    if(switchalg>29 || switchalg<0){
        showinfo(4);
        exit(1);
    }
    title=title+"  演算法:"+a_name[switchalg];
    TESTNAME=TESTNAME+"_"+a_name[switchalg];
    int sss=0;
    if(switchalg>=14 && switchalg<=29 ){
        if(n<6){
            showinfo(4);
            exit(1);
        }else{
            sss=1;
            title=title+string(argv[5]);
            TESTNAME=TESTNAME+string(argv[5]);
            randnum=atoi(argv[5]);
        }
    }
    unsigned int seed = 0;
    string name="";
    if(n>5+sss){
        name=string(argv[5+sss]);
        title=title+"  檔案備註:"+name;
    }
    if(n>6+sss){
        seed = atoi(argv[6+sss]);
        title=title+"  亂數種子:"+string(argv[6+sss]);
    }
    cout<<title<<endl;



    std::default_random_engine RAND_generator = std::default_random_engine( seed );
    srand( seed );
    //CPolynomial GXP=generateGX(T,1);
    int errorbit[2]={0,0};
    int errorsym[2]={0,0};
    int errorcod[2]={0,0};
    unsigned long long de_count=0;
    unsigned long long allbit = count*MASK*_m_;
    unsigned long long allsym = count*MASK;

    if(switchalg==0){
        for(unsigned long long i=0;i<count;i++){
            system("clear");
            cout<<title<<endl;
            cout<<Noise<<":"<<i<<"/"<<count<<endl;
            CPolynomial M=getMessage(MASK-T*2);//rand
            CPolynomial C = RSencode(M,CGLOBAL::GXP);
            CSignal *bpsk_c=new CSignal(C,type);
            CSignal *noise=CSignal::new_gen_noise(type,Noise,RAND_generator);//rand
            CSignal *bpsk_r = bpsk_c->new_add(noise);
            bpsk_r->computing_probability(Noise);
            CPolynomial TEST[2];
            CPolynomial R = bpsk_r->get_data();
            TEST[0] = R;
            TEST[1] = Berlekamp_Massey(R); 
            //if(TEST[1]!=C){
            //    CPolynomial SSD = algorithm_SSCA_CHU3(R,*bpsk_r,de_count,NULL);
            //    cout<<de_count<<endl;  
            //}
            for(int i=0;i<2;i++){
                errorbit[i]+=computing_BER(C,TEST[i]);
                errorsym[i]+=computing_SER(C,TEST[i]);
                errorcod[i]+=computing_CER(C,TEST[i]);
            }
            delete bpsk_c;
            delete noise;
            delete bpsk_r;
        }
    }else if(switchalg<14){
        CPolynomial (*funcPtr1)(CPolynomial&,CSignal&,unsigned long long &,bool*);
        switch(switchalg){
            case 1:funcPtr1=&algorithm_EC2;break;
            case 2:funcPtr1=&algorithm_CHU;break;
            case 3:funcPtr1=&algorithm_CHU2;break;
            case 4:funcPtr1=&algorithm_CHU3;break;
            case 5:funcPtr1=&algorithm_CHU_and_EC2;break;
            case 6:funcPtr1=&algorithm_CHU2_and_EC2;break;
            case 7:funcPtr1=&algorithm_CHU3_and_EC2;break;
            case 8:funcPtr1=&algorithm_CHU_SORT;break;
            case 9:funcPtr1=&algorithm_CHU2_SORT;break;
            case 10:funcPtr1=&algorithm_CHU3_SORT;break;
            case 11:funcPtr1=&algorithm_CHU_and_EC2_SORT;break;
            case 12:funcPtr1=&algorithm_CHU2_and_EC2_SORT;break;
            case 13:funcPtr1=&algorithm_CHU3_and_EC2_SORT;break;
        }
        for(unsigned long long i=0;i<count;i++){
            system("clear");
            cout<<title<<endl;
            cout<<Noise<<":"<<i<<"/"<<count<<endl;
            CPolynomial M=getMessage(MASK-T*2);//rand
            CPolynomial C = RSencode(M,CGLOBAL::GXP);
            CSignal *bpsk_c=new CSignal(C,type);
            CSignal *noise=CSignal::new_gen_noise(type,Noise,RAND_generator);//rand
            CSignal *bpsk_r = bpsk_c->new_add(noise);
            bpsk_r->computing_probability(Noise);
            CPolynomial TEST[2];
            CPolynomial R = bpsk_r->get_data();
            TEST[0] = R;
            TEST[1] = (*funcPtr1)(R,*bpsk_r,de_count,NULL);
            //CPolynomial HDD=Berlekamp_Massey(R);
            //if(HDD!=TEST[1] || HDD!=C){
            //    cout<<"C:"<<C<<endl;
            //    cout<<"R:"<<R<<endl;
            //    cout<<"H:"<<HDD<<endl;
            //    cout<<"S:"<<TEST[1]<<endl;
            //    exit(1);
            //}
            for(int i=0;i<2;i++){
                errorbit[i]+=computing_BER(C,TEST[i]);
                errorsym[i]+=computing_SER(C,TEST[i]);
                errorcod[i]+=computing_CER(C,TEST[i]);
            }
            delete bpsk_c;
            delete noise;
            delete bpsk_r;
        }
    }else if(switchalg<30){
        CPolynomial (*funcPtr1)(Code,CPolynomial&,CSignal&,unsigned long long &,bool*);
        switch(switchalg){
            case 14:funcPtr1=&algorithm_SCA_EC2;break;
            case 15:funcPtr1=&algorithm_SCA_CHU;break;
            case 16:funcPtr1=&algorithm_SCA_CHU2;break;
            case 17:funcPtr1=&algorithm_SCA_CHU3;break;
            case 18:funcPtr1=&algorithm_SCA_CHU_and_EC2;break;
            case 19:funcPtr1=&algorithm_SCA_CHU2_and_EC2;break;
            case 20:funcPtr1=&algorithm_SCA_CHU3_and_EC2;break;
            case 21:funcPtr1=&algorithm_SSCA;break;
            case 22:funcPtr1=&algorithm_SSCA_EC2;break;
            case 23:funcPtr1=&algorithm_SSCA_CHU;break;
            case 24:funcPtr1=&algorithm_SSCA_CHU2;break;
            case 25:funcPtr1=&algorithm_SSCA_CHU3;break;
            case 26:funcPtr1=&algorithm_SSCA_CHU_and_EC2;break;
            case 27:funcPtr1=&algorithm_SSCA_CHU2_and_EC2;break;
            case 28:funcPtr1=&algorithm_SSCA_CHU3_and_EC2;break;
            case 29:funcPtr1=&algorithm_ENCODE_CHU3;break;
        }
        for(unsigned long long i=0;i<count;i++){
            system("clear");
            cout<<title<<endl;
            cout<<Noise<<":"<<i<<"/"<<count<<endl;
            CPolynomial M=getMessage(MASK-T*2);//rand
            CPolynomial C = RSencode(M,CGLOBAL::GXP);
            CSignal *bpsk_c=new CSignal(C,type);
            CSignal *noise=CSignal::new_gen_noise(type,Noise,RAND_generator);//rand
            CSignal *bpsk_r = bpsk_c->new_add(noise);
            bpsk_r->computing_probability(Noise);
            CPolynomial TEST[2];
            CPolynomial R = bpsk_r->get_data();
            TEST[0] = R;
            CGLOBAL::RAND_generator=std::default_random_engine(i);
            unsigned long long delta=de_count;
            TEST[1] = (*funcPtr1)(randnum,R,*bpsk_r,de_count,NULL);
            //CPolynomial HDD=Berlekamp_Massey(R);
            //cout<< setprecision(20)<<(double)de_count/(i+1)<<endl;
            //if(delta==de_count){
            //    exit(1);
            //}
            /*if(HDD!=TEST[1]){
               cout<<"C:"<<C<<endl;
               cout<<"R:"<<R<<endl;
               cout<<"H:"<<HDD<<endl;
               cout<<"S:"<<TEST[1]<<endl;
               //exit(1);s
            }*/
            for(int i=0;i<2;i++){
                errorbit[i]+=computing_BER(C,TEST[i]);
                errorsym[i]+=computing_SER(C,TEST[i]);
                errorcod[i]+=computing_CER(C,TEST[i]);
            }
            delete bpsk_c;
            delete noise;
            delete bpsk_r;
        }
    }
    string filename="DATA/RS(";
    filename=filename+to_string(MASK)+","+to_string(MASK-T*2)+")"
    +TESTNAME+"_N"+to_string(Noise)+"("+name+")";

    ofstream file(filename);
    file<<"R    BER:"<<(double)errorbit[0]/allbit<<endl;
    file<<"     BER:"<<(double)errorbit[1]/allbit<<endl;
    file<<"R    SER:"<<(double)errorsym[0]/allsym<<endl;
    file<<"     SER:"<<(double)errorsym[1]/allsym<<endl;
    file<<"R    CER:"<<(double)errorcod[0]/count<<endl;
    file<<"     CER:"<<(double)errorcod[1]/count<<endl;
    file<<"     CNT:"<< setprecision(20)<<(double)de_count/count<<endl;
    file.close();
/*
int errorbit[2]={0,0};
    int errorsym[2]={0,0};
    int errorcod[2]={0,0};
    unsigned long long de_count=0;
*/

    /*double Noise = strtod(argv[1],NULL);
    unsigned long long count = atol(argv[2]);
    unsigned int seed = 0;//time(NULL);
    string name="";
    if(n>=4){
        name=string(argv[3]);
    }
    if(n==5){
        seed = atoi(argv[4]);
    }
    std::default_random_engine RAND_generator = std::default_random_engine( seed );
    srand( seed );
    CPolynomial GXP=generateGX(T,1);//b=1

    //for(double Noise=7.0;Noise<=7.0;Noise+=1.0){
    modulation type=psk16;
    //unsigned long long count=1000000;//--------------
    int errorbit[12]={0,0,0,0,0,0,0,0,0,0,0,0};
    int errorsym[12]={0,0,0,0,0,0,0,0,0,0,0,0};
    int errorcod[12]={0,0,0,0,0,0,0,0,0,0,0,0};

    unsigned long long de_count[12]={0,0,0,0,0,0,0,0,0,0,0,0};
    unsigned long long allbit = count*MASK*_m_;
    unsigned long long allsym = count*MASK;
    for(unsigned long long i=0;i<count;i++){
        system("clear");
        cout<<Noise<<":"<<i<<"/"<<count<<endl;
        CPolynomial M=getMessage(MASK-T*2);//rand
        CPolynomial C = RSencode(M,GXP);
        
        CSignal *bpsk_c=new CSignal(C,type);
        CSignal *noise=CSignal::new_gen_noise(type,Noise,RAND_generator);//rand
        
        CSignal *bpsk_r = bpsk_c->new_add(noise);
        bpsk_r->computing_probability(Noise);
        CPolynomial TEST[12];
        CPolynomial R = bpsk_r->get_data();
        TEST[0] = R;
        TEST[1] = Berlekamp_Massey(R); 

        //TEST[2]=algorithm_EC2(R,*bpsk_r,de_count[2]);
        TEST[3]=algorithm_CHU(R,*bpsk_r,de_count[3]);
        //TEST[4]=algorithm_CHU_and_EC2_SORT(R,*bpsk_r,de_count[4]);
        TEST[5]=algorithm_CHU2(R,*bpsk_r,de_count[5]);
        TEST[6]=algorithm_CHU22(R,*bpsk_r,de_count[6]);
        if(TEST[3]!=TEST[5]){
            for(int i=0;i<MASK;i++){
                cout<<(*bpsk_r)[i]->reliability()<<endl;
            }
            cout<<R<<endl;
            cout<<C<<endl;
            cout<<TEST[3]<<endl;
            cout<<TEST[5]<<endl;
            exit(77);
        }
        for(int i=0;i<7;i++){
            errorbit[i]+=computing_BER(C,TEST[i]);
            errorsym[i]+=computing_SER(C,TEST[i]);
            errorcod[i]+=computing_CER(C,TEST[i]);
        }
        delete bpsk_c;
        delete noise;
        delete bpsk_r;
    }
    string filename="DATA/RS(";
    filename=filename+to_string(MASK)+","+to_string(MASK-T*2)+")"+TESTNAME+"_N"+to_string(Noise)+"("+name+")";
    ofstream file(filename);
    string NAME[12]={
        "R        ",
        "HDD      ",
        "EC2      ",
        "CHU      ",
        "CHUEC2   ",
        "CHU2     ",
        "CHU2EC2  ",
        "SCA 32          ",
        "CHUSCA 32       ",
        "CHUEC2SCA 32    ",
        "CHU2SCA 32      ",
        "CHU2EC2SCA 32   "};
    for(int i=0;i<7;i++)
        file<<NAME[i]<<"BER:"<<(double)errorbit[i]/allbit<<endl;
    for(int i=0;i<7;i++)
        file<<NAME[i]<<"SER:"<<(double)errorsym[i]/allsym<<endl;
    for(int i=0;i<7;i++)
        file<<NAME[i]<<"CER:"<<(double)errorcod[i]/count<<endl;
    for(int i=2;i<7;i++)    
        file<<NAME[i]<<"CNT:"<<(double)de_count[i]/count<<endl;
    for(int i=0;i<7;i++)
        file<<(double)errorbit[i]/allbit<<endl;
    for(int i=0;i<7;i++)
        file<<(double)errorsym[i]/allsym<<endl;
    for(int i=0;i<7;i++)
        file<<(double)errorcod[i]/count<<endl;
    for(int i=2;i<7;i++)    
        file<<(double)de_count[i]/count<<endl;
    //}
C:24X^30+22X^29+12X^28+12X^27+25X^26+29X^25+14X^24+21X^23+16X^22+9X^21+12X^20+9X^19+27X^18+21X^17+24X^16+30X^15+5X^14+6X^13+28X^12+22X^11+1X^10+22X^9+30X^8+24X^7+14X^6+12X^5+15X^4+2X^3+18X^2+25X^1+27
R:24X^30+22X^29+12X^28+12X^27+27X^26+31X^25+10X^24+21X^23+16X^22+9X^21+12X^20+9X^19+27X^18+21X^17+24X^16+22X^15+5X^14+6X^13+28X^12+22X^11+1X^10+22X^9+30X^8+24X^7+14X^6+12X^5+15X^4+2X^3+18X^2+25X^1+27
H:24X^30+22X^29+12X^28+12X^27+27X^26+31X^25+10X^24+21X^23+16X^22+9X^21+12X^20+9X^19+27X^18+21X^17+24X^16+22X^15+5X^14+6X^13+28X^12+22X^11+1X^10+22X^9+30X^8+24X^7+14X^6+12X^5+15X^4+2X^3+18X^2+25X^1+27
S:24X^30+22X^29+12X^28+12X^27+25X^26+29X^25+14X^24+21X^23+16X^22+9X^21+12X^20+9X^19+27X^18+21X^17+24X^16+30X^15+5X^14+6X^13+28X^12+22X^11+1X^10+22X^9+30X^8+24X^7+14X^6+12X^5+15X^4+2X^3+18X^2+25X^1+27

    */
    return 0;
}
