/*
   This code can be compiled and run ok.
     防疫模擬問題

   usage (how to run):
     ./hw6

   input file:
     data_4_6_1.txt

   output file:
     none

   compile (how to compile):
     g++ -o hw6 hw6.cpp

   pseudocode:
   主要分成有投藥還是沒投藥，每一格都會受成長、擴散、原生影響，詳見老師說明檔


   coded by 江冠駒, ID: H24031354, email: kevin040208@gmail.com
   date: 2018.06.29
*/


#include<iostream>
#include<iomanip>
#include<fstream>
#include<cstdlib>
#include<cstdio>
#include<ctime>
#include<vector>
#include<math.h>

using namespace std;

struct Zone{      //此為各網格之 struct
    int a;     //儲存該網格之縱軸座標，即由上而下第 a 列, a=0,…,M
    int b;     //儲存該網格之橫軸座標，即由左而右第 b 行, b=0,…,N
    int id;    //儲存其編號，id=0,…,MN-1
    float r;  //儲存該網格本期之殘存率，初始化為 1，
    //若該網格有投藥，則將有 X 比例被撲殺，亦即殘留了 1-X 比例。
    //倘若另有一網格 k’其所投放的藥效範圍會波及本網格，則會將網格
    //k’對本網格的殘餘量（每多一單位則殺傷效果減半）拿來乘以
    //本網格原先的殘餘量，以此類推。
    bool isM;    //儲存某網格於本期期初是否被投藥，是則為 1，否則為 0
    Zone() : r(1.) {} //將存殘存率預設為1
};
struct Virus{
    int F;   //原生病毒出現之強度(或可理解為數量)
    int p;     //原生病毒該強度之發生機率(%)
};


void readFile(char *filename,int &M,int &N,int &T,int &s,int &q, int &n_Vir,
              Zone *&V,Virus *&Vir,float *&Q0,float *&Q1,float *&Q2,float *&Q3);
void print_colored_2Darray(float *Q0,int t,int M,int N);
void print_mono_2Darray(float *Q1,float *Q2,float *Q3,int t,int M,int N);
void cal_Q1(float *&Q1,float *Q0,int s,int M,int N,Zone *V);
void cal_Q2(float *&Q2,float *Q0,int M,int N,int q,Zone *V);
void cal_Q3(float *&Q3,int M,int N,int n_Vir,Virus *Vir);
void cal_r(int *row,int *col,int l,int *X,int *R,int M,int N,Zone *V);

int main(){
    //---begin--- PART 1: Declaration and read file-----------
    int M,      //區域縱軸座標範圍(1,2,…,M)
        N,      //區域橫軸座標範圍(1,2,…,N)
        n_Vir,  //原生病毒出現之種類個數，譬如圖二(b)之 n_Vir=3 代表 3 種
            //病毒強度，分別為 60,40,0 機率 10%,20%,70%
        T,      //模擬之總期數
        s,      //給定的病毒每期成長比率(%)
        q,      //給定的病毒每期擴散比率(%)
        i,j,k,l;    //for loop 用
    Zone *V;    //動態宣告 MN 長度之網格 struct 陣列;V[k],k=0,1,…,MN-1
    Virus *Vir; //動態宣告 n_Vir 長度之原生病毒陣列;
    float   *Q0,//Q0[k],k=0,…,MN-1 存網格 k 於本期期末總病毒數量 Q0=(Q1+Q2+Q3)*r
            *Q1,//Q1[k]存 k 之上期 Q0 於本期成長所衍生的總病毒數量 Q1=Q0*(1+s)
            *Q2,//Q2[k]存 k 之相鄰網格上期 Q0 於本期擴散移入之病毒總量 Q2+=鄰格 Q0*q
            *Q3;//Q3[k]儲存網格 k 經過本期而新隨機出現的原生病毒數量，為某個強度 F
    int n_M,  //記錄當期投藥次數，若為 0 則代表不投藥
        *row,  //記錄當下投藥之網格縱軸座標
        *col,  //記錄當下投藥之網格橫軸座標
        *X,    //記錄當下投藥之撲殺率(%)
        *R;    //記錄當下投藥之影響範圍，以曼哈頓距離表示，若為 0 則代表僅影響那一格
    char filename[50];//記錄輸入之測試檔名
    
    cout << "Enter filename: ";
    cin >> filename;
    readFile(filename,M,N,T,s,q,n_Vir,V,Vir,Q0,Q1,Q2,Q3);
    cout <<"M="<<M<<",N="<<N<<", T="<<T<<", s="<<s<<"%, q="<<q<<"%, n_Vir="<<n_Vir<<"; F={";
    for(i=0;i<n_Vir;i++){
        if(i==(n_Vir-1)){cout<<Vir[i].F<<"}";}
        else{cout << Vir[i].F << ",";}
    }
    cout<<", p={";
    for(i=0;i<n_Vir;i++){
        if(i==(n_Vir-1)){cout<<Vir[i].p<<"};";}
        else{cout << Vir[i].p << ",";}
    }
    cout << endl;
    print_colored_2Darray(Q0,0,M,N);
    //--- end --- PART 1: Declaration and read file-------------


    //---begin--- PART 2,3: virus grow,spread,generate and kill--------------------
    for(i=0;i<(M*N);i++){   //初始化網格參數
            V[i].b = i%N;
            V[i].a = (i-V[i].b)/N;
            V[i].id = i;
    }
    for(i=1;i<=T;i++){
        for(j=0;j<(M*N);j++){   //預設每隔都沒投藥
            V[j].isM = 0;
        }
        cout << "Enter number of medicine: ";
        cin >> n_M;
        if(n_M){
            row = new int[n_M];
            col = new int[n_M];
            X = new int[n_M];
            R = new int[n_M];
            cout << "Enter medicine location(row,column), kill rate(%), influence area: ";
            for(l=0;l<n_M;l++){
                cin >> row[l] >> col[l] >> X[l] >> R[l];
                for(j=0;j<(M*N);j++){
                    if(V[j].a==row[l] && V[j].b==col[l]){
                        V[j].isM = 1;
                    }
                }
                cal_r(row,col,l,X,R,M,N,V);
            }
        }
        cal_Q1(Q1,Q0,s,M,N,V);
        cal_Q2(Q2,Q0,M,N,q,V);
        cal_Q3(Q3,M,N,n_Vir,Vir);
        for(j=0;j<(M*N);j++){
            Q0[j] = (Q1[j]+Q2[j]+Q3[j])*V[j].r;
        }
        print_mono_2Darray(Q1,Q2,Q3,i,M,N);
        print_colored_2Darray(Q0,i,M,N);
        cout << endl;
    }
    //---end--- PART 2,3: virus grow,spread,generate and kill----------------------
    delete [] V;
    delete [] Vir;
    delete [] Q0;
    delete [] Q1;
    delete [] Q2;
    delete [] Q3;
    delete [] row;
    delete [] col;
    delete [] X;
    delete [] R;
    return 0;
}



void readFile(char *filename,int &M,int &N,int &T,int &s,int &q, int &n_Vir
              ,Zone *&V,Virus *&Vir,float *&Q0,float *&Q1,float *&Q2,float *&Q3){
    int i;
    fstream infile;
    infile.open(filename,ios::in);
    if(!infile){
        cout << "Error in opening the file\n";
    }
    else{
        infile >> M >> N >> T >> s >> q >> n_Vir ;
        Vir = new Virus [n_Vir];
        V = new Zone [M*N];
        for(i=0;i<n_Vir;i++){
            infile >> Vir[i].F >>Vir[i].p;
        }
        Q0 = new float [M*N];
        for(i=0;i<(M*N);i++){
            infile >> Q0[i];
        }
    }
}

void print_colored_2Darray(float *Q0,int t,int M,int N){
    int i,j;
    string color[4]={"\x1b[;32;1m","\x1b[;33;1m","\x1b[;31;1m","\x1b[;35;1m"}; //存文字顏色
    string reset = "\x1b[0m"; //恢復成系統顏色
    printf("Q0[k,%d]\n",t);
    printf("%6s","V");
    for(i=0;i<N;i++){printf("%6d",i);}
    cout << endl;
    for(i=0;i<M;i++){
        printf("%6d",i);
        for(j=0;j<N;j++){
            if(Q0[i*N+j] >= 0 && Q0[i*N+j] < 33){
                printf("%s%6.1f%s",color[0].c_str(),Q0[i*N+j],reset.c_str());
            }
            else if(Q0[i*N+j] >= 33 && Q0[i*N+j] < 66){
                printf("%s%6.1f%s",color[1].c_str(),Q0[i*N+j],reset.c_str());
            }
            else if(Q0[i*N+j] >= 66 && Q0[i*N+j] < 100){
                printf("%s%6.1f%s",color[2].c_str(),Q0[i*N+j],reset.c_str());
            }
            else{
                printf("%s%6.1f%s",color[3].c_str(),Q0[i*N+j],reset.c_str());
            }
            fflush(stdout);
        }
        cout << endl;
    }
}

void print_mono_2Darray(float *Q1,float *Q2,float *Q3,int t,int M,int N){
    int i,j;
    printf("Q1[k,%d]\n",t);
    printf("%6s","V");
    for(i=0;i<N;i++){printf("%6d",i);}
    cout << endl;
    for(i=0;i<M;i++){
        printf("%6d",i);
        for(j=0;j<N;j++){
            printf("%6.1f",Q1[i*N+j]);
        }
        cout << endl;
    }
    cout << endl;
    printf("Q2[k,%d]\n",t);
    printf("%6s","V");
    for(i=0;i<N;i++){printf("%6d",i);}
    cout << endl;
    for(i=0;i<M;i++){
        printf("%6d",i);
        for(j=0;j<N;j++){
            printf("%6.1f",Q2[i*N+j]);
        }
        cout << endl;
    }
    cout << endl;
    printf("Q3[k,%d]\n",t);
    printf("%6s","V");
    for(i=0;i<N;i++){printf("%6d",i);}
    cout << endl;
    for(i=0;i<M;i++){
        printf("%6d",i);
        for(j=0;j<N;j++){
            printf("%6.1f",Q3[i*N+j]);
        }
        cout << endl;
    }
    cout << endl;
}

void cal_Q1(float *&Q1,float *Q0,int s,int M,int N,Zone *V){
    int i;
    float srate;
    srate = float(s)/100;
    Q1 = new float[M*N];
    for(i=0;i<(M*N);i++){
        Q1[i] = 0;
    }
    for(i=0;i<(M*N);i++){
        if(V[i].isM){Q1[i] = Q0[i];}
        else{Q1[i] = Q0[i]*(1+srate);}
    }
}

void cal_Q2(float *&Q2,float *Q0,int M,int N,int q,Zone *V){
    int i,j;
    float qrate;
    qrate = float(q)/100;
    Q2 = new float[M*N];
    for(i=0;i<(M*N);i++){
        Q2[i] = 0;
    }
    for(i=0;i<(M*N);i++){
        for(j=0;j<(M*N);j++){
            if((abs(V[i].a-V[j].a)+abs(V[i].b-V[j].b))==1){
                if(V[j].isM){Q2[i] = Q2[i];}
                else{Q2[i] = Q2[i] + Q0[j]*qrate;}
            }
        }
    }
}

void cal_Q3(float *&Q3,int M,int N,int n_Vir,Virus *Vir){
    int i,j,
        rn,     //random number
        lb,
        ub;
    Q3 = new float[M*N];
    for(i=0;i<(M*N);i++){
        Q3[i] = 0;
    }
    for(i=0;i<(M*N);i++){
        rn = rand()%100;
        lb = 0;
        ub = Vir[0].p;
        for(j=0;j<n_Vir;j++){
            if(rn>=lb && rn<ub){Q3[i] = Vir[j].F;}
            lb = lb + Vir[j].p;
            ub = ub + Vir[j+1].p;
        }
    }
}

void cal_r(int *row,int *col,int l,int *X,int *R,int M,int N,Zone *V){
    int j,k;
    for(j=0;j<=R[l];j++){
        for(k=0;k<(M*N);k++){
            if((abs(V[k].a-row[l])+abs(V[k].b-col[l]))==j){
                V[k].r = V[k].r*(1-((float(X[l])/100)/pow(2,j)));
            }
        }
    }
}



