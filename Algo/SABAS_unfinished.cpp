/*
 * Author: Miao Shenhao
 * Date: 9.19
 * Motto: Face to the weakness, expect for the strength.
*/
#pragma GCC optimize(3, "Ofast", "inline")
#include<cstdio>
#include<cstring>
#include<algorithm>
#include<iostream>
#include<cstdlib>
#include<map>
#include<ctime>
#include<utility>
#include<functional>
#include<cmath>
#include<vector>
#include<assert.h>
using namespace std;
#define reg register
#define endfile fclose(stdin);fclose(stdout);
typedef long long ll;
typedef unsigned long long ull;
typedef double db;
typedef std::pair<int,int> pii;
typedef std::pair<ll,ll> pll;
namespace IO{
    char buf[1<<15],*S,*T;
    inline char gc(){
        if (S==T){
            T=(S=buf)+fread(buf,1,1<<15,stdin);
            if (S==T)return EOF;
        }
        return *S++;
    }
    inline int read(){
        reg int x;reg bool f;reg char c;
        for(f=0;(c=gc())<'0'||c>'9';f=c=='-');
        for(x=c^'0';(c=gc())>='0'&&c<='9';x=(x<<3)+(x<<1)+(c^'0'));
        return f?-x:x;
    }
    inline ll readll(){
        reg ll x;reg bool f;reg char c;
        for(f=0;(c=gc())<'0'||c>'9';f=c=='-');
        for(x=c^'0';(c=gc())>='0'&&c<='9';x=(x<<3)+(x<<1)+(c^'0'));
        return f?-x:x;
    }
}
using namespace IO;
const long long llINF = 9223372036854775807;
const int INF = 2147483647;
const int MAXN = 2 + 1;//维数
const int K = 10;//最大环境承载力
const int N = 88;//总迭代次数
const int P = 1;//初始种群数
const db C = 2.0;//两须距离与步长比
const db EPS = 1e-9;
const int maxn = 2e3 + 10;
const db alpha = 0.88;//步长因子
const db beta = 1.05;//繁殖因子
db STEP = 14140;
struct Beetle{
    db stat[MAXN];
    db delta, step;
    int age;
};
Beetle deque[K+10];
int frnt, bck, n, nn, m;
int AOI;
db sol[MAXN], fminv, sl[MAXN], sr[MAXN], birthnum;
db fleft, fright, Tempr;
db X[maxn], Y[maxn], len[maxn];
db temp[MAXN];
int sign(db x){
    if(x < 0.0) return -1;
    else return 1;
}
int rdm(int a, int b){
    return rand() % (b-a+1) + a;
}
db F(db x[]){
    //cerr << "$$$" << endl;
    for(int i = 1; i <= nn; i++){
        len[i] = (X[i] - x[1]) * (X[i] - x[1]) + (Y[i] - x[2]) * (Y[i] - x[2]);
        //cerr << len[i] << " ";
    }
    sort(len+1, len+1+nn);
    //cerr << endl;
    //cerr << sqrt(len[m]) << endl;
    return sqrt(len[m]);
}
void RAND(db A[]){
    db len = 0;
    for(int i = 1; i <= n; i++)
        *(A+i) = rdm(-100, 100), len += abs(*(A+i));
    for(int i = 1; i <= n; i++)
        *(A+i) /= len;
}
void BAS(){
    //memset(deque, 0, sizeof(deque));
    //memset(sl, 0, sizeof(sl)), memset(sr, 0, sizeof(sr));
    //memset(temp, 0, sizeof(temp));
    fleft = fright = 0;
    frnt = 1, bck = 0, birthnum = 1.0;
    for(int i = 1; i <= P; i++){
        deque[++bck].delta = alpha, deque[bck].step = STEP;
        deque[bck].age = 1;
        for(int j = 1; j <= n; j++){
            deque[bck].stat[j] = (db)rdm(-100, 100)/200 + sol[j];
                //cerr << deque[bck].stat[j] << " ";
        }
            //cerr << endl;
    }
    //cerr << "I'm in!" << endl;
    for(int i = 1; i <= N; i++){
        for(int j = frnt; j <= bck; j++){
            int p = j % (K+10);
            if(deque[p].step <= EPS) continue;
            db b[MAXN];
            RAND(b);
            deque[p].step = deque[p].delta * deque[p].step;
            db step = deque[p].step;
            db d = step/C;
            bool BeyondDim = false;
            RAND(b);
            for(int i1 = 1; i1 <= n; i1++){
                sl[i1] = deque[p].stat[i1] + d * b[i1] / 2.0;
                sr[i1] = deque[p].stat[i1] - d * b[i1] / 2.0;
                if(abs(sr[i1]) > 10000.0 || abs(sl[i1]) > 10000.0){//修改函数定义域
                    BeyondDim = true;
                    break;
                }
                //cerr << sl[i1] << " ";
            }
            //cerr << endl;
            if(BeyondDim) continue;
            fleft = F(sl);
            fright = F(sr);
            //cerr << fleft << " " << fright << endl;
            for(int i1 = 1; i1 <= n; i1++){
                temp[i1] = deque[p].stat[i1] - step * b[i1] * sign(fleft - fright);
                if(abs(temp[i1]) > 10000.0){//修改函数定义域
                    BeyondDim = true;
                    break;
                }
            }
            //cerr << endl;
            if(BeyondDim) continue;
            //db key = F(temp);
            //cerr << "$: " << key << endl;
            for(int i1 = 1; i1 <= n; i1++){
                deque[p].stat[i1] = temp[i1];
                //cerr << deque[p].stat[i1] << " ";
            }
            //cerr << endl;
            fleft = min(fleft, fright);
            if(fleft < fminv){
                fminv = fleft, memcpy(sol, deque[p].stat, sizeof(sol));
                deque[p].delta = alpha;//自适应因子
            }else{
                deque[p].delta = alpha - 0.2*((1.0+deque[p].age)/(5.0*N) + 0.5);//自适应因子--寻优状况不佳
            }
            //for(int )
            deque[p].age++;
        }
    }
}
double Time(){return (double)clock()/CLOCKS_PER_SEC;}
int main(){
    //freopen("fagttering17.in", "r", stdin);
    srand(123456);
    n = 2;
    scanf("%d%d", &nn, &m);
    for(int i = 1; i <= nn; i++){
        scanf("%lf%lf", &X[i], &Y[i]);
        //cerr << X[i] << "  " << Y[i] << endl;
    }
    fminv = 2147483647;
    int L = 1;
    Tempr = 200;
    for(L = 1; L <= 100; L++){
        
    }
    if(nn <= 50){
        while(Time() < 0.9){
            BAS();
            STEP *= 0.92;
        }
        printf("%f", fminv);
        return 0;
    }
    while(Time() < 2.9){
        BAS();
        STEP *= 0.94;
    }
    printf("%.8f", fminv);
    //for(int i = 1; i <= n; i++)
    //    cout << sol[i] << " ";
    return 0;
}
