/*
 * Author: Miao Shenhao
 * Date: 9.19
 * Motto: Face to the weakness, expect for the strength.
*/
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
const int MAXN = 100 + 1;//维数
const int K = 10;//最大环境承载力
const int N = 200;//总迭代次数
const int P = 10;//初始种群数
const db C = 2.0;//两须距离与步长比
const db STEP = 100;//初始步长
const db EPS = 1e-6;
const db alpha = 0.95;//步长因子
const db beta = 1.05;//繁殖因子
struct Beetle{
    db stat[MAXN];
    db delta, step;
    int age;
};
Beetle deque[K+10];
int frnt, bck, n;
int AOI;
db sol[MAXN], fminv, sl[MAXN], sr[MAXN], birthnum;
db fleft, fright;
db temp[MAXN];
int sign(db x){
    if(x < 0.0) return -1;
    else return 1;
}
int rdm(int a, int b){
    return rand() % (b-a+1) + a;
}
db F(db x[]){
    db res = 0.0;
    for(int i = 1; i <= n; i++){
        res += i*x[i]*x[i];
        //cerr << x[i] << " ";
    }
    //cerr << endl;
    return res;
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
            deque[bck].stat[j] = (db)rdm(-100, 100)/200;
                //cerr << deque[bck].stat[j] << " ";
        }
            //cerr << endl;
    }
    //cerr << "I'm in!" << endl;
    for(int i = 1; i <= N; i++){
        int cnt = bck - frnt + 1;
        db rate = birthnum / cnt;
        for(int j = frnt; j <= bck; j++){
            int p = j % (K+10);
            //cerr << "$" << endl;
            //*
            if(deque[p].step <= EPS){
                deque[p].age++;//走不动了，饿死了
                if(rate > 0.000001){
                    if(rdm(1, 1000) <= round(rate * 1000)){
                        int p1 = ++bck; p1 %= K+10;
                        deque[p1].delta = alpha;
                        deque[p1].step = STEP;
                        deque[p1].age = 0;
                        memcpy(deque[p1].stat, deque[p].stat, sizeof(sol));
                    }
                }else{
                    if(rdm(1, 1000) <= round(-rate * 1000))
                        frnt++;
                }
                frnt++;
                continue;
            }
            //*/
            AOI++;
            cnt = bck - frnt + 1;
            if(cnt <= K / 2){
                birthnum = abs(birthnum) * beta;
            }else if(cnt > K/2 && cnt <= K){
                birthnum = abs(birthnum) / beta;
            }else{
                birthnum = abs(birthnum) * (-beta);
            }
            rate = birthnum / cnt;
            db b[MAXN];
            RAND(b);
            deque[p].step = alpha * deque[p].step;
            db step = deque[p].step;
            //cerr << step << endl;
            db d = step/C;
            //cerr << d << endl;
            bool BeyondDim = false;
            RAND(b);
            //cerr << b[1] << endl;
            //cerr << "$$$" << endl;
            for(int i1 = 1; i1 <= n; i1++){
            //    cerr << sl[i1] << " " ;
                //cerr << b[i1] << " ";
                sl[i1] = deque[p].stat[i1] + d * b[i1] / 2.0;
                sr[i1] = deque[p].stat[i1] - d * b[i1] / 2.0;
                if(abs(sr[i1]) > 100.0 || abs(sl[i1]) > 100.0){//修改函数定义域
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
                if(abs(temp[i1]) > 100.0){//修改函数定义域
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
            //*
            if(rate > 0.000001){
                if(rdm(1, 1000) <= round(rate * 1000)){
                    int p1 = ++bck; p1 %= K+10;
                    deque[p1].delta = alpha;
                    deque[p1].step = STEP;
                    deque[p1].age = 0;
                    memcpy(deque[p1].stat, deque[p].stat, sizeof(sol));
                }
            }else{
                if(rdm(1, 1000) <= round(-rate * 1000))
                    frnt++;
            }
            // */
        }
    }
}
int main(){
    srand(time(0));
    n = read();
    db AOS = 0;
    db kminv = 2147483647.0;
    db kmaxv = 0.0;
    fminv = 2147483647.0;
    for(int i = 1; i <= 50; i++){
        fminv = 2147483647.0;
        BAS();
        AOS += fminv;
        //cout << "$$: " << fminv << endl;
        kminv = min(kminv, fminv);
        kmaxv = max(kmaxv, fminv);
    }
    AOI /= 50;
    AOS /= 50;
    //cerr << n << endl;
    //for(int ww = 1; ww <= 99; ww++)
    //    BAS();
    cout << "ASOI: " << AOI << endl;
    cout << "AOS: " << AOS << endl;
    cout << "MOS: " << kminv << endl;
    cout << kmaxv << endl;
    //for(int i = 1; i <= n; i++)
    //    cout << sol[i] << " ";
    return 0;
}
