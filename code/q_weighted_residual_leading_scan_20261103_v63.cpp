#include <array>
#include <algorithm>
#include <cassert>
#include <cstdio>
#include <fstream>
#include <numeric>
#include <vector>
using namespace std;
struct Scanner {
    FILE*f;char buf[1<<20];size_t size=0,pos=0;
    int get(){if(pos==size){size=fread(buf,1,sizeof(buf),f);pos=0;if(!size)return EOF;}return buf[pos++];}
    unsigned long number(){int c;do{c=get();assert(c!=EOF);}while(c<=32);unsigned long n=0;while(c>32&&c!=EOF){n=n*10+c-'0';c=get();}return n;}
};
int main(int argc,char**argv){
    assert(argc==4);ifstream weightsfile(argv[2]);array<int,34>weights;for(auto&w:weights)weightsfile>>w;
    vector<array<int,34>>orders(4);for(auto&o:orders)iota(o.begin(),o.end(),0);
    sort(orders[0].begin(),orders[0].end(),[&](int i,int j){return weights[i]!=weights[j]?weights[i]>weights[j]:i<j;});
    sort(orders[1].begin(),orders[1].end(),[&](int i,int j){return weights[i]!=weights[j]?weights[i]>weights[j]:i>j;});
    reverse(orders[3].begin(),orders[3].end());
    Scanner input;input.f=fopen(argv[1],"r");assert(input.f);ofstream out(argv[3]);
    out<<4<<'\n';for(auto&o:orders){for(int i:o)out<<i<<' ';out<<'\n';}
    for(int rownumber=0;rownumber<350;rownumber++){
        int row=input.number();long count=input.number();array<array<unsigned long,34>,4>best{};array<unsigned long,4>coeff{};
        for(long term=0;term<count;term++){
            unsigned long c=input.number();array<unsigned long,34>exponents;for(auto&e:exponents)e=input.number();
            for(int order=0;order<4;order++){
                bool better=term==0;
                if(term)for(int j:orders[order])if(exponents[j]!=best[order][j]){better=exponents[j]>best[order][j];break;}
                if(better){best[order]=exponents;coeff[order]=c;}
            }
        }
        for(int order=0;order<4;order++){out<<order<<' '<<row<<' '<<coeff[order];for(auto e:best[order])out<<' '<<e;out<<'\n';}
    }
    fclose(input.f);
}
