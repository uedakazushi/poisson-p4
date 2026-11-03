#include <flint/nmod_mpoly.h>
#include <flint/ulong_extras.h>
#include <cassert>
#include <chrono>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <sys/resource.h>
using namespace std;

struct Node { char op; ulong constant=0; int input=-1; vector<int> refs; };

int main(int argc,char**argv) {
    assert(argc==3);
    rlimit memory{6000000000ULL,6000000000ULL};assert(setrlimit(RLIMIT_AS,&memory)==0);
    rlimit cpu{300,301};assert(setrlimit(RLIMIT_CPU,&cpu)==0);
    ifstream in(argv[1]);int N,V,R;ulong prime,a;in>>N>>V>>R>>prime>>a;
    assert(V==34&&R==350&&prime==32003);
    vector<int> weights(V),targets(R),outputs(R);for(auto&x:weights)in>>x;
    for(auto&x:targets)in>>x;for(auto&x:outputs)in>>x;
    vector<Node> nodes(N);vector<long> uses(N);
    vector<vector<int>> emit(N);
    for(int i=0;i<R;i++){emit[outputs[i]].push_back(i);uses[outputs[i]]++;}
    for(int i=0;i<N;i++) {
        auto&n=nodes[i];in>>n.op;
        if(n.op=='q')in>>n.constant;
        else if(n.op=='v')in>>n.input;
        else {int count;in>>count;n.refs.resize(count);for(auto&j:n.refs){in>>j;assert(j<i);uses[j]++;}}
    }
    assert(in.good());string prefix=argv[2];
    FILE* sparse=fopen((prefix+".sparse").c_str(),"w");
    FILE* m2=fopen((prefix+".m2").c_str(),"w");assert(sparse&&m2);
    fprintf(m2,"R=GF(%lu)[",prime);for(int j=0;j<V;j++)fprintf(m2,"%sz%d",j?",":"",j);
    fprintf(m2,",Degrees=>{");for(int j=0;j<V;j++)fprintf(m2,"%s%d",j?",":"",weights[j]);
    fprintf(m2,"},MonomialOrder=>GRevLex];\n");
    nmod_mpoly_ctx_t ctx;nmod_mpoly_ctx_init(ctx,V,ORD_DEGREVLEX,prime);
    vector<nmod_mpoly_struct*> values(N,nullptr);
    vector<string> names;vector<const char*> nameptr;
    for(int j=0;j<V;j++)names.push_back("z"+to_string(j));for(auto&n:names)nameptr.push_back(n.c_str());
    nmod_mpoly_t temp;nmod_mpoly_init(temp,ctx);
    vector<long> lengths(R,-1);vector<ulong> evaluation(R);ulong point[34];
    for(int j=0;j<V;j++)point[j]=(37*j*j+19*j+11)%prime;
    long peak=0,live=0,maxnode=0,total=0;int completed=0;
    auto start=chrono::steady_clock::now();
    auto release=[&](int j){assert(uses[j]>0);if(--uses[j]==0&&values[j]){live-=values[j]->length;nmod_mpoly_clear(values[j],ctx);delete values[j];values[j]=nullptr;}};
    for(int index=0;index<N;index++) {
        auto&n=nodes[index];auto*p=new nmod_mpoly_struct;nmod_mpoly_init(p,ctx);values[index]=p;
        if(n.op=='q')nmod_mpoly_set_ui(p,n.constant,ctx);
        else if(n.op=='v') {if(n.input==0)nmod_mpoly_set_ui(p,a,ctx);else nmod_mpoly_gen(p,n.input-1,ctx);}
        else if(n.op=='i') {assert(n.refs.size()==1&&nmod_mpoly_is_ui(values[n.refs[0]],ctx));ulong c=nmod_mpoly_get_ui(values[n.refs[0]],ctx);assert(c);nmod_mpoly_set_ui(p,n_invmod(c,prime),ctx);}
        else {
            assert(n.op=='+'||n.op=='*');if(n.op=='*')nmod_mpoly_one(p,ctx);
            for(int j:n.refs) {
                assert(values[j]);
                if(n.op=='+')nmod_mpoly_add(temp,p,values[j],ctx);
                else nmod_mpoly_mul(temp,p,values[j],ctx);
                nmod_mpoly_swap(p,temp,ctx);
                assert(nmod_mpoly_length(p,ctx)<=2000000);
            }
        }
        long len=nmod_mpoly_length(p,ctx);live+=len;maxnode=max(maxnode,len);peak=max(peak,live);
        for(int row:emit[index]) {
            lengths[row]=len;total+=len;
            evaluation[row]=nmod_mpoly_evaluate_all_ui(p,point,ctx);
            fprintf(sparse,"%d %ld\n",row,len);
            ulong exponents[34];
            for(long term=0;term<len;term++) {
                nmod_mpoly_get_term_exp_ui(exponents,p,term,ctx);long degree=0;
                for(int j=0;j<V;j++)degree+=weights[j]*exponents[j];assert(degree==targets[row]);
                fprintf(sparse,"%lu",nmod_mpoly_get_term_coeff_ui(p,term,ctx));
                for(auto exponent:exponents)fprintf(sparse," %lu",exponent);fputc('\n',sparse);
            }
            fprintf(m2,"q%d=",row);nmod_mpoly_fprint_pretty(m2,p,nameptr.data(),ctx);fprintf(m2,";\n");
            completed++;release(index);
            if(completed%25==0) {fflush(sparse);fflush(m2);cerr<<"ROWS "<<completed<<" TERMS "<<total<<" MAXNODE "<<maxnode<<" LIVETERMS "<<live<<" SECONDS "<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<endl;}
        }
        for(int j:n.refs)release(j);
        if(uses[index]==0&&values[index]){live-=values[index]->length;nmod_mpoly_clear(values[index],ctx);delete values[index];values[index]=nullptr;}
        if(index%500==0)cerr<<"NODE "<<index<<" MAXNODE "<<maxnode<<" LIVETERMS "<<live<<endl;
    }
    assert(completed==R&&live==0);
    fprintf(m2,"I=ideal{");for(int i=0;i<R;i++)fprintf(m2,"%sq%d",i?",":"",i);fprintf(m2,"};\n");
    fclose(sparse);fclose(m2);ofstream summary(prefix+".summary");
    summary<<"TERMS";for(auto len:lengths){assert(len>=0);summary<<' '<<len;}summary<<'\n';
    summary<<"EVALUATIONS";for(auto value:evaluation)summary<<' '<<value;summary<<'\n';
    rusage usage;getrusage(RUSAGE_SELF,&usage);
    summary<<"TOTAL "<<total<<" MAXNODE "<<maxnode<<" PEAKLIVETERMS "<<peak<<" MAXRSS_KB "<<usage.ru_maxrss<<" SECONDS "<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<'\n';
    nmod_mpoly_clear(temp,ctx);nmod_mpoly_ctx_clear(ctx);flint_cleanup();
    cerr<<"EXPANSION COMPLETE "<<total<<" TERMS MAXRSS_KB "<<usage.ru_maxrss<<endl;
}
