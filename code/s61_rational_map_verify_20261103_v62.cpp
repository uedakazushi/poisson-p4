// Exact general-parameter S61 identities with weighted common denominators.
#include <cassert>
#include <chrono>
#include <flint/fmpq_mpoly.h>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <string>
#include <vector>
using namespace std;
fmpq_mpoly_ctx_t ctx;
int variables;
struct Poly {
    fmpq_mpoly_struct v[1];
    Poly(){fmpq_mpoly_init(v,ctx);}
    Poly(const Poly& p):Poly(){fmpq_mpoly_set(v,p.v,ctx);}
    Poly(Poly&& p):Poly(){fmpq_mpoly_swap(v,p.v,ctx);}
    Poly& operator=(Poly p){fmpq_mpoly_swap(v,p.v,ctx);return *this;}
    ~Poly(){fmpq_mpoly_clear(v,ctx);}
    bool zero()const{return fmpq_mpoly_is_zero(v,ctx);}
};
Poly multiply(const Poly& a,const Poly& b){Poly p;fmpq_mpoly_mul(p.v,a.v,b.v,ctx);return p;}
Poly read_poly(istream& in){
    long size;in>>size;assert(in&&size>=0);Poly p;
    fmpq_t q;fmpq_init(q);vector<ulong> exponents(variables);
    for(long i=0;i<size;i++){
        string value;in>>value;for(auto& e:exponents)in>>e;assert(in);
        assert(fmpq_set_str(q,value.c_str(),10)==0);fmpq_canonicalise(q);
        fmpq_mpoly_push_term_fmpq_ui(p.v,q,exponents.data(),ctx);
    }
    fmpq_mpoly_sort_terms(p.v,ctx);fmpq_mpoly_combine_like_terms(p.v,ctx);fmpq_clear(q);return p;
}
struct Term{vector<int> monomial;string scalar;};
using Equation=vector<Term>;
Equation read_equation(istream& in){
    int size;in>>size;assert(in&&size>=0);Equation result(size);
    for(auto& t:result){int n;in>>n>>t.scalar;t.monomial.resize(n);for(int& j:t.monomial)in>>j;assert(in);}
    return result;
}
int main(int argc,char**argv){
    assert(argc==2);ifstream in(argv[1]);in>>variables;assert(in&&variables>0);
    fmpq_mpoly_ctx_init(ctx,variables,ORD_DEGLEX);
    auto started=chrono::steady_clock::now();
    {
        vector<int> weights(61);for(int& w:weights)in>>w;
        Poly denominator=read_poly(in);vector<Poly> numerator(61),powers(6);
        fmpq_mpoly_one(powers[0].v,ctx);
        for(int i=1;i<4;i++)powers[i]=multiply(powers[i-1],denominator);
        for(auto& p:numerator)p=read_poly(in);
        long nodes=0,max_terms=0;
        function<Poly(const Equation&)> horner=[&](const Equation& equation)->Poly{
            ++nodes;if(equation.empty())return Poly();
            vector<int> frequency(61,0);
            for(const auto& t:equation)for(int j:t.monomial){assert(j>=0&&j<61);++frequency[j];}
            int variable=0;for(int j=1;j<61;j++)if(frequency[j]>frequency[variable])variable=j;
            if(!frequency[variable]){
                Poly p;fmpq_t q,sum;fmpq_init(q);fmpq_init(sum);
                for(const auto& t:equation){assert(fmpq_set_str(q,t.scalar.c_str(),10)==0);fmpq_canonicalise(q);fmpq_add(sum,sum,q);}
                fmpq_mpoly_set_fmpq(p.v,sum,ctx);fmpq_clear(sum);fmpq_clear(q);return p;
            }
            map<int,Equation> groups;int maximum=0;
            for(const auto& t:equation){
                Term reduced;reduced.scalar=t.scalar;int power=0;
                for(int j:t.monomial){if(j==variable)++power;else reduced.monomial.push_back(j);}
                groups[power].push_back(std::move(reduced));maximum=max(maximum,power);
            }
            Poly result;
            for(int power=maximum;power>=0;--power){
                if(power!=maximum)result=multiply(result,numerator[variable]);
                auto group=groups.find(power);
                if(group!=groups.end()){Poly coefficient=horner(group->second);fmpq_mpoly_add(result.v,result.v,coefficient.v,ctx);}
            }
            max_terms=max(max_terms,fmpq_mpoly_length(result.v,ctx));return result;
        };
        auto evaluate=[&](const Equation& equation){
            int weight=0;for(const auto& t:equation){int n=0;for(int i:t.monomial)n+=weights[i];weight=max(weight,n);}
            assert(weight<6);
            for(const auto& t:equation){
                int w=0;for(int i:t.monomial)w+=weights[i];
                assert(w==weight); // all encoded S61 rows are weighted homogeneous
            }
            return pair<Poly,int>(horner(equation),weight);
        };
        vector<Poly> reconstructed(150),reconstruction_denominator(150);
        for(int i=0;i<150;i++){
            string value;in>>value;fmpq_t q;fmpq_init(q);assert(fmpq_set_str(q,value.c_str(),10)==0);
            fmpq_canonicalise(q);fmpq_mpoly_set_fmpq(reconstructed[i].v,q,ctx);fmpq_clear(q);
            reconstruction_denominator[i]=powers[0];
        }
        for(int i=0;i<61;i++){int j;in>>j;assert(j>=0&&j<150);reconstructed[j]=numerator[i];reconstruction_denominator[j]=powers[weights[i]];}
        int eliminated;in>>eliminated;assert(eliminated==55);
        for(int i=0;i<eliminated;i++){
            int j;in>>j;auto p=evaluate(read_equation(in));
            reconstructed[j]=std::move(p.first);reconstruction_denominator[j]=powers[p.second];
        }
        int equations;in>>equations;assert(equations==134);
        for(int i=0;i<equations;i++){
            auto p=evaluate(read_equation(in));assert(p.first.zero());
            if((i+1)%20==0)cout<<"EQUATIONS "<<i+1<<" PASS\n"<<flush;
        }
        for(int i=0;i<150;i++){
            Poly expected_numerator=read_poly(in),expected_denominator=read_poly(in);
            Poly left=multiply(reconstructed[i],expected_denominator),right=multiply(expected_numerator,reconstruction_denominator[i]);
            fmpq_mpoly_sub(left.v,left.v,right.v,ctx);assert(left.zero());
        }
        assert(in);cout<<"PASS EQUATIONS 134 RECONSTRUCTION 150 HORNER_NODES "<<nodes<<" MAX_TERMS "<<max_terms
                        <<" SECONDS "<<chrono::duration<double>(chrono::steady_clock::now()-started).count()<<"\n";
    }
    fmpq_mpoly_ctx_clear(ctx);
}
