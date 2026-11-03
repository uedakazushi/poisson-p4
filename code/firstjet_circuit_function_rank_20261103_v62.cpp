// Exact differentiation and rank of rational S61 circuits over QQ(z)[r]/(r²-3).
#include <cassert>
#include <chrono>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <gmp.h>
#include <flint/fmpq.h>
#include <flint/fmpz_mpoly_q.h>
using namespace std;
fmpz_mpoly_ctx_t ctx;
struct Q {
    fmpz_mpoly_q_struct v[1];
    Q(){fmpz_mpoly_q_init(v,ctx);fmpz_mpoly_q_zero(v,ctx);}
    Q(const Q& p):Q(){fmpz_mpoly_q_set(v,p.v,ctx);}
    Q(Q&& p):Q(){fmpz_mpoly_q_swap(v,p.v,ctx);}
    Q& operator=(Q p){fmpz_mpoly_q_swap(v,p.v,ctx);return *this;}
    ~Q(){fmpz_mpoly_q_clear(v,ctx);}
    bool zero()const{return fmpz_mpoly_q_is_zero(v,ctx);}
};
Q integer(long n){Q q;fmpz_mpoly_q_set_si(q.v,n,ctx);return q;}
Q add(const Q&a,const Q&b){Q q;fmpz_mpoly_q_add(q.v,a.v,b.v,ctx);return q;}
Q neg(const Q&a){Q q;fmpz_mpoly_q_neg(q.v,a.v,ctx);return q;}
Q sub(const Q&a,const Q&b){return add(a,neg(b));}
Q mul(const Q&a,const Q&b){if(a.zero()||b.zero())return Q();Q q;fmpz_mpoly_q_mul(q.v,a.v,b.v,ctx);return q;}
Q divq(const Q&a,const Q&b){assert(!b.zero());Q q;fmpz_mpoly_q_div(q.v,a.v,b.v,ctx);return q;}
struct E{Q a,b;bool zero()const{return a.zero()&&b.zero();}};
E scalar(long n){E e;e.a=integer(n);return e;}
E plusE(const E&a,const E&b){return E{add(a.a,b.a),add(a.b,b.b)};}
E negative(const E&a){return E{neg(a.a),neg(a.b)};}
E minusE(const E&a,const E&b){return plusE(a,negative(b));}
E times(const E&a,const E&b){
    if(a.zero()||b.zero())return E();
    if(a.b.zero()&&b.b.zero())return E{mul(a.a,b.a),Q()};
    return E{add(mul(a.a,b.a),mul(integer(3),mul(a.b,b.b))),add(mul(a.a,b.b),mul(a.b,b.a))};
}
E inverse(const E&a){
    assert(!a.zero());
    if(a.b.zero())return E{divq(integer(1),a.a),Q()};
    Q norm=sub(mul(a.a,a.a),mul(integer(3),mul(a.b,a.b)));
    return E{divq(a.a,norm),divq(neg(a.b),norm)};
}
struct Node{char op;string constant;int arg=-1;vector<int> children;};
int main(int argc,char**argv){
    assert(argc==2);ifstream in(argv[1]);int variables,inputs,columns,count,output_count;
    in>>variables>>inputs>>columns>>count>>output_count;assert(in&&variables>0&&output_count==16);
    fmpz_mpoly_ctx_init(ctx,variables,ORD_DEGLEX);auto started=chrono::steady_clock::now();
    {
        vector<E> source(inputs);vector<int> column(inputs);
        for(int i=0;i<inputs;i++){
            char kind;int index;in>>kind>>index>>column[i];assert(in);
            if(kind=='v'){assert(index>=0&&index<variables);fmpz_mpoly_q_gen(source[i].a.v,index,ctx);}
            else if(kind=='c')source[i]=scalar(index);
            else{assert(kind=='r');source[i].b=integer(1);}
        }
        vector<Node> nodes(count);
        for(auto& n:nodes){
            in>>n.op;assert(in);
            if(n.op=='q')in>>n.constant;
            else if(n.op=='v'||n.op=='i')in>>n.arg;
            else{assert(n.op=='+'||n.op=='*');int size;in>>size;n.children.resize(size);for(int& j:n.children)in>>j;}
        }
        vector<int> outputs(output_count);for(int& j:outputs)in>>j;assert(in);
        vector<bool> needed(count,false);vector<int> pending=outputs;
        while(!pending.empty()){
            int j=pending.back();pending.pop_back();assert(j>=0&&j<count);if(needed[j])continue;needed[j]=true;
            if(nodes[j].op=='i')pending.push_back(nodes[j].arg);
            else for(int k:nodes[j].children)pending.push_back(k);
        }
        vector<E> values(count);vector<vector<E>> derivative(count);
        for(int i=0;i<count;i++)if(needed[i]){
            auto& n=nodes[i];auto& gradient=derivative[i];gradient.resize(columns);
            if(n.op=='q'){
                fmpq_t q;fmpq_init(q);assert(fmpq_set_str(q,n.constant.c_str(),10)==0);fmpq_canonicalise(q);
                fmpz_mpoly_q_set_fmpq(values[i].a.v,q,ctx);fmpq_clear(q);
            }else if(n.op=='v'){
                values[i]=source[n.arg];if(column[n.arg]>=0)gradient[column[n.arg]]=scalar(1);
            }else if(n.op=='+'){
                for(int j:n.children){values[i]=plusE(values[i],values[j]);
                    for(int c=0;c<columns;c++)if(!derivative[j][c].zero())gradient[c]=plusE(gradient[c],derivative[j][c]);}
            }else if(n.op=='*'){
                vector<E> prefix(n.children.size()+1);prefix[0]=scalar(1);
                for(int j=0;j<(int)n.children.size();j++)prefix[j+1]=times(prefix[j],values[n.children[j]]);
                values[i]=prefix.back();E suffix=scalar(1);
                for(int j=(int)n.children.size()-1;j>=0;j--){
                    int child=n.children[j];E factor=times(prefix[j],suffix);
                    if(!factor.zero())for(int c=0;c<columns;c++)if(!derivative[child][c].zero())
                        gradient[c]=plusE(gradient[c],times(factor,derivative[child][c]));
                    suffix=times(suffix,values[child]);
                }
            }else{
                assert(n.op=='i');values[i]=inverse(values[n.arg]);E factor=negative(times(values[i],values[i]));
                for(int c=0;c<columns;c++)if(!derivative[n.arg][c].zero())gradient[c]=times(factor,derivative[n.arg][c]);
            }
        }
        cout<<"EXACT_DIFFERENTIAL_READY SECONDS "<<chrono::duration<double>(chrono::steady_clock::now()-started).count()<<endl;
        vector<vector<E>> matrix;for(int j:outputs)matrix.push_back(derivative[j]);int rank=0;
        for(int c=0;c<columns&&rank<16;c++){
            int pivot=rank;while(pivot<16&&matrix[pivot][c].zero())pivot++;if(pivot==16)continue;
            if(pivot!=rank)swap(matrix[pivot],matrix[rank]);
            for(int row=rank+1;row<16;row++)if(!matrix[row][c].zero()){
                E factor=times(matrix[row][c],inverse(matrix[rank][c]));
                for(int j=c+1;j<columns;j++)matrix[row][j]=minusE(matrix[row][j],times(factor,matrix[rank][j]));
                matrix[row][c]=E();
            }
            rank++;cout<<"PIVOT "<<c<<" RANK "<<rank<<endl;
        }
        for(int row=rank;row<16;row++)for(int c=0;c<columns;c++)assert(matrix[row][c].zero());
        cout<<"EXACT_CIRCUIT_FUNCTION_FIELD_RANK "<<rank<<endl;
        cout<<"SECONDS "<<chrono::duration<double>(chrono::steady_clock::now()-started).count()<<endl;
    }
    fmpz_mpoly_ctx_clear(ctx);flint_cleanup();
}
