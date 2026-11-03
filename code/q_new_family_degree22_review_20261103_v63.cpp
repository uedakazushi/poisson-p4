// Independent QQ degree22 verification at the SAME21 selected old points.
// Trace gradients use forward dual matrix products, not cyclic K contractions.
#include <flint/fmpq.h>
#include <flint/fmpq_mat.h>
#include <cassert>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>
struct Mat {
    fmpq_mat_t x;
    Mat(long n=16,long m=16){fmpq_mat_init(x,n,m);}
    Mat(const Mat& a):Mat(fmpq_mat_nrows(a.x),fmpq_mat_ncols(a.x)){fmpq_mat_set(x,a.x);}
    Mat(Mat&& a):Mat(){fmpq_mat_swap(x,a.x);}
    Mat& operator=(Mat a){fmpq_mat_swap(x,a.x);return *this;}
    ~Mat(){fmpq_mat_clear(x);}
};
Mat product(const Mat& a,const Mat& b){Mat c(fmpq_mat_nrows(a.x),fmpq_mat_ncols(b.x));fmpq_mat_mul(c.x,a.x,b.x);return c;}
Mat transpose(const Mat& a){Mat c(fmpq_mat_ncols(a.x),fmpq_mat_nrows(a.x));fmpq_mat_transpose(c.x,a.x);return c;}
Mat sum(const Mat& a,const Mat& b){Mat c(a);fmpq_mat_add(c.x,a.x,b.x);return c;}
Mat input(std::istream& in){long n,m;in>>n>>m;assert(n>0&&m>0);Mat A(n,m);std::string q;
    for(long i=0;i<n;++i)for(long j=0;j<m;++j){in>>q;assert(in);assert(fmpq_set_str(fmpq_mat_entry(A.x,i,j),q.c_str(),10)==0);fmpq_canonicalise(fmpq_mat_entry(A.x,i,j));}return A;}
void trace(fmpq_t out,const Mat& a){fmpq_zero(out);for(long i=0;i<fmpq_mat_nrows(a.x);++i)fmpq_add(out,out,fmpq_mat_entry(a.x,i,i));}
void trace_product(fmpq_t out,const Mat& a,const Mat& b){fmpq_zero(out);fmpq_t term;fmpq_init(term);
    for(int i=0;i<16;++i)for(int j=0;j<16;++j){fmpq_mul(term,fmpq_mat_entry(a.x,i,j),fmpq_mat_entry(b.x,j,i));fmpq_add(out,out,term);}fmpq_clear(term);}
std::string text(const fmpq_t x){char* p=fmpq_get_str(nullptr,10,x);std::string s(p);flint_free(p);return s;}
void json_row(std::ostream& out,const Mat& a,int row,int n){out<<'[';for(int j=0;j<n;++j){if(j)out<<',';out<<'"'<<text(fmpq_mat_entry(a.x,row,j))<<'"';}out<<']';}
Mat combine(const std::vector<Mat>& rep,const Mat& coordinates){Mat A;fmpq_t term;fmpq_init(term);
    for(int k=0;k<10;++k)for(int i=0;i<16;++i)for(int j=0;j<16;++j){fmpq_mul(term,fmpq_mat_entry(coordinates.x,k,0),fmpq_mat_entry(rep[k].x,i,j));fmpq_add(fmpq_mat_entry(A.x,i,j),fmpq_mat_entry(A.x,i,j),term);}fmpq_clear(term);return A;}
int main(int argc,char** argv){
    assert(argc==3);std::ifstream in(argv[1]);std::ofstream out(argv[2]);int count,width;in>>count>>width;assert(count==21&&width==24);
    Mat omega=input(in),oi=input(in),gi=input(in);std::vector<Mat> rep,orho;
    for(int i=0;i<10;++i){rep.push_back(input(in));orho.push_back(product(omega,rep.back()));assert(fmpq_mat_equal(orho.back().x,transpose(orho.back()).x));}
    std::vector<std::string> words8(5),words18(9);for(auto& w:words8){in>>w;assert(w.size()==4);}for(auto& w:words18){in>>w;assert(w.size()==9);}
    std::vector<std::pair<int,int>> pairs(3);for(auto& p:pairs){in>>p.first>>p.second;assert(p.first>=1&&p.second<=3&&p.first<p.second);}
    std::vector<int> selected(21);for(auto& c:selected){in>>c;assert(c>=0&&c<24);}
    std::vector<Mat> dV,dU;
    for(int k=0;k<16;++k){Mat A(16,10);for(int i=0;i<16;++i)for(int j=0;j<10;++j)fmpq_set(fmpq_mat_entry(A.x,i,j),fmpq_mat_entry(orho[j].x,i,k));dV.push_back(A);dU.push_back(product(A,gi));}
    Mat values(21,24),factors(21,18),minor(21,21);std::vector<Mat> all_gradients;
    fmpq_t half,work,term;fmpq_init(half);fmpq_set_si(half,1,2);fmpq_init(work);fmpq_init(term);
    for(int point=0;point<count;++point){
        Mat b=input(in);assert(fmpq_mat_nrows(b.x)==16&&fmpq_mat_ncols(b.x)==1);Mat V(16,10);
        for(int j=0;j<10;++j){Mat col=product(orho[j],b);for(int i=0;i<16;++i)fmpq_set(fmpq_mat_entry(V.x,i,j),fmpq_mat_entry(col.x,i,0));}
        Mat mu=transpose(product(transpose(b),V));fmpq_mat_scalar_mul_fmpq(mu.x,mu.x,half);
        Mat coefficients=product(gi,mu),theta=combine(rep,coefficients),U=product(V,gi);
        Mat H=sum(theta,product(product(oi,U),transpose(V)));
        std::map<std::string,Mat> words;Mat unit;fmpq_mat_one(unit.x);words.emplace("",unit);
        for(int length=1;length<=5;++length)for(int bits=0;bits<(1<<length);++bits){std::string w;for(int j=length-1;j>=0;--j)w+=char('0'+((bits>>j)&1));words.emplace(w,product(words.at(w.substr(0,w.size()-1)),w.back()=='0'?theta:H));}
        Mat quartic=product(transpose(mu),coefficients);fmpq_set(fmpq_mat_entry(factors.x,point,0),fmpq_mat_entry(quartic.x,0,0));
        for(int i=0;i<5;++i)trace(fmpq_mat_entry(factors.x,point,1+i),words.at(words8[i]));
        for(int i=0;i<9;++i)trace_product(fmpq_mat_entry(factors.x,point,9+i),words.at(words18[i].substr(0,4)),words.at(words18[i].substr(4)));
        Mat gradients(16,3);
        for(int k=0;k<16;++k){
            Mat dmu(10,1);
            for(int j=0;j<10;++j){fmpq_set(work,fmpq_mat_entry(V.x,k,j));for(int i=0;i<16;++i){fmpq_mul(term,fmpq_mat_entry(b.x,i,0),fmpq_mat_entry(orho[j].x,i,k));fmpq_add(work,work,term);}fmpq_mul(fmpq_mat_entry(dmu.x,j,0),work,half);}
            Mat dtheta=combine(rep,product(gi,dmu));
            Mat dH=sum(dtheta,product(oi,sum(product(dU[k],transpose(V)),product(U,transpose(dV[k])))));
            for(int q=0;q<3;++q){
                const std::string& word=words8[q+1];Mat jet;
                for(int r=0;r<4;++r){const Mat& value=(word[r]=='0'?theta:H);const Mat& derivative=(word[r]=='0'?dtheta:dH);
                    jet=sum(product(jet,value),product(words.at(word.substr(0,r)),derivative));}
                trace(fmpq_mat_entry(gradients.x,k,q),jet);
            }
        }
        Mat bracket=product(product(transpose(gradients),oi),gradients);
        for(int i=0;i<3;++i)fmpq_set(fmpq_mat_entry(factors.x,point,6+i),fmpq_mat_entry(bracket.x,pairs[i].first-1,pairs[i].second-1));
        for(int i=0;i<9;++i)fmpq_mul(fmpq_mat_entry(values.x,point,i),fmpq_mat_entry(factors.x,point,0),fmpq_mat_entry(factors.x,point,9+i));
        for(int i=0;i<5;++i)for(int j=0;j<3;++j)fmpq_mul(fmpq_mat_entry(values.x,point,9+3*i+j),fmpq_mat_entry(factors.x,point,1+i),fmpq_mat_entry(factors.x,point,6+j));
        for(int i=0;i<21;++i)fmpq_set(fmpq_mat_entry(minor.x,point,i),fmpq_mat_entry(values.x,point,selected[i]));
        all_gradients.push_back(gradients);
    }
    std::string extra;assert(!(in>>extra));fmpq_mat_det(work,minor.x);assert(!fmpq_is_zero(work));
    out<<"{\"all_passed\":true,\"factors\":[";for(int i=0;i<21;++i){if(i)out<<',';json_row(out,factors,i,18);}out<<"],\"evaluation\":[";
    for(int i=0;i<21;++i){if(i)out<<',';json_row(out,values,i,24);}out<<"],\"gradients_of_Inv8_indices1_2_3\":[";
    for(int i=0;i<21;++i){if(i)out<<',';out<<'[';for(int k=0;k<16;++k){if(k)out<<',';json_row(out,all_gradients[i],k,3);}out<<']';}
    out<<"],\"exact_selected21_minor\":\""<<text(work)<<"\"}\n";
    std::cout<<"OWN_FORWARD_DUAL_QQ_ALL378_FACTORS_ALL504_PRODUCTS_21_MINOR_PASS\n";
    fmpq_clear(half);fmpq_clear(work);fmpq_clear(term);
}
