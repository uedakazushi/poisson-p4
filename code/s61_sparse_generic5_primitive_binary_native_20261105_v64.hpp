// Exact uncompressed primitive integer codec. Requires P, ctx and need().
#include <climits>
void binary_uint(string&out,ulong n){while(n>=128){out.push_back(char((n&127)|128));n>>=7;}out.push_back(char(n));}
void binary_z(string&out,const fmpz_t z){
    char*raw=fmpz_get_str(nullptr,16,z);string h(raw);flint_free(raw);bool neg=!h.empty()&&h[0]=='-';if(neg)h.erase(0,1);
    if(h=="0")h.clear();if(h.size()%2)h="0"+h;binary_uint(out,2*(h.size()/2)+neg);
    auto digit=[](char c){return c<='9'?c-'0':c-'a'+10;};
    for(size_t j=0;j<h.size();j+=2)out.push_back(char(16*digit(h[j])+digit(h[j+1])));
}
string binary_encode(const P&p){
    string out;binary_z(out,fmpq_numref(p.v->content));binary_z(out,fmpq_denref(p.v->content));
    for(long s:p.shift){ulong u=s>=0?2*ulong(s):2*ulong(-(s+1))+1;binary_uint(out,u);}binary_uint(out,p.n);
    fmpz_t z;fmpz_init(z);ulong e[5];
    for(long j=0;j<p.n;j++){fmpz_mpoly_get_term_coeff_fmpz(z,p.v->zpoly,j,ctx->zctx);binary_z(out,z);fmpq_mpoly_get_term_exp_ui(e,p.v,j,ctx);for(ulong v:e)binary_uint(out,v);}
    fmpz_clear(z);return out;
}
struct BinaryReader{
    const string&s;size_t at=0;
    ulong u(){
        ulong n=0;unsigned shift=0;size_t start=at;
        while(true){need(at<s.size()&&shift<8*sizeof(ulong),"binary varint bounds");unsigned b=(unsigned char)s[at++];ulong part=b&127;
            need(part<=(ULONG_MAX>>shift),"binary integer overflow");n|=part<<shift;if(!(b&128))break;shift+=7;}
        string canonical;binary_uint(canonical,n);need(s.compare(start,at-start,canonical)==0,"canonical varint");return n;
    }
    void z(fmpz_t out){
        ulong tag=u(),n=tag/2;need(n<=s.size()-at,"binary magnitude bounds");bool neg=tag&1;
        if(!n){need(!neg,"no negative zero");fmpz_zero(out);return;}
        need((unsigned char)s[at]!=0,"canonical magnitude");static const char*hex="0123456789abcdef";string h;h.reserve(2*n);
        for(ulong k=0;k<n;k++){unsigned b=(unsigned char)s[at++];h.push_back(hex[b>>4]);h.push_back(hex[b&15]);}
        need(fmpz_set_str(out,h.c_str(),16)==0,"binary FLINT integer");if(neg)fmpz_neg(out,out);
    }
};
P binary_decode(const string&s,long expected){
    BinaryReader x{s};P p;fmpz_t num,den,z;fmpz_init(num);fmpz_init(den);fmpz_init(z);x.z(num);x.z(den);
    need(!fmpz_is_zero(num)&&fmpz_sgn(den)>0,"binary content");fmpq_t content,q;fmpq_init(content);fmpq_init(q);fmpq_set_fmpz_frac(content,num,den);
    for(int j=0;j<5;j++){ulong v=x.u(),a=v>>1;need(a<=ulong(LONG_MAX),"binary signed exponent");p.shift[j]=(v&1)?-long(a)-1:long(a);}
    ulong n=x.u();need(n==ulong(expected)&&n>0,"binary term count");ulong e[5];
    for(ulong k=0;k<n;k++){x.z(z);need(!fmpz_is_zero(z),"binary nonzero term");for(auto&v:e)v=x.u();fmpq_mul_fmpz(q,content,z);fmpq_mpoly_set_coeff_fmpq_ui(p.v,q,e,ctx);p.sync();need(p.n==long(k+1),"binary duplicate term");}
    need(x.at==s.size(),"binary full payload consumed");
    need(binary_encode(p)==s,"native value-to-binary roundtrip");
    fmpq_clear(content);fmpq_clear(q);fmpz_clear(num);fmpz_clear(den);fmpz_clear(z);return p;
}
