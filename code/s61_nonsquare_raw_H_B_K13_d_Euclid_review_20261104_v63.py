#!/usr/bin/env python3
"""Independent fixed quotient-ring identities; no d-Euclid/inverse/GCD search."""
import hashlib,json,math,resource,signal,sys,time
from fractions import Fraction as Q
from pathlib import Path
sys.set_int_max_str_digits(0)
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_RAW_H_B_K13_D_EUCLID_20261104_v63_REVIEW'
def load(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def ip(p):return {int(i):int(v)for i,v in p if int(v)}
def plus(p,q,c=1):
    r=p.copy()
    for i,v in q.items():
        z=r.get(i,0)+c*v
        if z:r[i]=z
        else:r.pop(i,None)
    return r
def times(p,q):
    r={}
    for i,c in p.items():
        for j,v in q.items():r[i+j]=r.get(i+j,0)+c*v
    return {i:v for i,v in r.items()if v}
def exactdiv(p,g):
    r=p.copy();q={};n=max(g);lc=g[n]
    while r and max(r)>=n:
        j=max(r)-n;c=r[max(r)];assert c%lc==0;c//=lc;q[j]=c
        r=plus(r,{i+j:v for i,v in g.items()},-c)
    assert not r and times(g,q)==p
    return q

class OwnAlgebra:
    """Dense13 numerator plus a rational denominator, with fixed monomial table.
    The table a^0..a^24 is certified once by exact integer division. Arithmetic
    does not choose coefficient pivots or compute inverses.
    """
    def __init__(self,g):
        self.g=g;self.L=g[13];self.C=self.L**12;self.ops=0;self.red=[]
        for k in range(25):
            p={k:1};den=1
            while p and max(p)>=13:
                j=max(p)-13;c=p[max(p)];p=plus({i:v*self.L for i,v in p.items()},{i+j:c*v for i,v in g.items()},-1);den*=self.L
            assert self.C%den==0;r={i:v*(self.C//den)for i,v in p.items()}
            exactdiv(plus({k:self.C},r,-1),g);self.red.append([r.get(i,0)for i in range(13)])
        self.z=self.make({},1);self.o=self.make({0:1},1)
    def make(self,p,d=1):
        assert d and max(p,default=-1)<13
        if not p:return ((0,)*13,1)
        if d<0:p={i:-v for i,v in p.items()};d=-d
        c=math.gcd(d,*(abs(v)for v in p.values()));return(tuple(p.get(i,0)//c for i in range(13)),d//c)
    def dec(self,r):return self.make(ip(r['numerator']),int(r['denominator']))
    def uni(self,p):
        q={int(i):Q(v)for i,v in p if Q(v)};d=math.lcm(*(v.denominator for v in q.values()))
        return self.make({i:int(v*d)for i,v in q.items()},d)
    def c(self,v):
        v=Q(v);return self.make({0:v.numerator},v.denominator)if v else self.z
    def neg(self,x):return(tuple(-v for v in x[0]),x[1])
    def add(self,x,y):
        self.ops+=1;den=math.lcm(x[1],y[1]);return self.make({i:x[0][i]*(den//x[1])+y[0][i]*(den//y[1])for i in range(13)if x[0][i]or y[0][i]},den)
    def sub(self,x,y):return self.add(x,self.neg(y))
    def mul(self,x,y):
        self.ops+=1;conv=[0]*25
        for i,c in enumerate(x[0]):
            if c:
                for j,v in enumerate(y[0]):
                    if v:conv[i+j]+=c*v
        ans=[sum(v*self.red[j][i]for j,v in enumerate(conv)if v)for i in range(13)]
        return self.make({i:v for i,v in enumerate(ans)if v},x[1]*y[1]*self.C)
    def poly(self,p):
        a=[self.dec(r)for r in p]
        return self.trim(a)
    def trim(self,p):
        p=list(p)
        while p and p[-1]==self.z:p.pop()
        return p
    def pa(self,p,q,sign=1):
        return self.trim([self.add(p[i]if i<len(p)else self.z,(q[i]if sign==1 else self.neg(q[i]))if i<len(q)else self.z)for i in range(max(len(p),len(q)))])
    def ps(self,p,c):return self.trim([self.mul(v,c)for v in p])
    def pm(self,p,q):
        if not p or not q:return []
        r=[self.z for _ in range(len(p)+len(q)-1)]
        for i,c in enumerate(p):
            for j,v in enumerate(q):r[i+j]=self.add(r[i+j],self.mul(c,v))
        return self.trim(r)

def main(out):
    started=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));out.mkdir(exist_ok=False,parents=True)
    pb=ROOT/'reports/S61_NONSQUARE_RAW_H_B_K13_D_EUCLID_20261104_v63';r=load(pb/'report.json');assert r['all_passed']
    paths={Path(__file__),ROOT/'research/notes'/(ID+'.md'),pb/'report.json'}
    for rel,h in r['input_sha256'].items():assert sha(ROOT/rel)==h,rel;paths.add(ROOT/rel)
    for rel,h in r['evidence_sha256'].items():assert sha(pb/rel)==h,rel;paths.add(pb/rel)
    c=load(pb/'checkpoint.json');result=load(pb/'result.json');g=ip(c['g13']);assert max(g)==13 and g[0] and math.gcd(*map(abs,g.values()))==1
    A=OwnAlgebra(g)
    sb=ROOT/'reports/S61_NONSQUARE_RAW_H_B_K13_GK_INPUT_20261104_v63';sr=load(sb/'report.json');assert sr['all_passed'];s=load(sb/'source-input.json');assert ip(s['g13'])==g
    for rel,h in sr['evidence_sha256'].items():assert sha(sb/rel)==h,rel;paths.add(sb/rel)
    paths.add(sb/'report.json')
    hb=ROOT/'reports/S61_NONSQUARE_RAW_H_B_K13_H4_INPUT_20261104_v63/H4-monic-field-coefficients.json';hf=load(hb);paths.add(hb)
    bridge=ROOT/'reports/S61_NONSQUARE_RAW_H_B_K13_SOURCE_BRIDGE_20261104_v63_REVIEW/report.json';assert load(bridge)['all_passed'];paths.add(bridge)
    H=A.poly(s['H4']);assert H==[A.uni(p)for p in hf['H4_coefficients_d0_to4']] and len(H)==5 and H[-1]==A.o
    master=ROOT/'reports/S61_NONSQUARE_RAW_H_B_35_MINORS_20261104_v63/certificate.json';allPi=load(master)['localized_ideal_generators_monic_QQ_after_namedunit_and_QQ_only_scaling'];paths.add(master)
    aa=A.make({1:1});tt=A.make(ip(hf['e_numerator']),int(hf['e_denominator']))
    sourcepolys=[allPi[15],allPi[20]];maxa=max(m[0]for p in sourcepolys for m,v in p);maxe=max(m[1]for p in sourcepolys for m,v in p)
    ap=[A.o];tp=[A.o]
    for _ in range(maxa):ap.append(A.mul(ap[-1],aa))
    for _ in range(maxe):tp.append(A.mul(tp[-1],tt))
    # Source provenance: all147 actual minor terms, then both exact H4 reductions.
    gk=[]
    for k,rec in enumerate(s['gk']):
        p=sourcepolys[k];assert {tuple(m):Q(v)for m,v in p}=={tuple(m):Q(v)for m,v in rec['source_polynomial']}
        raw=[A.z for _ in range(5)]
        for m,v in p:raw[m[2]]=A.add(raw[m[2]],A.mul(A.c(v),A.mul(ap[m[0]],tp[m[1]])))
        raw=A.trim(raw);assert raw==A.poly(rec['raw_field_coefficients'])
        quotient=A.poly(rec['H4_quotient']);remainder=A.poly(rec['H4_remainder']);assert A.pa(A.pm(H,quotient),remainder)==raw and len(remainder)<=4
        gk.append(remainder)
    inputs=[H]+gk;assert list(map(len,inputs))==[5,4,4]
    units=[]
    for u in c['units']:
        x=A.dec(u['input']);p=ip(u['primitive_numerator']);content=int(u['numerator_content']);assert content
        assert {i:v*content for i,v in p.items()}=={i:v for i,v in enumerate(x[0])if v}
        cert=u['native_or_constant_certificate'];S,T=ip(cert['S']),ip(cert['T']);R=int(cert['R']);assert R
        assert plus(times(S,p),times(T,g))=={0:R}
        inv=A.dec(u['inverse']);assert inv==A.make({i:v*x[1]for i,v in S.items()},content*R) and A.mul(x,inv)==A.o
        units.append((x,inv))
    assert len(units)==4
    nodes=c['polynomial_nodes'];slp=c['source_cofactor_slp'];assert len(nodes)==len(slp)==15
    values=[];cofs=[]
    for j,(n,operator)in enumerate(zip(nodes,slp)):
        assert n['id']==operator['id']==j;p=A.poly(n['coefficients']);assert len(p)-1==n['degree'];kind=n['kind']
        if kind=='source':
            i=['H4','gR','kR'].index(n['name']);assert p==inputs[i] and operator=={'id':j,'op':'source_basis','index':i}
            cof=[[A.o]if i==k else []for k in range(3)]
        elif kind=='unit_normalization':
            i=n['source'];uid=n['unit'];assert i<j and uid<len(units)
            x,inv=units[uid];assert values[i][-1]==A.dec(n['old_LC'])==x and p==A.ps(values[i],inv) and p[-1]==A.o
            assert operator['op']=='scale' and operator['source']==i and A.dec(operator['scalar'])==inv
            cof=[A.ps(z,inv)for z in cofs[i]]
        else:
            assert kind=='monic_division_remainder';i,k=n['dividend'],n['divisor'];assert i<j and k<j and values[k][-1]==A.o
            q=A.poly(n['quotient']);assert A.pa(A.pm(q,values[k]),p)==values[i] and len(p)<len(values[k])
            assert operator['op']=='subtract_product' and operator['left']==i and operator['right']==k and A.poly(operator['polynomial'])==q
            cof=[A.pa(cofs[i][z],A.pm(q,cofs[k][z]),-1)for z in range(3)]
        combined=[]
        for z in range(3):combined=A.pa(combined,A.pm(cof[z],inputs[z]))
        assert combined==p;values.append(p);cofs.append(cof)
    j=result['gcd_node'];assert j==10 and result['gcd_degree']==1 and A.poly(result['gcd_coefficients'])==values[j]
    assert values[j][-1]==A.o and len(values[j])==2 and A.dec(result['root_s'])==A.neg(values[j][0])
    finalcofs=[A.poly(p)for p in result['source_Bezout_cofactors']];assert finalcofs==cofs[j]
    combined=[]
    for p,q in zip(finalcofs,inputs):combined=A.pa(combined,A.pm(p,q))
    assert combined==values[j]
    for k,row in enumerate(result['all_three_source_divisions']):
        n=nodes[row['zero_remainder_node']];assert row['source']==k and n['dividend']==k and n['divisor']==j and not values[n['id']]
    proof={'all_passed':True,'arithmetic':'Own dense13 rational algebra with a fixed certified monomial reduction table; no author arithmetic import or coefficient inverse/GCD search.',
      'all147_original_g15_k20_source_terms_and_H4_quotient_identities':True,'all4_LC_Bezout_units_and_inverse_products':True,
      'all15_polynomial_nodes_and_source_cofactor_SLP':True,'all3_final_zero_divisions_and_literal_source_Bezout':True,
      'actual_gcd_node':j,'gcd_degree':1,'exact_ideal_identity':'(H4,gR,kR)=(d-s) in (QQ[a]/g13)[d]',
      'original_source_gk_transport':'g=H4*qg+gR and k=H4*qk+kR checked against the147 original source terms. Thus d-s=A H4+B gR+C kR=(A-Bqg-Cqk)H4+B g+C k, without another inverse.',
      'whole_Ka13_no_field_irreducibility_reducedness_assumption':True,
      'scope':'Equality of the fixed necessary common-root cover with Ka13. Remaining actual B conditions are not evaluated by this audit. It proves neither actual Bbad nonemptiness nor emptiness, nor any full cubic/quartic assertion.',
      'quotient_arithmetic_operations':A.ops,'new_solve_polynomial_GCD_inverse_search_GB_point_prime_jobs':0}
    (out/'independent-proof.json').write_text(json.dumps(proof,indent=2)+'\n')
    report={'id':ID,'all_passed':True,'all_source_units_nodes_divisions_Bezout_checked':True,'whole_Ka13_linear_uppercover_identity':True,'actual_Bbad_existence_or_empty_not_claimed':True,
      'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in sorted(paths)},'evidence_sha256':{'independent-proof.json':sha(out/'independent-proof.json')},
      'seconds':time.monotonic()-started,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'frozen':True}
    assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024;signal.alarm(0)
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2),flush=True)
if __name__=='__main__':
    import argparse
    pa=argparse.ArgumentParser();pa.add_argument('--output',type=Path,default=ROOT/'reports'/ID);main(pa.parse_args().output.resolve())
