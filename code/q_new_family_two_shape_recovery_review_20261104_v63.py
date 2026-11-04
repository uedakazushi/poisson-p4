#!/usr/bin/env python3
"""Independent sparse coefficient audit of dense q/B recovery and two-shape constraints."""
import argparse
from fractions import Fraction as Q
import gzip
import hashlib
import json
from pathlib import Path
import resource
import signal
import time
import sympy as s

ROOT=Path(__file__).resolve().parents[1]
STEM='Q_NEW_FAMILY_TWO_SHAPE_RECOVERY_20261104_v63_REVIEW'

def read(p):
    with (gzip.open(p,'rt') if p.suffix=='.gz' else p.open()) as f:return json.load(f)
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()
def add(*polys):
    ans={}
    for p in polys:
        for e,c in p.items():ans[e]=ans.get(e,Q(0))+c
    return {e:c for e,c in ans.items() if c}
def scale(p,c):return {e:c*a for e,a in p.items() if c*a}
def mul(p,q):
    out={}
    for e,c in p.items():
        for f,d in q.items():
            g=tuple(a+b for a,b in zip(e,f));out[g]=out.get(g,Q(0))+c*d
    return {e:c for e,c in out.items() if c}
def power(p,n):
    assert n>=0
    one={(0,)*len(next(iter(p))):Q(1)}
    ans=one
    while n:
        if n&1:ans=mul(ans,p)
        n//=2
        if n:p=mul(p,p)
    return ans
def shift(p,ex):return {tuple(a+b for a,b in zip(e,ex)):c for e,c in p.items()}
def var(n,i):return {tuple(int(j==i) for j in range(n)):Q(1)}
def enc(p):return [[str(c.numerator),str(c.denominator),*e] for e,c in sorted(p.items(),reverse=True)]
def from_record(d):return {tuple(e):Q(int(a),int(b)) for a,b,*e in d['coefficients']}

def main(out):
    start=time.monotonic();signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
    out.mkdir(parents=True,exist_ok=False)
    base=ROOT/'reports/Q_NEW_FAMILY_LOW_INVARIANT_BIRATIONAL_COMPRESSION_20261104_v63'
    recovery=ROOT/'reports/Q_NEW_FAMILY_MOBIUS_SHAPE_RECOVERY_20261104_v63'
    hom=ROOT/'reports/Q_NEW_FAMILY_TWO_SHAPE_HOMOGENIZATION_20261104_v63'
    crp=ROOT/'reports/Q_NEW_FAMILY_LOW_INVARIANT_BIRATIONAL_COMPRESSION_20261104_v63_REVIEW/report.json'
    assert read(crp)['all_passed']
    inputs=[crp,Path(__file__),ROOT/'research/notes'/(STEM+'.md')]
    ar,hr=read(recovery/'report.json'),read(hom/'report.json')
    assert ar['all_passed'] and hr['all_passed']
    for folder,rp in [(recovery,ar),(hom,hr)]:
        for rel,d in rp['input_sha256'].items():assert sha(ROOT/rel)==d,rel
        for rel,d in rp['evidence_sha256'].items():assert sha(folder/rel)==d,rel
        inputs.extend([folder/'report.json',folder/'certificate.json'])
    rc,hc=read(recovery/'certificate.json'),read(hom/'certificate.json')
    source={};meta={}
    for name in ['F4','p0','p1','p2','p3','p4','N12']:
        path=base/(name+'.json.gz');d=read(path);inputs.append(path);meta[name]=d
        source[name]={tuple(e):Q(c)/Q(d['denominator_integer']) for c,*e in d['terms']}
    # The source numerator ring is QQ[r,s,q,B]. These identities are checked
    # by coefficient cross-products only; no rational-expression gcd/cancel.
    r,ss,q,B=[var(4,i) for i in range(4)];one={(0,0,0,0):Q(1)}
    rm1=add(r,scale(one,-1))
    F=add(scale(mul(mul(r,rm1),power(q,2)),100),scale(mul(r,ss),-4),power(ss,2),scale(r,4))
    hh=add(power(ss,2),scale(mul(r,power(q,2)),-25));alpha=Q(2,405)
    linear=add(scale(q,5),scale(ss,-1),one)
    DD=scale(mul(mul(power(r,4),rm1),mul(power(q,2),linear)),2359296000)
    assert source['F4']==scale(F,alpha)
    assert source['p0']==scale(mul(B,DD),alpha**2)
    assert meta['F4']['denominator_h_power']==0
    assert all(meta['p'+str(i)]['denominator_h_power']==1 for i in range(5))
    assert meta['N12']['denominator_h_power']==2
    aa=scale(power(r,2),80)
    bb=add(scale(mul(power(r,2),ss),-6),scale(power(r,2),6),scale(mul(r,ss),-11),scale(r,10),ss)
    cc=scale(mul(r,rm1),40);dd=scale(mul(mul(r,rm1),add(one,scale(ss,-1))),8)
    nR=add(mul(aa,q),bb);dR=add(mul(cc,q),dd)
    assert mul(add(source['p3'],scale(source['p4'],4)),dR)==mul(source['p0'],nR)
    J=add(mul(add(scale(r,10),scale(one,-1)),ss),scale(r,-10))
    determinant=add(mul(aa,dd),scale(mul(bb,cc),-1))
    assert determinant==scale(mul(mul(r,power(rm1,2)),J),-40) and determinant
    # Constant B^0 relations and exact rank, independently from the sparse
    # coefficient columns; this does not claim a new p1,p2 combination.
    constant=[{e[:3]:c for e,c in source['p'+str(i)].items() if e[3]==0} for i in range(1,5)]
    extra=mul(hh,mul(source['F4'],source['F4']))
    constant.append({e[:3]:c for e,c in extra.items()})
    monoms=sorted(set().union(*(set(v) for v in constant)))
    matrix=s.Matrix([[s.Rational(v.get(e,0).numerator,v.get(e,Q(0)).denominator)
                      if e in v else 0 for v in constant] for e in monoms])
    assert matrix.rank()==3
    assert [list(map(str,v)) for v in matrix.nullspace()]==rc['constant_B0_kernel']
    # Work now in QQ[r,s,R], using independent sparse coefficient arithmetic.
    r,ss,R=[var(3,i) for i in range(3)];one={(0,0,0):Q(1)};rm1=add(r,scale(one,-1))
    T=add(mul(rm1,R),scale(r,-2))
    L=add(scale(mul(mul(R,r),rm1),8),scale(power(r,2),-6),scale(r,-11),one)
    Nq=add(mul(add(ss,scale(one,-1)),L),scale(rm1,-1))
    dq=scale(mul(r,T),40)
    J=add(mul(add(scale(r,10),scale(one,-1)),ss),scale(r,-10))
    Hn=add(scale(mul(mul(r,power(T,2)),power(ss,2)),64),scale(power(Nq,2),-1))
    Fn=add(mul(rm1,power(Nq,2)),scale(mul(mul(r,power(T,2)),add(power(ss,2),scale(mul(r,ss),-4),scale(r,4))),16))
    factors={'r':r,'rm1':rm1,'T':T,'Nq':Nq,'J':J,'Hn':Hn,'Fn':Fn}
    assert all(factors.values())
    for name,p in factors.items():assert from_record(hc['factors'][name])==p
    # Forward and inverse Mobius compositions, through affine coefficients
    # rather than a generic substitution/cancellation routine.
    aa=scale(power(r,2),80)
    bb=add(scale(mul(power(r,2),ss),-6),scale(power(r,2),6),scale(mul(r,ss),-11),scale(r,10),ss)
    cc=scale(mul(r,rm1),40);dd=scale(mul(mul(r,rm1),add(one,scale(ss,-1))),8)
    assert add(mul(aa,Nq),mul(bb,dq))==mul(R,add(mul(cc,Nq),mul(dd,dq)))
    # Inverse composed with source R(q), using q as the third variable of a
    # separate polynomial ring; the second affine-coefficient identity
    # certifies inverse composition in the rational function field.
    q=var(3,2)
    nA={e[:2]+(0,):c for e,c in Nq.items() if e[2]==1}
    nB={e:c for e,c in Nq.items() if e[2]==0}
    dA={e[:2]+(0,):c for e,c in dq.items() if e[2]==1}
    dB={e:c for e,c in dq.items() if e[2]==0}
    nr=add(mul(aa,q),bb);dr=add(mul(cc,q),dd)
    assert add(mul(nA,nr),mul(nB,dr))==mul(q,add(mul(dA,nr),mul(dB,dr)))
    # All substitutions of the named open factors have exact polynomial
    # numerator identities. Nonvanishing follows from polynomial nonzero.
    hhom=add(mul(power(ss,2),power(dq,2)),scale(mul(r,power(Nq,2)),-25))
    fhom=add(scale(mul(mul(r,rm1),power(Nq,2)),100),mul(add(scale(mul(r,ss),-4),power(ss,2),scale(r,4)),power(dq,2)))
    assert scale(mul(mul(r,power(T,2)),hhom),64)==mul(Hn,power(dq,2))
    assert scale(mul(mul(r,power(T,2)),fhom),16)==mul(Fn,power(dq,2))
    linearhom=add(scale(Nq,5),mul(add(one,scale(ss,-1)),dq))
    assert scale(mul(mul(r,T),linearhom),8)==mul(mul(rm1,J),dq)
    dhom=scale(mul(mul(mul(power(r,4),rm1),power(Nq,2)),linearhom),2359296000)
    dback=scale(mul(mul(mul(mul(r,power(rm1,2)),power(Nq,2)),J),power(dq,3)),184320)
    assert mul(dhom,power(T,3))==dback
    # Every univariate-q numerator is homogenized independently. This checks
    # all 2310 coefficients, not only term counts or samples.
    split={};homogeneous={};degrees={}
    for label,src,bpower in [('A1','p1',0),('L1','p1',1),('A2','p2',0),('L2','p2',1),('M0','N12',1),('M1','N12',2)]:
        original={e[:3]:c for e,c in source[src].items() if e[3]==bpower}
        assert original
        degree=max(e[2] for e in original);acc={}
        for (i,j,k),c in original.items():
            term=mul(power(Nq,k),power(dq,degree-k))
            acc=add(acc,scale(shift(term,(i,j,0)),c))
        assert acc==from_record(hc['homogenized_input_polynomials'][label])
        assert degree==hc['homogenized_input_polynomials'][label]['q_degree']
        split[label]=original;homogeneous[label]=acc;degrees[label]=degree
    # A factored rational value is (constant, polynomial numerator,
    # signed factor exponents). Formula construction is separate from the
    # author's tracked division algorithm; polynomial cross-products check
    # the complete normalized output values.
    def exadd(*ee):
        d={}
        for e in ee:
            for k,v in e.items():d[k]=d.get(k,0)+v
        return {k:v for k,v in d.items() if v}
    def exscale(e,n):return {k:n*v for k,v in e.items() if n*v}
    hvalue=(Q(1,64),{'Hn':1,'r':-1,'T':-2})
    fvalue=(Q(1,16),{'Fn':1,'r':-1,'T':-2})
    dvalue=(Q(184320),{'r':1,'rm1':2,'Nq':2,'J':1,'T':-3})
    specs={'A1':('A1',2,[-1,-2,0]),'C1':('L1',2,[0,0,-1]),
           'A2':('A2',2,[-1,-2,0]),'C2':('L2',2,[0,0,-1]),
           'E1':('M0',3,[-1,-1,-1]),'E2':('M1',3,[0,1,-2])}
    norm={};pair_checks=[]
    for name,(inp,alphapower,powers) in specs.items():
        m=degrees[inp];co=Q(40)**(-m)/alpha**alphapower;ee={'r':-m,'T':-m}
        for (scalar,ex),exponent in zip([hvalue,fvalue,dvalue],powers):
            co*=scalar**exponent;ee=exadd(ee,exscale(ex,exponent))
        raw=hc['rational_coefficient_pairs'][name]
        pn=from_record(raw['numerator']);pc=Q(raw['constant']);ed={k:-v for k,v in raw['denominator_factors'].items()}
        left=scale(homogeneous[inp],co);right=scale(pn,pc)
        for k in factors:
            mn=min(ee.get(k,0),ed.get(k,0))
            left=mul(left,power(factors[k],ee.get(k,0)-mn))
            right=mul(right,power(factors[k],ed.get(k,0)-mn))
        assert left==right,name
        assert all(v>0 for v in raw['denominator_factors'].values())
        norm[name]=(pn,pc,raw['denominator_factors'])
        pair_checks.append({'name':name,'normalized_terms':len(pn),'full_cross_product_equal':True})
    own_constraints={}
    for name,names,zpowers in [('Z1',['A1','C1'],[0,1]),('Z2',['A2','C2'],[0,1]),('ZN',['E1','E2'],[1,2])]:
        aa,bb=[norm[a] for a in names]
        lcm={k:max(aa[2].get(k,0),bb[2].get(k,0)) for k in factors};lcm={k:v for k,v in lcm.items() if v}
        den=one
        for k,e in lcm.items():den=mul(den,power(factors[k],e))
        nums=[]
        for pn,pc,ed in [aa,bb]:
            nn=scale(pn,pc)
            for k,e in lcm.items():nn=mul(nn,power(factors[k],e-ed.get(k,0)))
            nums.append(nn)
        rr=hc['localized_constraints'][name]
        assert rr['Z0_powers']==zpowers and rr['denominator_factors']==lcm
        assert from_record(rr['coefficient1'])==nums[0]
        assert from_record(rr['coefficient2'])==nums[1]
        assert from_record(rr['denominator'])==den
        allterms=set(nums[0])|set(nums[1])|set(den)
        degrees=[max(e[i] for e in allterms) for i in range(3)]
        assert degrees[:2]==rr['shape_bidegree'] and degrees[2]==rr['R_degree']
        assert sum(map(len,nums+[den]))==rr['terms_with_independent_target_parameter']
        own_constraints[name]={'shape_bidegree':degrees[:2],'R_degree':degrees[2],
                               'literal_all_coefficients_equal':True,'terms':sum(map(len,nums+[den])),
                               'denominator_factors':lcm}
    evidence={'r_s_R_factors':{k:enc(v) for k,v in factors.items()},
              'B_recovery':'B=Z0*h*F^2/[2359296000*r^4*(r-1)*q^2*(5*q-s+1)]',
              'q_recovery':'q=((s-1)*(8*R*r*(r-1)-6*r^2-11*r+1)-(r-1))/(40*r*((r-1)*R-2*r))',
              'six_rational_pairs':pair_checks,'three_constraints':own_constraints,
              'full_homogenized_terms':{k:len(p) for k,p in homogeneous.items()},
              'source_field_identity':'QQ(xi,r,s,q,B)=QQ(xi,r,s,R,Z0)',
              'dense_chart_units':['xi','r','s','Z0','r-1','T','Nq','J','Hn','Fn'],
              'normalizations':{'Z0':'p0/F4^2','R':'(p3+4*p4)/p0','Z1':'p1/F4^2','Z2':'p2/F4^2','ZN':'N12/F4^3'},
              'actual_NEW_target_field_index_one_or_equation_not_asserted':True}
    (out/'independent-certificate.json').write_text(json.dumps(evidence,indent=2)+'\n')
    # Notes are finalized separately before this bounded review executes.
    for st in ['Q_NEW_FAMILY_MOBIUS_SHAPE_RECOVERY_20261104_v63','Q_NEW_FAMILY_TWO_SHAPE_HOMOGENIZATION_20261104_v63']:
        inputs.append(ROOT/'research/notes'/(st+'.md'))
    inputs.append(ROOT/'reports/Q_NEW_FAMILY_TWO_SHAPE_RECOVERY_20261104_v63/report.json')
    report={'id':STEM,'all_passed':True,
            'independent_sparse_coefficient_arithmetic_not_author_producer':True,
            'exact_B_free_R_source_identity_and_both_q_B_inverse_compositions':True,
            'Mobius_determinant_nonzero_dense_chart_and_whole_known_source_field_generation':True,
            'all_six_homogenized_input_polynomials_and_all_six_rational_pairs_verified':True,
            'all_three_full_cleared_constraint_coefficient_tables_verified':True,
            'localized_equivalence_only_no_boundary_extension_or_properness':True,
            'no_actual_NEW_equation_target_field_index_one_or_degree_claim':True,
            'failed_60s_generic_cancellation_preserved_not_retried':True,
            'no_GB_elimination_resultant_new_point_source26_or_generic_rational_gcd':True,
            'constraint_summaries':own_constraints,
            'input_sha256':{str(p.relative_to(ROOT)):sha(p) for p in dict.fromkeys(inputs)},
            'evidence_sha256':{'independent-certificate.json':sha(out/'independent-certificate.json')},
            'seconds':time.monotonic()-start,'frozen':True}
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    signal.alarm(0)
    print(json.dumps({k:report[k] for k in ['all_passed','seconds','constraint_summaries']},indent=2))

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=ROOT/'reports'/STEM)
    args=parser.parse_args();main(args.output.resolve())
