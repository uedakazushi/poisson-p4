#!/usr/bin/env python3
"""Own sparse Fraction polynomial audit; no producer, CAS cancellation or elimination."""
import hashlib,itertools,json,math,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='AMBIENT_B18_PROPER_RATIONAL_ELIMINATION_20261104_v63_REVIEW'
COUNT={'products':0,'full_coefficient_equalities':0}
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def zero(n):return (0,)*n
def one(n):return {zero(n):Q(1)}
def var(n,i):
    m=[0]*n;m[i]=1;return {tuple(m):Q(1)}
def add(a,b,c=1):
    r=a.copy()
    for m,v in b.items():
        v=r.get(m,Q())+c*v
        if v:r[m]=v
        else:r.pop(m,None)
    return r
def scale(a,c):return {m:c*v for m,v in a.items()if c*v}
def mul(a,b):
    COUNT['products']+=1;r={}
    for m,c in a.items():
        for n,v in b.items():
            k=tuple(i+j for i,j in zip(m,n));r[k]=r.get(k,Q())+c*v
    r={m:c for m,c in r.items()if c};assert len(r)<150000;return r
def power(a,n):
    r=one(len(next(iter(a))))
    for _ in range(n):r=mul(r,a)
    return r
def same(a,b):COUNT['full_coefficient_equalities']+=1;assert a==b
def decode(p):return {tuple(m):Q(n,d)for m,n,d in p if n}
def det(A):
    r=Q();n=len(A)
    for p in itertools.permutations(range(n)):
        term=Q(-1 if sum(p[i]>p[j]for i in range(n)for j in range(i+1,n))%2 else 1)
        for i,j in enumerate(p):term*=A[i][j]
        r+=term
    return r
def rank(A):
    a=[list(map(Q,row))for row in A];r=0
    for j in range(len(a[0])):
        k=next((i for i in range(r,len(a))if a[i][j]),None)
        if k is None:continue
        a[r],a[k]=a[k],a[r];z=a[r][j];a[r]=[c/z for c in a[r]]
        for i in range(r+1,len(a)):
            z=a[i][j]
            if z:a[i]=[c-z*d for c,d in zip(a[i],a[r])]
        r+=1
    return r
def main(out):
    started=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));out.mkdir(exist_ok=False,parents=True)
    pb=ROOT/'reports/AMBIENT_B18_PROPER_RATIONAL_ELIMINATION_20261104_v63';pr=read(pb/'report.json');assert pr['all_passed']
    paths={Path(__file__),ROOT/'research/notes'/(ID+'.md'),pb/'report.json'}
    for rel,h in pr['input_sha256'].items():assert sha(ROOT/rel)==h,rel;paths.add(ROOT/rel)
    for rel,h in pr['evidence_sha256'].items():assert sha(pb/rel)==h,rel;paths.add(pb/rel)
    parent=ROOT/'reports/AMBIENT_B18_PROPER_POLYNOMIAL_MODEL_20261104_v63_REVIEW/report.json';assert read(parent)['all_passed'];paths.add(parent)
    source=read(ROOT/'reports/AMBIENT_B18_PROPER_POLYNOMIAL_MODEL_20261104_v63/actual-proper-source14-equations17.json')
    cert=read(pb/'actual12-flag-divided-polynomials.json');guard=read(pb/'expansion-guard.json')
    names=['v'+str(i)for i in range(9)]+['a','b','c','d','e'];assert cert['original14_variable_order']==names
    Rindices=[9,12,2,3,4,5,6,7,8];Rnames=[names[i]for i in Rindices];assert cert['remaining9_variable_order']==Rnames
    def embed(p):
        out={}
        for m,c in p.items():
            z=[0]*14
            for i,k in zip(Rindices,m):z[i]=k
            out[tuple(z)]=c
        return out
    def restrict(p):
        assert all(all(not m[i]for i in range(14)if i not in Rindices)for m in p)
        return {tuple(m[i]for i in Rindices):c for m,c in p.items()}
    variables=[var(14,i)for i in range(14)];v=variables[:9];a,b,c,d,e=variables[9:]
    D=add(one(14),mul(d,d),-1);assert D[zero(14)]==1
    B=scale(mul(d,add(add(a,d),v[2],-1)),Q(1,2));V0=add(scale(d,2),a,-1)
    EN=add(add(add(add(mul(B,B),scale(mul(mul(d,d),B),Q(1,2))),mul(mul(a,mul(d,d)),add(d,v[2],-1)),-1),mul(mul(d,d),v[4]),-1),mul(d,v[6]))
    VN=add(scale(mul(B,D),Q(1,2)),scale(EN,2));CN=add(scale(mul(EN,add(add(a,d),v[2])),Q(1,2)),mul(mul(d,v[3]),D),-1)
    final={0:(V0,0),1:(VN,1),10:(B,0),11:(CN,1),13:(EN,1)}
    for i,(num,m)in final.items():
        saved=cert['five_explicit_rational_roots'][names[i]];assert saved['D_denominator_power']==m
        same(num,embed(decode(saved['numerator'])))
    equations=list(map(decode,source['equations14variable_sparse']));assert len(equations)==17
    Dpow=[power(D,k)for k in range(7)]
    def substitute(p,roots):
        clearing=max((sum(ex[i]*roots[i][1]for i in roots)for ex in p),default=0);answer={};cache={}
        for ex,coef in p.items():
            base=list(ex);den=0
            for i in roots:base[i]=0
            term={tuple(base):coef}
            for i,(num,powd)in roots.items():
                if ex[i]:
                    key=(i,ex[i])
                    if key not in cache:cache[key]=power(num,ex[i])
                    term=mul(term,cache[key]);den+=powd*ex[i]
            term=mul(term,Dpow[clearing-den]);answer=add(answer,term)
        return answer,clearing
    # Literal triangular ideal equalities, in the actual source14 ring.
    stage_v1=add(scale(B,Q(1,2)),scale(e,2))
    stages=[(2,0,V0,0,{zero(14):Q(1,3)}),(0,10,B,0,{zero(14):Q(2)}),
      (4,1,stage_v1,0,{zero(14):Q(2)}),(3,13,EN,1,scale(D,-1)),(1,11,CN,1,{zero(14):Q(2)})]
    prior={}
    for k,(row,pivot,num,den,lead)in enumerate(stages):
        rec=cert['five_triangular_unit_identities'][k];assert rec['source_equation_row']==row and rec['pivot']==names[pivot]
        left,lp=substitute(equations[row],prior);right=mul(lead,add(mul(Dpow[den],variables[pivot]),num,-1))
        same(mul(left,Dpow[den]),mul(right,Dpow[lp]));prior[pivot]=(num,den)
    for row in range(5):same(substitute(equations[row],final)[0],{})
    H=[];recordchecks=[]
    for k,rec in enumerate(cert['all12_actual_remaining_and_divided_records']):
        row=rec['original_equation_row'];assert row==k+5
        numerator,m=substitute(equations[row],final);assert m==rec['D_denominator_power']==guard['D_clearing_powers12'][k]
        numerator=restrict(numerator);C0=rec['constant_denominator'];content=rec['primitive_content_unit'];assert C0 and content
        integer=decode(rec['cleared_numerator']);Hi=decode(rec['divided_primitive_polynomial'])
        same(scale(numerator,C0),integer);same(integer,scale(mul(var(9,1),Hi),content))
        assert all(c.denominator==1 for c in integer.values()) and math.lcm(*(c.denominator for c in numerator.values()))==C0
        assert math.gcd(*(abs(int(z))for z in Hi.values()))==1 and not Hi.get(zero(9),Q())
        # An independent rational source cross-product, without cancellation.
        same(scale(numerator,C0),scale(mul(var(9,1),Hi),content))
        H.append(Hi);recordchecks.append({'source_equation':row,'D_power':m,'QQ_denominator':C0,'primitive_content':content,'H_terms':len(Hi)})
    J=[[p.get(tuple(int(k==j)for k in range(9)),Q())for j in range(9)]for p in H]
    sameJ=[[Q(c)for c in row]for row in cert['origin_Jacobian12_by9']];assert J==sameJ and rank(J)==5
    rows=cert['candidate_unit_pivot_rows'];cols=cert['candidate_unit_pivot_columns'];minor=det([[J[i][j]for j in cols]for i in rows])
    assert minor==Q(cert['candidate_unit_pivot_minor'])==19660800
    assert [Rnames[i]for i in cols]==cert['candidate_pivot_names']==['v3','v4','v5','v7','v8']
    actualdegrees=[max(sum(m[j]for j in cols)for m in H[i])for i in rows];assert actualdegrees==cert['candidate_equation_degrees_in_joint_pivots']
    proof={'all_passed':True,'arithmetic':'Own sparse QQ dictionaries and denominator cross-products; no producer import or symbolic cancel.',
      'all5_literal_triangular_ideal_equalities':True,'only_variable_denominator':'D=1-d^2, D(0)=1',
      'all5_final_roots_bound_to_exact14_source':True,'all12_source_rational_numerators_and_exact_d_divisions':True,
      'records':recordchecks,'own_divided_Jacobian_rank':5,'own_fixed_pivot_minor':str(minor),
      'whole_E4':'QQ[[a,d,v2,...,v8]]/(d*H0,...,d*H11)',
      'core_only':'The separately divided ideal (H0,...,H11) defines an auxiliary core. Its selected H pivots may be eliminated there; setting them to zero in whole E4 would change the scheme and is not done.',
      'proof':'The five source Fi are triangular unit multiples of the five graph generators in QQ[variables,1/D]. Consequently the selected5 quotient equals the literal rational graph, not merely its point set. The twelve independently checked source cross-products have Fi(root)=content*d*Hi/(constant*D^mi), with all prefactors units at the origin. The previously accepted full source-gauge E4 model therefore has exactly the displayed completed ideal.',
      'no_D0_CM_depth_socle_or_next_elimination_claim':True,'new_GB_sourcejet_point_prime_or_CAS_jobs':0,'counts':COUNT}
    (out/'independent-proof.json').write_text(json.dumps(proof,indent=2)+'\n')
    report={'id':ID,'all_passed':True,'exact5_triangular_and12_flag_division_identities':True,'E4_not_divided_core_not_D0':True,
      'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in sorted(paths)},'evidence_sha256':{'independent-proof.json':sha(out/'independent-proof.json')},
      'seconds':time.monotonic()-started,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'frozen':True}
    assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024;signal.alarm(0)
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2),flush=True)
if __name__=='__main__':
    import argparse
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,default=ROOT/'reports'/ID)
    main(ap.parse_args().output.resolve())
