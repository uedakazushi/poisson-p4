#!/usr/bin/env python3
"""Independent sparse Fraction verification of normalized E4 nil generators."""
import argparse,hashlib,itertools,json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='AMBIENT_B18_NORMALIZER_NIL_GENERATORS_20261104_v63_REVIEW';Z=(0,0,0,0)
counts={'products':0,'pair_products':0,'identities':0,'max_terms_one_poly':0};start=0
def read(p):
    with p.open()as f:return json.load(f)
def sha(p):
    with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def dump(p,d):p.write_text(json.dumps(d,indent=2)+'\n')
def ck(p):
    counts['max_terms_one_poly']=max(counts['max_terms_one_poly'],len(p));assert len(p)<150000 and time.monotonic()-start<52;return p
def dec(v):
    p={}
    for m,c in v:
        m=tuple(m);assert len(m)==4 and all(x>=0 for x in m);c=Q(c);assert c and m not in p;p[m]=c
    return ck(p)
def const(c):return {Z:Q(c)}if c else{}
def var(j):m=list(Z);m[j]=1;return{tuple(m):Q(1)}
def add(a,b):
    p=dict(a)
    for m,c in b.items():
        p[m]=p.get(m,Q(0))+c
        if not p[m]:del p[m]
    return ck(p)
def neg(a):return{m:-c for m,c in a.items()}
def sub(a,b):return add(a,neg(b))
def mul(a,b):
    counts['products']+=1;counts['pair_products']+=len(a)*len(b);p={}
    for m,c in a.items():
        for n,d in b.items():
            k=tuple(x+y for x,y in zip(m,n));p[k]=p.get(k,Q(0))+c*d
    return ck({m:c for m,c in p.items()if c})
def sm(c,a):return{m:Q(c)*v for m,v in a.items()if c*v}
def sump(xs):
    p={}
    for x in xs:p=add(p,x)
    return p
def comp(p,images):
    out={}
    for m,c in p.items():
        t=const(c)
        for j,k in enumerate(m):
            for _ in range(k):t=mul(t,images[j])
        out=add(out,t)
    return out
def eq(a,b):counts['identities']+=1;assert a==b
def mat(v):return[[dec(p)for p in row]for row in v]
def trans(A):return[list(c)for c in zip(*A)]
def mm(A,B):
    assert A and B and len(A[0])==len(B)
    return[[sump(mul(A[i][k],B[k][j])for k in range(len(B)))for j in range(len(B[0]))]for i in range(len(A))]
def meq(A,B):
    assert len(A)==len(B)and len(A[0])==len(B[0])
    for rowa,rowb in zip(A,B):
        for a,b in zip(rowa,rowb):eq(a,b)
def scalarI(p,n):return[[p if i==j else{} for j in range(n)]for i in range(n)]
def origin(p):return p.get(Z,Q(0))
def det_constant(A):
    n=len(A);return sum((-1)**sum(perm[i]>perm[j]for i in range(n)for j in range(i+1,n))*Q(1)*__import__('functools').reduce(lambda x,y:x*y,(A[i][perm[i]]for i in range(n)),Q(1))for perm in itertools.permutations(range(n)))
def line(p,v):
    out={}
    for m,c in p.items():
        q=c
        for x,e in zip(v,m):q*=x**e
        k=(sum(m),0,0,0);out[k]=out.get(k,Q(0))+q
    return{k:c for k,c in out.items()if c}
def minmons(ms):
    ms=set(ms);return sorted(m for m in ms if not any(n!=m and all(a<=b for a,b in zip(n,m))for n in ms))
def main(out):
    global start
    start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(56);out.mkdir(parents=True,exist_ok=False);R=ROOT/'reports';b=R/'AMBIENT_B18_NORMALIZER_NIL_GENERATORS_20261104_v63';r=read(b/'report.json')
    assert sha(b/'report.json')=='9dccf13ca56ed2f5d8676680d9e023b093ebfa740de8bbaa1dcafc5fde7a1af3'and r['all_passed'];ep=b/'actual-E4-nil-two-generators-and-relations.json';assert sha(ep)=='3ec81afe22c28db5dea27eb1f8834b7383d8668fd3dea63d2acb65fa8b6fbc0e';e=read(ep)
    for p,h in r['input_sha256'].items():assert sha(ROOT/p)==h
    source=read(R/'AMBIENT_B18_PROPER_CORE_FIVE_LINES_20261104_v63/five-line-and-generic-models.json')
    for dep in ['AMBIENT_B18_PROPER_CORE_FIVE_LINES_20261104_v63_REVIEW','AMBIENT_B18_PROPER_CORE_CI_MEMBERSHIP_20261104_v63_REVIEW','AMBIENT_B18_PROPER_CORE_ELIMINATION_20261104_v63_REVIEW']:assert read(R/dep/'report.json')['all_passed']
    old=[var(j)for j in range(4)];a,d,x,t=old;new=[var(j)for j in range(4)];c,rr,ll,h=new
    fwd=[sm(Q(1,2),add(d,x)),sm(Q(1,6),sub(sub(d,x),sm(2,a))),sm(Q(1,3),add(sub(d,x),a)),sub(add(sub(a,d),x),sm(2,t))]
    back=[add(sm(-2,rr),ll),add(add(c,rr),ll),sub(sub(c,rr),ll),sm(Q(1,2),sub(sub(sm(-4,rr),ll),h))]
    for j,k in enumerate(['c','r','l','h']):eq(fwd[j],dec(e['forward_GL4'][k]));eq(comp(fwd[j],back),new[j])
    for j,k in enumerate(['a','d','x','t']):eq(back[j],dec(e['inverse_GL4'][k]));eq(comp(back[j],fwd),old[j])
    J=[[fwd[i].get(tuple(int(k==j)for k in range(4)),Q(0))for j in range(4)]for i in range(4)];assert det_constant(J)==Q(1,3)==Q(e['GL4_determinant'])
    Qold=[dec(source[k])for k in ['source_Q0','source_Q1','source_f']]
    for i,k in enumerate(['old_source_Q0','old_source_Q1','old_source_f']):eq(Qold[i],dec(e[k]))
    eq(Qold[2],mul(add(a,x),mul(sub(sub(x,a),d),sub(sub(x,a),d))))
    core=[comp(p,back)for p in Qold]
    for p,q in zip(core,e['actual_core_generators3']):eq(p,dec(q))
    g=[h,mul(c,rr),mul(c,ll),mul(rr,ll)]
    for p,q in zip(g,e['I3_generators4']):eq(p,dec(q))
    V=mat(e['source_factorization_coefficients4_by3']);meq(mm([g],V),[core])
    C=[[origin(p)for p in row]for row in V];assert C==[[Q(0),Q(0),Q(0)],[Q(32),Q(0),Q(0)],[Q(11),Q(0),Q(0)],[Q(-9),Q(-90),Q(0)]]
    assert [[Q(v)for v in row]for row in e['minimal_generator_class_matrix4_by3']]==C
    assert e['quad_pivots_in_I3_generators']==[1,3]and e['retained_quad_index']==2
    Delta=dec(e['unit_Delta']);assert len(Delta)==24 and origin(Delta)==-2880==Q(e['unit_Delta_origin'])
    A=[[V[1][0],V[3][0]],[V[1][1],V[3][1]]];eq(Delta,sub(mul(A[0][0],A[1][1]),mul(A[0][1],A[1][0])))
    E=mat(e['old_I3_generators_Delta_reduction_E4_by2']);T=mat(e['old_I3_generators_Delta_source_Q_cofactors_T4_by2']);N=mat(e['unit_change_N4_by4']);ET=[er+tr for er,tr in zip(E,T)]
    Nexpected=[[const(1),{},{},{}],[{},{},const(1),{}],*[list(col)for col in trans(V)[:2]]];meq(N,Nexpected)
    meq(mm(N,ET),scalarI(Delta,4));meq(mm(ET,N),scalarI(Delta,4));assert det_constant([[origin(p)for p in row]for row in N])!=0
    small=[h,g[2]];rhs=mm(E,[[p]for p in small]);QT=mm(T,[[p]for p in core[:2]])
    for j in range(4):eq(mul(Delta,g[j]),add(rhs[j][0],QT[j][0]))
    sy=mat(e['I3_complete_syzygies4_by5']);meq(mm([g],sy),[[{}for _ in range(5)]])
    # All six monomial pair-lcm syzygies reduce to the saved five.
    exps=[next(iter(p))for p in g];pair_checks=0
    for i,j in itertools.combinations(range(4),2):
        lc=tuple(max(v,w)for v,w in zip(exps[i],exps[j]));v=[{}for _ in range(4)];v[i]={tuple(x-y for x,y in zip(lc,exps[i])):Q(1)};v[j]={tuple(x-y for x,y in zip(lc,exps[j])):Q(-1)}
        if i==0:w=[sy[k][j-1]for k in range(4)]
        elif(i,j)==(1,2):w=[sy[k][3]for k in range(4)]
        elif(i,j)==(1,3):w=[sy[k][4]for k in range(4)]
        else:w=[sub(sy[k][4],sy[k][3])for k in range(4)]
        for p,q in zip(v,w):eq(p,q)
        pair_checks+=1
    W=[sy[i]+[V[i][2]]for i in range(4)];Rel=mat(e['nil_complete_core_S4_relations2_by6']);wit=mat(e['all6_relation_source_cofactors3_by6']);meq(Rel,mm(trans(E),W))
    eqw=mm(trans(T),W)
    for i in range(2):
        for j in range(6):eq(wit[i][j],neg(eqw[i][j]))
    for j in range(6):eq(wit[2][j],Delta if j==5 else{})
    meq(mm([small],Rel),mm([core],wit));assert all(origin(p)==0 for row in Rel for p in row)
    ng=[mul(back[1],p)for p in small]
    for p,q in zip(ng,e['nil_generators_in_new_core']):eq(p,dec(q))
    original=[comp(p,fwd)for p in ng];want=[mul(d,sub(add(sub(a,d),x),sm(2,t))),sm(Q(1,6),mul(mul(d,add(d,x)),sub(add(a,d),x)))]
    for i in range(2):eq(original[i],want[i]);eq(original[i],dec(e['nil_generators_in_original_core'][i]));assert min(sum(m)for m in original[i])>=2
    lines={'N':[Q(0),Q(0),Q(0),Q(1)],'C':[Q(0),Q(1),Q(1),Q(0)],'R':[Q(-2),Q(1),Q(-1),Q(-2)],'P':[Q(2),Q(0),Q(-2),Q(1)],'L':[Q(1),Q(1),Q(-1),Q(-1,2)]}
    axis={}
    for key,v in lines.items():
        for pol in Qold:eq(line(pol,v),{})
        prime=[dec(pol)for pol in source['support']['primes'][key]]
        for pol in prime:eq(line(pol,v),{})
        fixedcols={'N':[0,1,2],'C':[0,3,2],'R':[0,2,3],'P':[1,0,2],'L':[0,2,3]}[key]
        coeffs=[[pol.get(tuple(int(k==j)for k in range(4)),Q(0))for j in fixedcols]for pol in prime];assert det_constant(coeffs)!=0
        if key in ['C','R','L']:
            values=[sum(coef*__import__('functools').reduce(lambda q,z:q*z,(u**k for u,k in zip(v,m)),Q(1))for m,coef in pol.items())for pol in fwd];expected=[Q(int(j==['C','R','L'].index(key)))for j in range(4)];assert values==expected;axis[key]=list(map(str,values))
    assert {k for k,v in lines.items()if not v[1]}=={'N','P'}
    # Monomial intersection of the three axis primes, without a GB.
    axis_primes=[[tuple(int(k==j)for k in range(4))for j in [1,2,3]],[tuple(int(k==j)for k in range(4))for j in [0,2,3]],[tuple(int(k==j)for k in range(4))for j in [0,1,3]]]
    inter=axis_primes[0]
    for prime in axis_primes[1:]:inter=minmons(tuple(max(x,y)for x,y in zip(a,b))for a in inter for b in prime)
    assert set(inter)==set(exps)
    assert e['full_E4_nil_module_presentation_shape']==[2,16]and e['nil_minimal_generator_count']==2 and e['no_D0_nil_generators_claimed']
    proof={'radical':'sqrt(dJ)=(d) intersect(u5+I3)=d(u5+I3): N/P contain d, C/R/L transform to the three axes, and d=c+r+l is regular on their reduced intersection.','nil_module':'Ambient S9 has d regular, so multiplication by d gives d(u5+I3)/d(u5+Jcore) ≅ I3/Jcore; allu5 annihilate this module, not the wholeE4.','complete_relations':'The six pair-lcm syzygies form the Taylor first presentation of monomial(h,cr,cl,rl); the five saved columns generate allsix. Completion is flat. Add the three actualJcore columns. UnitN inverse=[E,T]/Delta removes Q0,Q1 and leaves exactly six presentation columns; append ten u_j*e_k over fullS9/E4.','minimality':'I3/mI3 has4 generators; the actual3Jcore classes have rank2. Two remain, equivalently every entry of the complete2x16 relation matrix vanishes at the maximal ideal. Ambient nilgenerators lieinmE4^2 but that is not mE4*Nil.','unit_scope':'Delta is the actual24-term formal unit with constant-2880; no constant inverse is substituted for Delta.','scope':'NormalizedE4 only, no claim about originalD0 nil generators, depth, or classification.','relation_count_not_claimed_minimal':True}
    dump(out/'proof.json',proof);audit={'all_literal_source_and_matrix_identities_passed':True,'GL4_determinant':'1/3','Delta_term_count':24,'Delta_constant':'-2880','both_N_inverse_numerator_products_verified':True,'all6_relation_source_cofactors_verified':True,'all6_pair_syzygies_generated_by_saved5':pair_checks,'actual_Q_class_matrix_rank':2,'complete_presentation_core':[2,6],'complete_presentation_S9_E4':[2,16],'minimal_nil_generator_count':2,'axes':axis,'new_source_jet_point_prime_or_GB':False,**counts}
    dump(out/'arithmetic-audit.json',audit);note=ROOT/f'research/notes/{ID}.md';assert note.exists()
    result={'id':ID,'all_passed':True,'frozen':True,'normalized_E4_only':True,'minimal_nil_generators':2,'complete_core_presentation':[2,6],'complete_full_presentation':[2,16],'arithmetic':audit,'input_report_sha256':sha(b/'report.json'),'input_evidence_sha256':sha(ep),'source_sha256':{str(p.relative_to(ROOT)):sha(p)for p in [Path(__file__),note]},'evidence_sha256':{p.name:sha(p)for p in out.iterdir()if p.is_file()},'no_original_D0_nil_claim':True,'seconds':time.monotonic()-start}
    assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024;signal.alarm(0);dump(out/'report.json',result);print(json.dumps({'id':ID,'all_passed':True,'seconds':result['seconds'],'report_sha256':sha(out/'report.json'),'arithmetic':audit},indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output.resolve())
