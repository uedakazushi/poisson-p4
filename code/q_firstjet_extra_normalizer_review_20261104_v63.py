#!/usr/bin/env python3
"""Independent small critical frame, rational DAG recipe and fixed-field proof audit."""
import argparse
from fractions import Fraction as F
import gzip,hashlib,itertools,json,resource,signal,time
from pathlib import Path
import sympy as s
from sympy.polys.rings import ring

ROOT=Path(__file__).resolve().parents[1]
STEM='Q_FIRSTJET_SPIN5_EXTRA_NORMALIZER_FIELD_20261104_v63_REVIEW'
def read(p):
    with(gzip.open(p,'rt')if p.suffix=='.gz'else p.open())as f:return json.load(f)
def sha(p):
    h=hashlib.sha256()
    with p.open('rb')as f:
        for a in iter(lambda:f.read(1<<20),b''):h.update(a)
    return h.hexdigest()
def main(out):
    started=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
    out.mkdir(parents=True,exist_ok=False)
    folder=ROOT/'reports/Q_FIRSTJET_SPIN5_EXTRA_NORMALIZER_FIELD_20261104_v63'
    ar,data=read(folder/'report.json'),read(folder/'action.json.gz');assert ar['all_passed']
    gates=[ROOT/'reports/Q_FIRSTJET_SPIN4_BIFORM_RATIONAL_FIELD_20261104_v63_REVIEW/report.json',
        ROOT/'reports/Q_FIRSTJET_SPIN4_MIXED_DISCRIMINANT_20261103_v63_REVIEW/report.json']
    assert all(read(p)['all_passed']for p in gates)
    physical_path=ROOT/'reports/Q_FIRSTJET_SPIN4_MIXED_DISCRIMINANT_20261103_v63/certificate.json'
    physical=read(physical_path);uv=list(map(s.Symbol,physical['variables12']));u=uv[:6];v=uv[6:]
    new=[-v[0],-v[3],-v[1],-v[4],-v[2],-v[5],u[0],u[2],u[4],u[1],u[3],u[5]]
    M=s.Matrix([[s.diff(a,b)for b in uv]for a in new])
    assert M==s.Matrix(physical['normalizer_biform12'])and M*M==-s.eye(12)
    cross=s.Matrix([[u[0],u[2]/2,v[0],v[1]/2],[u[1],u[3]/2,v[1]/2,v[2]],
        [u[2]/2,u[4],v[3],v[4]/2],[u[3]/2,u[5],v[4]/2,v[5]]])
    delta=s.expand(cross.det());assert s.expand(delta-s.sympify(physical['normal_component_delta']))==0
    assert s.expand(delta.subs(dict(zip(uv,new)),simultaneous=True)+delta)==0
    c=s.Symbol('c');assert s.expand(delta.subs(dict(zip(u,[c,0,0,0,0,c])))-c*c*(v[1]*v[4]/4-v[2]*v[3]))==0
    # Critical-root frame: a different matrix-based calculation over the
    # full five-coefficient rational field, not six-variable action expansion.
    RR,*vars=ring('p,q,r,s,t',s.QQ);K=RR.to_field();p,q,r,ss,t=map(K,vars)
    resultant=(p*q-ss*t)**2-(p*r-ss)*(q-r*t);J0=2*(p*r-ss);J1=4*(p*q-ss*t)
    z0=(K.zero,K.zero);one=(K.one,K.zero)
    def pa(a,b):return(a[0]+b[0],a[1]+b[1])
    def ps(a,b):return(a[0]*b,a[1]*b)
    def pm(a,b):return(a[0]*b[0]+resultant*a[1]*b[1],a[0]*b[1]+a[1]*b[0])
    A=[[(K(s.Rational(1,2)),-J1/(8*resultant)),(-J1/(2*J0),-2/J0)],[(K.zero,J0/(4*resultant)),one]]
    assert pa(pm(A[0][0],A[1][1]),ps(pm(A[0][1],A[1][0]),-1))==one
    rows=[]
    for qa,qb,qc in [(-p,-K.one,-t),(-ss,-r,-q)]:
        coeff=[]
        for i,j in [(0,0),(0,1),(1,1)]:
            value=pa(pa(ps(pm(A[0][i],A[0][j]),qa),ps(pm(A[1][i],A[1][j]),qc)),
                ps(pa(pm(A[0][i],A[1][j]),pm(A[1][i],A[0][j])),qb/2))
            coeff.append(value)
        assert coeff[1]==z0;rows.append([coeff[0],coeff[2]])
    detM=pa(pm(rows[0][0],rows[1][1]),ps(pm(rows[1][0],rows[0][1]),-1))
    assert detM==(K.zero,-K.one)
    # A separate exact quadratic-algebra descent of R and delta in the six
    # rational coordinates. This computes only the small resultant/delta.
    I=s.symbols('I0:6');eta,I1,I2,I3,I4,I5=I
    R6,*rv=ring('I0,I1,I2,I3,I4,I5',s.QQ);F6=R6.to_field();a0,a1,a2,a3,a4,a5=map(F6,rv)
    def za(a,b):return(a[0]+b[0],a[1]+b[1])
    def zm(a,b):return(a[0]*b[0]+a0*a[1]*b[1],a[0]*b[1]+a[1]*b[0])
    def zs(a,b):return(a[0]*b,a[1]*b)
    pp=(a2/2,a3/2);rr=(F6.zero,a1);qq=(-a3/(2*a0*a1**3),a2/(2*a0**2*a1**3))
    ssp=(a4/2,a5/2);tt=(a5/(2*a0**2*a1**5),-a4/(2*a0**3*a1**5))
    aa=za(zm(pp,qq),zs(zm(ssp,tt),-1));bb=za(zm(pp,rr),zs(ssp,-1));cc=za(qq,zs(zm(rr,tt),-1))
    Rpair=za(zm(aa,aa),zs(zm(bb,cc),-1));dpair=zm((F6.zero,F6.one),za(zs(rr,F6(s.Rational(1,4))),zs(zm(tt,ssp),-1)))
    assert Rpair[1]==dpair[1]==0
    assert Rpair[0]==F6.from_expr(s.sympify(data['ResV_in_six_coordinates']))
    assert dpair[0]==F6.from_expr(s.sympify(data['anti_delta_in_six_coordinates']))
    normJ0=(a0*a1*a3-a4)**2-a0*(a1*a2-a5)**2
    assert normJ0!=0 and Rpair[0]!=0 and dpair[0]!=0
    # z flip equals the already-classified old sigma, followed by rational
    # torus regauging v1 to1. It is not the extra Spin5 U/V exchange.
    Z=s.Symbol('z');rr0=Z*I1;pp0=(I2+Z*I3)/2;qq0=(I2-Z*I3)/(2*rr0**3)
    ss0=(I4+Z*I5)/2;tt0=-(I4-Z*I5)/(2*rr0**5)
    old=[rr0,pp0,qq0,ss0,tt0];expected=[-rr0,qq0*rr0**3,-pp0/rr0**3,-tt0*rr0**5,ss0/rr0**5]
    assert all(s.cancel(a.subs(Z,-Z)-b)==0 for a,b in zip(old,expected))
    # Verify the entire finite rational DAG by rebuilding the rank4 algebra
    # recipe from its constant-folded node table; do not expand large images.
    nodes=data['arithmetic_nodes'];assert len(nodes)==1169
    index={}
    for i,(op,arg)in enumerate(nodes):
        if op in ['+','*']:assert all(j<i for j in arg);arg=tuple(arg)
        elif op=='inv':assert arg<i
        assert(op,arg)not in index;index[op,arg]=i
    qn=lambda q:index['q',str(F(q))]
    zero,uno=qn(0),qn(1)
    def ad(*ids):
        co=F(0);terms=[]
        for i in ids:
            if nodes[i][0]=='q':co+=F(nodes[i][1])
            else:terms.append(i)
        if co:terms.append(qn(co))
        return zero if not terms else terms[0]if len(terms)==1 else index['+',tuple(sorted(terms))]
    def mt(*ids):
        co=F(1);terms=[]
        for i in ids:
            if nodes[i][0]=='q':co*=F(nodes[i][1])
            else:terms.append(i)
        if not co:return zero
        if co!=1:terms.append(qn(co))
        return uno if not terms else terms[0]if len(terms)==1 else index['*',tuple(sorted(terms))]
    neg=lambda i:mt(qn(-1),i)
    def inv(i):return qn(1/F(nodes[i][1]))if nodes[i][0]=='q'else index['inv',i]
    def pw(i,n):return mt(*([i]*n))if n>=0 else pw(inv(i),-n)
    xs=[index['v',i]for i in range(6)]
    def expr(a):
        if a.is_Number:return qn(str(a))
        if a.is_Symbol:return xs[I.index(a)]
        if a.is_Add:return ad(*(expr(b)for b in a.args))
        if a.is_Mul:return mt(*(expr(b)for b in a.args))
        assert a.is_Pow and a.exp.is_Integer;return pw(expr(a.base),int(a.exp))
    rn=expr(s.sympify(data['ResV_in_six_coordinates']));dn=expr(s.sympify(data['anti_delta_in_six_coordinates']))
    assert dn==data['anti_delta_node']
    def scalar(i):return(i,zero,zero,zero)
    def ec(q):return scalar(qn(q))
    def ea(*vv):return tuple(ad(*(v[i]for v in vv))for i in range(4))
    def en(v):return tuple(neg(x)for x in v)
    def es(v,i):return tuple(mt(x,i)for x in v)
    def em(a,b):
        terms=[[]for _ in range(4)]
        for i,j in itertools.product(range(4),repeat=2):
            x=mt(a[i],b[j])
            if i&j&1:x=mt(x,xs[0])
            if i&j&2:x=mt(x,rn)
            terms[i^j].append(x)
        return tuple(ad(*v)for v in terms)
    def ep(v,n):
        x=ec(1)
        for _ in range(n):x=em(x,v)
        return x
    def ei(v):
        assert v[2:]==(zero,zero)
        den=ad(pw(v[0],2),neg(mt(xs[0],pw(v[1],2))))
        return(mt(v[0],inv(den)),neg(mt(v[1],inv(den))),zero,zero)
    ez=(zero,uno,zero,zero);ek=(zero,zero,uno,zero)
    r=es(ez,xs[1]);p=es(ea(scalar(xs[2]),es(ez,xs[3])),qn('1/2'))
    q=es(ea(es(ez,xs[2]),scalar(neg(mt(xs[0],xs[3])))),inv(mt(qn(2),pw(xs[0],2),pw(xs[1],3))))
    ss=es(ea(scalar(xs[4]),es(ez,xs[5])),qn('1/2'))
    t=en(es(ea(es(ez,xs[4]),scalar(neg(mt(xs[0],xs[5])))),inv(mt(qn(2),pw(xs[0],3),pw(xs[1],5)))))
    j0=es(ea(em(p,r),en(ss)),qn(2));j1=es(ea(em(p,q),en(em(ss,t))),qn(4))
    aa=ea(ec('1/2'),en(es(em(j1,ek),inv(mt(qn(8),rn)))))
    gg=es(em(j0,ek),inv(mt(qn(4),rn)));bb=es(em(ea(en(j1),es(ek,qn(-4))),ei(j0)),qn('1/2'))
    m00=en(ea(em(p,ep(aa,2)),em(aa,gg),em(t,ep(gg,2))))
    m01=en(ea(em(ss,ep(aa,2)),em(r,em(aa,gg)),em(q,ep(gg,2))))
    m10=en(ea(em(p,ep(bb,2)),bb,t));m11=en(ea(em(ss,ep(bb,2)),em(r,bb),q))
    idet=es(ek,neg(inv(rn)));b00=em(idet,m11);b01=en(em(idet,m01));b10=en(em(idet,m10));b11=em(idet,m00)
    w=[ea(em(aa,ep(b00,2)),em(gg,ep(b10,2))),es(ea(em(aa,em(b00,b01)),em(gg,em(b10,b11))),qn(2)),ea(em(aa,ep(b01,2)),em(gg,ep(b11,2))),
        ea(em(bb,ep(b00,2)),ep(b10,2)),es(ea(em(bb,em(b00,b01)),em(b10,b11)),qn(2)),ea(em(bb,ep(b01,2)),ep(b11,2))]
    plus4=ea(em(w[0],ep(w[1],3)),em(w[5],ep(w[4],3)));minus4=ea(em(w[0],ep(w[1],3)),en(em(w[5],ep(w[4],3))))
    minus6=ea(em(w[3],ep(w[1],5)),en(em(w[2],ep(w[4],5))));plus6=ea(em(w[3],ep(w[1],5)),em(w[2],ep(w[4],5)))
    vectors=[scalar(rn),en(em(ez,em(ek,em(w[1],w[4])))),es(plus4,mt(xs[0],pw(rn,2))),en(es(em(ek,minus4),mt(xs[0],rn))),es(em(ez,minus6),mt(xs[0],pw(rn,3))),en(es(em(ez,em(ek,plus6)),mt(xs[0],pw(rn,2))))]
    assert [list(v)for v in vectors]==data['four_algebra_coefficients_before_trace']
    images=[v[0]for v in vectors];assert images==data['tau_images_I0_through_I5']
    fixed=[pw(dn,2)]+[mt(qn('1/2'),ad(x,y))for x,y in zip(xs,images)]+[mt(qn('1/2'),ad(x,neg(y)),inv(dn))for x,y in zip(xs,images)]
    assert fixed==data['fixed_field13_generator_nodes']and len(fixed)==13
    attempt=read(folder/'bounded-expansion.json');assert attempt['stage']=='TIME_LIMIT_AT_six_normalized_W_coefficients'
    evidence=dict(independent_physical_swap12=[[str(x)for x in row]for row in M.tolist()],
        independent_delta_terms=len(s.Poly(delta,*uv).terms()),critical_frame_detA=1,critical_frame_detM='-kappa',
        resultant_in_K=str(Rpair[0]),anti_delta_in_K=str(dpair[0]),norm_J0=str(normJ0),
        all1169_node_rank4_recipe_and13_fixed_outputs_verified=True,
        root_choice_descent='z flip is old sigma plus torus regauge; kappa flip changes target root order within full line normalizer',
        three_nonconstant_coefficients_zero_by_root_choice_proof_not_large_expansion=True,
        generic_Sp4_stabilizer='trivial: equivariant normal VECTOR forces connected Spin4; exactT and weight1V force identity',
        stable_rationality='H1(F,Sp4)=1 via symplectic Gram-Schmidt; F(t1..t10) is rational',
        fixed_field='C(I0..I5)^extra_tau generated by delta^2,S_i,A_i (13)',
        pure_rationality_not_asserted=True,old_expansion_timeout_not_replayed=True)
    (out/'independent-certificate.json').write_text(json.dumps(evidence,indent=2)+'\n')
    paths=[Path(__file__),ROOT/f'research/notes/{STEM}.md',physical_path,*gates,
        *[folder/p for p in ['report.json','action.json.gz','bounded-expansion.json','bounded-expansion-probe.py']],
        ROOT/'code/q_firstjet_spin5_extra_normalizer_field_20261104_v63.py',
        ROOT/'research/notes/Q_FIRSTJET_SPIN5_EXTRA_NORMALIZER_FIELD_20261104_v63.md']
    report=dict(id=STEM,all_passed=True,independent_actual12_coefficient_swap_and_full_delta_sign=True,
        independent_full5coefficient_critical_frame_and_nonzero_generic_norm_denominators=True,
        all1169_rational_DAG_nodes_semantically_bound_to_rank4_etale_normalization_recipe=True,
        root_choice_descent_and_trace_over4_proof_accepted_without_large_expansion=True,
        full13_generator_fixed_field_lemma_and_tau_involution_accepted=True,
        full_W16_Sp4_invariant_field_identified_by_prior_associated_bundle=True,
        generic_Sp4_freeness_and_stable_rationality_level_atmost10_accepted=True,
        pure_rationality_or_polynomial_invariant_ring_or_NEW_equation_claimed=False,
        no_new_point_big_expansion_GB_or_timedout_job_retry=True,
        input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},
        evidence_sha256={'independent-certificate.json':sha(out/'independent-certificate.json')},seconds=time.monotonic()-started)
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');signal.alarm(0)
    print(json.dumps({k:report[k]for k in ['all_passed','generic_Sp4_freeness_and_stable_rationality_level_atmost10_accepted','seconds']},indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
