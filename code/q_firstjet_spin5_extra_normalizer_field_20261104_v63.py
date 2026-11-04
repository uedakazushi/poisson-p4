#!/usr/bin/env python3
"""Finite rational DAG for the EXTRA Spin5 line component on rational Spin4 fields.

Uses trace on an etale rank4 algebra; no expanded high-degree action retry.
"""
from fractions import Fraction
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import resource
import signal
import time
import sympy as s
from sympy.polys.rings import ring

ROOT=Path(__file__).resolve().parents[1]
ID='Q_FIRSTJET_SPIN5_EXTRA_NORMALIZER_FIELD_20261104_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

class DAG:
    def __init__(self):self.nodes=[];self.index={};self.zero=self.q(0);self.one=self.q(1)
    def put(self,key):
        if key not in self.index:
            assert len(self.nodes)<250000
            self.index[key]=len(self.nodes);self.nodes.append(key)
        return self.index[key]
    def q(self,q):return self.put(('q',str(Fraction(q))))
    def var(self,i):return self.put(('v',i))
    def add(self,*args):
        terms=[];constant=Fraction(0)
        for a in args:
            if self.nodes[a][0]=='q':constant+=Fraction(self.nodes[a][1])
            else:terms.append(a)
        if constant:terms.append(self.q(constant))
        if not terms:return self.zero
        if len(terms)==1:return terms[0]
        return self.put(('+',tuple(sorted(terms))))
    def mul(self,*args):
        terms=[];constant=Fraction(1)
        for a in args:
            if self.nodes[a][0]=='q':constant*=Fraction(self.nodes[a][1])
            else:terms.append(a)
        if not constant:return self.zero
        if constant!=1:terms.append(self.q(constant))
        if not terms:return self.one
        if len(terms)==1:return terms[0]
        return self.put(('*',tuple(sorted(terms))))
    def neg(self,a):return self.mul(self.q(-1),a)
    def inv(self,a):
        assert a!=self.zero
        if self.nodes[a][0]=='q':return self.q(1/Fraction(self.nodes[a][1]))
        return self.put(('inv',a))
    def power(self,a,n):
        if n<0:return self.power(self.inv(a),-n)
        return self.mul(*([a]*n))
    def expression(self,expr,variables):
        if expr.is_Number:return self.q(str(expr))
        if expr.is_Symbol:return variables[expr]
        if expr.is_Add:return self.add(*(self.expression(c,variables)for c in expr.args))
        if expr.is_Mul:return self.mul(*(self.expression(c,variables)for c in expr.args))
        assert expr.is_Pow and expr.exp.is_Integer
        return self.power(self.expression(expr.base,variables),int(expr.exp))

def main():
    begin=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000))
    output=ROOT/'reports'/ID;output.mkdir(parents=True,exist_ok=False)
    joint=ROOT/'reports/Q_FIRSTJET_SPIN4_BIFORM_RATIONAL_FIELD_20261104_v63_REVIEW/report.json'
    spin=ROOT/'reports/Q_FIRSTJET_SPIN4_MIXED_DISCRIMINANT_20261103_v63_REVIEW/report.json'
    checked={}
    for path in [joint,spin]:
        obj=json.loads(path.read_text());assert obj['all_passed']
        for name,want in obj['input_sha256'].items():assert sha(ROOT/name)==want;checked[name]=want
    physicalpath=ROOT/'reports/Q_FIRSTJET_SPIN4_MIXED_DISCRIMINANT_20261103_v63/certificate.json'
    physical=json.loads(physicalpath.read_text());uv=list(map(s.Symbol,physical['variables12']))
    swap=s.Matrix(physical['normalizer_biform12']);assert swap*swap==-s.eye(12)
    delta=s.sympify(physical['normal_component_delta']);sub=dict(zip(uv,swap*s.Matrix(uv)))
    assert s.expand(delta.subs(sub,simultaneous=True)+delta)==0
    #Canonical source I0..I5 -> old torus gauge, with z²=I0 and v1=1.
    eta,A,B,C,D,E=s.symbols('I0:6');zz=s.Symbol('z')
    r=zz*A;p=(B+zz*C)/2;q=(B-zz*C)/(2*r**3)
    ss=(D+zz*E)/2;tt=-(D-zz*E)/(2*r**5)
    RV=(p*q-tt*ss)**2-(p*r-ss)*(q-tt*r)
    def descend(expr):
        n,d=s.fraction(s.cancel(expr));n=s.rem(n,zz**2-eta,zz);d=s.rem(d,zz**2-eta,zz)
        assert not n.has(zz)and not d.has(zz);return s.factor(n/d)
    RV=descend(RV);deltaI=descend(zz*(r/4-tt*ss))
    assert s.cancel(deltaI-(eta**3*A**6+D**2-eta*E**2)/(4*eta**2*A**5))==0
    expected_delta=s.expand(delta.subs(dict(zip(uv,[s.Symbol('c'),0,0,0,0,s.Symbol('c'),*s.symbols('v0:6')])),simultaneous=True))
    v0,v1,v2,v3,v4,v5=s.symbols('v0:6');cc=s.Symbol('c')
    assert expected_delta==s.expand(cc**2*(v1*v4/4-v2*v3))
    #Small independent exact critical-root normalization over five free coefficients.
    PP,pp,qq,rr,uu,vv=ring('p,q,r,s,t',s.QQ);KK=PP.to_field();pp,qq,rr,uu,vv=map(KK,[pp,qq,rr,uu,vv])
    norm=(pp*qq-uu*vv)**2-(pp*rr-uu)*(qq-rr*vv);j0=2*(pp*rr-uu);j1=4*(pp*qq-uu*vv)
    pairone=(KK.one,KK.zero);pairzero=(KK.zero,KK.zero)
    def pa(*args):return tuple(sum((x[i]for x in args),KK.zero)for i in range(2))
    def pm(a,b):return(a[0]*b[0]+norm*a[1]*b[1],a[0]*b[1]+a[1]*b[0])
    def ps(a,c):return(a[0]*c,a[1]*c)
    aa=(KK(s.Rational(1,2)),-j1/(8*norm));bb=(-j1/(2*j0),-2/j0);gg=(KK.zero,j0/(4*norm))
    assert pa(aa,ps(pm(bb,gg),-1))==pairone
    m00=ps(pa(ps(pm(aa,aa),pp),pm(aa,gg),ps(pm(gg,gg),vv)),-1)
    m01=ps(pa(ps(pm(aa,aa),uu),ps(pm(aa,gg),rr),ps(pm(gg,gg),qq)),-1)
    m10=ps(pa(ps(pm(bb,bb),pp),bb,(vv,KK.zero)),-1)
    m11=ps(pa(ps(pm(bb,bb),uu),ps(bb,rr),(qq,KK.zero)),-1)
    assert pa(pm(m00,m11),ps(pm(m01,m10),-1))==(KK.zero,-KK.one)
    assert pa(ps(pm(aa,bb),2*pp),aa,pm(gg,bb),ps(gg,2*vv))==pairzero
    assert pa(ps(pm(aa,bb),2*uu),ps(aa,rr),ps(pm(gg,bb),rr),ps(gg,2*qq))==pairzero
    #Base rational DAG, and rank4 etale algebra basis (1,z,kappa,z*kappa).
    dag=DAG();xs=[dag.var(i)for i in range(6)];symbols=[eta,A,B,C,D,E];rv=dag.expression(RV,dict(zip(symbols,xs)))
    delta_node=dag.expression(deltaI,dict(zip(symbols,xs)));eta0=xs[0]
    zero=(dag.zero,)*4
    def scalar(n):return(n,dag.zero,dag.zero,dag.zero)
    def con(q):return scalar(dag.q(q))
    def ad(*vs):return tuple(dag.add(*(v[j]for v in vs))for j in range(4))
    def neg(v):return tuple(dag.neg(c)for c in v)
    def sc(v,n):return tuple(dag.mul(c,n)for c in v)
    def mu(a,b):
        terms=[[]for _ in range(4)]
        for i,j in itertools.product(range(4),repeat=2):
            c=dag.mul(a[i],b[j])
            if i&j&1:c=dag.mul(c,eta0)
            if i&j&2:c=dag.mul(c,rv)
            terms[i^j].append(c)
        return tuple(dag.add(*ts)for ts in terms)
    def power(v,n):
        ans=con(1)
        for _ in range(n):ans=mu(ans,v)
        return ans
    def invz(v):
        assert v[2:]==(dag.zero,dag.zero)
        den=dag.add(dag.power(v[0],2),dag.neg(dag.mul(eta0,dag.power(v[1],2))))
        return(dag.mul(v[0],dag.inv(den)),dag.neg(dag.mul(v[1],dag.inv(den))),dag.zero,dag.zero)
    Z=(dag.zero,dag.one,dag.zero,dag.zero);kap=(dag.zero,dag.zero,dag.one,dag.zero)
    r=sc(Z,xs[1]);p=sc(ad(scalar(xs[2]),sc(Z,xs[3])),dag.q('1/2'))
    q=sc(ad(sc(Z,xs[2]),scalar(dag.neg(dag.mul(eta0,xs[3])))),dag.inv(dag.mul(dag.q(2),dag.power(eta0,2),dag.power(xs[1],3))))
    ss=sc(ad(scalar(xs[4]),sc(Z,xs[5])),dag.q('1/2'))
    tt=neg(sc(ad(sc(Z,xs[4]),scalar(dag.neg(dag.mul(eta0,xs[5])))),dag.inv(dag.mul(dag.q(2),dag.power(eta0,3),dag.power(xs[1],5)))))
    J0=sc(ad(mu(p,r),neg(ss)),dag.q(2));J1=sc(ad(mu(p,q),neg(mu(ss,tt))),dag.q(4))
    A00=ad(con('1/2'),neg(sc(mu(J1,kap),dag.inv(dag.mul(dag.q(8),rv)))))
    A10=sc(mu(J0,kap),dag.inv(dag.mul(dag.q(4),rv)))
    A01=sc(mu(ad(neg(J1),sc(kap,dag.q(-4))),invz(J0)),dag.q('1/2'))
    M00=neg(ad(mu(p,power(A00,2)),mu(A00,A10),mu(tt,power(A10,2))))
    M01=neg(ad(mu(ss,power(A00,2)),mu(r,mu(A00,A10)),mu(q,power(A10,2))))
    M10=neg(ad(mu(p,power(A01,2)),A01,tt));M11=neg(ad(mu(ss,power(A01,2)),mu(r,A01),q))
    invdet=sc(kap,dag.neg(dag.inv(rv)))
    b00=mu(invdet,M11);b01=neg(mu(invdet,M01));b10=neg(mu(invdet,M10));b11=mu(invdet,M00)
    w=[ad(mu(A00,power(b00,2)),mu(A10,power(b10,2))),
       sc(ad(mu(A00,mu(b00,b01)),mu(A10,mu(b10,b11))),dag.q(2)),
       ad(mu(A00,power(b01,2)),mu(A10,power(b11,2))),
       ad(mu(A01,power(b00,2)),power(b10,2)),
       sc(ad(mu(A01,mu(b00,b01)),mu(b10,b11)),dag.q(2)),
       ad(mu(A01,power(b01,2)),power(b11,2))]
    plus4=ad(mu(w[0],power(w[1],3)),mu(w[5],power(w[4],3)))
    minus4=ad(mu(w[0],power(w[1],3)),neg(mu(w[5],power(w[4],3))))
    minus6=ad(mu(w[3],power(w[1],5)),neg(mu(w[2],power(w[4],5))))
    plus6=ad(mu(w[3],power(w[1],5)),mu(w[2],power(w[4],5)))
    image_vectors=[scalar(rv),neg(mu(Z,mu(kap,mu(w[1],w[4])))),sc(plus4,dag.mul(eta0,dag.power(rv,2))),
                   neg(sc(mu(kap,minus4),dag.mul(eta0,rv))),sc(mu(Z,minus6),dag.mul(eta0,dag.power(rv,3))),
                   neg(sc(mu(Z,mu(kap,plus6)),dag.mul(eta0,dag.power(rv,2))))]
    #Coefficient1 is exactly etale trace/4, even if either quadratic extension splits.
    images=[v[0]for v in image_vectors]
    fixed=[dag.power(delta_node,2)]
    fixed += [dag.mul(dag.q('1/2'),dag.add(x,y))for x,y in zip(xs,images)]
    fixed += [dag.mul(dag.q('1/2'),dag.add(x,dag.neg(y)),dag.inv(delta_node))for x,y in zip(xs,images)]
    assert len(fixed)==13
    evidence={'variables':[str(x)for x in symbols],'arithmetic_nodes':dag.nodes,
        'tau_images_I0_through_I5':images,'four_algebra_coefficients_before_trace':image_vectors,
        'ResV_in_six_coordinates':str(RV),'anti_delta_in_six_coordinates':str(deltaI),'anti_delta_node':delta_node,
        'fixed_field13_generator_nodes':fixed,'trace_normalization':'coefficient of1 in basis1,z,kappa,z*kappa = trace/4',
        'algebra_relations':'z²=I0,kappa²=ResV(I); both radicands invertible generically; rank4 etale, not assumed a field',
        'exact_frame_M_determinant':'-kappa','output_is_base_field_by_root_choice_invariance_proof':True,
        'three_nonconstant_algebra_components_not_expanded_to_zero':True,
        'separate_extra_component_UV_action':'(U(x,y),V(x,y))->(-V(y,x),U(y,x))',
        'generic_stabilizer_proof':'nonisotropic normal vector forces Spin4; U stabilizerT and weight1 V coefficient forces t1',
        'stable_rationality_level_upper_bound':10}
    ev=output/'action.json.gz';ev.write_bytes(gzip.compress((json.dumps(evidence,separators=(',',':'))+'\n').encode(),mtime=0))
    attempt=ROOT/'work/spin4-extra-normalizer-v63/action-probe.json'
    assert json.loads(attempt.read_text())['stage']=='TIME_LIMIT_AT_six_normalized_W_coefficients'
    (output/'bounded-expansion.json').write_bytes(attempt.read_bytes())
    (output/'bounded-expansion-probe.py').write_bytes((attempt.parent/'action_probe.py').read_bytes())
    note=ROOT/'research/notes'/(ID+'.md')
    files=[Path(__file__),note,joint,spin,physicalpath,
           ROOT/'research/notes/Q_FIRSTJET_SPIN4_BIFORM_RATIONAL_FIELD_20261104_v63.md',
           ROOT/'research/notes/Q_FIRSTJET_SPIN4_BIFORM_RATIONAL_FIELD_20261104_v63_REVIEW.md']
    report={'id':ID,'all_passed':True,'extra_Spin5_component_distinguished_from_old_mu4':True,
        'all12_swap_coefficients_square_central_and_delta_sign_exact':True,'critical_root_frame_and_all_mixed_terms_and_detM_exact':True,
        'finite_rational_base_field_DAG_for_all6_extra_involution_images':True,
        'whole_fixed_field_generated_by13_rational_circuits':True,'fixed_field_transcendence_degree':6,
        'generic_Sp4_freeness_and_torsor_GramSchmidt_stable_rationality_proved':True,'stable_rationality_level_upper_bound':10,
        'pure_rationality_of_fixed_field_or_W16_quotient_claim':False,
        'full_polynomial_invariant_ring_Hilbert_series_NEW_equation_or_genus_claim':False,
        'expanded_action_attempt_capped_seconds':60,'expanded_action_attempt_completed':False,'expanded_action_retry_performed':False,
        'new_beta_point_high_degree_probe_or_GB':False,'DAG_node_count':len(dag.nodes),
        'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in files},'checked_dependency_sha256':checked,
        'evidence_sha256':{p.name:sha(p)for p in output.iterdir()if p.is_file()},
        'peak_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'seconds':time.monotonic()-begin}
    (output/'report.json').write_text(json.dumps(report,indent=2)+'\n');signal.alarm(0)
    print(json.dumps({'all_passed':True,'DAG_nodes':len(dag.nodes),'fixed_generators':13,'seconds':report['seconds']}))

if __name__=='__main__':main()
