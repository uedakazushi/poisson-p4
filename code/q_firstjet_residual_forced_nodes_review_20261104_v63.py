#!/usr/bin/env python3
"""Own sparse universal arithmetic and conceptual two-node audit; no CAS job."""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import resource
import signal
import time

ROOT=Path(__file__).resolve().parents[1]
BASE='Q_FIRSTJET_RESIDUAL_QUARTIC_FORCED_NODES_20261104_v63'
STEM=BASE+'_REVIEW'
NAMES=['c','d','f','u','v','r','t','D','a','x','y','k','B']
ZERO=(0,)*len(NAMES)


class Poly:
    def __init__(self, terms):
        if isinstance(terms,(int,F)):
            terms={ZERO:F(terms)}
        self.terms={e:F(c) for e,c in terms.items() if c}
    def __add__(self,other):
        other=other if isinstance(other,Poly) else Poly(other)
        out=dict(self.terms)
        for e,c in other.terms.items():out[e]=out.get(e,F(0))+c
        return Poly(out)
    __radd__=__add__
    def __neg__(self):return Poly({e:-c for e,c in self.terms.items()})
    def __sub__(self,other):return self+-coerce(other)
    def __rsub__(self,other):return coerce(other)+-self
    def __mul__(self,other):
        other=coerce(other);out={}
        for e,c in self.terms.items():
            for f,d in other.terms.items():
                g=tuple(a+b for a,b in zip(e,f))
                out[g]=out.get(g,F(0))+c*d
        return Poly(out)
    __rmul__=__mul__
    def __pow__(self,n):
        out=Poly(1)
        for _ in range(n):out=out*self
        return out
    def coeff(self,name,k):
        i=NAMES.index(name);out={}
        for e,c in self.terms.items():
            if e[i]==k:
                f=list(e);f[i]=0;out[tuple(f)]=c
        return Poly(out)
    def subst(self,name,value):
        i=NAMES.index(name);out=Poly(0)
        for e,c in self.terms.items():
            f=list(e);f[i]=0
            out+=Poly({tuple(f):c})*value**e[i]
        return out
    def __eq__(self,other):return self.terms==coerce(other).terms


def coerce(p):return p if isinstance(p,Poly) else Poly(p)
def var(name):
    e=list(ZERO);e[NAMES.index(name)]=1;return Poly({tuple(e):1})
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())


def main(out):
    start=time.monotonic()
    resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60)
    primary=ROOT/'reports'/BASE
    report=read(primary/'report.json')
    assert report['all_passed']
    assert report['generic_two_distinct_ordinary_nodes_on_m0_proved']
    assert report['all_geometric_quotient_components_normalization_genus_upper']==1
    assert not report['only_these_nodes_irreducibility_or_exact_genus_claimed']
    quotient=ROOT/'reports/Q_FIRSTJET_RESIDUAL_RECIPROCAL_QUOTIENT_20261104_v63'
    qreport=read(quotient/'report.json');qcert=read(quotient/'certificate.json')
    assert qreport['all_passed'] and qreport['actual_plane_total_degree_r_m']==4
    groups=qcert['actual_primitive_plane_coefficients_r_m']
    assert len(groups)==15 and sum(map(len,groups.values()))==1592
    assert max(sum(map(int,k.split(','))) for k in groups)==4
    assert groups['4,0']==[[[4,8,0,0,0],'16']]
    failed=ROOT/'reports/Q_FIRSTJET_RESIDUAL_QUARTIC_SMOOTH_SPECIALIZATION_20261104_v63'
    old=read(failed/'report.json')
    assert old['all_passed'] is False and old['run_status']=='completed'
    oldlog=(failed/'M2.log').read_text()
    assert 'PROJECTIVE_SINGULAR_STOP' in oldlog
    assert 'PROJECTIVE_SMOOTH_PASS' not in oldlog
    checked={}
    for path in [primary/'report.json',quotient/'report.json',failed/'report.json']:
        data=read(path)
        for name,want in data.get('input_sha256',{}).items():
            target=ROOT/name;assert sha(target)==want,str(target)
            checked[str(target.relative_to(ROOT))]=want
        for name,want in data.get('evidence_sha256',{}).items():
            target=path.parent/name;assert sha(target)==want,str(target)
            checked[str(target.relative_to(ROOT))]=want

    c,d,f,u,v,r,t,D,a,x,y,k,B=[var(n) for n in NAMES]
    # Reconstruct the universal identity from sparse multiplication, not the producer.
    qR=r*r+4*f*u*v;qT=t*t+4*c*u*v
    W=2*u*v*D+(u+v)*r*t
    psi=r*r*t*t+(u+v)*r*t*D+u*v*D*D-(u-v)**2*(c*r*r+f*t*t+4*c*f*u*v)
    assert W*W-(u-v)**2*qR*qT-4*u*v*psi==0
    D0=B+2*a*r+2*d*t
    P0=r*t+u*D0
    assert psi.subst('v',u).subst('D',D0)==P0*P0
    assert W.subst('v',u).subst('D',D0)==2*u*P0
    # Clear only y or x in the linear substitution, keeping K10 as independent k.
    A=-x
    Br=k+2*u*a*y-2*u*d*x
    Cr=u*B*y+2*u*d*k
    yP0=A*r*r+Br*r+Cr
    assert yP0.coeff('r',2)==-x
    disc=Br*Br-4*A*Cr
    assert disc.coeff('k',2)==1
    assert yP0.coeff('k',1)==r+2*u*d
    xP0=-y*t*t+(k+2*u*d*x-2*u*a*y)*t+u*B*x+2*u*a*k
    assert xP0.coeff('k',1)==t+2*u*a
    # Actual m=0 coefficient identities, after clearing the single y^2 denominator.
    ell=(a+x)**2-y*y
    dnum=-(a*x+x*x+y*y)
    assert dnum*dnum-x*x*ell==y*y*(2*a*x+3*x*x+y*y)
    assert a*a-ell==-2*a*x-x*x+y*y
    assert (r*r+4*f*u*u).subst('r',-2*u*d)==4*u*u*(d*d+f)
    assert (t*t+4*c*u*u).subst('t',-2*u*a)==4*u*u*(a*a+c)

    proof={
        'own_sparse_universal_W_identity_and_m0_square':True,
        'P0_r_leading':'-x/y',
        'P0_discriminant_K10_squared_coefficient':'1/y^2',
        'P0_K10_coefficient':'(r+2*u0*d)/y',
        'symmetric_P0_K10_coefficient':'(t+2*u0*a)/x',
        'qR_obstruction':'d^2+f0=2*a*x+3*x^2+y^2 !=0',
        'qT_obstruction':'a^2+c0=-2*a*x-x^2+y^2 !=0',
        'ordinary_node_argument':'dW/dr=2*u0*P0prime !=0; quadratic tangent W^2-qR0*qT0*m^2 nondegenerate',
        'normalization_genus_bound_for_geometric_plane_components':1,
        'exact_genus_or_irreducibility_or_only_two_singularities':False,
        'original_quadratic_cover_genus_claimed':False,
        'failed_specialization_retested_or_promoted_to_generic':False,
    }
    out.mkdir(parents=True,exist_ok=False)
    evidence=out/'independent-proof.json'
    evidence.write_text(json.dumps(proof,indent=2)+'\n')
    paths=[Path(__file__),ROOT/'research/notes'/(STEM+'.md'),primary/'report.json',
           ROOT/'research/notes'/(BASE+'.md'),quotient/'report.json',quotient/'certificate.json',
           failed/'report.json',failed/'M2.log']
    for p in paths:checked[str(p.relative_to(ROOT))]=sha(p)
    result={
        'id':STEM,'all_passed':True,
        'own_sparse_arithmetic_W_identity_square_discriminant_and_unit_obstructions':True,
        'generic_two_distinct_ordinary_nodes_over_geometric_coefficient_field':True,
        'all_geometric_plane_quotient_components_normalization_genus_upper':1,
        'transcendental_K10_root_unit_argument_reviewed':True,
        'exact_genus1_irreducibility_or_no_additional_singularities_claimed':False,
        'original_quadratic_cover_or_full_residual_curve_genus_claimed':False,
        'old_failed_smoothness_probe_preserved_without_rerun_or_generic_negative_inference':True,
        'new_CAS_GB_point_prime_source_or_invariant_job':False,
        'input_sha256':checked,'evidence_sha256':{'independent-proof.json':sha(evidence)},
        'seconds':time.monotonic()-start,'frozen':True,
    }
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ['all_passed','seconds']}))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True)
    main(p.parse_args().output)
