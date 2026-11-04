#!/usr/bin/env python3
"""Independent fixed pair/monic/norm/finite112 exact arithmetic audit.

Integer numerators with one common denominator per Laurent polynomial avoid
per-term Fraction products. No author executable, GB, point or prime is run.
"""
import argparse,hashlib,json,math,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
import sympy as sp
from sympy.polys.rings import ring

ROOT=Path(__file__).resolve().parents[1]
PREFIX='S61_NONSQUARE_RAW_H_B_'
ID=PREFIX+'PAIR_NORM_FINITE112_20261104_v63_REVIEW'
MAXTERMS=0
class P:
    def __init__(self,n=None,d=1):
        global MAXTERMS
        n={tuple(m):int(c) for m,c in (n or {}).items() if c};d=int(d);assert d
        if d<0:n={m:-c for m,c in n.items()};d=-d
        if not n:d=1
        elif d!=1:
            g=math.gcd(d,*(abs(c) for c in n.values()))
            if g>1:n={m:c//g for m,c in n.items()};d//=g
        self.n=n;self.d=d;MAXTERMS=max(MAXTERMS,len(n));assert len(n)<=150000
    def __eq__(self,b):return isinstance(b,P) and self.n==b.n and self.d==b.d
    def __neg__(self):return P({m:-c for m,c in self.n.items()},self.d)
    def __add__(self,b):
        d=math.lcm(self.d,b.d);a={m:c*(d//self.d) for m,c in self.n.items()}
        for m,c in b.n.items():a[m]=a.get(m,0)+c*(d//b.d)
        return P(a,d)
    def __sub__(self,b):return self+-b
    def __mul__(self,b):
        a={}
        for m,c in self.n.items():
            for n,e in b.n.items():
                k=(m[0]+n[0],m[1]+n[1],m[2]+n[2]);a[k]=a.get(k,0)+c*e
        return P(a,self.d*b.d)
    def scale(self,q):
        q=Q(q);return P({m:c*q.numerator for m,c in self.n.items()},self.d*q.denominator)
    def shift(self,ex):return P({tuple(a+b for a,b in zip(m,ex)):c for m,c in self.n.items()},self.d)
    def power(self,n):
        v=ONE
        for _ in range(n):v=v*self
        return v
    def coefficients(self):return {m:Q(c,self.d) for m,c in self.n.items()}
    def degree(self,k):return max((m[k] for m in self.n),default=-1)
    def lc(self,k):
        degree=self.degree(k);out={}
        for m,c in self.n.items():
            if m[k]==degree:
                n=list(m);n[k]=0;out[tuple(n)]=c
        return P(out,self.d)
ZERO=P();ONE=P({(0,0,0):1})
def decode(a):
    q={tuple(m):Q(c) for m,c in a};assert len(q)==len(a) and all(c and len(m)==3 and m[1]>=0 and m[2]>=0 for m,c in q.items())
    d=math.lcm(*(c.denominator for c in q.values())) if q else 1
    return P({m:c.numerator*(d//c.denominator) for m,c in q.items()},d)
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def total(a):
    out=ZERO
    for p in a:out=out+p
    return out
def det3(A):
    return total(A[0][j]*(A[1][(j+1)%3]*A[2][(j+2)%3]-A[1][(j+2)%3]*A[2][(j+1)%3]) for j in range(3))

def main(out):
    start=time.monotonic();signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
    out=out.resolve();out.mkdir(parents=True,exist_ok=False)
    names=['COPRIME_PAIR_PREFLIGHT','PAIR_MONIC_COMPLETION','PAIR_NORM_INPUT','PAIR_NORM','FINITE_CURVE_COVER']
    folders={n:ROOT/'reports'/f'{PREFIX}{n}_20261104_v63' for n in names}
    hashes={}
    for base in folders.values():
        r=read(base/'report.json')
        for p,h in r['input_sha256'].items():assert sha(ROOT/p)==h;hashes[p]=h
        for p,h in r['evidence_sha256'].items():assert sha(base/p)==h;hashes[str((base/p).relative_to(ROOT))]=h
        hashes[str((base/'report.json').relative_to(ROOT))]=sha(base/'report.json')
    gate=ROOT/'reports/S61_NONSQUARE_RAW_H_B_MINORS_GCD_20261104_v63_REVIEW/report.json'
    assert read(gate)['all_passed'];hashes[str(gate.relative_to(ROOT))]=sha(gate)
    source=read(ROOT/'reports/S61_NONSQUARE_RAW_H_B_35_MINORS_20261104_v63/certificate.json')
    pair=read(folders['COPRIME_PAIR_PREFLIGHT']/'certificate.json')
    monic=read(folders['PAIR_MONIC_COMPLETION']/'certificate.json')
    inp=read(folders['PAIR_NORM_INPUT']/'certificate.json')
    norm=read(folders['PAIR_NORM']/'native-certificate.json')
    cover=read(folders['FINITE_CURVE_COVER']/'certificate.json')
    gens=source['localized_ideal_generators_monic_QQ_after_namedunit_and_QQ_only_scaling']
    assert pair['selected_pair']==pair['selected_pair_source']['indices']==monic['source_pair_indices']==[10,15]
    f,g=map(decode,pair['selected_pair_source']['polynomials'])
    assert f==decode(gens[10]) and g==decode(gens[15])
    assert [len(f.n),len(g.n)]==[68,59] and f.degree(2)==g.degree(2)==4
    rr,aa,ee,dd=ring('a,e,d',sp.QQ)
    def rp(p):return rr.from_dict({m:sp.QQ(c,p.d) for m,c in p.n.items()})
    assert rp(f).gcd(rp(g))==rr.one
    # Preserve the old metadata's TOP d^4 piece, and separately verify the
    # new TRUE leading coefficients with that d^4 exponent removed.
    for p,profile in zip((f,g),pair['selected_pair_source']['profiles']):
        top=P({m:c for m,c in p.n.items() if m[2]==4},p.d)
        assert top==decode(profile['d_leading_coefficient'])
    lf,lg=f.lc(2),g.lc(2)
    assert [lf,lg]==list(map(decode,monic['actual_d_leading_coefficients']))
    l0,l1=lf.shift((-3,0,0)),lg.shift((-4,0,0))
    assert [l0,l1]==list(map(decode,monic['normalized_LC_by_known_a3_a4']))
    a=ONE.shift((1,0,0));e=ONE.shift((0,1,0));d=ONE.shift((0,0,1))
    q=a+e.scale(Q(13,50))+ONE.scale(Q(-13,75));assert q==decode(monic['q'])
    assert l0-q*l1==a.scale(Q(-11,9))
    H=decode(monic['monic_H_laurent']);alpha=decode(monic['H_membership_f_cofactor']);beta=decode(monic['H_membership_g_cofactor'])
    assert alpha==ONE.shift((-4,0,0)).scale(Q(-9,11))
    assert beta==q.shift((-5,0,0)).scale(Q(9,11))
    assert H==alpha*f+beta*g and H.shift((5,0,0)).scale(-11)==(a*f-q*g).scale(9)
    assert H.degree(2)==4 and H.lc(2)==ONE and min(m[0] for m in H.n)==-5
    assert len(H.n)==74 and monic['unknown_localizers_inverted']==[] and monic['known_localizer_inverted']==['a']
    coeff=list(map(decode,monic['monic_H_coefficients_by_d_degree']))
    assert total(p.shift((0,0,j)) for j,p in enumerate(coeff))==H
    assert H==decode(inp['H']) and g==decode(inp['g'])
    M=[[decode(p) for p in row] for row in inp['matrix_rows']]
    quotient=list(map(decode,inp['column_quotients_by_H']))
    assert len(M)==4 and all(len(row)==4 for row in M) and len(quotient)==4
    assert sum(len(p.n) for row in M for p in row)==1284
    for j in range(4):
        remainder=total(M[i][j].shift((0,0,i)) for i in range(4))
        assert all(p.degree(2)==0 for p in (M[i][j] for i in range(4)))
        assert g.shift((0,0,j))==H*quotient[j]+remainder
    bounds=inp['determinant_cost_bounds_NOT_value']
    D=Q(bounds['matrix_uniform_QQ_denominator']);assert D==2035768318176269531250000
    assert bounds['matrix_uniform_a_clear_power']==10
    s=ONE.shift((10,0,0)).scale(D)
    N=[[s*p for p in row] for row in M]
    assert all(p.d==1 and all(min(m)>=0 and m[2]==0 for m in p.n) for row in N for p in row)
    # Independently compute all four cofactors and the norm determinant by
    # first-row Laplace expansion. No determinant-degree guess or Nv-only test.
    adj=[]
    for j in range(4):
        small=[[N[i][k] for k in range(4) if k!=j] for i in range(1,4)]
        p=det3(small)
        adj.append(p if j%2==0 else -p)
    assert adj==list(map(decode,norm['adjugate_first_column']))
    R=total(N[0][j]*adj[j] for j in range(4))
    assert R==decode(norm['det_cleared_matrix_R']) and len(R.n)==361 and R.d==1
    for i in range(4):assert total(N[i][j]*adj[j] for j in range(4))==(R if i==0 else ZERO)
    V=total(adj[j].shift((0,0,j)) for j in range(4))
    Qsum=total(adj[j]*quotient[j] for j in range(4))
    assert V==decode(norm['monic_V']) and Qsum==decode(norm['monic_Qsum'])
    assert s*g*V-s*H*Qsum==R
    uf=-s*Qsum*alpha;vg=s*V-s*Qsum*beta
    assert uf==decode(norm['literal_membership_f_cofactor']) and vg==decode(norm['literal_membership_g_cofactor'])
    assert uf*f+vg*g==R
    content=Q(norm['R_QQ_content']);primitive=decode(norm['R_primitive'])
    assert primitive.scale(content)==R
    assert content==math.gcd(*(abs(c) for c in R.n.values()))
    U=decode(source['FULLU_product']);C=a.scale(100)+e.power(2).scale(18)+e.scale(-9)+ONE.scale(-2)
    current=primitive;removed=[]
    assert len(norm['localizer_removal_steps'])==31
    for i,step in enumerate(norm['localizer_removal_steps']):
        factor=decode(step['removed_gcd_with_U']);after=decode(step['remaining_after']);uq=decode(step['U_quotient'])
        assert factor*after==current and factor*uq==U
        assert factor==(a*C).scale(Q(1,100)) if i<4 else factor==a
        removed.append(factor);current=after
    h=decode(norm['remaining_h_coprime_U']);assert current==h
    assert decode(norm['final_gcd_h_U'])==ONE
    assert rp(U).gcd(rp(h))==rr.one
    assert len(h.n)==225 and h.degree(0)==14 and h.degree(1)==28 and h.degree(2)==0
    hprim=decode(cover['h_primitive_integer'])
    assert h.scale(Q(1,10**8))==hprim and hprim.d==1
    assert R==hprim*a.power(31)*C.power(4).scale(content)
    assert Q(cover['R_QQ_content'])==content and Q(cover['h_QQ_content'])==10**8
    E=decode(cover['monic_E_e28']);assert E==h.scale(Q(cover['E_from_h_QQ_unit']))
    assert E.degree(1)==28 and E.lc(1)==ONE and len(E.n)==225
    assert decode(cover['monic_H_d4'])==H
    assert cover['basis_e_i_d_j']==[[i,j] for j in range(4) for i in range(28)]
    assert cover['FULLU_inverse_NOT_adjoined_to_cover_yet'] is True
    proof=dict(selected22_generators10_15_source_bound_and_coprime_QQ=True,
      old_TOP_d4_metadata_and_TRUE_LC_separately_verified=True,
      full_monic_H_identity74terms_min_a_minus5_unknown_inverses_none=True,
      all4_monic_division_identities_and16_norm_entries1284terms=True,
      uniform_clear_D_a10_produces_integer_matrix=True,
      all4_cofactors_independently_computed_and_Laplace_detR361terms=True,
      all4_Nv_equals_R_e0_and_full_V_Q_uf_vg_memberships=True,
      all31_removed_factors_unit_on_FULLU_by_literal_U_quotients=True,
      R_factorization_content_a31_C4_hprim_without_unknown_factorization=True,
      h225terms_degree_e28_a14_and_exact_gcd_h_U1=True,
      E_monic_e28_and_H_monic_d4_give_finite_free112_cover=True,
      all_FULLU_B_bad_points_contained_not_claimed_nonempty=True,
      no_pair_quotient_freeness_or_FULLU_localization_finite_claim=True,
      no_author_worker_GB_newpoint_prime_or_degree_guess=True,
      peak_integer_polynomial_terms=MAXTERMS)
    evidence=out/'independent-pair-norm-finite112-proof.json';evidence.write_text(json.dumps(proof,indent=2)+'\n')
    for p in (Path(__file__),ROOT/f'research/notes/{ID}.md'):hashes[str(p.relative_to(ROOT))]=sha(p)
    report=dict(id=ID,all_passed=True,input_sha256=hashes,evidence_sha256={evidence.name:sha(evidence)},
      source_bound_pair_monic_norm_and_finite112_cover_accepted=True,
      no_bad_nonempty_fullU_finite_or_global_H7_claim=True,
      seconds=time.monotonic()-start,max_rss_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
    print(json.dumps(report))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True)
    main(p.parse_args().output)
