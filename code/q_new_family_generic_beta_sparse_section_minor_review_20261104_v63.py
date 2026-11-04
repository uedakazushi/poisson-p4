#!/usr/bin/env python3
"""Independent tensor Lie action, exact Schur binding and Bareiss transcript replay."""
import argparse
from fractions import Fraction as F
import gzip
import itertools
import json
import math
from pathlib import Path
import resource
import signal
import time
import sympy as s
import q_new_family_degree28_exact_equation_review_20261104_v63 as q

ROOT=Path(__file__).resolve().parents[1]
BASE='Q_NEW_FAMILY_GENERIC_BETA_SPARSE_SECTION_MINOR_20261104_v63'
STEM=BASE+'_REVIEW'

def kernel_rref(A):
    a=[[F(v) for v in row] for row in A];r=0;piv=[]
    for j in range(len(a[0])):
        choices=[i for i in range(r,len(a)) if a[i][j]]
        if not choices:continue
        k=choices[0];a[r],a[k]=a[k],a[r];v=a[r][j];a[r]=[c/v for c in a[r]]
        for i in range(len(a)):
            if i!=r and a[i][j]:v=a[i][j];a[i]=[x-v*y for x,y in zip(a[i],a[r])]
        piv.append(j);r+=1
        if r==len(a):break
    free=[j for j in range(len(a[0])) if j not in piv]
    B=[[F(0) for _ in free] for _ in range(len(a[0]))]
    for c,j in enumerate(free):
        B[j][c]=1
        for i,k in enumerate(piv):B[k][c]=-a[i][j]
    return B,free

def enc(p):return q.record(p)
def decode_rows(rows):return q.from_records([(row[1:],row[0]) for row in rows])
def derivative(p,j):
    out={}
    for m,c in p.a.items():
        e=list(q.old.unpack(m))
        if e[j]:v=e[j];e[j]-=1;out[q.old.pack(e)]=c*v
    return q.QP(out,p.d)

def main(out):
    started=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
    def timeout(*_):raise TimeoutError('independent sparse section review60s cap')
    signal.signal(signal.SIGALRM,timeout);signal.alarm(60)
    folder=ROOT/'reports'/BASE;pre=ROOT/'reports/Q_NEW_FAMILY_GENERIC_BETA_SPARSE_SECTION_PREFLIGHT_20261104_v63'
    report=q.old.read(folder/'report.json');pc=q.old.read(pre/'certificate.json')
    transcript=q.old.read(folder/'native-certificate.json.gz');section=q.old.read(folder/'section6.json')
    assert report['all_passed'] and transcript['all_passed'] and report['nonzero_deleted_row']==4
    checked={}
    gates=[ROOT/'reports/Q_NEW_FAMILY_V_ZERO_BETA_SECTION_20261104_v63_REVIEW/report.json',ROOT/'reports/Q_NEW_FAMILY_DEGREE28_EXACT_EQUATION_20261104_v63_REVIEW/report.json']
    for rp in [folder/'report.json',pre/'report.json',*gates]:
        d=q.old.read(rp);assert d['all_passed']
        for name,h in d['input_sha256'].items():assert q.old.sha(ROOT/name)==h,name;checked[name]=h
        for name,h in d.get('evidence_sha256',{}).items():
            p=rp.parent/name;assert q.old.sha(p)==h,str(p);checked[str(p.relative_to(ROOT))]=h
    ip=ROOT/'reports/Q_NEW_FAMILY_HYPERSURFACE_RECOGNITION_20261103_v63/invariant-certificates.json'
    info=q.old.read(ip)
    # Independent primitive basis from d(beta)=0 and Lambda(beta)=0 in24 tensors.
    pairs=list(itertools.combinations(range(4),2));triples=list(itertools.combinations(range(4),3))
    coords=[(i,j,k) for i,j in pairs for k in range(4)];ci={v:i for i,v in enumerate(coords)}
    A=[[0]*24 for _ in range(8)]
    for col,(i,j,k) in enumerate(coords):
        if k not in [i,j]:
            order=(k,i,j);sgn=(-1)**sum(order[r]>order[t] for r in range(3) for t in range(r+1,3))
            A[triples.index(tuple(sorted(order)))][col]=sgn
        if (i,j) in [(0,1),(2,3)]:A[4+k][col]=1
    B,free=kernel_rref(A)
    assert len(free)==16 and B==[list(map(F,row)) for row in info['primitive24_by16_basis']]
    J=[[0,1,0,0],[-1,0,0,0],[0,0,0,1],[0,0,-1,0]]
    generators=[]
    for i in range(4):
        for j in range(i,4):
            S=[[0]*4 for _ in range(4)];S[i][j]=S[j][i]=1
            G=[[sum(J[r][k]*S[k][t] for k in range(4)) for t in range(4)] for r in range(4)]
            assert all(sum(G[k][r]*J[k][t]+J[r][k]*G[k][t] for k in range(4))==0 for r in range(4) for t in range(4))
            generators.append(G)
    def tensor_entry(v,i,j,k):return 0 if i==j else v[ci[i,j,k]] if i<j else -v[ci[j,i,k]]
    rho=[]
    for G in generators:
        columns=[]
        for c in range(16):
            v=[row[c] for row in B]
            w=[sum(G[l][i]*tensor_entry(v,l,j,k)+G[l][j]*tensor_entry(v,i,l,k)+G[l][k]*tensor_entry(v,i,j,l) for l in range(4)) for i,j,k in coords]
            z=[w[r] for r in free]
            assert all(sum(B[r][j]*z[j] for j in range(16))==w[r] for r in range(24))
            columns.append(z)
        rho.append([list(row) for row in zip(*columns)])
    assert rho==[[list(map(F,row)) for row in G] for G in info['all_ten_representation_matrices']]
    weights=info['Cartan_weights']
    for row,(i,j,k) in enumerate(coords):
        w1=sum([1,-1,0,0][v] for v in [i,j,k]);w2=sum([0,0,1,-1][v] for v in [i,j,k])
        for c in range(16):assert not B[row][c] or weights[c]==[w1,w2]
    H,L,ss,qq,lam,xi,r,BB=s.symbols('H L s0 q0 lambda xi0 r0 B0');vs=(H,L,ss,qq)
    def parse(z):return s.sympify(z.replace('lambda','lam_parameter'),locals={'lam_parameter':lam})
    def expr(p):return sum((s.Rational(c,p.d)*s.prod(v**e for v,e in zip(vs,q.old.unpack(m))) for m,c in p.a.items()),s.S.Zero)
    def poly(z):
        z=s.Poly(s.cancel(z),*vs)
        return q.from_records([(e,str(c)) for e,c in z.terms() if c])
    bp=ROOT/'reports/Q_NEW_FAMILY_V_ZERO_BETA_SECTION_20261104_v63/certificate.json'
    oldbeta=q.old.read(bp)['beta_bar0_16']
    pp=ROOT/'reports/Q_NEW_FAMILY_V_ZERO_DEGREE26_SUPPORT_PREFLIGHT_20261104_v63/certificate.json'
    beta=[q.from_records(z) for z in q.old.read(pp)['beta16_sparse']]
    gauge=[lam*L/H**2,H**2/(lam*L),lam*L**2/H**2,H**2/(lam*L**2)]
    assert s.cancel(gauge[0]*gauge[1]-1)==s.cancel(gauge[2]*gauge[3]-1)==0
    for i,(w1,w2) in enumerate(weights):
        pulled=s.sympify(oldbeta[i]).subs({xi:lam**2,r:H**2,BB:L**2/H**3})*gauge[0]**w1*gauge[2]**w2
        assert s.cancel(pulled-lam*expr(beta[i]))==0
        assert s.cancel(parse(pc['new_beta16'][i])-lam*expr(beta[i]))==0
    assert s.Matrix([[2,0,0],[0,2,0],[0,-3,2]]).det()==8
    orbit=[[q.total(q.scale(b,c) for b,c in zip(beta,row)) for row in G] for G in rho]
    db=[[derivative(p,j) for p in beta] for j in range(4)]
    jaccols=orbit+[beta]+db
    Jac=[list(row) for row in zip(*jaccols)]
    assert Jac==[[decode_rows(z) for z in row] for row in pc['stripped_Jacobian16x15']]
    assert pc['actual_lambda_column_factors']==[1]*10+[0]+[1]*4
    piv=[2,8,0,1,5];remaining=[i for i in range(16) if i not in piv]
    assert piv==pc['source_pivot_rows'] and remaining==pc['remaining_rows']
    # Solve the fixed tiny source-coordinate block independently, not with author formulas.
    S=s.Matrix([[expr(Jac[i][j]) for j in range(10,15)] for i in piv])
    detS=s.factor(S.det());assert detS==45000*H
    Sinv=S.inv();ownK=[[q.ZERO for _ in range(10)] for _ in remaining]
    identities=0
    for c in range(10):
        factors=Sinv*s.Matrix([expr(orbit[c][i]) for i in piv])
        stored=pc['exact_elimination_coefficients'][c]
        expected=[expr(decode_rows(z)) for z in stored[:4]]+[expr(decode_rows(stored[4]))/H]
        assert all(s.cancel(a-b)==0 for a,b in zip(factors,expected))
        for i in range(16):
            residual=s.cancel(expr(orbit[c][i])-sum(factors[j]*expr(Jac[i][10+j]) for j in range(5)))
            if i in piv:assert residual==0
            else:ownK[remaining.index(i)][c]=poly(residual)
            # Multiplication by H proves the saved elimination over the polynomial ring.
            assert s.cancel(H*(expr(orbit[c][i])-sum(expected[j]*expr(Jac[i][10+j]) for j in range(5))-residual))==0
            identities+=1
    assert identities==160 and ownK==[[decode_rows(z) for z in row] for row in pc['Schur11x10']]
    # Replay supplied fraction-free trace only: no polynomial division or new pivot search.
    stepcount=0;attempt_summary=[]
    for attempt in transcript['attempts']:
        deleted=attempt['deleted_original_row'];selected=[i for i in remaining if i!=deleted]
        assert selected==attempt['selected_original_rows']
        M=[ownK[remaining.index(i)].copy() for i in selected]
        swaps=attempt['row_swap_pairs_flat'];assert len(swaps)%2==0
        swapmap=dict(zip(swaps[::2],swaps[1::2]));assert len(swapmap)*2==len(swaps)
        steps=attempt['exact_division_steps'];pos=0;previous=q.ONE;sign=1;computed=None
        for k in range(9):
            if pos==len(steps) or steps[pos]['stage']!=k:
                assert all(not M[i][k] for i in range(k,10));computed=q.ZERO;break
            if k in swapmap:
                pivot=swapmap[k];assert k<pivot<10 and not M[k][k] and M[pivot][k]
                M[k],M[pivot]=M[pivot],M[k];sign=-sign
            pivot=M[k][k];assert pivot and previous
            for i in range(k+1,10):
                left=M[i][k]
                for j in range(k+1,10):
                    t=steps[pos];pos+=1;stepcount+=1
                    assert (t['stage'],t['row'],t['column'])==(k,i,j)
                    assert q.saved_poly(t['pivot'])==pivot and q.saved_poly(t['old_ij'])==M[i][j]
                    assert q.saved_poly(t['old_ik'])==left and q.saved_poly(t['old_kj'])==M[k][j]
                    assert q.saved_poly(t['divisor'])==previous
                    answer=q.saved_poly(t['quotient'])
                    assert q.mul(answer,previous)==q.add(q.mul(pivot,M[i][j]),q.scale(q.mul(left,M[k][j]),-1))
                    M[i][j]=answer
                M[i][k]=q.ZERO
            previous=pivot
        assert pos==len(steps)
        if computed is None:computed=q.scale(M[9][9],sign)
        assert computed==q.saved_poly(attempt['determinant'])
        attempt_summary.append({'deleted_row':deleted,'divisions':len(steps),'nonzero':bool(computed)})
    assert stepcount==854 and [x['deleted_row'] for x in attempt_summary]==[7,3,4]
    assert [x['nonzero'] for x in attempt_summary]==[False,False,True]
    detK=computed;Gamma=poly(1+5*qq+8*ss+15*H**2+115*H**2*qq)
    expected=poly(s.Rational(2**30,5605041796875)*H**2*L**2*expr(Gamma));assert detK==expected
    found=4;sortedrows=[i for i in range(16) if i!=found]
    reordered=piv+[i for i in remaining if i!=found]
    perm=[sortedrows.index(i) for i in reordered]
    row_sign=(-1)**sum(perm[i]>perm[j] for i in range(15) for j in range(i+1,15))
    column_sign=(-1)**(5*10);laplace_sign=(-1)**(found+15)
    assert row_sign==section['row_permutation_sign']==-1 and column_sign==1 and laplace_sign==-1
    det16=q.scale(q.mul(poly(45000*H),detK),row_sign*column_sign*laplace_sign)
    wanted=poly(s.Rational(2**33,996451875)*H**3*L**2*expr(Gamma));assert det16==wanted
    assert section['full16_Jacobian_determinant']['lambda_power']==14
    assert det16==decode_rows(section['full16_Jacobian_determinant']['HLsq_polynomial_terms'])
    assert len(det16.a)==5 and Gamma.a.get(0)==Gamma.d
    assert section['deleted_row']==4 and not beta[4]
    for i,z in enumerate(section['sparse_beta6']):
        eps=s.Symbol('epsilon');assert s.cancel(parse(z)-lam*expr(beta[i])-(eps if i==4 else 0))==0
    assert sum(bool(b) for b in beta)+1==10 and sum(len(b.a) for b in beta)+1==13
    out.mkdir(parents=True,exist_ok=False)
    ep=out/'independent-proof.json'
    proof={'primitive24_by16_basis_QQ':[[str(c) for c in row] for row in B],
        'independent_original10_rho16':[[[str(c) for c in row] for row in G] for G in rho],
        'Lie_source_column_order':['orbit0..9','d_lambda','d_H','d_L','d_s0','d_q0'],
        'source_block_determinant':'45000*H','Schur11x10_exactQQ':[[enc(p) for p in row] for row in ownK],
        'all160_original_elimination_identities_verified':True,'all854_Bareiss_remultiplication_identities_verified':True,
        'attempts':attempt_summary,'source5_orbit10_column_permutation_sign':column_sign,
        'row_permutation_sign':row_sign,'epsilon_laplace_sign':laplace_sign,
        'full16_Jacobian_lambda_power':14,'full16_Jacobian_HLsq':enc(det16),
        'nonzero_constant_of_Gamma':1,'finite_torus_cover_degree':8,
        'dominant_generically_finite_actual_Sp4_times_A6_to_W16':True,
        'birational_global_finite_or_generic_beta_exhaustion_claim':False}
    ep.write_text(json.dumps(proof,indent=2)+'\n')
    for p in [Path(__file__),ROOT/'code/q_new_family_degree28_exact_equation_review_20261104_v63.py',ROOT/'research/notes'/(STEM+'.md'),ROOT/'research/notes'/(BASE+'.md'),folder/'report.json',folder/'section6.json',folder/'native-certificate.json.gz',pre/'report.json',pre/'certificate.json',ip,bp,pp,*gates]:checked[str(p.relative_to(ROOT))]=q.old.sha(p)
    result={'id':STEM,'all_passed':True,'primitive_basis_and_all10_Lie16_matrices_from_actual24_tensors':True,
        'all16_actual_finite_gauge_and_lambda_derivatives_bound':True,
        'all16by15_entries_and160_source_elimination_identities_verified':True,
        'all854_saved_exact_divisions_replayed_without_author_solver':True,
        'full16_Jacobian_positive_2pow33_over996451875_lambda14_H3_L2_Gamma':True,
        'actual_six_parameter_cover_dominant_generically_finite':True,
        'faithful_generic_beta_field_basechange_accepted':True,
        'birationality_global_finiteness_or_generic_beta_exhaustion_claim':False,
        'input_sha256':checked,'evidence_sha256':{ep.name:q.old.sha(ep)},
        'seconds':time.monotonic()-started,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'frozen':True}
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');signal.alarm(0)
    print(json.dumps({k:result[k] for k in ['all_passed','seconds','max_rss_kib']}))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
