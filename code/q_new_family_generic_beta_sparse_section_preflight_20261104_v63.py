#!/usr/bin/env python3
"""Small exact Jacobian entries/support only; no transversality rank/minor job."""
import argparse,json,math,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
import sympy as s
from q_new_family_v_zero_degree26_support_preflight_20261104_v63 import read,sha,add,mul,scale,l1,ceilq
ROOT=Path(__file__).resolve().parents[1]
ID='Q_NEW_FAMILY_GENERIC_BETA_SPARSE_SECTION_PREFLIGHT_20261104_v63'
Z=(0,0,0,0)
def dec(p):return {tuple(m):Q(c)for m,c in p}
def enc(p):return [[str(c),*m]for m,c in sorted(p.items())]
def deriv(p,j):
    ans={}
    for m,c in p.items():
        if m[j]:e=list(m);e[j]-=1;ans[tuple(e)]=c*m[j]
    return ans
def expr(p,vs):return sum((s.Rational(c.numerator,c.denominator)*s.prod(v**e for v,e in zip(vs,m))for m,c in p.items()),s.Integer(0))
def stats(ps):
    return dict(nonzero=sum(bool(p)for p in ps),terms=sum(map(len,ps)),maxterms=max(map(len,ps)),max_total_degree=max((sum(m)for p in ps for m in p),default=0),max_HL_degree=max((m[0]+m[1]for p in ps for m in p),default=0),max_sq_degree=max((m[2]+m[3]for p in ps for m in p),default=0))
def main(out):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False)
    ip=ROOT/'reports/Q_NEW_FAMILY_HYPERSURFACE_RECOGNITION_20261103_v63/invariant-certificates.json'
    pp=ROOT/'reports/Q_NEW_FAMILY_V_ZERO_DEGREE26_SUPPORT_PREFLIGHT_20261104_v63/certificate.json'
    bp=ROOT/'reports/Q_NEW_FAMILY_V_ZERO_BETA_SECTION_20261104_v63/certificate.json'
    gates=[ROOT/'reports/Q_NEW_FAMILY_V_ZERO_BETA_SECTION_20261104_v63_REVIEW/report.json',ROOT/'reports/Q_NEW_FAMILY_DEGREE28_EXACT_EQUATION_20261104_v63_REVIEW/report.json']
    assert all(read(p)['all_passed']for p in gates)
    info,pre,old=map(read,[ip,pp,bp]);beta=[dec(p)for p in pre['beta16_sparse']];db=[[deriv(p,j)for p in beta]for j in range(4)];rho=info['all_ten_representation_matrices'];weights=info['Cartan_weights']
    H,L,ss,q,lam,xi,r,B=s.symbols('H L s0 q0 lambda xi0 r0 B0');vs=[H,L,ss,q]
    newg=[lam*L/H**2,H**2/(lam*L),lam*L**2/H**2,H**2/(lam*L**2)]
    for i,(text,w)in enumerate(zip(old['beta_bar0_16'],weights)):
        original=s.sympify(text).subs({xi:lam**2,r:H**2,B:L**2/H**3})
        assert s.cancel(original*newg[0]**w[0]*newg[2]**w[1]-lam*expr(beta[i],vs))==0
    # All Lie entries are actual original16 representation actions, not matrix-pair substitutes.
    orbit=[]
    for A in rho:
        col=[]
        for row in A:
            p={}
            for c,b in zip(row,beta):p=add(p,scale(b,Q(str(c))))
            col.append(p)
        orbit.append(col)
    jaccols=orbit+[beta]+db
    J=[[jaccols[j][i]for j in range(15)]for i in range(16)]
    # Five source-coordinate pivots only; the15-column transversality rank is NOT computed.
    pivrows=[2,8,0,1,5];block=s.Matrix([[expr(J[i][j],vs)for j in range(10,15)]for i in pivrows]);dblock=s.factor(block.det());assert dblock==45000*H
    remaining=[i for i in range(16)if i not in pivrows];assert len(remaining)==11
    h={ (1,0,0,0):Q(1)};spoly={(0,0,1,0):Q(1)};lpoly={(0,1,0,0):Q(1)}
    schur=[];cof=[]
    for v in orbit:
        a=scale(v[2],Q(1,2400));bc=add(scale(v[8],Q(-3,32)),scale(mul(h,a),-1));sc=add(scale(v[0],Q(-16,45)),scale(mul(spoly,a),-1));lc=add(scale(v[1],Q(-1,30)),scale(mul(lpoly,a),-1))
        rem=[add(add(add(add(v[i],scale(mul(a,beta[i]),-1)),scale(mul(bc,db[0][i]),-1)),scale(mul(sc,db[2][i]),-1)),scale(mul(lc,db[1][i]),-1))for i in range(16)]
        qc_numerator=scale(rem[5],48)
        schur.append([rem[i]for i in remaining]);cof.append([a,bc,lc,sc,qc_numerator])
        for i in range(16):
            lhs=add(mul(h,rem[i]),scale(mul(qc_numerator,db[3][i]),-1));rhs=mul(h,rem[i])if i in remaining else{}
            assert lhs==rhs
    K=[[schur[j][i]for j in range(10)]for i in range(11)]
    constant_entries=[]
    for i,row in enumerate(K):
        for j,p in enumerate(row):
            if p and set(p)=={Z}:constant_entries.append(dict(original_row=remaining[i],Lie_column=j,value=str(p[Z])))
    rowstats=[stats(row)for row in K];bounds=[]
    for deleted in range(11):
        rr=[z for i,z in enumerate(rowstats)if i!=deleted];dh=sum(z['max_HL_degree']for z in rr);ds=sum(z['max_sq_degree']for z in rr);dt=sum(z['max_total_degree']for z in rr)
        count=sum(1 for a in range(dh+1)for b in range(dh-a+1)for c in range(ds+1)for e in range(ds-c+1)if a+b+c+e<=dt)
        bounds.append(dict(deleted_original_row=remaining[deleted],HL_degree=dh,sq_degree=ds,total_degree=dt,monomial_bound=count))
    flat=[p for row in K for p in row];den=math.lcm(*(c.denominator for p in flat for c in p.values()));M=max(map(l1,flat));bit=ceilq(math.factorial(10)*(den*M)**10).bit_length()
    cert=dict(source_parameters=['lambda','H','L','s0','q0'],old_source_parameters=['xi0','r0','s0','q0','B0'],
        finite_cover=dict(xi0='lambda^2',r0='H^2',B0='L^2/H^3',degree8_on_torus=True,localizers=['lambda','H','L']),
        positive_pullback_gauge=list(map(str,newg)),all16_gauge_identities=True,
        new_beta16=[str(lam*expr(p,vs))for p in beta],polynomial_beta_nonzeros9_terms12=True,
        Lie_column_convention='rho_j*beta; original10 exact16x16 primitive Lie matrices',
        Jacobian_column_order=['orbit0..9','d_lambda','d_H','d_L','d_s0','d_q0'],actual_lambda_column_factors=[1]*10+[0]+[1]*4,all15_minors_lambda_factor14=True,
        stripped_Jacobian16x15=[[enc(p)for p in row]for row in J],Jacobian_statistics=stats([p for row in J for p in row]),
        source_pivot_rows=pivrows,source_pivot_determinant='45000*H',source_block_only_not_full_transversality=True,
        remaining_rows=remaining,Schur11x10=[[enc(p)for p in row]for row in K],Schur_row_statistics=rowstats,
        exact_elimination_coefficients=[[enc(p)for p in row]for row in cof],elimination_coefficient_order=['a','H_coefficient','L_coefficient','s_coefficient','H*q_coefficient'],all160_polynomial_elimination_identities=True,
        nonzero_constant_entries=constant_entries,all11_determinant_support_bounds=bounds,
        determinant10_common_constant_denominator=str(den**10),determinant10_integer_coefficient_bit_bound=bit,
        full_Jacobian_deleted_minor='plus/minus lambda^14*(45000*H)*det(Schur deleting same remaining row)',
        candidate_transverse_coordinates=remaining,zero_beta_coordinate_candidates=[i for i in remaining if not beta[i]],
        no_transversality_minor_rank_point_or_prime_job=True,
        proposed_next_job='ONE symbolic fractionfree exact determinant job on these10x10 Schur candidates, begin deletedrow7 then finite remaining rows only if needed within60s/2GB/32MiB; all divisions exact verified; save first nonzero polynomial/coefficient witness and determinant identity. No point probe.',
        geometric_inference='Nonzero deletedrowj minor implies map Sp4 x (lambda,H,L,s0,q0,epsilon) -> W16, g*(lambda*beta+epsilon*e_j), is dominant and generically finite inchar0. This supplies faithful field basechange for generic-beta exhaustion, not birationality or an already-verified rank.')
    (out/'certificate.json').write_text(json.dumps(cert,indent=2)+'\n')
    paths=[Path(__file__),ROOT/f'research/notes/{ID}.md',ip,pp,bp,*gates]
    report=dict(id=ID,all_passed=True,only_low_entry_support_and_finite_gauge_preflight=True,no_full15_minor_rank_point_prime=True,
        input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},evidence_sha256={'certificate.json':sha(out/'certificate.json')},seconds=time.monotonic()-start)
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(dict(all_passed=True,seconds=report['seconds'],Schur_rows=remaining,row_stats=rowstats,constant_entries=constant_entries,max_minor_terms=max(z['monomial_bound']for z in bounds),minor_integer_bits=bit)))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args();main(a.output)
