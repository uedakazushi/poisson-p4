#!/usr/bin/env python3
"""Exact primitive/S61 adapter and bounded five-parameter input preflight."""
import argparse,gzip,itertools,json,math,resource,signal,time
from collections import Counter,defaultdict
from fractions import Fraction as Q
from functools import lru_cache
from pathlib import Path
import sympy as s
from s61_model import S61,X,J as J5,BASIS,coefficients,bivector_from_coefficients
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from s61_infinity_flags import old_beta_basis
from elliptic_weight1_dominance import projective_affine
from weight1_structure import COORDS,J
from q_new_family_v_zero_degree26_support_preflight_20261104_v63 import read,sha
ROOT=Path(__file__).resolve().parents[1]
ID='S61_SPARSE_GENERIC5_INPUT_PREFLIGHT_20261104_v63'
ZERO=(0,0,0,0,0)
def add(a,b):
    c=dict(a)
    for m,v in b.items():
        c[m]=c.get(m,Q(0))+v
        if not c[m]:del c[m]
    return c
def scale(p,c):return {m:v*c for m,v in p.items()if v*c}
def mul(a,b):
    c={}
    for m,x in a.items():
        for n,y in b.items():
            k=tuple(i+j for i,j in zip(m,n));c[k]=c.get(k,Q(0))+x*y
    return {m:v for m,v in c.items()if v}
def enc(p):return [[list(m),str(v)]for m,v in sorted(p.items())]
def summary(ps):
    ps=list(ps);return dict(entries=len(ps),nonzero=sum(bool(p)for p in ps),terms=sum(map(len,ps)),maxterms=max(map(len,ps),default=0),degree=max((sum(m)for p in ps for m in p),default=0))
def matching(M):
    chosen={}
    def dfs(i,seen):
        for j,p in enumerate(M[i]):
            if p and j not in seen:
                seen.add(j)
                if j not in chosen or dfs(chosen[j],seen):chosen[j]=i;return True
        return False
    return sum(dfs(i,set())for i in range(len(M)))
def plane_cost(Ai,weights,U,plane,chosen):
    uw=[]
    for i in range(16):
        ww={tuple(weights[j])for j in range(16)if Ai[i,j]};assert len(ww)==1;uw.append(ww.pop())
    std=[(1,0),(-1,0),(0,1),(0,-1)];cw=[tuple(sum(std[i][h]for i in t)for h in range(2))for t in itertools.combinations_with_replacement(range(4),3)];rw=[]
    for row in plane['all51_homogeneous_linear_consequences']:
        ww=set()
        for j,f in enumerate(row):
            for powers,c in s.Poly(s.sympify(f),*U).terms():
                if c:ww.add(tuple(cw[j][h]+sum(p*w[h]for p,w in zip(powers,uw))for h in range(2)))
        assert len(ww)<=1;rw.append(next(iter(ww))if ww else None)
    selector=plane['generic_plane_recipe'];dw=tuple(sum(rw[i][h]for i in selector['L_rows'])-sum(cw[j][h]for j in selector['L_pivot_columns'])for h in range(2));assert dw==(-1,-3)
    degrees_rows=[max((sum(e)for p in row for e in p),default=0)for row in chosen];degrees_columns=[max((sum(e)for row in chosen for e in row[j]),default=0)for j in range(18)];totaldegree=min(sum(degrees_rows),sum(degrees_columns));support=0
    for b in range(19):
        for nt in range(19-b):
            if (18+dw[1]-3*b+nt)%2:continue
            np=(18+dw[1]-3*b+nt)//2;nm=(18-dw[1]+b-3*nt)//2
            if min(np,nm)<0 or (b+nt-dw[0])%2:continue
            for a in range(19-b-nt):
                if (2*a+3*b+3*nt+dw[0]-2*dw[1])%4:continue
                support+=sum(a+b+nt+es+eq<=totaldegree for es in range(np+1)for eq in range(np-es+1))
    den=math.lcm(*(c.denominator for row in chosen for p in row for c in p.values()));M=max(sum(map(abs,p.values()),Q(0))for row in chosen for p in row);bound=math.factorial(18)*(M*den)**18;bits=((bound.numerator+bound.denominator-1)//bound.denominator).bit_length()
    return dict(historical_u_Cartan_weights=uw,cubic20_weights=cw,L51_row_weights=rw,Delta_ambient_beta_degree=18,Delta_Cartan_weight=dw,
        selected18_row_degree_bounds=degrees_rows,selected18_column_degree_bounds=degrees_columns,Delta_total_parameter_degree_bound=totaldegree,Delta_support_bound=support,
        source_support_rule='Let b=Lexp,nt=t0exp,a=Hexp; np=(d+w2-3b+nt)/2,nm=(d-w2+b-3nt)/2 nonnegative integers; a<=d-b-nt,s+q<=np,b+nt=w1 mod2,2a+3b+3nt=-w1+2w2 mod4',
        Delta_constant_denominator=str(den**18),Delta_integer_coefficient_bit_bound=bits,selected18_nonzero_constant_entries=sum(bool(p)and all(not any(e)for e in p)for row in chosen for p in row),
        plane_direction_numerator_beta_degree_bound19=True,plane3_minor_numerator_beta_degree_bound40=True,no18_determinant_or_inverse=True)
def main(out):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False)
    model=S61();weights=model.reduced['parameter_weights'];one=[i for i,w in enumerate(weights)if w==1];two=[i for i,w in enumerate(weights)if w==2];index={old:new for new,old in enumerate(two)}
    ip=ROOT/'reports/Q_NEW_FAMILY_HYPERSURFACE_RECOGNITION_20261103_v63/invariant-certificates.json';info=read(ip);K=s.Matrix(info['primitive24_by16_basis'])
    z,U,beta=primitive_beta(model);M=s.Matrix([s.diff(beta[i,j],z[k])for (i,j),k in COORDS]).jacobian(U);basisrows=list(K.T.rref()[1]);A=K[basisrows,:].inv()*M[basisrows,:];Ai=A.inv();assert K*A==M and A.det()==s.Rational(36,625)
    # Independently bind the first-jet columns to the original constant-J reconstruction.
    raw=s.zeros(150,16)
    for col,ui in enumerate(one):raw[model.certificate['raw_free_indices'][ui],col]=1
    for row in model.data['eliminated_raw_variables']:
        for mon,n,d in row['polynomial']:
            if len(mon)==1 and mon[0]in one:raw[row['raw_index'],one.index(mon[0])]+=s.Rational(n,d)
    rawlinear=list(raw*s.Matrix(U));P1=bivector_from_coefficients(rawlinear);yp=s.symbols('y0:4');aff=projective_affine(P1,X,yp);origin=dict.fromkeys(yp,0);LL=sum((aff.diff(y).subs(origin)*y for y in yp),s.zeros(4))
    zz,bases=old_beta_basis(model);old=sum((u*b for u,b in zip(U,bases)),s.zeros(4));oldL=-J*old*J
    assert (LL-oldL.subs(dict(zip(zz,yp)),simultaneous=True)).applyfunc(s.expand)==s.zeros(4)
    tau=LL[0,1]+LL[2,3];alpha=s.Matrix([s.diff(tau,y)/5 for y in yp]);alphaM=alpha.jacobian(U)
    H,L,ss,qq,tt,lam=s.symbols('H L s0 q0 t0 lambda');params=[H,L,ss,qq,tt]
    sectionpath=ROOT/'reports/Q_NEW_FAMILY_GENERIC_BETA_SPARSE_SECTION_MINOR_20261104_v63/section6.json';section=read(sectionpath);b16=s.Matrix([s.expand(s.sympify(v.replace('lambda','lam_src'),locals={'lam_src':lam}).subs({lam:1,s.Symbol('epsilon'):tt}))for v in section['sparse_beta6']]);us=(Ai*b16).applyfunc(s.expand)
    up=[]
    for e in us:up.append({m:Q(str(c))for m,c in s.Poly(e,*params).terms()if c})
    assert sum(bool(p)for p in up)==9 and sum(map(len,up))==10
    allu={i:{}for i in range(61)}
    for i,p in zip(one,up):allu[i]=p
    # Weighted covariance of all original equations and reconstruction is exact.
    eqweights=[]
    for eq in model.reduced['equations']:
        ww={sum(weights[i]for i in mon)for mon,n,d in eq};assert len(ww)==1;eqweights.append(ww.pop())
    assert Counter(eqweights)=={4:70,5:64}
    rawweights=[]
    for (i,j),mon in BASIS:
        ex=s.Poly(mon,*X).monoms()[0];rawweights.append(2-int(i>0)-int(j>0)+sum(ex[1:]))
    for i,ui in zip(model.certificate['raw_free_indices'],range(61)):assert rawweights[i]==weights[ui]
    for e in model.data['eliminated_raw_variables']:
        assert all(sum(weights[i]for i in mon)==rawweights[e['raw_index']]for mon,n,d in e['polynomial'])
    @lru_cache(None)
    def mono(mon):
        if not mon:return {ZERO:Q(1)}
        p=allu[mon[0]]
        return mul(p,mono(mon[1:]))if p else{}
    def eval_u(text):
        p={}
        for powers,c in s.Poly(s.sympify(text),*U).terms():
            mon=tuple(i for i,e in zip(one,powers)for _ in range(e));p=add(p,scale(mono(mon),Q(str(c))))
        return p
    equations=[];stats=[];linear=[];leadingM=[];leadmons=[]
    for eq in model.reduced['equations']:
        ff={}
        for mon,n,d in eq:
            qm=tuple(index[i]for i in mon if i in index);bm=tuple(i for i in mon if i not in index)
            for e,c in mono(bm).items():
                key=(qm,e);ff[key]=ff.get(key,Q(0))+Q(n,d)*c
        ff={k:c for k,c in ff.items()if c};equations.append(ff);byq=Counter(len(m)for (m,e),c in ff.items());degrees={str(k):max((sum(e)for (m,e),c in ff.items()if len(m)==k),default=0)for k in range(3)}
        stats.append(dict(terms=len(ff),qdegree_counts=dict(byq),parameter_degree_by_qdegree=degrees))
        row=[{}for _ in range(45)]
        for (m,e),c in ff.items():
            if len(m)==1:row[m[0]][e]=c
        linear.append(row);top=max((e for p in row for e in p),default=None);leadmons.append(top);leadingM.append([s.Rational(p.get(top,0))if top is not None else 0 for p in row])
    MM=s.Matrix(leadingM);lr=list(MM.T.rref()[1]);lc=list(MM[lr,:].rref()[1])if lr else[];det=MM.extract(lr,lc).det(method='domain-ge')if lr else s.Integer(1);assert det
    ranklower=len(lr)
    # Global known graph formulas and the single universal affine3-plane.
    path51=ROOT/'reports/Q51_UNIVERSAL_SECTION_20261101.json';path52=ROOT/'reports/Q52_UNIVERSAL_SECTION_20261101.json';pathplane=ROOT/'reports/TORIC_UNIVERSAL_THREE_PLANE_20261101.json';q51,q52,plane=map(read,[path51,path52,pathplane]);assert q51['all_passed']and q52['all_passed']and plane['all_passed']
    assert plane['beta_variables']==list(map(str,U))
    for data in[q51,q52]:
        assert data['weight1_u_indices_zero_based']==one and data['weight2_u_indices_zero_based']==two
        assert all(s.Poly(s.sympify(f),*U).is_zero or set(sum(m)for m,c in s.Poly(s.sympify(f),*U).terms())=={2}for f in data['quadratic_section_q'])
    for name in['all51_homogeneous_linear_consequences','q_linear_matrix_N']:
        assert all(not s.sympify(f)or set(sum(m)for m,c in s.Poly(s.sympify(f),*U).terms())=={1}for row in plane[name]for f in row)
    Q51=list(map(eval_u,q51['quadratic_section_q']));Q52=list(map(eval_u,q52['quadratic_section_q']));center=list(map(eval_u,plane['q_center']));LLin=[[eval_u(f)for f in row]for row in plane['all51_homogeneous_linear_consequences']];NN=[[eval_u(f)for f in row]for row in plane['q_linear_matrix_N']]
    assert all(scale(c,3)==add(scale(a,2),b)for c,a,b in zip(center,Q51,Q52));direction=[add(a,scale(b,-1))for a,b in zip(Q52,center)]
    selector=plane['generic_plane_recipe'];chosen=[[LLin[i][j]for j in selector['L_pivot_columns']]for i in selector['L_rows']];chosenstat=summary(p for row in chosen for p in row);chosenmatch=matching(chosen);chartcost=plane_cost(Ai,info['Cartan_weights'],U,plane,chosen)
    pulled=dict(parameter_variables=list(map(str,params)),q_to_u_indices_zero_based=two,equations=[[[list(m),list(e),str(c)]for (m,e),c in sorted(f.items())]for f in equations],linear_q_matrix134x45=[[enc(p)for p in row]for row in linear])
    with gzip.open(out/'pulled134-input.json.gz','wt')as f:json.dump(pulled,f,separators=(',',':'))
    known=dict(Q51=[enc(p)for p in Q51],Q52=[enc(p)for p in Q52],q_center=[enc(p)for p in center],Q52_minus_center=[enc(p)for p in direction],L51x20=[[enc(p)for p in row]for row in LLin],N45x20=[[enc(p)for p in row]for row in NN],saved_selector=selector)
    with gzip.open(out/'known-five-input.json.gz','wt')as f:json.dump(known,f,separators=(',',':'))
    cert=dict(parameter_variables=list(map(str,params)),primitive_beta16=list(map(str,b16)),historical_weight1_names=list(map(str,U)),historical_weight1_indices=one,historical_weight2_indices=two,
        adapter_b16_equals_A_u16=[[str(c)for c in row]for row in A.tolist()],adapter_u16_equals_Ainv_b16=[[str(c)for c in row]for row in Ai.tolist()],adapter_determinant='36/625',primitive24_basis=info['primitive24_by16_basis'],primitive_projection24x16=[[str(c)for c in row]for row in M.tolist()],original_raw_weight1_columns150x16=[[str(c)for c in row]for row in raw.tolist()],
        Lee_alpha4x16=[[str(c)for c in row]for row in alphaM.tolist()],source_Lee_alpha=list(map(str,(alphaM*us).applyfunc(s.expand))),source_u16=list(map(str,us)),source_nonzero9_terms10=True,
        full150_reconstruction_weight_covariance=True,weighted134_distribution=dict(Counter(eqweights)),
        physical_Gm='lambda^2 times pushforward diag(1,lambda^-1,lambda^-1,lambda^-1,lambda^-1); raw affine P(y) -> P(lambda*y)',
        scalar_strip='epsilon=lambda*t0, q45=lambda^2*Q45; F_i(lambda*u,lambda^2*Q)=lambda^weight_i*F_i(u,Q)',
        Sp4_alone_lambda1_not_dominant_asserted=True,dominant_cover_uses_Sp4_and_Gm=True,
        equation_statistics=stats,zero_equations=[i for i,f in enumerate(equations)if not f],total_pulled_terms=sum(map(len,equations)),
        linear_q_part=summary(p for row in linear for p in row),linear_q_rank_lower_bound=ranklower,linear_q_rank_upper_bound=45,linear_q_exact_rank=45 if ranklower==45 else None,
        linear_q_rank_certificate=dict(method='row-wise leading monomial, lex H,L,s0,q0,t0; product of row leading monomials bounds determinant and selected coefficient is nonzero',rows=lr,columns=lc,leading_monomials=[list(leadmons[i])for i in lr],nonzero_leading_determinant=str(det),matrix=[[str(c)for c in row]for row in MM.extract(lr,lc).tolist()]),
        Q51_statistics=summary(Q51),Q52_statistics=summary(Q52),center_statistics=summary(center),direction_statistics=summary(direction),L51x20_statistics=summary(p for row in LLin for p in row),N45x20_statistics=summary(p for row in NN for p in row),chosen18_block_statistics=chosenstat,chosen18_support_matching=chosenmatch,
        chosen18_minor_or_inverse_NOT_computed=True,chosen18_denominator_name='Delta=det(L_saved18x18)',chosen18_and_plane_cost=chartcost,
        plane_numerator_recipe='V=Delta*N_free-N_pivot*adj(L_pivot)*L_free; twoNTdirections=V/Delta; third=Q52-center',
        one_affine3plane_not_three_planes=True,chosen3direction_rows=[0,1,2],chosen3direction_minor_NOT_computed=True,
        fixed_chart_caution='Dominance is for Sp4 times source. A nonzero non-invariant beta-only chart minor need not stay nonzero at g=identity on this source. Fixed18/3 selectors require their own polynomial check or a different adapted chart; retainDelta as named condition.',
        known_five_Gm='L,N scale lambda; c scales lambda; center,Q51,Q52 scale lambda^2; TnormalizedkerL scales0; affineplane coordinates(cfree2,sigma) weights(1,1,0). Known3+1+1 persist generically via actual equivariant models, not all special beta.',
        no_42normal_J0_substitution_GB_point_prime_or_18_inverse=True,
        next_cost='134equations already available overQQ[H,L,s0,q0,t0][Q45]. First check fixed18minor/plane chart by bounded exact elimination before constructing42normal complement. No generic field5 GB authorized here.')
    (out/'certificate.json').write_text(json.dumps(cert,indent=2)+'\n')
    paths=[Path(__file__),ROOT/f'research/notes/{ID}.md',ip,sectionpath,path51,path52,pathplane,ROOT/'code/s61_model.py',ROOT/'code/primitive_jet_elliptic_curve_20261101.py',ROOT/'code/s61_infinity_flags.py',ROOT/'legacy/v02/data/symplectic_slice_certificate.json',ROOT/'legacy/v02/data/symplectic_slice_equations.json',ROOT/'legacy/v02/data/symplectic_slice_reduced_equations.json',ROOT/'reports/Q_NEW_FAMILY_GENERIC_BETA_SPARSE_SECTION_MINOR_20261104_v63_REVIEW/report.json']
    assert sum(p.stat().st_size for p in out.iterdir()if p.is_file())<32*1024*1024
    report=dict(id=ID,all_passed=True,exact_adapter_and_weighted_input_preflight=True,primitive_to_historical_det='36/625',pulled134_equations=len(equations),pulled_terms=cert['total_pulled_terms'],linear_q_rank_lower=ranklower,linear_q_rank_upper=45,
        fixed18_minor_not_run=True,no_GB_J0_or42normal_substitution_point_prime=True,
        input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},evidence_sha256={p.name:sha(p)for p in out.iterdir()if p.is_file()},seconds=time.monotonic()-start)
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(dict(all_passed=True,seconds=report['seconds'],terms=cert['total_pulled_terms'],zeroeq=cert['zero_equations'],linear_rank_lower=ranklower,L=cert['L51x20_statistics'],N=cert['N45x20_statistics'],block18=chosenstat,block18_matching=chosenmatch)))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args();main(a.output)
