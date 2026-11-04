#!/usr/bin/env python3
"""Literal cached-source2/adjoint98 input for seven core quadrics.

Only degree<=2 cache composition and constant QQ linear algebra are done.
The proposed cubic evaluation has not been performed by this preflight.
"""
import argparse,hashlib,itertools,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
from ambient_B18_top_graph_telescoping_20261104_v63 import forward,linear_solve,ancestors
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_CORE4_QUADRATIC_PREFLIGHT_20261104_v63'
M9=[()]+[(i,) for i in range(9)]+list(itertools.combinations_with_replacement(range(9),2))
M4=[()]+[(i,) for i in range(4)]+list(itertools.combinations_with_replacement(range(4),2))
P4={m:i for i,m in enumerate(M4)}
Z4=(F(0),)*15
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def enc(x):
    if isinstance(x,F):return str(x)
    if isinstance(x,(list,tuple)):return [enc(v) for v in x]
    if isinstance(x,dict):return {str(k):enc(v) for k,v in x.items()}
    return x
def mul(a,b):
    r=[F(0)]*15
    for i,m in enumerate(M4):
        if not a[i]:continue
        for j,n in enumerate(M4):
            if b[j] and len(m)+len(n)<=2:r[P4[tuple(sorted(m+n))]]+=a[i]*b[j]
    return tuple(r)
def compose(a,sub):
    r=[F(a[0])]+[F(0)]*14
    for c,m in zip(a[1:],M9[1:]):
        c=F(c)
        if not c:continue
        v=sub[m[0]] if len(m)==1 else mul(sub[m[0]],sub[m[1]])
        for j,x in enumerate(v):r[j]+=c*x
    return tuple(r)
def main(out):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False)
    cp=ROOT/'reports/AMBIENT_B18_NORMALIZER_CORE4_REDUCTION_20261104_v63/typed-core4-input.json'
    ep=ROOT/'reports/AMBIENT_B18_ZERO_NORMALIZER_QUADRATIC_20261104_v63/cached-source2-and-exact-E4-quadrics.json'
    jp=ROOT/'reports/AMBIENT_B18_T3_FULL9_FIRSTJET_20261104_v63/full9-firstjet-certificate.json'
    C,E,J=map(read,[cp,ep,jp]);nodes=C['DAG_scalar_nodes'];unknown=C['unit_schur_unknown_nodes98'];eq=C['unit_schur_equations98']
    assert len(eq)==len(unknown)==98
    # Keep x112,x117,tauA,tauB; all five selected physical tangent X are zero.
    free9_keep=[2,4,7,8]
    linear4=[]
    for j in range(9):
        a=[F(0)]*15
        if j in free9_keep:a[1+free9_keep.index(j)]=1
        linear4.append(tuple(a))
    schur11=[compose(a,linear4) for a in E['all11_variables_as_remaining9_through2']]
    physical_rho=schur11[:9];bind={};sourceproof={}
    for name in ['CA','CB']:
        f=E['families'][name];cache=J['both57_source2'][name]
        family_F9=[compose(a,physical_rho) for a in f['inverse_F9_in_common_R0_through2']]
        data={int(i):compose(a,family_F9) for i,a in cache['all57_source2jets']}
        for i,a in data.items():
            if i in bind:assert bind[i]==a
            bind[i]=a
        sourceproof[name]={'cached57_unknown_order':cache['unknown57'],
          'family_F9_in_core4_through2':family_F9,'all57_source2_in_core4':list(data.items())}
    # The source57 data contains the exact independent physical x variables.
    param=C['actual9_parameter_nodes']
    for i,node in enumerate(param):
        if node in bind:assert bind[node]==linear4[i]
        else:bind[node]=linear4[i]
    assert set(unknown)<=set(bind)
    for node in unknown:assert bind[node][0]==0
    # Seven residuals: row72-row60 is the sole nontrivial constant row operation.
    residual_indices=[5,6,7,8,9,10,11]
    expressions=[]
    for i in residual_indices:
        expressions.append([(C['all12_actual_difference_nodes'][i],F(1))]+([(C['all12_actual_difference_nodes'][4],F(-1))] if i==6 else []))
    for i in residual_indices:
        a=list(map(F,E['all14_equations_after_exact_linear_Schur_through2'][i+2]))
        if i==6:a=[x-F(y) for x,y in zip(a,E['all14_equations_after_exact_linear_Schur_through2'][6])]
        assert not any(a),('residual not zero through2',i)
    # Actual98 origin Jacobian. The gradient is on the original scalar DAG,
    # before its coefficient and monic equations are imposed.
    values,grad=forward(nodes,set(unknown))
    A=[[grad[r].get(v,F(0)) for v in unknown] for r in eq]
    B=[[sum(c*grad[r].get(v,F(0)) for r,c in expr) for v in unknown] for expr in expressions]
    lamT,det=linear_solve([list(c) for c in zip(*A)],[list(c) for c in zip(*B)])
    lam=[list(c) for c in zip(*lamT)];assert len(lam)==7 and all(len(r)==98 for r in lam)
    assert all(sum(lam[i][k]*A[k][j] for k in range(98))==B[i][j] for i in range(7) for j in range(98))
    roots=eq+[r for e in expressions for r,c in e]
    need=sorted(ancestors(nodes,roots));unbound=[i for i in need if nodes[i][0]=='v' and i not in bind];assert not unbound,unbound
    allmons=[()]+[(i,) for i in range(4)]+list(itertools.combinations_with_replacement(range(4),2))+list(itertools.combinations_with_replacement(range(4),3))
    pairs=sum(1 for a in allmons for b in allmons if len(a)+len(b)<=3)
    mult=sum(nodes[i][0]=='m' for i in need)
    evidence={'original_scalar_DAG_file':str(cp.relative_to(ROOT)),
      'actual98_unknown_nodes':unknown,'actual98_equation_nodes':eq,'actual98_origin_Jacobian':A,
      'seven_residual_linear_combinations':expressions,'residuals_physical_rows':[67,72,63,64,68,73,74],
      'row72_replaced_by_row72_minus_row60':True,'seven_raw_unknown_gradients':B,
      'adjoint98_rows_lambda':lam,'all686_lambda_J_equals_dF_verified':True,'unit98_determinant':det,
      'source2_cache_reuse':sourceproof,'four_parameter_order':['x112','x117','tau_CA','tau_CB'],
      'monomials4_through2':M4,'all_needed_source_and_parameter2_bindings':list(bind.items()),
      'source_DAG_ancestor_nodes':need,'four_parameter_monomials_through3':allmons,
      'known_L_linear_coefficients':[5,-20,24,6],
      'no_source3_correction_needed':'F3=lambda*E3+true-composed-F3; dF/du-lambda*dE/du=0 cancels all98 unknown third jets',
      'no_Npair_L2_or_L3_needed':'all seven modified F have zero2jet, so F3=L1*G2',
      'no_five_IFT2_needed':'all seven modified G have zero1jet, so restriction to tangent U=0 preserves G2',
      'source_multiplication_nodes':mult,'monomial_convolution_pairs35':pairs,
      'worst_dense_multiplication_pairs_upper':mult*pairs,
      'source3_solve_or_cubic_coefficient_evaluation_done':False}
    p=out/'literal-adjoint98-core4-preflight.json';p.write_text(json.dumps(enc(evidence),separators=(',',':'))+'\n');assert p.stat().st_size<32*1024*1024
    paths=[Path(__file__),ROOT/'research/notes'/f'{STEM}.md',cp,ep,jp,
      ROOT/'reports/AMBIENT_B18_N_NORMALIZER_PAIR_PREFLIGHT_20261104_v63_REVIEW/report.json',
      ROOT/'reports/AMBIENT_B18_NORMALIZER_CORE4_REDUCTION_20261104_v63_REVIEW/report.json']
    rep={'id':STEM,'all_passed':True,'actual_unit98_input_and_lambda_identity_verified':True,
      'all98_unknown_source3_coefficients_cancel':True,'four_parameters':4,'monomials_through3':35,
      'seven_modified_residuals_have_zero2jet':True,'Npair_order3_or_source3_solve_required':False,
      'proposed_target':'seven exact core homogeneous quadrics + all35-coefficient L1 division identities + QQ rank',
      'coefficient_computation_done':False,'no_auto_core_cubic_or_source4_job':True,
      'ancestor_nodes':len(need),'worst_dense_multiplication_pairs_upper':mult*pairs,
      'input_sha256':{str(p.relative_to(ROOT)):sha(p) for p in paths},'evidence_sha256':{p.name:sha(p)},
      'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    (out/'report.json').write_text(json.dumps(rep,indent=2)+'\n');print(json.dumps(rep))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
