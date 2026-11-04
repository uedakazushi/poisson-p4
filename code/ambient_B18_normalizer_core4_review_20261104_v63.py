#!/usr/bin/env python3
"""Independent typed-source/linear-division/H0 audit of the E4 core4 recipe."""
import argparse,hashlib,itertools,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
BASE='AMBIENT_B18_NORMALIZER_CORE4_REDUCTION_20261104_v63';STEM=BASE+'_REVIEW'
def read(p):return json.loads(p.read_text())
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()
def enc(a):
    if isinstance(a,F):return str(a)
    if isinstance(a,(list,tuple)):return [enc(x) for x in a]
    if isinstance(a,dict):return {str(k):enc(v) for k,v in a.items()}
    return a
def main(out):
    start=time.monotonic();signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False)
    cp=ROOT/'reports'/BASE/'typed-core4-input.json';cr=cp.parent/'report.json'
    sp=ROOT/'reports/AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63/actual-source-coefficients.json'
    ep=ROOT/'reports/AMBIENT_B18_ZERO_NORMALIZER_QUADRATIC_20261104_v63/cached-source2-and-exact-E4-quadrics.json'
    np=ROOT/'reports/AMBIENT_B18_N_NORMALIZER_PAIR_PREFLIGHT_20261104_v63/actual-N-pair-preflight-evidence.json'
    nr=ROOT/'reports/AMBIENT_B18_N_NORMALIZER_PAIR_PREFLIGHT_20261104_v63_REVIEW/report.json'
    ne=nr.parent/'independent-N-pair-and-root-proof.json'
    gates=[cr,nr,ROOT/'reports/AMBIENT_B18_ZERO_NORMALIZER_QUADRATIC_20261104_v63_REVIEW/report.json',
      ROOT/'reports/AMBIENT_B18_ZERO_NORMALIZER_TENSOR_MODEL_20261104_v63_REVIEW/report.json']
    assert all(read(p)['all_passed'] for p in gates)
    C,S,E,N,NP=map(read,[cp,sp,ep,np,ne]);nodes=[tuple(n) for n in C['DAG_scalar_nodes']]
    assert nodes[:len(S['DAG_scalar_nodes'])]==[tuple(n) for n in S['DAG_scalar_nodes']]
    lookup={n:i for i,n in enumerate(nodes)};assert len(lookup)==len(nodes)
    zero=lookup[('c',0,1)];one=lookup[('c',1,1)]
    values=[]
    for i,n in enumerate(nodes):
        op=n[0]
        if op=='c':v=F(n[1],n[2])
        elif op=='v':v=F(0)
        elif op=='n':v=-values[n[1]]
        elif op=='a':v=values[n[1]]+values[n[2]]
        elif op=='m':v=values[n[1]]*values[n[2]]
        elif op=='i':assert values[n[1]];v=1/values[n[1]]
        else:raise ValueError(op)
        values.append(v)
    def neg(i):
        n=nodes[i]
        if i==zero:return zero
        if n[0]=='n':return n[1]
        if n[0]=='c':return lookup[('c',(-values[i]).numerator,(-values[i]).denominator)]
        return lookup[('n',i)]
    def add(i,j):
        if i==zero:return j
        if j==zero:return i
        if nodes[i]==('n',j) or nodes[j]==('n',i):return zero
        if nodes[i][0]==nodes[j][0]=='c':
            q=values[i]+values[j];return lookup[('c',q.numerator,q.denominator)]
        return lookup[('a',*sorted((i,j)))]
    def mul(i,j):
        if i==zero or j==zero:return zero
        if i==one:return j
        if j==one:return i
        if nodes[i][0]==nodes[j][0]=='c':
            q=values[i]*values[j];return lookup[('c',q.numerator,q.denominator)]
        return lookup[('m',*sorted((i,j)))]
    x=[lookup[('v',f'x{r}')] for r in S['physical_zero_base_rows9']]
    tau=C['tau_nodes'];assert nodes[tau['CA']]==('v','E4_CORE_CA_tau') and nodes[tau['CB']]==('v','E4_CORE_CB_tau')
    us=[];eq=[];squares=[]
    for family in ['CA','CB']:
        rec=S['families'][family]
        us += [i for p in rec['coefficient_unknown_pairs'].values() for i in p]+[rec['p_quadratic_relation'][k] for k in ['s','t']]
        eq+=rec['all48_zero_IFT_quadratic_coefficient_equations']
        squares.append(add(mul(tau[family],tau[family]),neg(rec['actual_tracefree_nu'])))
    assert C['unit_schur_unknown_nodes98']==us+[x[0],x[8]]
    assert C['unit_schur_equations98']==eq+squares
    assert C['actual9_parameter_nodes']==x[1:8]+[tau['CA'],tau['CB']]
    targets=[18,19,28,29,60,67,72,63,64,68,73,74];diff=[]
    for row in targets:
        q=[]
        for family in ['CA','CB']:
            a,b=S['families'][family]['complete_normalized_physical_w_coefficients'][0][row]
            q.append(add(a,mul(b,tau[family])))
        diff.append(add(q[0],neg(q[1])))
    assert diff==C['all12_actual_difference_nodes']
    assert all(values[i]==0 for i in eq+squares+diff)
    linear=[[F(v) for v in r] for r in E['linear_Jacobian14_by11']]
    schur=linear[0][0]*linear[1][8]-linear[0][8]*linear[1][0];assert schur==F(-20,9)
    ell=list(map(F,NP['own_normal_linear_equation']));assert ell[-1]==6
    mons=list(itertools.combinations_with_replacement(range(9),2));pos={m:i for i,m in enumerate(mons)}
    G=[]
    for q in E['all14_equations_after_exact_linear_Schur_through2'][2:]:
        assert not any(map(F,q[:10]));q=list(map(F,q[10:]));g=[F(0)]*9
        g[8]=q[pos[(8,8)]]/ell[8]
        for i in range(8):g[i]=(q[pos[(i,8)]]-ell[i]*g[8])/ell[8]
        assert all(v==(ell[i]*g[i] if i==j else ell[i]*g[j]+ell[j]*g[i]) for (i,j),v in zip(mons,q))
        G.append(g)
    assert G==[[F(v) for v in g] for g in C['G12_firstjets_from_existing_quadratic_cache']]
    selected=[targets.index(r) for r in [18,19,60,28,29]];corevars=[2,4,7,8];uvars=[0,1,3,5,6]
    A=[[G[i][j] for j in uvars] for i in selected]
    assert A==[[F(i==j,24 if i==2 else 12) for j in range(5)] for i in range(5)]
    assert all(G[i][j]==0 for i in selected for j in corevars)
    det5=F(1)
    for i in range(5):det5*=A[i][i]
    assert det5==F(1,497664)==F(C['five_minor'])
    remaining=[i for i in range(12) if i not in selected]
    assert [targets[i] for i in remaining]==[67,72,63,64,68,73,74]
    assert all(G[i][j]==0 for i in remaining for j in corevars)
    ops=C['typed_operations'];assert [p['id'] for p in ops]==['SCHUR98','N_BETA','L','G12','CORE5','K7']
    assert ops[0]['unknown_nodes']==us+[x[0],x[8]] and ops[0]['equation_nodes']==eq+squares
    assert ops[0]['output12_difference_nodes']==diff
    assert ops[1]['finite_function_inputs']==N['actual_finite_regular_function_input']
    assert ops[1]['no_vanishing_denominator_unit_inverse'] and ops[1]['minor']=='-1/2'
    assert ops[2]['firstjet']==N['unique_linear_equation_coefficients'] and ops[2]['formula']=='6*(tau_CB-N_BETA(X,tau_CA))'
    assert not ops[3]['localization_at_L'] and ops[3]['number_of_outputs']==12
    assert ops[4]['equations_G12_indices']==selected and ops[4]['remaining_coordinates']==[112,117,'tau_CA','tau_CB']
    assert ops[5]['functions_G12_indices']==remaining
    evidence={'all_original_source_nodes_and96_equations_retained':True,'all12_output_nodes_independently_rebuilt':diff,
      'typed98_unknowns_and_equations_checked':True,'actual_effective_monic_Schur_minor':schur,
      'own_all12_quotient_firstjets_by_exact_quadric_linear_division':G,
      'own_selected5_unit_Jacobian':A,'own_selected5_unit_minor':det5,
      'core4_coordinates':[112,117,'tau_CA','tau_CB'],'remaining7_physical_rows':[targets[i] for i in remaining],
      'all6_typed_operator_interfaces_checked':True,
      'coordinate_clarification':'after CORE5, dLprime/dtau_CB has origin6 and is a unit; it is not asserted to be the constant function6',
      'exact_sequence_proof':'S regular domain and L nonzero; La in LJ implies a in J. Quotient is S/(L), regular8.',
      'H0_E4_equals_H0_core4':True,'L_times_selected5_equations_retained':True,
      'full_core_coefficients_socle_or_D0_depth_claimed':False}
    p=out/'independent-typed-core4-proof.json';p.write_text(json.dumps(enc(evidence),separators=(',',':'))+'\n')
    paths=[Path(__file__),ROOT/'research/notes'/f'{STEM}.md',cp,cr,sp,ep,np,nr,ne,*gates,
      ROOT/'code/ambient_B18_normalizer_core4_reduction_20261104_v63.py',ROOT/'research/notes'/f'{BASE}.md']
    rep={'id':STEM,'all_passed':True,'actual_typed_source98_all12_outputs_and_operator_interfaces_independently_checked':True,
      'own_quotient_linear5_minor':'1/497664','H0_E4_equals_H0_four_variable_core_accepted':True,
      'whole_E4_identified_with_core4_or_D0_socle_depth_claimed':False,
      'no_producer_rerun_higher_source_jets_GB_or_denominator_unit_assumption':True,
      'input_sha256':{str(p.relative_to(ROOT)):sha(p) for p in paths},'evidence_sha256':{p.name:sha(p)},
      'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    (out/'report.json').write_text(json.dumps(rep,indent=2)+'\n');print(json.dumps({k:rep[k] for k in ['all_passed','own_quotient_linear5_minor','seconds','max_rss_kib']}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
