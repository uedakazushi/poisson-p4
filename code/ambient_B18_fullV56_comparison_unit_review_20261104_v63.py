#!/usr/bin/env python3
"""Independent actual-H20, fixed24minor and typed regular comparison audit.

No Artin105 solve, source evaluation, new inverse/rank search or author worker.
"""
import argparse,hashlib,itertools,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
PRIMARY='AMBIENT_B18_FULLV56_COMPARISON_UNIT_20261104_v63';STEM=PRIMARY+'_REVIEW'
DIVISION='AMBIENT_B18_CBAR_COMPARISON_DIVISION_PREFLIGHT_20261104_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def rec(p):return {i:F(c) for i,c in p if F(c)}
def plus(*ps):
    r={}
    for p in ps:
        for i,c in p.items():r[i]=r.get(i,F(0))+c
    return {i:c for i,c in r.items() if c}
def scale(p,c):return {i:c*v for i,v in p.items() if c*v}
def combine(columns,v):return plus(*(scale(columns[i],c) for i,c in v.items()))
def r5mul(p,q):
    r={}
    for i,c in p.items():
        for j,d in q.items():
            if i+j<5:r[i+j]=r.get(i+j,F(0))+c*d
    return {i:c for i,c in r.items() if c}
def mat(M):return [[F(x) for x in row] for row in M]
def determinant(M):
    a=[r[:] for r in M];det=F(1)
    for j in range(len(a)):
        p=next(i for i in range(j,len(a)) if a[i][j])
        if p!=j:a[p],a[j]=a[j],a[p];det=-det
        t=a[j][j];det*=t
        for i in range(j+1,len(a)):
            if not a[i][j]:continue
            c=a[i][j]/t
            for k in range(j+1,len(a)):a[i][k]-=c*a[j][k]
            a[i][j]=0
    return det
def mm(A,B):
    out=[[F(0)]*len(B[0]) for _ in A]
    for i,row in enumerate(A):
        for k,c in enumerate(row):
            if c:
                for j,d in enumerate(B[k]):
                    if d:out[i][j]+=c*d
    return out

def main(out):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False)
    folder=ROOT/'reports'/PRIMARY;rp=folder/'report.json';A=read(rp)
    assert sha(rp)=='fd35519e085bea0afee317b83f0a58f886bd66cb22929bd608cc65f3d973b8ce' and A['all_passed']
    ep=ROOT/A['evidence']['path'];assert sha(ep)==A['evidence']['sha256'];E=read(ep)
    D={}
    for k,r in A['input_hashes'].items():
        p=ROOT/r['path'];assert sha(p)==r['sha256'];D[k]=read(p)
    assert E['input_hashes']==A['input_hashes']
    divrp=ROOT/'reports'/DIVISION/'report.json';div=read(divrp)
    assert sha(divrp)=='7e9518f84c67ddd0619b5d017af412de7bc73355a9560943fdce0b86b8fea929' and div['all_passed']
    divep=ROOT/div['evidence']['path'];assert sha(divep)==div['evidence']['sha256'];dv=read(divep)
    for r in div['input_hashes'].values():assert sha(ROOT/r['path'])==r['sha256']
    gates=[ROOT/'reports'/s/'report.json' for s in [
      'AMBIENT_B18_RAW_COKERNEL_ARTIN20_EVALUATION_20261104_v63_REVIEW',
      'AMBIENT_B18_CBAR_LITERAL_RESOLUTION_20261104_v63_REVIEW',
      'AMBIENT_B18_NORMALIZER_CHAIN_LIFT_640_20261104_v63_REVIEW',
      'AMBIENT_B18_SPECIALIZED_COMPARISON_T_PREFLIGHT_20261104_v63_REVIEW']]
    assert all(read(p)['all_passed'] for p in gates)
    (out/'phase.json').write_text(json.dumps({'phase':'all_source_hashes_and_prior_peers_bound'})+'\n')
    H,S,B,R=[D[k] for k in ['H20','source105','chainlift','resolution']]
    assert S['all105_equations_exactly_zero'] and B['all640_UP_equals_QB_verified']
    basis=[tuple(m) for m in H['quotient_basis_monomials']];assert len(basis)==20
    actions=[[rec(c) for c in M] for M in H['multiplication8_column_sparse']]
    # A different multiplication construction from the author's495 reductions:
    # act with each literal basis monomial on every other basis vector.
    table={}
    for i,m in enumerate(basis):
        for j in range(20):
            p={j:F(1)}
            for k,n in enumerate(m):
                for _ in range(n):p=combine(actions[k],p)
            table[i,j]=p
    def amul(p,q):return plus(*(scale(table[i,j],c*d) for i,c in p.items() for j,d in q.items()))
    cs=list(map(F,R['CB_shift_coefficients7']));maps=[]
    for shifted in [False,True]:
        columns=[]
        for m in basis:
            degree=m[7]+2*sum(m[:7]);c=F(1)
            for i,n in enumerate(m[:7]):c*=(-cs[i] if shifted else F(0))**n
            columns.append({degree:c} if c and degree<5 else {})
        maps.append(columns)
    assert [[rec(v) for v in M] for M in E['H20_to_r5_projections20_column_sparse']]==maps
    for columns in maps:
        assert columns[0]=={0:F(1)}
        for i,j in itertools.product(range(20),repeat=2):
            assert combine(columns,table[i,j])==r5mul(columns[i],columns[j])
    # Bind every finite field value to the accepted SAME original source roots.
    sourcepath=ROOT/B['scalar_source_artifacts']['source']['path'];source=read(sourcepath)
    assert sha(sourcepath)==B['scalar_source_artifacts']['source']['sha256']
    roots={i:rec(p) for i,p in S['all_requested_original_source_root_values']}
    nus=[source['families'][f]['actual_tracefree_nu'] for f in ['CA','CB']]
    assert S['nu_roots_CA_CB']==nus
    nuA,nuB=[rec(p) for p in S['nu_CA_CB_values']]
    assert [roots[i] for i in nus]==[nuA,nuB]
    def cmul(v,w):
        result=[{} for _ in range(4)]
        for i,x in enumerate(v):
            for j,y in enumerate(w):
                p=amul(x,y)
                if i&j&1:p=amul(p,nuA)
                if i&j&2:p=amul(p,nuB)
                result[i^j]=plus(result[i^j],p)
        return result
    V=[];labels=[]
    for j,r in enumerate(S['all12_physical_coefficients']):
        row=r['physical_row'];a,b=source['families']['CA']['complete_normalized_physical_w_coefficients'][0][row]
        c,d=source['families']['CB']['complete_normalized_physical_w_coefficients'][0][row]
        assert r['roots']=={'aA_node':a,'aB_node':c,'bA_node':b,'bB_node':d}
        assert [rec(r[k]) for k in ['aA','aB','bA','bB']]==[roots[i] for i in [a,c,b,d]]
        delta=plus(roots[a],scale(roots[c],-1));assert rec(r['delta'])==delta
        f=[delta,roots[b],scale(roots[d],-1),{}]
        for k in range(4):
            e=[{} for _ in range(4)];e[k]={0:F(1)};V.append(cmul(f,e)[1:])
            labels.append({'Q48_column':4*j+k,'physical_row':row,'tensor_basis':k,'exact_H20_value':True})
    for i in range(7):
        assert basis[i+1]==tuple(int(k==i) for k in range(8))
        V.append([{i+1:F(-1)},{},{}]);labels.append({'negative_U8_column':i,'actual_function':'-zA'+str(i),'exact_H20_value':True})
    r3=basis.index((0,0,0,0,0,0,0,3));assert H['gA_after_ell']==[[list(basis[r3]),'125/24']]
    V.append([{r3:F(-125,24)},{},{}]);labels.append({'negative_U8_column':7,'actual_function':'-gA_plane(r,r)',
      'exact_H20_value':False,'m0_relevant_representative':'-(125/24)r3','discarded_order4_tail_does_not_affect_any_m0_slot':True})
    assert E['V56_column_order']==labels
    assert [[rec(p) for p in c] for c in E['V56_H20_values_or_explicit_m0_representative']]==V
    mcols=[];curve=[]
    for v in V:
        c=[];rr=[]
        for s,n in enumerate([3,3,2]):
            p=combine(maps[int(s==1)],v[s]);assert all(not p.get(k,F(0)) for k in range(n))
            c.extend(v[s].get(i,F(0)) for i in range(1,8));c.append(p.get(n,F(0)));rr.append(p)
        mcols.append(c);curve.append(rr)
    assert [[rec(p) for p in c] for c in E['all56_sector_curve_restrictions']]==curve
    m0=[list(c) for c in zip(*mcols)];assert m0==mat(E['all24_by56_m0'])
    selected=E['fixed_selected24_columns'];assert len(selected)==len(set(selected))==24
    M=[[row[j] for j in selected] for row in m0];assert M==mat(E['fixed24_minor_matrix'])
    det=determinant(M);assert det==F(-397267578125,50096498540544)==F(E['fixed24_minor_determinant'])
    inverse=mat(E['fixed24_minor_inverse']);identity=[[F(i==j) for j in range(24)] for i in range(24)]
    assert mm(M,inverse)==mm(inverse,M)==identity
    factors=mat(E['all56_column_factorization']);assert mm(M,factors)==m0
    # Check every typed map argument and node signature against actual chain refs;
    # no scalar Taylor coefficient other than the proven center is invented.
    ops=E['typed_regular_function_operators'];assert len(ops)==2502
    zero=ops.index(['zero']);functions=E['actual_V56_function_operator_refs'];expected=[]
    for c in B['exact_existing_Q48_C4']:
        expected.append([['actual_chain_scalar_restrict_ell0',i] if i else ['zero'] for i in c[1:]])
    for c in B['exact_existing_U9_C4'][1:]:
        expected.append([['negative_actual_chain_scalar_restrict_ell0',i] if i else ['zero'] for i in c[1:]])
    assert [[[ops[i] for i in c] for c in functions]][0]==expected
    typedm=E['typed_m24_by56'];assert len(typedm)==24 and all(len(row)==56 for row in typedm)
    checkedslots=0
    for j,f in enumerate(functions):
        for s,root in enumerate(f):
            for k in range(8):
                node=ops[typedm[8*s+k][j]]
                if root==zero:assert node==['zero'];continue
                if k<7:
                    tag,up,low,ss,kk=node
                    assert tag=='regular_coordinate_Hadamard' and (ss,kk)==(s,k)
                    assert ops[up]==['adapted_sector_prefix',root,s,k+1]
                    assert ops[low]==['adapted_sector_prefix',root,s,k]
                else:
                    tag,powerroot,n,why=node
                    assert tag=='exact_monomial_r_division' and n==[3,3,2][s]
                    assert ops[powerroot]==['adapted_sector_prefix',root,s,0]
                    assert why=='actual_sector_annihilator_membership_not_numerical_zeroing'
                checkedslots+=1
    u=E['typed_unit24_inverse_node'];selected_typed=[[row[j] for j in selected] for row in typedm]
    assert ops[u]==['formal_unit_matrix_inverse',selected_typed,'actual_center_matrix_inverse',E['fixed24_minor_inverse']]
    n=E['typed_n56_by24'];assert len(n)==56
    position={j:i for i,j in enumerate(selected)}
    for j,row in enumerate(n):
        assert len(row)==24
        for k,i in enumerate(row):
            assert ops[i]==(['unit_matrix_inverse_entry',u,position[j],k] if j in position else ['zero'])
    assert E['conditional_T_source_shape_before_freezero16_cancellation']==[84+(56-24),168+768]==[116,936]
    assert dv['surjective_m_and_exact_cancellation']['T_complete_conditional_shape']==[100,920]
    assert 116-16==100 and 936-16==920
    assert not E['scalar_m_n_coefficients_or_348_T_or_K_computed']
    ev={'all800_H20_curve_map_products_verified_from_independent_coordinate_action_table':True,
      'all48_Q_values_bound_to_actual105_original_source_IDs':True,
      'all56_V_columns_m0_and_sector_orders_verified':True,
      'last_U_gA_unknown_order4_tail_retained':True,'last_U_only_representative_for_center':True,
      'fixed24_minor_determinant':str(det),'both576_unit_inverse_products_and_all1344_factor_entries_verified':True,
      'typed_operator_nodes':len(ops),'all_typed_m_n_nodes_and_actual_chain_refs_verified':True,
      'nonzero_m_operator_slots_checked':checkedslots,
      'regular_operator_proof':'adapted prefix differences telescope in each sector; last prefix is divisible by r^n from exact actual Cbar annihilation; all are regular compositions before restriction',
      'unit_comparison':'M has the displayed nonzero actual center determinant; n=E_S M^-1, mn=I24, Vn=d1,h=0',
      'cycle_contraction':'remove last component via seven coordinate-Hadamards after domain cancellation of r^n; remaining ordinary coordinate Koszul relation has its canonical contraction',
      'conditional_full_T_cokernel':[116,936],'conditional_proper_T_cokernel':[100,920],
      'actual768_cycle_lifts_Tlength_Ktorsion_not_computed':True}
    out_ev=out/'independent-H20-m0-and-regular-operator-checks.json';out_ev.write_text(json.dumps(ev,indent=2)+'\n')
    note=ROOT/'research/notes'/f'{STEM}.md';files=[rp,ep,divrp,divep,Path(__file__),note,sourcepath]+gates
    files += [ROOT/r['path'] for r in A['input_hashes'].values()]+[ROOT/r['path'] for r in div['input_hashes'].values()]
    report={'id':STEM,'all_passed':True,'actual_fullV56_m0_rank':24,'fixed24_minor':str(det),
      'independent_actual_source_table_and_fixed_minor_verification':True,'regular_typed_m_n_and_h_zero_accepted':True,
      'division_preflight_accepted':True,'gA_order4_tail_not_zeroed':True,
      'actual_cycle_lifts_Tlength_Ktorsion_D0_depth_computed':False,
      'no_authorworker_source105_solve_newpoint_jet_ranksearch_or_scalar_inverse_search':True,
      'input_sha256':{str(p.relative_to(ROOT)):sha(p) for p in files},'evidence_sha256':{out_ev.name:sha(out_ev),'phase.json':sha(out/'phase.json')},
      'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:report[k] for k in ['id','all_passed','seconds','max_rss_kib']}))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
