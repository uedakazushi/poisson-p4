#!/usr/bin/env python3
"""Actual first640 X columns modulo m² from frozen scalar-jet caches only.
No new source/implicit solve, division, point, source jet order, or Gaussian.
"""
import argparse,hashlib,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_COMPARISON_X_FIRST640_LINEAR_20261104_v63'
ZERO=(F(0),)*9
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def add(*ps):return tuple(sum(xs,F(0)) for xs in zip(*ps)) if ps else ZERO
def scale(p,c):return tuple(c*x for x in p)
def mul(p,q):return (p[0]*q[0],)+tuple(p[0]*q[i]+q[0]*p[i] for i in range(1,9))
def const(q):return (F(q),)+(F(0),)*8
def full9_to_bar(p):
    assert len(p)==10
    return tuple(F(x) for x in p[:8])+(F(p[8])+F(p[9]),)
def encode(p):return [[i,str(c)] for i,c in enumerate(p) if c]

def main(out):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out=out.resolve();out.mkdir(parents=True,exist_ok=False)
    paths={
      'chain':'reports/AMBIENT_B18_NORMALIZER_CHAIN_LIFT_640_NEGATIVE_BOUNDARY_20261104_v63/actual-UP-QB-640-typed-operator.json',
      'C0':'reports/AMBIENT_B18_T0_HONEST_C0_20261104_v63/honest-C0-D0-D4-DAG.json',
      'firstjets':'reports/AMBIENT_B18_T3_FULL9_FIRSTJET_20261104_v63/full9-firstjet-certificate.json',
      'quartic':'reports/AMBIENT_B18_T0_QUOTIENT_FIRSTS_QUARTIC_20261104_v63/quotient-firsts-and-actual-D0-firstjets.json',
      'source105':'reports/AMBIENT_B18_RAW_COKERNEL_ARTIN20_EVALUATION_20261104_v63/actual-source105-Artin20-values.json',
      'origin':'reports/AMBIENT_B18_COMPARISON_X_ORIGIN_20261104_v63/actual-X0-and-unit-cancellation.json',
      'H20':'reports/AMBIENT_B18_CONDUCTOR_PRODUCT_SPECIAL_FIBER20_20261104_v63/actual-special-fiber20-coefficients.json'}
    d={k:read(ROOT/p) for k,p in paths.items()};B,C0,J,Q,S,X0,H=[d[k] for k in paths]
    ops=B['operator_scalar_nodes'];columns=B['B48_by640_sparse_columns'];roots=[j for col in columns for i,j in col]
    needed=set();todo=roots[:]
    while todo:
        i=todo.pop()
        if i in needed:continue
        needed.add(i);n=ops[i]
        if n[0]=='sum':todo.extend(j for j,a,b in n[1])
        elif n[0]=='product':todo.extend(n[3])
    atomids=sorted({ops[i][1] for i in needed if ops[i][0]=='actual_C0_scalar_function'})
    Praw=Q['all640_exact_D0_firstjet81_columns'];assert len(Praw)==640
    Pnz=sum(x!='0' for col in Praw for x in col)
    bound=9*(len(needed)+len(roots)+len(atomids))+Pnz
    assert len(needed)==1302 and len(roots)==2808 and len(atomids)==79
    assert bound<150000
    preflight={'typed_needed_nodes':len(needed),'needed_atom_count':len(atomids),'existing_B_scalar_roots':len(roots),
      'max_typed_all9slots':9*len(needed),'max_B_output_all9slots':9*len(roots),
      'existing_P_full9_derivative_nonzeros':Pnz,'known_stage_total_slot_upper_bound':bound,
      'original_source_translation_dynamic_total_term_guard':150000,
      'no_dense_56_by640_by9_array':True,'no_native128_derivative_queries':True}
    (out/'preflight.json').write_text(json.dumps(preflight,indent=2)+'\n')
    (out/'phase.json').write_text(json.dumps({'phase':'sparse_size_guard_passed_before_dual_allocation'})+'\n')
    # Scalar derivative cache, always in ACTUAL CA F9, then ell: s=r.
    values={int(i):full9_to_bar(p) for key in ['all792_CA_h_firstjets','all784_CB_h_firstjets_in_YA'] for i,p in J[key].items()}
    basis=H['quotient_basis_monomials'];assert S['A20_basis']==basis
    def artin_first(p):
        data={i:F(c) for i,c in p};v=[data.get(0,F(0))]
        for j in range(8):
            mon=[0]*8;mon[j]=1;v.append(data.get(basis.index(mon),F(0)))
        return tuple(v)
    values.update({i:artin_first(p) for i,p in S['all_requested_original_source_root_values']})
    for i,p in enumerate(S['actual_physical_R0_readback'],2):values[i]=artin_first(p)
    for j,i in enumerate(C0['actual_C0_basis8_scalar_nodes']):
        v=[F(0)]*9
        if j<7:v[j+1]=F(1)
        values[i]=tuple(v)
    monindex={tuple(m):i for i,m in enumerate(J['source2_monomials'])}
    prefixes=[]
    for rec in C0['all35_regular_coordinate_divisions']:
        f=rec['original_pure_CA_source_function'];p=J['both57_source2']['CA']['pure_function2jets'][str(f)]
        for j,r in enumerate(rec['all7_regular_coordinate_divisions']):
            full=[F(p[monindex[(j,)]])]
            for k in range(9):
                full.append(F(p[monindex[tuple(sorted((j,k)))]])) if (k<=j or k>=7) else full.append(F(0))
            v=full9_to_bar(full);assert v[0]==F(r['origin']);values[r['node']]=v
            prefixes.append({'node':r['node'],'actual_source_function':f,'prefix_index':j,'first9_bar':encode(v),
              'retained_z_indices':list(range(j+1)),'always_retained_r_s':True})
    # Independent source4 packet proved all8 quotient origins0 and gradients;
    # check its numerator starts degree4, then reuse the whole9 derivatives.
    gnodes=Q['original_regular_function_nodes_unmodified']['old_CB']+Q['original_regular_function_nodes_unmodified']['new_CA']
    for family in ['CB','CA']:
        for p in Q['families'][family]['four_numerator_all15_coefficients']:assert all(F(c)==0 for c in p[:10])
    for i,gradient in zip(gnodes,Q['all8_regular_function_first_gradients_in_FULL9_CA_F9']):values[i]=full9_to_bar([F(0)]+list(map(F,gradient)))
    gb=next(r['gB'] for r in B['accepted_actual_regular_function_identities'] if r['kind']=='accepted_actual_CB_f63_tail');values[gb]=ZERO
    full=C0['DAG_scalar_nodes'];memo=dict(values);initial=set(memo)
    def original(i):
        if i in memo:return memo[i]
        n=full[i]
        if n[0]=='c':p=const(F(n[1],n[2]))
        elif n[0]=='a':p=add(original(n[1]),original(n[2]))
        elif n[0]=='n':p=scale(original(n[1]),-1)
        elif n[0]=='m':
            a=original(n[1]);p=ZERO if a==ZERO else mul(a,original(n[2]))
        else:raise AssertionError(('MISSING_ACTUAL_REGULAR_FIRST',i,n))
        memo[i]=p;assert 9*len(memo)<150000;return p
    atom_values={i:original(i) for i in atomids}
    origins={i:F(c) for i,c in X0['source_function_constant_refs']}
    assert all(v[0]==origins[i] for i,v in atom_values.items())
    (out/'phase1-actual-atom-firsts.json').write_text(json.dumps({'actual_variable_order':['zA0','zA1','zA2','zA3','zA4','zA5','zA6','rA'],
      'all79_atom_values':[[i,encode(v)] for i,v in sorted(atom_values.items())],
      'all35_ordered_CA_prefix_firsts':prefixes,'all8_g_quotient_overrides':[[i,encode(values[i])] for i in gnodes],
      'all79_origins_match_previous_X0':True},separators=(',',':'))+'\n')
    typed={}
    for i in sorted(needed):
        n=ops[i]
        if n[0]=='c':v=const(F(n[1],n[2]))
        elif n[0]=='actual_C0_scalar_function':v=atom_values[n[1]]
        elif n[0]=='sum':v=add(*(scale(typed[j],F(a,b)) for j,a,b in n[1]))
        elif n[0]=='product':
            v=const(F(n[1],n[2]))
            for j in n[3]:v=mul(v,typed[j])
        else:raise AssertionError(n)
        typed[i]=v
    outcols=[];nonzeros=0
    for j,saved in enumerate(columns):
        entries={i:typed[k] for i,k in saved if typed[k]!=ZERO}
        # P's row0 is removed by the exact unit-target cancellation; retain8.
        for row in range(8):
            q=Praw[j][9*(row+1):9*(row+2)];assert len(q)==9
            v=full9_to_bar([F(0)]+list(map(F,q)))
            if v!=ZERO:entries[48+row]=v
        for row in range(56):assert entries.get(row,ZERO)[0]==F(X0['full_X0_56_by768'][row][j])
        outcols.append([[i,encode(v)] for i,v in sorted(entries.items())]);nonzeros+=sum(sum(bool(x) for x in v) for v in entries.values())
    typednz=sum(sum(bool(x) for x in v) for v in typed.values())
    sourcenz=sum(sum(bool(x) for x in memo[i]) for i in memo if i not in initial)
    atomnz=sum(sum(bool(x) for x in v) for v in atom_values.values())
    total=typednz+sourcenz+atomnz+nonzeros;assert total<150000
    ev={'coordinate_order':['constant','zA0','zA1','zA2','zA3','zA4','zA5','zA6','rA'],
      'source_columns640_rows56_sparse_dual9':outcols,'all35840_origins_match_saved_X0':True,
      'all_existing_P8_derivatives_restricted_s_equals_r':True,'native128_X1_not_evaluated':True,
      'no_fresh_r_g_division_or_source_solve':True,'actual_nonzero_QQ_coefficients':nonzeros,
      'expanded_live_stage_coefficient_bound_checked':total,'all_source_cache_refs_pinned':True}
    ep=out/'actual-first640-X-mod-m2.json';ep.write_text(json.dumps(ev,separators=(',',':'))+'\n')
    (out/'phase.json').write_text(json.dumps({'phase':'FINAL_all640_source_firsts_completed_native128_pending','coefficient_nonzeros':nonzeros})+'\n')
    note=ROOT/'research/notes'/f'{STEM}.md'
    report={'id':STEM,'all_passed':True,'actual_first640_X1_source_bound':True,'actual640_X0_columns_match':True,
      'all35_ordered_prefix_firsts_and8_g_source4_overrides_bound':True,
      'no_new_source_implicit_solve_jetorder_point_GB_or_Gaussian':True,
      'native128_X1_T_mu_length_K_torsion_D0_depth_computed':False,
      'source_atom_count':len(atomids),'typed_closure_nodes':len(needed),'nonzero_output_coefficients':nonzeros,
      'expanded_live_coefficient_total':total,'input_sha256':{p:sha(ROOT/p) for p in paths.values()},
      'code_sha256':sha(Path(__file__)),'note_sha256':sha(note),
      'evidence_sha256':{p.name:sha(p) for p in out.iterdir()},
      'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:report[k] for k in ['id','all_passed','nonzero_output_coefficients','expanded_live_coefficient_total','seconds','max_rss_kib']}))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
