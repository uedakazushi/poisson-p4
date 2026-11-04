#!/usr/bin/env python3
"""Add the actual new cubic initial relation; homogeneous length upper bound."""
import argparse,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
from ambient_B18_T_known_initial_module_length_upper24_20261104_v63 import sha,enc,add,order,leading,divides,exponent,reduce
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_T_KNOWN_INITIAL_LENGTH_UPPER17_20261104_v63'
def dump(out,name,x):(out/name).write_text(json.dumps(x,separators=(',',':'))+'\n')
def main(out):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
    out=out.resolve();out.mkdir(parents=True,exist_ok=False)
    paths={
      'old':'reports/AMBIENT_B18_T_KNOWN_INITIAL_MODULE_LENGTH_UPPER24_20261104_v63/known-initial-module-and-length-bound.json',
      'oldpeer':'reports/AMBIENT_B18_T_KNOWN_INITIAL_MODULE_LENGTH_UPPER24_20261104_v63_REVIEW/report.json',
      'values':'reports/AMBIENT_B18_T_LINEAR_SYZYGY_CUBIC_VALUES_20261104_v63/actual-linear-syzygy-cubic-source-values.json',
      'valuesreport':'reports/AMBIENT_B18_T_LINEAR_SYZYGY_CUBIC_VALUES_20261104_v63/report.json',
      'Pi2':'reports/AMBIENT_B18_T_LINEAR_SYZYGY_CUBIC_ADJOINT_20261104_v63/phase1-cubic-Pi2-and-source-adjoint.json',
      'closure':'reports/AMBIENT_B18_T_KNOWN_QUADRATIC_CUBIC_CLOSURE_20261104_v63/actual-known-quadratic-cubic-closure-certificate.json',
      'rankK':'reports/AMBIENT_B18_CRITICAL_CONTACT_AND_COMPARISON_RANK4_20261104_v63/report.json'}
    hashes={p:sha(ROOT/p) for p in paths.values()};dump(out,'phase.json',{'phase':'frozen_source_inputs_before_small_module_check','input_sha256':hashes})
    D={k:json.loads((ROOT/p).read_text()) for k,p in paths.items()};O,OP,V,VR,PI,C,K=[D[k] for k in paths]
    assert OP['all_passed'] and VR['all_passed'] and K['rank_D_K']==4
    assert OP['actual_I_squared_ann_m10_and_maximal_ideal_binding_accepted']
    assert V['actual_new_cubic_relation_channel_rank']==1
    basis=V['complete_fixed_rank_certificate'];j=basis['selected_original_columns'][0];a=F(basis['selected_minor_determinant'])
    assert j==71 and a==F(-251099375,1051776)
    actual=list(map(F,V['all119_two_actual_cubic_relation_values'][j]['two_actual_cubic_relation_values']));assert actual==[a,F(0)]
    qC={(4,exponent([5,7,7])):F(1)}
    # Cubic standard row3 is exactly e4*z5*r², with quotient coordinates(1,0).
    dual=[{i:F(c) for i,c in row} for row in PI['two_quotient_dual_rows600']]
    mons3=PI['monomials3_120'];index=4*120+mons3.index([5,7,7]);assert [d.get(index,F(0)) for d in dual]==[F(1),F(0)]
    gs=[{(t,tuple(m)):F(c) for t,m,c in g} for g in O['homogeneous_generators39']]+[qC]
    assert len(gs)==40
    dump(out,'preflight.json',{'module_generators40_QQ_terms':sum(map(len,gs)),
      'all_pair_upper_bound':40*39//2,'degree_upper_bound_for_new_S_pairs':4,'coefficient_trace_cap':150000,
      'no_new_source_X3_or_jet':True})
    pairs=[];trace_support=0
    for i in range(40):
        for j0 in range(i):
            p,q=leading(gs[i]),leading(gs[j0])
            if p[0]!=q[0]:continue
            lcm=tuple(max(u,v) for u,v in zip(p[1],q[1]));s={}
            add(s,gs[i],1/gs[i][p],tuple(u-v for u,v in zip(lcm,p[1])))
            add(s,gs[j0],-1/gs[j0][q],tuple(u-v for u,v in zip(lcm,q[1])))
            rem,tr=reduce(s,gs);assert not rem,(i,j0,enc(rem))
            trace_support+=len(tr);assert trace_support<150000
            pairs.append({'i':i,'j':j0,'lcm':list(lcm),'zero_reduction_trace':tr})
    # The leading terms give all-degree standard monomials directly.
    assert {leading(g) for g in gs[:37]}=={(t,exponent([v])) for t in range(5) for v in range(8)}-{(3,exponent([7])),(4,exponent([5])),(4,exponent([7]))}
    assert {leading(g) for g in gs[37:]}=={(4,exponent([5,5])),(3,exponent([7,7])),(4,exponent([5,7,7]))}
    hf=[5,3,2]+[1]*7;assert sum(hf)==17
    proof={'homogeneous_module_generators40':list(map(enc,gs)),
      'new_actual_cubic_source':{'linear_syzygy119_index':j,'nonzero_Pi2_value':str(a),
       'normalization':str(1/a),'canonical_initial_relation':enc(qC),
       'source_existence':'actual syzygy cubic/a minus e4*z5*r² lies in the known quadratic-relation-times-variables space, whose kernel is Pi2'},
      'all_same_position_S_pairs_exact_zero':pairs,'S_pair_count':len(pairs),'trace_QQ_support':trace_support,
      'standard_degree0_count':5,'standard_degree1':['e3*r','e4*z5','e4*r'],
      'standard_degree2':['e4*z5*r','e4*r²'],'standard_degree_at_least3':['e4*r^d'],
      'known_initial_submodule_quotient_Hilbert_series':'5+3t+2t²+t³/(1-t)',
      'known_initial_submodule_NOT_claimed_full_initial_module':True,
      'actual_m10_T_zero_from_accepted_source_base':True,'actual_T_layer_upper_bounds_degrees0_through9':hf,
      'actual_T_length_interval':[10,17],'actual_mu_D_K_interval':[6,13],'actual_mu_D_torsion_K_interval':[2,9],
      'actual_T_length_exact':None,'actual_HF3_exact':None,'actual_H0_length_or_nilcube_computed':False,
      'actual_new_cubic_source_independent_review_status':'PENDING'}
    dump(out,'known-initial40-and-length-upper17.json',proof)
    dump(out,'phase.json',{'phase':'FINAL','length_T_bounds':[10,17],'source_cubic_peer':'PENDING'})
    note=ROOT/'research/notes'/f'{STEM}.md'
    report={'id':STEM,'all_passed':True,'known_homogeneous_module_Groebner_basis_size':40,
      'all_S_pair_count':len(pairs),'actual_T_length_lower_bound':10,'actual_T_length_upper_bound':17,
      'actual_T_length_exact':None,'actual_HF3_or_nilcube_claimed':False,
      'actual_cubic_source_peer_pending':True,'input_sha256':hashes,
      'code_sha256':sha(Path(__file__)),'module_arithmetic_helper_sha256':sha(ROOT/'code/ambient_B18_T_known_initial_module_length_upper24_20261104_v63.py'),
      'note_sha256':sha(note),'evidence_sha256':{p.name:sha(p) for p in out.iterdir()},
      'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
    dump(out,'report.json',report);print(json.dumps({k:report[k] for k in ['id','all_passed','all_S_pair_count','actual_T_length_upper_bound','seconds']}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);main(p.parse_args().output)
