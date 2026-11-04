#!/usr/bin/env python3
"""Sufficient cubic annihilation test using only frozen true quadratic forms.

Linear rewrites use a compatible module order. All overlap remainders are
retained via the full known degree-two space; the37 forms are not assumed GB.
"""
import argparse,itertools,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
from ambient_B18_comparison_T_constant_rank_20261104_v63 import read,sha,rec,encode,addscaled,sparse_rank
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_T_KNOWN_QUADRATIC_CUBIC_CLOSURE_20261104_v63'
def dump(out,name,x):(out/name).write_text(json.dumps(x,separators=(',',':'))+'\n')
def main(out):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
    out=out.resolve();out.mkdir(parents=True,exist_ok=False)
    paths={
      'low':'reports/AMBIENT_B18_COMPARISON_T_LINEAR_QUADRATIC_20261104_v63/actual-linear-quadratic-T-certificate.json',
      'low_peer':'reports/AMBIENT_B18_COMPARISON_T_LINEAR_QUADRATIC_20261104_v63_REVIEW/report.json',
      'rank1':'reports/AMBIENT_B18_T_ACTUAL_PROJECTED_QUADRATIC_RANK_20261104_v63/actual-three-by620-quadratic-image-certificate.json',
      'rank1_report':'reports/AMBIENT_B18_T_ACTUAL_PROJECTED_QUADRATIC_RANK_20261104_v63/report.json',
      'duals':'reports/AMBIENT_B18_COMPARISON_T_NEXT3_SOURCE_PREFLIGHT_20261104_v63/three-quadratic-duals.json',
      'complete5':'reports/AMBIENT_B18_T_MINIMAL5_SCHUR_OPERATOR_20261104_v63_REVIEW/report.json'}
    hashes={p:sha(ROOT/p) for p in paths.values()};dump(out,'phase.json',{'phase':'frozen_inputs_bound_before_closure','input_sha256':hashes})
    D={k:read(ROOT/p) for k,p in paths.items()};L,LP,R,RP,P,C=[D[k] for k in paths]
    assert LP['all_passed'] and C['all_passed'] and RP['all_passed']
    assert C['actual_complete_minimal_generator_presentation']==[5,825]
    assert L['complete_linear_rank_certificate']['rank']==37
    assert L['quadratic_certificate']['rank']==177 and R['actual_projected_quadratic_rank']==1
    assert R['actual_HF_through_degree2']==[5,3,2] and R['actual_length_T_mod_m3']==10
    linear=list(map(rec,L['complete_linear_rank_certificate']['basis_normalized_columns']))
    lead={min(v):j for j,v in enumerate(linear)};assert len(lead)==37
    assert all(v[i]==1 and all(k>=i for k in v) for i,j in lead.items() for v in [linear[j]])
    mons2=list(itertools.combinations_with_replacement(range(8),2));m2={m:i for i,m in enumerate(mons2)}
    mons3=list(itertools.combinations_with_replacement(range(8),3));m3={m:i for i,m in enumerate(mons3)};assert len(mons3)==120
    free1=[i for i in range(40) if i not in lead]
    standard=[120*t+j for t in range(5) for j,m in enumerate(mons3) if all(8*t+v in free1 for v in m)]
    assert len(free1)==3 and len(standard)<=10
    sindex={i:j for j,i in enumerate(standard)}
    quads=list(map(rec,L['quadratic_certificate']['basis_normalized_columns']))
    # The actual rank-one relation has known image in the three-dimensional
    # quotient. The canonical representative is an actual initial relation:
    # subtract its difference in the true177 relation span. No full unknown
    # source quadratic coefficient is assigned a value by this argument.
    actual_normal=rec(R['complete_fixed_rank_certificate']['basis_normalized_columns'][0]);representatives=P['quotient_representative_rows']
    qnew={representatives[i]:c for i,c in actual_normal.items()};quads.append(qnew)
    duals=list(map(rec,P['projection3_by180_rows_sparse']))
    assert [sum((a*qnew.get(j,F(0)) for j,a in d.items()),F(0)) for d in duals]==[actual_normal.get(i,F(0)) for i in range(3)]
    pre={'linear_free_monomials':[[i//8,i%8] for i in free1],'cubic_standard_target_count':len(standard),
      'all600_normal_form_support_upper':600*len(standard),'known_quadratic_basis_count':178,
      'known_quadratic_input_QQ_nonzeros':sum(map(len,quads)),
      'cubic_relation_columns':178*8,'all_reduced_cubic_column_slots_upper':178*8*len(standard),
      'linear_overlap_check_columns':37*36,'cap':150000}
    pre['TOTAL_coefficient_bound']=600*len(standard)+178*8*len(standard)*3+sum(map(len,quads))*8+37*36*len(standard)+15000
    assert pre['TOTAL_coefficient_bound']<150000
    dump(out,'preflight.json',pre)
    memo={};steps={}
    def order(t,m):return (-t,)+tuple(m.count(i) for i in range(8))
    def nf(t,m):
        key=120*t+m3[m]
        if key in memo:return memo[key]
        candidates=[v for v in m if 8*t+v in lead]
        if not candidates:
            ans={sindex[key]:F(1)};step={'standard':True}
        else:
            v=min(candidates);p=lead[8*t+v];rest=list(m);rest.remove(v);ans={};children=[]
            for index,c in linear[p].items():
                if index==8*t+v:continue
                target,var=divmod(index,8);child=tuple(sorted(rest+[var]))
                assert order(target,child)<order(t,m),'module term order fails termination'
                addscaled(ans,nf(target,child),-c);children.append([120*target+m3[child],str(-c)])
            step={'linear_relation_basis_index':p,'leading_variable':v,'children':children}
        memo[key]=ans;steps[key]=step;return ans
    for t in range(5):
        for m in mons3:nf(t,m)
    assert len(memo)==600
    def reduce_cubic(col):
        v={}
        for i,c in col.items():addscaled(v,memo[i],c)
        return v
    cubic=[];labels=[];raw_support=0
    for j,q in enumerate(quads):
        for k in range(8):
            raw={}
            for i,c in q.items():
                t,r=divmod(i,36);idx=120*t+m3[tuple(sorted(mons2[r]+(k,)))];raw[idx]=raw.get(idx,F(0))+c
            raw_support+=len(raw);cubic.append(reduce_cubic(raw));labels.append([j,k])
    rank=sparse_rank(cubic,len(standard));rho=rank['rank'];full=rho==len(standard)
    # Do not assume the37 linear forms constitute a module Groebner basis.
    # Explicitly retain/check every linear overlap after the fixed rewrite.
    rb=list(map(rec,rank['basis_normalized_columns']));rlead={min(v):j for j,v in enumerate(rb)}
    overlaps=[]
    for p,l in enumerate(linear):
        for m in mons2:
            raw={}
            for i,c in l.items():
                t,v=divmod(i,8);idx=120*t+m3[tuple(sorted(m+(v,)))];raw[idx]=raw.get(idx,F(0))+c
            v=reduce_cubic(raw);original=dict(v);coeff={}
            while v:
                i=min(v);assert i in rlead,'linear overlap omitted from known quadratic closure'
                j=rlead[i];c=v[i];coeff[j]=coeff.get(j,F(0))+c;addscaled(v,rb[j],-c)
            overlaps.append({'linear_basis':p,'quadratic_multiplier':m,'normal_form':encode(original),'cubic_relation_span_coordinates':encode(coeff)})
    support=sum(map(len,memo.values()))+sum(map(len,cubic))+sum(map(len,rank['basis_normalized_columns']))+sum(map(len,rank['all_column_factorization_in_normalized_basis']))+raw_support
    assert support<150000
    dump(out,'linear-module-rewrite-certificate.json',{'module_order':'position-over-lex: smaller target is larger; z0>z1>...>r; fixed degree3',
      'monomials3_120':mons3,'standard_global600_positions':standard,
      'standard_target_monomials':[[i//120,mons3[i%120]] for i in standard],
      'linear_basis37':list(map(encode,linear)),
      'all600_rewrite_steps':[[i,steps[i],encode(memo[i])] for i in range(600)],
      'all1332_linear_overlap_remainders_retained_and_spanned':overlaps,
      'linear_basis_assumed_Groebner':False})
    dump(out,'actual-known-quadratic-cubic-closure-certificate.json',{'quadratic_basis178':list(map(encode,quads)),
      'actual_new_quad_canonical_lift_from_rank1_source':{'source620_column':R['complete_fixed_rank_certificate']['selected_original_columns'][0],
       'normalized_three_class_vector':encode(actual_normal),'canonical_quad_representative':encode(qnew),
       'source_existence_proof':'actual source quadratic minus canonical representative lies in the true177 initial relation space'},
      'cubic_column_labels':labels,'all1424_reduced_cubic_columns':list(map(encode,cubic)),
      'complete_reduced_rank_certificate':rank,'reduced_target_count':len(standard),'certified_cubic_rank':rho,
      'all600_cubic_monomials_in_certified_initial_relation_span':full,
      'm_cubed_T_zero_proved':full,'total_T_length_proved':10 if full else None,
      'actual_degree3_HF':0 if full else None,'degree3_HF_upper_bound_only':len(standard)-rho,
      'partial_span_does_not_give_a_T_length_lower_bound':True,
      'new_actual_X3_or_source_jet_computed':False,'expanded_live_QQ_support':support})
    dump(out,'phase.json',{'phase':'FINAL','standard_cubic_count':len(standard),'certified_rank':rho,'m3_T_zero':full,'total_T_length':10 if full else None})
    note=ROOT/'research/notes'/f'{STEM}.md'
    report={'id':STEM,'all_passed':True,'actual_known_linear_rank':37,'actual_known_quadratic_rank':178,
      'cubic_standard_target_count':len(standard),'known_cubic_closure_rank':rho,
      'm_cubed_T_zero_proved':full,'actual_total_T_length':10 if full else None,
      'actual_degree3_HF_if_proved':0 if full else None,'only_degree3_HF_upper_bound_if_partial':len(standard)-rho,
      'no_X3_new_source_jet_GB_point_or_348_matrix':True,'input_sha256':hashes,
      'code_sha256':sha(Path(__file__)),'note_sha256':sha(note),
      'evidence_sha256':{p.name:sha(p) for p in out.iterdir()},'seconds':time.monotonic()-start,
      'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
    dump(out,'report.json',report);print(json.dumps({k:report[k] for k in ['id','all_passed','cubic_standard_target_count','known_cubic_closure_rank','m_cubed_T_zero_proved','actual_total_T_length','seconds']}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
