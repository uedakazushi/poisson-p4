#!/usr/bin/env python3
"""Fixed seven monomial reductions and pure uniserial-tail consequence."""
import argparse
import hashlib
import json
from fractions import Fraction as Q
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_T_UNISERIAL_DEEP_TAIL_20261104_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def decode(g):return {(p,tuple(e)):Q(c)for p,e,c in g}
def leading(f):return max(f,key=lambda z:(-z[0],z[1]))
def add_term(f,p,e,c):
    key=(p,e);f[key]=f.get(key,Q(0))+c
    if not f[key]:del f[key]
def main(out):
    out.mkdir(parents=True,exist_ok=False)
    files={
      'known40':'reports/AMBIENT_B18_T_KNOWN_INITIAL_LENGTH_UPPER17_20261104_v63/known-initial40-and-length-upper17.json',
      'primary':'reports/AMBIENT_B18_T_KNOWN_INITIAL_LENGTH_UPPER17_20261104_v63/report.json',
      'peer':'reports/AMBIENT_B18_T_KNOWN_INITIAL_LENGTH_UPPER17_20261104_v63_REVIEW/report.json',
      'cubic_peer':'reports/AMBIENT_B18_T_LINEAR_SYZYGY_CUBIC_VALUES_20261104_v63_REVIEW/report.json',
      'HF2_peer':'reports/AMBIENT_B18_T_ACTUAL_PROJECTED_QUADRATIC_RANK_20261104_v63_REVIEW/report.json',
      'minimal5_peer':'reports/AMBIENT_B18_T_MINIMAL5_SCHUR_OPERATOR_20261104_v63_REVIEW/report.json',
      'm10_peer':'reports/AMBIENT_B18_CONDUCTOR_SQUARE_CLOSED348_20261104_v63_REVIEW/report.json',
    }
    expected={
      'primary':'616f82154894153a4d4926a9329202c4e0f68d2b327422033d3cdd6e75816e87',
      'peer':'5bedafd8efeed3dee33bc43321a78d4caf5a83869ee033230c38b2ed070ab0bd',
      'cubic_peer':'f830eaa73ce76134bcdd9bec9170bad1e3378550f86711ce320b68449ec4ed41',
      'HF2_peer':'143af8f9218d2e43626f58db9abd34827dee9e419f6ec4a9626b5ca5d36caae1',
      'minimal5_peer':'ca127c4dc0f3c12a4d4d80dfc1ca316f42dee8b90ecc74cb3e5e423752c906b7',
      'm10_peer':'b937f7f936645761081ffd4e5b8663d095806f19ce14a771b5ae0e6b204a4108',
    }
    for k,h in expected.items():assert sha(ROOT/files[k])==h,k
    for k in ['primary','peer','cubic_peer','HF2_peer','minimal5_peer','m10_peer']:assert read(ROOT/files[k])['all_passed']
    data=read(ROOT/files['known40']);gens=[decode(g)for g in data['homogeneous_module_generators40']]
    assert len(gens)==40 and data['actual_m10_T_zero_from_accepted_source_base']
    assert data['standard_degree_at_least3']==['e4*r^d']
    assert read(ROOT/files['HF2_peer'])['actual_HF_through_degree2']==[5,3,2]
    # Independent normal-form division, NO GB search or S-pair rerun.
    lt=[leading(g)for g in gens]
    def nf(f):
        f=dict(f);rem={};trace=[]
        while f:
            pos,e=leading(f);c=f[(pos,e)]
            for j,(gp,ge)in enumerate(lt):
                if gp==pos and all(a>=b for a,b in zip(e,ge)):
                    m=tuple(a-b for a,b in zip(e,ge));q=c/gens[j][(gp,ge)]
                    for (p,x),v in gens[j].items():add_term(f,p,tuple(a+b for a,b in zip(x,m)),-q*v)
                    trace.append([j,list(m),str(q)]);break
            else:
                rem[(pos,e)]=c;del f[(pos,e)]
        return rem,trace
    receipts=[]
    for i in range(7):
        e=[0]*8;e[i]=1;e[7]=3
        f={(4,tuple(e)):Q(1)};rem,trace=nf(f);assert not rem
        back={}
        for j,m,c in trace:
            for (p,x),v in gens[j].items():add_term(back,p,tuple(a+b for a,b in zip(x,m)),Q(c)*v)
        assert back==f
        receipts.append({'variable_index':i,'input':[4,e,'1'],'zero_normal_form':True,'exact_known40_cofactors':trace})
    pure=[]
    for degree in range(3,10):
        e=(0,)*7+(degree,);f={(4,e):Q(1)};rem,tr=nf(f)
        assert rem==f and not tr;pure.append([4,list(e),'1'])
    proof={
      'actual_coefficient_ring':'S=R0/(ell)=k[[zA0,...,zA6,rA]]',
      'actual_generator_e4':'last minimal five target, source V56-row55',
      'r_not_DVR_parameter_ell':True,
      'all7_zi_e4_r3_known_initial_zero_reductions':receipts,
      'all_pure_r3_through_r9_standard_in_known_upper_module':pure,
      'actual_low_HF_exact':[5,3,2],
      'known_upper_basis_not_claimed_whole_initial_module':True,
      'M_definition':'M=m^3*T, g=rA^3*e4',
      'M_generated_by_g_over_S_by_Nakayama':True,
      'mM_equals_rM_by_Nakayama_on_mM_mod_rM':True,
      'all_m_power_d_ge3':'m^d*T=rA^(d-3)*M',
      'M_is_cyclic_k_r_module':'k[[rA]]/(rA^s), 0<=s<=7, s=0 means M=0',
      'M_uniserial_as_S_module':True,
      'all_zi_actions_unique_mod_r_power_s':'multiplication by hi(rA) in rA^2*k[[rA]]/(rA^s)',
      'exact_Hilbert_series_form':'5+3t+2t^2+t^3*(1-t^s)/(1-t)',
      'exact_length_formula':'length(T)=10+s, 0<=s<=7',
      'remaining_actual_gate':'annihilator of actual g=rA^3*e4 under rA',
      'actual_s_hi_coefficients_H0_length_DVR_torsion_exponents_or_nilcube_computed':False,
      'new_source_jet_point_GB_search_or_matrix_expansion':False,
    }
    ep=out/'seven-reductions-and-uniserial-tail-proof.json';ep.write_text(json.dumps(proof,indent=2)+'\n')
    rels=list(files.values())+[f'research/notes/{STEM}.md',str(Path(__file__).relative_to(ROOT))]
    report={
      'id':STEM,'all_passed':True,'status':'PROVED_actual_T_deep_tail_uniserial_one_integer_remaining',
      'all7_fixed_monomial_reductions_source_bound':True,
      'actual_T_deep_tail':'m^3*T is k[[rA]]/(rA^s), 0<=s<=7',
      'actual_T_full_Hilbert_function_form':'(5,3,2, then s consecutive ones, then zero)',
      'actual_T_total_length':'10+s',
      'all_zi_deep_actions_in_rA_squared':True,
      'rA_distinguished_from_DVR_ell':True,
      's_or_H0_DVR_exponents_or_nilcube_decided':False,
      'input_sha256':{p:sha(ROOT/p)for p in rels},'evidence_sha256':{ep.name:sha(ep)},
    }
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:report[k]for k in ['id','all_passed','actual_T_deep_tail','actual_T_total_length']}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);main(p.parse_args().output.resolve())
