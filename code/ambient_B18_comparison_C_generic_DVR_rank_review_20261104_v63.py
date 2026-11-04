#!/usr/bin/env python3
"""Independent fixed source-polynomial and conceptual comparison-rank review.

No author worker, source solve, jet, root search, factorization, point or GB.
Sparse Fraction polynomial operations come only from an earlier independent peer.
"""
import argparse, gzip, hashlib, json, resource, signal, time
from fractions import Fraction as F
from pathlib import Path
from ambient_B18_proper_core_five_lines_review_20261104_v63 import (
    parse, substitute, variable, add, mul, scale, power, encode)

ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_COMPARISON_C_GENERIC_DVR_RANK_20261104_v63_REVIEW'
PRIMARY='AMBIENT_B18_COMPARISON_C_GENERIC_DVR_RANK_20261104_v63'
def read(p):
    b=p.read_bytes()
    return json.loads(gzip.decompress(b) if p.suffix=='.gz' else b)
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()

def main(out):
    start=time.monotonic();signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
    out.mkdir(parents=True,exist_ok=False)
    folder=ROOT/'reports'/PRIMARY;primary=folder/'report.json';A=read(primary)
    assert sha(primary)=='e5272f143977cfce50ed5d26df4d6910fa3e6a583f2874c4730c0bbe50e4818a'
    assert A['all_passed']
    for p,h in A['input_sha256'].items(): assert sha(ROOT/p)==h,p
    for p,h in A['evidence_sha256'].items(): assert sha(folder/p)==h,p
    proof=read(folder/'fixed-source-curve-and-rank-proof.json')
    inputs={}
    for p in A['input_sha256']:
        if p.startswith('reports/'):
            obj=read(ROOT/p);inputs[p]=obj
            if 'REVIEW' in p: assert obj['all_passed'],p
    def get(s): return inputs['reports/'+s]
    ca=get('AMBIENT_B18_CA_CRITICAL_SCHEME_20261103_v63/certificate.json.gz')
    cb=get('AMBIENT_B18_CB_CRITICAL_SCHEME_20261103_v63/certificate.json.gz')
    assert ca['affine_slice_free_indices']==cb['slice_free_indices']==[1,2,3,9,10,14,15,17,18,19]
    assert ca['center_fixed20'][:3]==cb['fiber_center19'][:3]==['-1','0','0']
    t=variable(1,0);z={}
    Ac={k:substitute(parse(p,['u','v','w']),[z,t,z],1) for k,p in ca['inverse_substitution'].items()}
    assert Ac['s0']==scale(t,F(1,10)) and not Ac['s1'] and not Ac['s2']
    assert all(not substitute(parse(p,[f's{i}' for i in range(10)]),[Ac[f's{i}'] for i in range(10)],1)
               for p in ca['seven_monic_eliminations'])
    assert {k:encode(v) for k,v in Ac.items()}==proof['actual_CA_V_source_coefficients']
    Bc={f's{i}':z for i in range(10)}
    Bc.update(s0=t,s3=scale(t,F(1,2)),s4=scale(power(t,2,1),F(-1,4)),
              s7=scale(t,F(-1,2)),s8=scale(power(t,2,1),F(1,4)),s9=scale(power(t,3,1),F(-1,24)))
    for p in cb['seven_monic_eliminations']:
        assert not substitute(parse(p,[f's{i}' for i in range(10)]),[Bc[f's{i}'] for i in range(10)],1)
    assert not substitute(parse(cb['u'],[f's{i}' for i in range(10)]),[Bc[f's{i}'] for i in range(10)],1)
    assert substitute(parse(cb['v'],[f's{i}' for i in range(10)]),[Bc[f's{i}'] for i in range(10)],1)==t
    assert {k:encode(v) for k,v in Bc.items()}==proof['actual_CB_V_source_coefficients']
    disc=add(mul(Bc['s0'],Bc['s0']),scale(Bc['s1'],4));assert disc==power(t,2,1)
    assert cb['full_ratio_source_formula']=='P=sum f_i F_i -q sum f_i G_i +sum H_i HH_i'
    f9=get('AMBIENT_B18_T3_FULL9_FIRSTJET_20261104_v63/full9-firstjet-certificate.json')
    coeff=get('AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63/actual-source-coefficients.json')
    bindings={}
    for family in ['CA','CB']:
        c=coeff['families'][family];q=f9['both57_source2'][family]
        expected=[p[1] for p in c['actual_zero_odd_I0_generators'][:7]]+[p[1] for p in c['actual_positive_odd_J0_generators']]
        assert q['F9_nodes']==expected and len(expected)==9
        assert F(q['determinant']) and q['all57_equations_all55_coefficients_verified']
        bindings[family]={'actual_F9_node_ids':expected,'accepted_unit_J57_determinant':q['determinant'],
                          'I0_selected_count':7,'J0_count':2}
    ga=get('AMBIENT_B18_T3_JOINT_CURVE_RESTRICTIONS_20261104_v63/joint-degree3-certificate.json')
    gb=get('AMBIENT_B18_CB_G_DIVISION_ORIGINS_20261104_v63/degree3-certificate.json')
    r,s=variable(2,0),variable(2,1)
    pa={tuple(m):F(c) for m,c in zip(ga['monomials'],ga['g_CA_all10_coefficients']) if F(c)}
    pb={tuple(m):F(c) for m,c in zip(gb['ordered_monomials'],gb['g_CB_all10_coefficients']) if F(c)}
    assert pa==scale(mul(power(add(scale(r,3),scale(s,2)),2,2),add(r,scale(s,4))),F(1,24))
    assert pb==scale(mul(add(r,scale(s,-1)),power(add(scale(r,2),scale(s,3)),2,2)),-6)
    source=get('AMBIENT_B18_CONDUCTOR_DVR_COMPARISON_PREFLIGHT_20261104_v63/DVR20-raw-cokernel26-and-comparison-preflight.json')
    assert list(map(F,source['CA_ell_on_own_plane_linear']))==[F(-5),F(5)]
    assert list(map(F,source['CB_ell_on_own_plane_linear']))==[F(35,3),F(-85,3)]
    tangents={'CAcritical':(4,-1),'CAself':(2,-3),'CBcritical':(1,1),'CBself':(3,-2)}
    values={};linearvalues={}
    for label,(i,j) in tangents.items():
        ell=add(scale(r,-5),scale(s,5)) if label.startswith('CA') else add(scale(r,F(35,3)),scale(s,F(-85,3)))
        val=substitute(ell,[scale(t,i),scale(t,j)],1)
        assert val==parse({'CAcritical':'-25*t','CAself':'-25*t','CBcritical':'-50*t/3','CBself':'275*t/3'}[label],['t'])
        values[label]=encode(val);linearvalues[label]=str(next(iter(val.values())))
    assert values==proof['exact_ell_linear_values_on_four_tangent_directions']
    artin=get('AMBIENT_B18_RAW_COKERNEL_ARTIN20_EVALUATION_20261104_v63_REVIEW/report.json')
    assert artin['actual_C_mod_ell_dimension']==8
    lengths={'CAcritical':1,'CAself':1,'CBcritical':0,'CBself':2}
    assert lengths==proof['generic_DVR_branch_contributions'] and sum(lengths.values())==4
    # This evidence states the separately checked mathematical gates, not an
    # automated substitute for their proofs in the dedicated review note.
    evidence={'fixed_source_CA_s10':{k:encode(v) for k,v in Ac.items()},
      'fixed_source_CB_s10':{k:encode(v) for k,v in Bc.items()},
      'CB_binary_discriminant':encode(disc),'actual_common_base_F9_bindings':bindings,
      'all_four_actual_tangent_ell_values':linearvalues,
      'canonical_conductor_embedding':'A=R+I*tau inside B=R+R*tau implies conductor IB and A/IB=R/I',
      'full_to_zero_gate':'I=I0+wJ0; I0+J0=m; generic J0 unit forces w*C=0; both tau generators weight0',
      'G19_cut_is_common_equivariant_action_product_not_arbitrary_basechange':True,
      'Pa_cokernel_after_right_exact_tensor':'meeting(A1,A2) tensor CB = S reduced, length1',
      'Qb_cokernel_after_right_exact_tensor':'meeting(B1,B2)=C2 is closed in CA, triplemeeting length2',
      'critical_CA_reduced_conductor_and_literal_N_incidence_length1':True,
      'critical_CB_unique_proper_flag_nonsquare_excludes_N_then_RminusN_smooth_excludes_CA':True,
      'all_generic_branch_reduced_fields_equal_fraction_QQ_ell_by_smooth_QQ_curves_and_ell_order1':True,
      'actual_generic_DVR_rank':4,'mu_from_separate_actual_Artin_peer':8,
      'torsion_minimal_generator_count':4,'no_torsion_exponents_or_length_K_or_D0_conclusion':True}
    ev=out/'independent-fixed-source-and-proof-gates.json';ev.write_text(json.dumps(evidence,indent=2)+'\n')
    note=ROOT/'research/notes'/f'{STEM}.md'
    additional_notes=[
      'AMBIENT_CA_SELF_PAIR_MEETING_DOUBLE_20261103_v63_REVIEW.md',
      'AMBIENT_R_COMPLEMENT_20261103_v63_REVIEW.md',
      'AMBIENT_B18_CA_N_CRITICAL_SCHEMES_20261103_v63_REVIEW.md',
      'AMBIENT_B18_CB_CRITICAL_SCHEME_20261103_v63_REVIEW.md']
    files=[primary,Path(__file__),ROOT/'code/ambient_B18_proper_core_five_lines_review_20261104_v63.py',note]
    files += [ROOT/p for p in A['input_sha256']]+[folder/p for p in A['evidence_sha256']]
    files += [ROOT/'research/notes'/p for p in additional_notes]
    report={'id':STEM,'all_passed':True,'primary_report_sha256':sha(primary),
      'independent_sparse_Fraction_source_checks':True,'independent_conceptual_all_four_branch_and_transport_audit':True,
      'actual_comparison_C_generic_DVR_rank':4,'actual_mu_D_C':8,'minimal_torsion_generator_count':4,
      'Pa_triple_meeting_length':1,'Qb_triple_meeting_length':2,'CAcritical_length':1,'CBcritical_length':0,
      'same_R0_canonical_conductor_embedding_and_G19_transport_accepted':True,
      'torsion_exponents_length_comparison_K_D0_depth_or_socle_claimed':False,
      'no_author_worker_new_source_jet_GB_point_root_or_factor_search':True,
      'input_sha256':{str(p.relative_to(ROOT)):sha(p) for p in files},'evidence_sha256':{ev.name:sha(ev)},
      'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:report[k] for k in ['id','all_passed','seconds','max_rss_kib']}))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
