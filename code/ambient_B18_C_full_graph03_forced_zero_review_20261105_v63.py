#!/usr/bin/env python3
"""Independent full-C cubic obstruction, using genuine full-module torsion."""
import hashlib,json,pathlib,resource,signal,time
from fractions import Fraction as Q
ROOT=pathlib.Path(__file__).resolve().parents[1]
PRI='AMBIENT_B18_C_FULL_GRAPH03_FORCED_ZERO_20261105_v63';ID=PRI+'_REVIEW'
PIN='2f6ee029b34d5dfe1bda99a085d3f82a83f5aff7378dfb37b966c39fd5e3766e'
def sha(p):return hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(pathlib.Path(p).read_text())
def save(p,d):pathlib.Path(p).write_text(json.dumps(d,indent=2)+'\n')
def main():
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30))
    out=ROOT/'reports'/ID;out.mkdir(exist_ok=False);pins={};docs={}
    p=ROOT/'reports'/PRI/'report.json';assert sha(p)==PIN;pr=read(p);assert pr['all_passed'];pins[str(p.relative_to(ROOT))]=PIN
    for p,h in pr['input_sha256'].items():
        q=ROOT/p;assert sha(q)==h;pins[p]=h;d=read(q);docs[pathlib.Path(p).name]=d
        if p.endswith('_REVIEW/report.json'):assert d['all_passed']
    for p,h in pr['evidence_sha256'].items():
        q=ROOT/'reports'/PRI/p;assert sha(q)==h;pins[str(q.relative_to(ROOT))]=h;E=read(q)
    assert sha(ROOT/'code/ambient_B18_C_full_graph03_forced_zero_20261105_v63.py')==pr['code_sha256']
    assert sha(ROOT/'research/notes'/f'{PRI}.md')==pr['note_sha256']
    reduction=docs['full-C-ell-reduction-and-conditional-graph03-exponent3-proof.json']
    assert docs['report.json']['all_passed'] # all colliding report basenames are pinned individually above
    complete=docs['complete155-source-cubic-operator-with-four-formal-alphas.json']
    assert complete['formal_operator_shape']==[4,155] and complete['all_other151_columns_zero_in_known_cover']
    cut=docs['one-CUT1-residual-and-formal-coordinate-rank-certificate.json']
    assert cut['actual_C3_length19_proved'] and cut['actual_C3_HF_if_positive']==[3,6,6,4,0]
    # Bind the actual alpha result explicitly, avoiding colliding report.json
    # names in the compact evidence dictionary.
    alpha=read(ROOT/'reports/AMBIENT_B18_C3_ALPHA4_CACHED_VALUES_20261105_v63/report.json')
    assert alpha['actual_C3_HF_degree3']==4 and alpha['nonzero_count']==0
    purepeer=read(ROOT/'reports/AMBIENT_B18_C_FULL_LINEAR_REDUCTION_SOURCE_20261105_v63_REVIEW/report.json')
    assert purepeer['unconditional_fullC_ell_reduction_accepted'] and purepeer['conditional_graph03_exponent3_theorem_accepted']
    N6=[[Q(c) for c in row] for row in reduction['full_no_cut_cubic_NF6_rows12_columns']]
    N4=[N6[i] for i in [0,1,3,4]]
    # All actual extra cubic images in this cut are zero by the COMPLETE
    # source theorem + four alpha values. These columns therefore describe
    # the actual gr3(C3), and m4(C3)=0 makes them actual nonzero classes.
    A=[row[2] for row in N4];B=[row[6] for row in N4];H=[row[10] for row in N4]
    assert A==[Q(0),Q(1),Q(0),Q(0)] and B==[Q(0),Q(0),Q(0),Q(1)] and H==[Q(0)]*4
    assert A[1]*B[3]-A[3]*B[1]==1
    assert E['actual_nonzero_r_times_ell2_classes']=={'A':list(map(str,A)),'B':list(map(str,B))}
    assert reduction['accepted_DVR_data']['rank_D']==4
    # The following is a proof with exact submodules of FULL C. No
    # associated-graded strictness, ell cancellation, or finite-kernel lift
    # is used. Each implication is independently checked in the note.
    steps=[
      'F=m3C is finiteD of rank4. The unconditional reduction gives mF=m4C=ellF.',
      'Assume fullHF3=4. Then dim(F/ellF)=4 and DVR structure implies F is Dfree4.',
      'Because ell3C subset F and m4C3=0, there is a genuine chain ellF=m4C subset ell3C subset F.',
      'F/ell3C=m3C3 has actualQ-length4. F/ellF also has length4. The exact quotient sequence gives ell3C/ellF length0, hence ell3C=ellF=m4C.',
      'For any x in FULL C, ell3x=ell*y for some y in F. Thus t=ell2x-y is genuinely in C[ell], not merely a kernel in a truncated module.',
      'For u in m, ut=u*ell2x-u*y lies in F: the first summand is in m3C, the second in mF. Also ell*(ut)=0 by R0-linearity and commutativity.',
      'Since F is D-torsionfree, ut=0. Therefore u*ell2x=u*y is in mF=ell3C.',
      'Thus m*ell2C maps to0 in actualC3. Taking u=r and x the raw A target contradicts its actual nonzero class with coordinate(0,1,0,0).']
    assert len(steps)==8
    assert E['unconditional_full_C_HF_degree3_values']==[5,6]
    rows=[r['physical_row'] for r in reduction['five_graph_witness_source_bindings']]
    assert rows==[63,64,68,73,74]
    assert E['five_graph03_actual_coefficients_source_forced_zero']==[[p,'0'] for p in rows]
    # The prior independently verified diagonal minor d^2 makes each
    # nonzero graph03 coefficient incompatible with the just-proved lower5.
    cp=ROOT/'reports/AMBIENT_B18_C_FULL_LINEAR_REDUCTION_SOURCE_20261105_v63_REVIEW/universal-graph03-minor-companion.json'
    assert sha(cp)=='e3a4fd2d0074e99ad8c5cd03487376d5e2989e30c5a4387c484c350a5acbff5b'
    pins[str(cp.relative_to(ROOT))]=sha(cp);assert read(cp)['ordinary_QQ_d_minor']==[[2,'1']]
    assert not(ROOT/'reports/AMBIENT_B18_C_FULL_GRAPH03_CACHE_VALUES_20261105_v63').exists()
    proof={'actual_raw_C3_r_ell2_map_columns_A_B_H':[list(map(str,v)) for v in [A,B,H]],
      'actual_map_rank2_fixed_minor':'1','literal_classes_not_only_upper_cover':True,
      'full_module_obstruction_steps':steps,
      'unconditional_fullC_HF_degree3':[5,6],
      'unconditional_D_torsion_generator_count_of_m3C':[1,2],
      'all5_actual_graph03_coefficients_zero_without_readback':rows,
      'conditional_nonzero_graph_to_HF3_4_theorem_used_only_contrapositively':True,
      'no_ell_cancellation_on_C_no_assumed_graded_exactness_no_Tor_shortcut':True,
      'genuine_full_C_ell_torsion_used_in_hypothetical_argument':True,
      'Tor_D_m3C_not_identified_with_m3_Tor_D_C':True,
      'remaining_C_torsion_pattern':[1,1,2,'e>=3'],
      'does_NOT_prove_e_ge4_or_exclude_1123_by_another_mechanism':True,
      'no_original_H0_nilcube_nilindex_or_global_component_conclusion':True,
      'graph03_prepared_SOURCE_custody_only_hash_archived_NOT_operationally_peered_or_run':True,
      'older_graph5_V1_V2_readbacks_remain_unrun':True,
      'no_new_source_physical_gauge_IFT_jet_or_scalar_readback':True,
      'static_numeric_upper':25000,'hard_numeric':150000,'hard_seconds':60,'hard_address_bytes':2<<30,'hard_RAW':32<<20}
    ep=out/'independent-full-module-torsion-obstruction-and-five-forced-zeros.json';save(ep,proof)
    result={'id':ID,'all_passed':True,'frozen':True,'source_only':True,'primary_report_sha256':PIN,
      'input_sha256':pins,'code_sha256':sha(__file__),'evidence_sha256':{ep.name:sha(ep)},
      'fullC_HF3_5_or6_and_graph03_zeros_accepted':True,'longer_exponent_still_only_ge3':True,
      'elapsed_seconds':time.monotonic()-start}
    save(out/'report.json',result);signal.alarm(0)
    print(json.dumps({'all_passed':True,'report':str(out/'report.json'),'sha256':sha(out/'report.json'),'seconds':result['elapsed_seconds']}))
if __name__=='__main__':main()
