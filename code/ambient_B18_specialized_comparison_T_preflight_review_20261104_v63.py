#!/usr/bin/env python3
"""Independent theorem/source review of the specialized comparison T preflight.
No actual comparison coefficients, source solves or large T matrix are computed.
"""
import argparse,hashlib,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
PRIMARY='AMBIENT_B18_SPECIALIZED_COMPARISON_T_PREFLIGHT_20261104_v63';STEM=PRIMARY+'_REVIEW'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def padd(a,b):
    r=dict(a)
    for w,c in b.items():r[w]=r.get(w,F(0))+c
    return {w:c for w,c in r.items() if c}
def neg(p):return {w:-c for w,c in p.items()}
def pmul(a,b):
    r={}
    for u,c in a.items():
        for v,d in b.items():r[u+v]=r.get(u+v,F(0))+c*d
    return {w:c for w,c in r.items() if c}
def reduce_words(p):
    # Typed identities only: d2*h=I24-m*n and d2*d3=0.
    r={};todo=list(p.items())
    while todo:
        w,c=todo.pop();changed=False
        for i in range(len(w)-1):
            if w[i:i+2]==('d2','d3'):changed=True;break
            if w[i:i+2]==('d2','h'):
                todo.append((w[:i]+w[i+2:],c));todo.append((w[:i]+('m','n')+w[i+2:],-c));changed=True;break
        if not changed:r[w]=r.get(w,F(0))+c
    return {w:c for w,c in r.items() if c}

def main(out):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False)
    folder=ROOT/'reports'/PRIMARY;rp=folder/'report.json';A=read(rp)
    assert sha(rp)=='4853d5e66d82b7b221495a4e55ca94667accc6a5d8e33e20ffec0c78f2f9e1b7' and A['all_passed']
    ep=ROOT/A['evidence']['path'];assert sha(ep)==A['evidence']['sha256'];E=read(ep)
    D={}
    for k,r in A['source_hashes'].items():
        p=ROOT/r['path'];assert sha(p)==r['sha256'];D[k]=read(p);assert D[k]['all_passed']
    assert E['source_hashes']==A['source_hashes']
    respeer=ROOT/'reports/AMBIENT_B18_CBAR_LITERAL_RESOLUTION_20261104_v63_REVIEW/report.json'
    assert read(respeer)['all_passed']
    assert D['literal_Cbar_resolution']['ranks']==[3,24,84,168,210,168,84,24,3]
    assert D['generic_C_rank_peer']['actual_comparison_C_generic_DVR_rank']==4
    assert D['generic_C_rank_peer']['actual_mu_D_C']==8
    low=E['source_low_complex'];adapter=E['complete_kernel_adapter'];finite=E['finite_support']
    assert low['V_shape']==[3,40] and low['X_shape']==[40,752]
    assert low['source_P']==[9,640] and low['source_B']==[48,640]
    assert low['Q_reduced']==[4,32] and low['SyzQ_reduced']==[32,112]
    assert 48-4*4==32 and 32+9-1==40 and 4-1==3 and 640+112==752
    assert not low['actual_X752_and_V40_scalar_matrices_serialized_here']
    assert adapter['m_shape']==[24,40] and adapter['n_shape']==[40,24] and adapter['h_shape']==[84,24]
    assert adapter['required_actual_identities']==['d1*m=Vbar','Vbar*n=d1','d2*h=I24-m*n']
    assert adapter['S_shape']==[40,124] and adapter['relation_shape']==[124,192]
    assert adapter['T_cokernel_shape']==[124,944]
    assert 84+40==124 and 168+24==192 and 192+752==944
    # Universal relation-block products in a free noncommutative alphabet:
    # [n*d2, I-n*m] * [[d3,-h],[0,n]] = 0, with exact signs.
    left={('n','d2','d3'):F(1)}
    right={('n','d2','h'):F(-1),('n',):F(1),('n','m','n'):F(-1)}
    assert not reduce_words(left) and not reduce_words(right)
    assert not adapter['actual_m_n_h_and752_cycle_lifts_constructed']
    assert E['DVR_count']['mu_D_C']==8 and E['DVR_count']['generic_DVR_rank_C']==4
    assert E['DVR_count']['dimension_QQ_C_ell']==8-4==4
    assert finite['I_squared_annihilates_T'] and finite['m10_T_zero'] and not finite['I_annihilates_T_claimed']
    assert not finite['H20_naive_kernel_replacement_permitted']
    assert E['conservative_raw_348_counts']['QQ_generators']==348*124==43152
    assert E['conservative_raw_348_counts']['QQ_relation_columns']==348*944==328512
    # Conditional cancellation: minimality d2(0)=0 implies m0*n0=I24.
    # The bottom n block in24 relation columns therefore has a unit24minor;
    # ordinary local presentation cancellation removes24 rows AND24 columns.
    conditional_reduced=[124-24,944-24];assert conditional_reduced==[100,920]
    ev={'actual_source_low_shapes':{'P':[9,640],'U':[4,9],'Q':[4,48],'Qred':[4,32],
         'syzQred':[32,112],'V':[3,40],'X':[40,752]},
      'source_unit_reductions_only_conditional_until_actual_coefficient_adapter_saved':True,
      'ell_torsionfree_J_reason':'J is a submodule of free R^3',
      'specialization_proof':['Tor1(J,Rbar)=0 embeds L/ellL into F40bar',
       'Tor1(C,Rbar)=C[ell] is ker(Jbar -> F3bar)',
       'preimage yields0 -> L/ellL -> kerVbar -> C[ell] ->0',
       'quotient by actual imageXbar and right exactness yield0 -> K/ellK -> T -> C[ell] ->0'],
      'C_ell_dimension':4,'formula_mu_D_K':'length_QQ(T)-4',
      'complete_S_generators':124,'complete_S_relations':192,'complete_T_cokernel_shape_conditional':[124,944],
      'free_noncommutative_relation_block_zero_checked':True,
      'completeness_proof':'For S(w,x)=0 put y=mx-d2w. Then x=ny and d2(w+hy)=0, so w=d3a-hy; all relations are the stated192columns.',
      'cycle_lift_proof':'For Vbar*v=0 exactness d1*m*v=0 gives d2*w=m*v and S(w,v)=v.',
      'conditional_24_unit_cancellations':{'premise':'actual m,n,h identities constructed; d2 minimal',
          'residue_identity':'m0*n0=I24','reduced_shape':[100,920],
          'actual_source_transports_or_cancellation_coefficients_constructed':False},
      'annihilator':'I kills both Kbar and C[ell], hence I^2*T=0; m^5 subset I gives m^10*T=0',
      'correct_finite_base':'Rbar/I^2; exact coefficient restrictions and completed T cokernel required first',
      'T_length_K_mu_or_generic_rank_torsion_D0_depth_not_computed':True}
    out_ev=out/'independent-specialization-and-comparison-proof.json';out_ev.write_text(json.dumps(ev,indent=2)+'\n')
    note=ROOT/'research/notes'/f'{STEM}.md'
    files=[rp,ep,Path(__file__),note,respeer]+[ROOT/r['path'] for r in A['source_hashes'].values()]
    files += [ROOT/'research/notes/AMBIENT_B18_COMPLETE_NORMALIZER_SYZYGY_OPERATOR_20261104_v63_REVIEW.md']
    report={'id':STEM,'all_passed':True,'proof_and_source_preflight_only':True,
      'exact_specialization_sequence_accepted':'0 -> K/ellK -> T -> C[ell] ->0',
      'C_ell_QQ_dimension':4,'mu_D_K_formula':'length_QQ(T)-4',
      'conditional_T_cokernel_shape':[124,944],'conditional_24_unit_reduced_shape':[100,920],
      'I_squared_T_zero_m10_T_zero_accepted':True,'I_T_zero_claimed':False,
      'actual_comparison_coefficients_T_length_K_rank_torsion_D0_depth_computed':False,
      'no_authorworker_GB_source_solve_newpoint_jet_rank_search':True,
      'input_sha256':{str(p.relative_to(ROOT)):sha(p) for p in files},'evidence_sha256':{out_ev.name:sha(out_ev)},
      'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:report[k] for k in ['id','all_passed','seconds','max_rss_kib']}))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
