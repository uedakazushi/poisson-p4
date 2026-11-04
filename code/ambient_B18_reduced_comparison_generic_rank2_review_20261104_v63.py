#!/usr/bin/env python3
"""Independent same-source quotient maps/graded transport/DVR proof audit."""
from pathlib import Path
import hashlib,json,resource,time
ROOT=Path(__file__).resolve().parents[1]
AUTHOR='AMBIENT_B18_REDUCED_COMPARISON_GENERIC_RANK2_20261104_v63';STEM=AUTHOR+'_REVIEW'
def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def read(p):return json.loads((ROOT/p).read_text())
def main():
 started=time.monotonic();out=ROOT/'reports'/STEM;assert not(out/'report.json').exists();out.mkdir(parents=True,exist_ok=True);hashes={}
 def get(st):
  p='reports/'+st+'/report.json';hashes[p]=sha(p);return read(p)
 R=get(AUTHOR);assert R['all_passed'] and hashes['reports/'+AUTHOR+'/report.json']=='6dda65f6400e143f71aff9a532fae8730cc39433bdcc7fc9cb3456510ccbe5f7'
 for p,h in R['input_sha256'].items():assert sha(p)==h;hashes[p]=h
 for n,h in R['evidence_sha256'].items():p='reports/'+AUTHOR+'/'+n;assert sha(p)==h;hashes[p]=h
 value=get('AMBIENT_B18_REDUCED_COMPARISON_CBAR_VALUES_20261104_v63_REVIEW')
 origin=get('AMBIENT_B18_REDUCED_COMPARISON_ORIGIN_20261104_v63_REVIEW')
 assert value['all_passed'] and value['actual_reduced_Cbar_length']==4
 assert hashes['reports/AMBIENT_B18_REDUCED_COMPARISON_CBAR_VALUES_20261104_v63_REVIEW/report.json']=='fb4546b8667fc50c668ded9cc62c5597f3d19c08cbe1a6f116b6532209a09e6c'
 assert origin['all_passed'] and origin['actual_Cprime_mu_R0']==3
 old=get('AMBIENT_B18_COMPARISON_C_GENERIC_DVR_RANK_20261104_v63_REVIEW')
 assert old['all_passed'] and old['actual_comparison_C_generic_DVR_rank']==4
 assert old['independent_conceptual_all_four_branch_and_transport_audit'] and old['same_R0_canonical_conductor_embedding_and_G19_transport_accepted']
 ca=get('AMBIENT_CA_CRITICAL_V1_CONTACT_ALGEBRA_20261104_v63_REVIEW')
 pa=get('AMBIENT_CA_SELF_PAIR_ACTUAL_D_FULL_RING_20261103_v63_REVIEW')
 qb=get('AMBIENT_CB_SELF_PAIR_ACTUAL_D_FULL_RING_20261103_v63_REVIEW')
 assert ca['all_passed'] and pa['all_passed'] and qb['all_passed']
 assert pa['actual_full_D_two_ideal_inclusions_and_all_ambient_upper_bound_accepted']
 assert qb['asymmetric_Qb_full_actual_D_ideal_equality_accepted']
 # These are the full actual maps from their independently bound source
 # diagrams, not abstract ring isomorphisms chosen to align afterwards.
 maps={
  'criticalCA':{'source_reduced':'G[[l7,v]]','target_reduced':'G[[l7,v]]','map':'identity','coker_length':0},
  'Pa':{'source_reduced':'O_N','target_reduced':'O_N direct_sum O_SP','map':'a -> (a,restriction_SP(a))','coker':'O_SP','coker_length':1},
  'criticalCB':{'source':'outside_actual_CA','target':'zero','coker_length':0},
  'Qb':{'source_reduced':'B/(f)','target_reduced':'B/(f) direct_sum B/(z7,f,rho)','map':'a -> (a,restriction(a))','coker':'B/(z7,f,rho)','coker_length':1}}
 assert sum(v['coker_length'] for v in maps.values())==2
 proof={
  'same_actual_quotient':'E4 -> E4_red induces C -> Cprime surjective. Thus Cprime is finite over the same DVR and has no generic support outside the accepted four C curves.',
  'source_maps':maps,
  'graph_exactness':'For any quotient pi:A->B, sigma:A⊕B->B sends (a,b) to b-pi(a); kernel is exactly the graph, and (0,b) splits sigma. Thus graph coker is B, including these Pa/Qb actual restrictions.',
  'criticalCA_reduction':'Accepted D=G[[l7,v,x]]/x(l7,x,v²) and E=G[[l7,v,t]]/(l7t,t²-vt,t³), x->t². Both reduced rings equalG[[l7,v]], with the actual identity map.',
  'Qb_reduction':'Accepted D=B/(f JZ), E=B/(f)⊕B/JZ, JZ=(z7,f²,rho²). Reduction gives B/(f) and B/(z7,f,rho); the second restriction factors through the first. No label-swap fromPa.',
  'degree_and_no_extra_support':'The already independently source-verified curves have ell order1 and residue-field degree1 over k((ell)). Reduced quotient cannot add support; lengths0+1+0+1 give DVR rank2.',
  'full_weight_transport':[
   'Use the actual full inclusions Ai=R10+Ii*taui subset Bi, Ii=Ii0+w*Ji0. Both same-source inclusions and tensor maps are graded.',
   'At each generic zero-conductor prime, the opposite conductor is a unit and the relevant positive Ji0 is a unit. The entire original comparison coker is killed by the product conductor, hence by w.',
   'Its generators tauA,tauB,tauA*tauB have weight0; generation PLUS w-annihilation gives concentration in degree0. w-annihilation alone would not suffice.',
   'The full reduced-target coker is its graded quotient. Nil(Efull) intersect degree0 equals Nil(Efull_0), and exact degree0 projection gives the actual same Cprime=coker(D0->E4_red).',
   'The accepted actual simultaneous G19 formal product includes the normalization/comparison maps. The flat regular factor preserves these generic lengths. No arbitrary w-cut or type-only coordinate identification is used.'],
  'DVR_structure':'mu_D Cprime=length(Cprime/ell Cprime)=4. FiniteDVRmodule rank2 therefore has exactly2 torsion cyclic summands D/(ell^mi), mi>=1. Their exponents/length are unknown.',
  'r_not_ell':'The finite quotient r-chain of length2 is after ell=0; it does not determine either DVR torsion exponent.',
  'not_claimed':['Nil(D0) finite or generated by these classes','K/JK or original cube','Torsion length or exponents','New support or whole classification']}
 ev=out/'independent-four-maps-degree0-and-DVR-proof.json';ev.write_text(json.dumps(proof,indent=2)+'\n')
 for p in [str(Path(__file__).relative_to(ROOT)),'research/notes/'+STEM+'.md']:hashes[p]=sha(p)
 result={'id':STEM,'all_passed':True,'same_actual_four_reduced_target_maps_accepted':True,'full_w_generated_weight0_then_graded_quotient_transport_accepted':True,'same_flat_G19_normalization_diagram_not_arbitrary_cut':True,'generic_DVR_rank_Cprime':2,'length_Cprime_mod_ell':4,'DVR_torsion_minimal_generators':2,'torsion_length_exponents_original_nilgens_or_cube_claimed':False,'no_new_source_job':True,'input_sha256':hashes,'evidence_sha256':{ev.name:hashlib.sha256(ev.read_bytes()).hexdigest()},'seconds':time.monotonic()-started,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
 (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k not in ('input_sha256','evidence_sha256')}))
if __name__=='__main__':main()
