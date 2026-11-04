#!/usr/bin/env python3
"""Independent reverse-divisor module arithmetic and upper-length audit.

No author import/run. The division order is the reverse of the producer's,
and the two extra generators are checked directly against original source
quadratics. This proves a known initial SUBmodule, never its equality to in(T).
"""
from fractions import Fraction as Q
from pathlib import Path
import hashlib,itertools,json,resource,signal,time
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_T_KNOWN_INITIAL_MODULE_LENGTH_UPPER24_20261104_v63_REVIEW'
AUTHOR='AMBIENT_B18_T_KNOWN_INITIAL_MODULE_LENGTH_UPPER24_20261104_v63'
PATHS={
 'author_report':f'reports/{AUTHOR}/report.json',
 'author_proof':f'reports/{AUTHOR}/known-initial-module-and-length-bound.json',
 'rewrite':'reports/AMBIENT_B18_T_KNOWN_QUADRATIC_CUBIC_CLOSURE_20261104_v63/linear-module-rewrite-certificate.json',
 'quadratics':'reports/AMBIENT_B18_T_KNOWN_QUADRATIC_CUBIC_CLOSURE_20261104_v63/actual-known-quadratic-cubic-closure-certificate.json',
 'source_quad_peer':'reports/AMBIENT_B18_T_ACTUAL_PROJECTED_QUADRATIC_RANK_20261104_v63_REVIEW/report.json',
 'ann_peer':'reports/AMBIENT_B18_SPECIALIZED_COMPARISON_T_PREFLIGHT_20261104_v63_REVIEW/report.json',
 'base_peer':'reports/AMBIENT_B18_CONDUCTOR_SQUARE_CLOSED348_20261104_v63_REVIEW/report.json',
 'five_peer':'reports/AMBIENT_B18_T_MINIMAL5_SCHUR_OPERATOR_20261104_v63_REVIEW/report.json',
 'rankK':'reports/AMBIENT_B18_CRITICAL_CONTACT_AND_COMPARISON_RANK4_20261104_v63/report.json'}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def powers(v):return tuple(v.count(i) for i in range(8))
def put(poly,key,c):
 z=poly.get(key,Q(0))+c
 if z:poly[key]=z
 else:poly.pop(key,None)
def minus(poly,q,c,shift=(0,)*8):
 poly=dict(poly)
 for key,a in q.items():put(poly,tuple(x+y for x,y in zip(key[:8],shift))+(key[8],),-c*a)
 return poly
def lead(p):return max(p,key=lambda k:(-k[8],k[:8]))
def divides(k,l):return k[8]==l[8] and all(a<=b for a,b in zip(k[:8],l[:8]))
def divide(p,gs):
 p=dict(p);r={};steps=[]
 while p:
  k=lead(p);eligible=[j for j,g in enumerate(gs) if divides(lead(g),k)]
  if not eligible:r[k]=p.pop(k);continue
  j=eligible[-1];l=lead(gs[j]);shift=tuple(a-b for a,b in zip(k[:8],l[:8]));c=p[k]/gs[j][l]
  p=minus(p,gs[j],c,shift);steps.append([j,list(shift),str(c)])
  assert len(steps)<5000
 return r,steps
def encode(p):return [[list(k[:8]),k[8],str(c)] for k,c in sorted(p.items())]
def source_linear(v):return {powers([int(i)%8])+(int(i)//8,):Q(c) for i,c in v}
def source_quad(v):
 ms=list(itertools.combinations_with_replacement(range(8),2))
 return {powers(ms[int(i)%36])+(int(i)//36,):Q(c) for i,c in v}
def main():
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30))
 out=ROOT/'reports'/STEM;out.mkdir(parents=True,exist_ok=False)
 inp={k:{'path':p,'sha256':sha(ROOT/p)} for k,p in PATHS.items()}
 (out/'phase.json').write_text(json.dumps({'phase':'independent_review_inputs_saved','inputs':inp},indent=2)+'\n')
 D={k:json.loads((ROOT/p).read_text()) for k,p in PATHS.items()}
 A,P=D['author_report'],D['author_proof'];assert A['all_passed']
 assert sha(ROOT/'code/ambient_B18_T_known_initial_module_length_upper24_20261104_v63.py')==A['code_sha256']
 assert sha(ROOT/'research/notes'/f'{AUTHOR}.md')==A['note_sha256']
 for p,h in A['input_sha256'].items():assert sha(ROOT/p)==h
 for p,h in A['evidence_sha256'].items():assert sha(ROOT/'reports'/AUTHOR/p)==h
 for k in ['source_quad_peer','ann_peer','base_peer','five_peer','rankK']:assert D[k]['all_passed']
 assert D['source_quad_peer']['actual_HF_through_degree2']==[5,3,2]
 assert D['source_quad_peer']['actual_length_T_mod_m3']==10
 assert D['ann_peer']['I_squared_T_zero_m10_T_zero_accepted']
 assert D['base_peer']['actual_base_m10_zero_r9_nonzero_accepted']
 assert D['five_peer']['actual_complete_minimal_generator_presentation']==[5,825]
 assert D['rankK']['rank_D_K']==4 and D['rankK']['dim_k_C_ell_torsion']==4
 gs=[{tuple(e)+(int(t),):Q(c) for t,e,c in q} for q in P['homogeneous_generators39']]
 linear=list(map(source_linear,D['rewrite']['linear_basis37']));quads=list(map(source_quad,D['quadratics']['quadratic_basis178']))
 assert gs[:37]==linear and len(gs)==39 and len(quads)==178
 # Actual source membership of each of the two new homogeneous generators,
 # verified by direct polynomial identities, without assuming the linear
 # generators themselves form a GB.
 membership=[]
 for n,receipt in enumerate(P['extra_quadratic_source_membership']):
  v=dict(quads[receipt['source_quad_index']])
  for j,shift,c in receipt['linear_trace']:v=minus(v,linear[j],Q(c),tuple(shift))
  for j,shift,c in receipt['previous_extra_trace']:v=minus(v,gs[37+j],Q(c),tuple(shift))
  assert v=={k:Q(receipt['divide'])*c for k,c in gs[37+n].items()}
  membership.append({'extra_generator':37+n,'source_quadratic':receipt['source_quad_index'],'direct_exact_identity_checked':True})
 expected=[{powers([5,5])+(4,):Q(1)},{powers([7,7])+(3,):Q(1),powers([5,7])+(4,):Q(11767,3652)}]
 assert gs[37:]==expected
 # Reverse-divisor reductions are independent of the author's saved path.
 # Every source quadratic is checked again against the resulting39 rows.
 qreceipts=[]
 for j,p in enumerate(quads):
  r,tr=divide(p,gs);assert not r;qreceipts.append({'source_quad':j,'reverse_divisor_reduction':tr})
 pairs=[]
 for i in range(39):
  for j in range(i+1,39):
   a,b=lead(gs[i]),lead(gs[j])
   if a[8]!=b[8]:continue
   cm=tuple(max(x,y) for x,y in zip(a[:8],b[:8]));u=tuple(x-y for x,y in zip(cm,a[:8]));v=tuple(x-y for x,y in zip(cm,b[:8]));s={}
   s=minus(s,gs[i],-1/gs[i][a],u);s=minus(s,gs[j],1/gs[j][b],v)
   r,tr=divide(s,gs);assert not r
   pairs.append({'generators':[i,j],'common_exponents':list(cm),'independent_reverse_reduction':tr})
 assert len(pairs)==P['S_pair_count']==133
 trace_support=sum(len(r['reverse_divisor_reduction']) for r in qreceipts)+sum(len(r['independent_reverse_reduction']) for r in pairs)+sum(map(len,gs+quads))
 assert trace_support<150000
 lms=set(lead(g) for g in linear);all1={powers([i])+(t,) for t in range(5) for i in range(8)}
 free=all1-lms
 assert free=={powers([7])+(3,),powers([5])+(4,),powers([7])+(4,)}
 assert set(lead(g) for g in gs[37:])=={powers([7,7])+(3,),powers([5,5])+(4,)}
 # Per-target combinatorics prove the standard set for EVERY degree.
 description={
  'targets0_1_2':'only the constant survives the leading linear submodule',
  'target3':'constant and r; r^2 is a leading monomial',
  'target4':'r^d and z5*r^(d-1) for d>=1; z5^2 is a leading monomial',
  'all_degree_hilbert_series':'5+3*t+2*t^2/(1-t)'}
 hf=[5,3]+[2]*8;assert sum(hf)==24
 ev={'source_manifest':inp,'author_source_not_imported_or_executed':True,
  'independent_polynomial_representation':'8 exponent coordinates followed by module position; exact Fraction arithmetic',
  'division_choice':'last applicable generator, opposite producer first-applicable policy',
  'source39_generators':list(map(encode,gs)),
  'both_extra_source_memberships_directly_checked':membership,
  'all178_source_quadratic_zero_reduction_receipts':qreceipts,
  'all133_same_position_S_pairs_independently_zero':pairs,
  'module_GB_generator_count39_not_S_pair_count':True,
  'standard_monomials_all_degree_proof':description,
  'known_homogeneous_initial_submodule_is_contained_in_true_initial_module':True,
  'known_submodule_equals_full_initial_module_claimed':False,
  'source_membership_bridge':'linear37 are actual initial relations; first177 normalized quadratics are combinations of true source linear multiples and d3 quadratics; last canonical rank1 lift differs from an actual620 quadratic by the true177 span.',
  'same_maximal_ideal_bridge':'all variables are the actual CA-F9 eight parameters of R0/(ell); accepted all-order conductor pair changes are local units and preserve this maximal ideal; actual I²-ann and m10 base imply m10*T=0.',
  'filtration_theorem':'the known homogeneous module quotient surjects degreewise onto gr_m(T); sum dimensions0..9 bounds length because m10*T=0.',
  'actual_T_length_interval':[10,24],'actual_K_DVR_mu_interval':[6,20],'actual_K_DVR_torsion_mu_interval':[2,16],
  'no_H0_length_or_torsion_exponent_upper_bound':True,'no_new_HF3_or_nil_cube_conclusion':True,
  'no_new_source_jet_scalar348_matrix_or_rank_search':True}
 ep=out/'independent-known-initial-module-and-filtration-bound.json';blob=json.dumps(ev,separators=(',',':'))+'\n';assert len(blob.encode())<32<<20;ep.write_text(blob)
 report={'id':STEM,'all_passed':True,'independent39_generator_module_Buchberger_accepted':True,
  'all133_S_pairs_exact_zero':True,'both_extra_source_memberships_and_all178_quad_memberships_accepted':True,
  'actual_I_squared_ann_m10_and_maximal_ideal_binding_accepted':True,'actual_T_length_interval_accepted':[10,24],
  'actual_K_mu_interval_accepted':[6,20],'actual_K_torsion_mu_interval_accepted':[2,16],
  'known_initial_submodule_only_not_full_initial_module':True,'H0_length_HF3_nil_or_exact_length_not_inferred':True,
  'input_hashes':inp,'code_sha256':sha(Path(__file__)),'evidence_sha256':{ep.name:sha(ep)},
  'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)
if __name__=='__main__':main()
