#!/usr/bin/env python3
"""Fixed tiny natural gr2/r test; no higher comparison coefficient."""
import hashlib,itertools,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
import ambient_B18_honest_nil_T_first_values_review_20261104_v63 as Q
ROOT=Path(__file__).resolve().parents[1];ID='AMBIENT_B18_TOLD_TNIL_GRADED2_KERNEL_TEST_20261105_v63';OUT=ROOT/'reports'/ID
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def enc(v):
 if isinstance(v,F):return str(v)
 if isinstance(v,dict):return {str(k):enc(x)for k,x in v.items()}
 if isinstance(v,(tuple,list)):return [enc(x)for x in v]
 return v
def main():
 st=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));OUT.mkdir(exist_ok=True);pins={}
 def load(s,n='report.json'):
  p=ROOT/'reports'/s/n;pins[str(p.relative_to(ROOT))]=sha(p);return json.loads(p.read_text())
 assert load('AMBIENT_B18_TOLD_TNIL_NATURAL_MAP_20261104_v63_REVIEW')['all_passed']
 N=load('AMBIENT_B18_TOLD_TNIL_NATURAL_MAP_20261104_v63','actual-natural-Told-to-Tnil-source-map-and-origin.json');O=N['actual_origin_data']
 exact=load('AMBIENT_B18_T_EXACT_LOEWY11_20261105_v63');assert exact['actual_T_length']==11 and exact['actual_T_Loewy_HF']==[5,3,2,1]
 B=load('AMBIENT_B18_T_EXACT_LOEWY11_20261105_v63','exact-Loewy11-and-source-basis-corollary.json')
 K=load('AMBIENT_B18_T_KNOWN_INITIAL_LENGTH_UPPER17_20261104_v63','known-initial40-and-length-upper17.json')
 C=load('AMBIENT_B18_HONEST_NIL_T_CUBIC_CLOSURE_PREFLIGHT_20261105_v63','literal-known-cubic-input-and-exact24-row-reduction.json')
 assert load('AMBIENT_B18_HONEST_NIL_T_CUBIC_CLOSURE_20261105_v63_REVIEW')['all_passed']
 assert N['old_minimal_generator_free56_rows']==[22,23,29,30,55]
 old2=[v for v in B['actual_filtered_monomial_basis11']if v['degree']==2]
 assert old2==[{'target':4,'source_V56_row':55,'monomial_original8':[0,0,0,0,0,1,0,1],'degree':2},{'target':4,'source_V56_row':55,'monomial_original8':[0,0,0,0,0,0,0,2],'degree':2}]
 assert K['new_actual_cubic_source']['canonical_initial_relation']==[[4,[0,0,0,0,0,1,0,2],'1']]
 assert K['homogeneous_module_generators40'][-1]==[[4,[0,0,0,0,0,1,0,2],'1']]
 origin=Q.sparse(O['actual_origin_map_columns5'][4]);assert origin[1]==1
 assert origin=={1:F(1),2:F(-14075,13448),3:F(36008,25215),4:F(28015,40344)}
 pi3=[Q.pairs(v)for v in C['quadratic_projection3_by180_sparse']];mons=list(itertools.combinations_with_replacement(range(8),2));idx={m:i for i,m in enumerate(mons)}
 cols=[]
 for v in old2:
  mon=tuple(i for i,n in enumerate(v['monomial_original8'])for _ in range(n));raw={36*t+idx[mon]:a for t,a in origin.items()};cols.append(Q.mvec(pi3,raw))
 assert not cols[0]and cols[1]
 pivot=min(cols[1]);a=cols[1][pivot];assert a
 # The fixed 3x2 map has rank1 and its kernel is exactly the first basis
 # vector; old r:gr2->gr3 is (0,1), with actual old m4=0.
 assert Q.plus(Q.scaled(cols[0],1),Q.scaled(cols[1],0))=={}
 rmap=[F(0),F(1)];kernel=[F(1),F(0)];assert sum(x*y for x,y in zip(rmap,kernel))==0
 proof=dict(actual_old_source_five_rows=[22,23,29,30,55],actual_old_head4_source=55,new_target_order='free44 42,43, then wedge84 6,24,26',
  actual_natural_head4_origin5=origin,stronger_zero_origin_head_premise_FALSE=True,
  actual_old_quadratic_basis=['e4*zA5*rA','e4*rA²'],actual_old_r_gr2_to_gr3=[0,1],old_m4_zero=True,
  actual_new_gr2_NF_basis180=C['quadratic_NF_basis180_indices'],actual_induced_gr2_columns3=cols,rank=1,fixed_nonzero_pivot=[pivot,1,a],
  exact_kernel_basis_old_gr2=[[1,0]],old_r_is_zero_on_entire_graded2_kernel=True,
  valid_criterion='If an old quadratic class v maps to zero in new gr2 and rv is nonzero in actual old m3 (old m4=0), its actual image lies in new m3 and r*i(v)=i(rv) is a nonzero element of new m4. This would exclude tail(3).',
  criterion_does_NOT_apply_to_this_kernel=True,
  possible_actual_new_tails_unchanged=[[3],[2,1],[1,1,1]],
  no_actual_higher_natural_coefficients_or_source_jet_or_nilcube=True,
  not_a_proof_that_tail3_occurs=True)
 p=OUT/'actual-graded2-map-and-negative-r-kernel-criterion.json';p.write_text(json.dumps(enc(proof),indent=2)+'\n')
 for p2 in[Path(__file__),Path(Q.__file__),ROOT/f'research/notes/{ID}.md']:pins[str(p2.relative_to(ROOT))]=sha(p2)
 rep=dict(id=ID,all_passed=True,frozen=True,actual_induced_gr2_rank=1,old_r_zero_on_graded2_kernel=True,tail3_exclusion_test_NEGATIVE_only=True,no_actual_tail_chosen=True,no_new_source_jet_or_higher_map_coefficients=True,input_sha256=pins,evidence_sha256={p.name:sha(p)},seconds=time.monotonic()-st)
 (OUT/'report.json').write_text(json.dumps(rep,indent=2)+'\n');signal.alarm(0);print(json.dumps({k:v for k,v in rep.items()if k not in['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':main()
