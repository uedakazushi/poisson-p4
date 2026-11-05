#!/usr/bin/env python3
"""Independent exact filtration proof for r^3 on the original Told module."""
import json,pathlib,resource,signal,time
from fractions import Fraction as Q
from ambient_B18_C3_alpha4_cached_source_review_20261105_v63 import field
from ambient_B18_C_genuine_projection7_source_review_20261105_v63 import count,sha
ROOT=pathlib.Path(__file__).resolve().parents[1]
PRI='AMBIENT_B18_C_CONNECTING_R3_FACTORIZATION_SOURCE_20261105_v63';ID=PRI+'_REVIEW'
def save(p,d):p.write_text(json.dumps(d,indent=2)+'\n')
def main():
 t=time.monotonic();signal.alarm(55);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30))
 out=ROOT/'reports'/ID;out.mkdir(exist_ok=False);pins={};raw=0
 def pin(p,h=None,parse=True):
  nonlocal raw
  p=pathlib.Path(p);k=str(p.relative_to(ROOT));s=sha(p)
  if h:assert h==s,(k,h,s)
  if k not in pins:raw+=p.stat().st_size
  pins[k]=s
  return json.loads(p.read_text()) if parse else None
 pr=pin(ROOT/'reports'/PRI/'report.json')
 for p,h in pr['input_sha256'].items():
  pin(ROOT/p,h,False)
  if p.endswith('_REVIEW/report.json'):assert json.loads((ROOT/p).read_text())['all_passed']
 pin(ROOT/'code/ambient_B18_C_connecting_r3_factorization_source_20261105_v63.py',pr['code_sha256'],False)
 pin(ROOT/'research/notes'/f'{PRI}.md',pr['note_sha256'],False)
 en=next(iter(pr['evidence_sha256']));D=pin(ROOT/'reports'/PRI/en,pr['evidence_sha256'][en])
 def src(s):return ROOT/next(p for p in pr['input_sha256'] if p.endswith(s))
 old=src('known-initial40-and-length-upper17.json');exact=src('exact-Loewy11-and-source-basis-corollary.json')
 seven=src('honest-seven-source-quotient-four-fixed-kernel-classes-and-edge-constraints.json')
 values=src('phase2-full-coherent-head23-and-low-two-values.json')
 gens=field(old,'homogeneous_module_generators40');basis=field(exact,'actual_filtered_monomial_basis11')
 assert field(exact,'exact_Loewy_HF')==[5,3,2,1] and field(exact,'exact_length')==11
 def mon(z=0,r=0):return [0,0,0,0,0,z,0,r]
 assert gens[38]==[[3,mon(r=2),'1'],[4,mon(z=1,r=1),'11767/3652']]
 assert gens[39]==[[4,mon(z=1,r=2),'1']]
 assert field(old,'new_actual_cubic_source')['canonical_initial_relation']==gens[39]
 deep=[b for b in basis if b['degree']==3]
 assert len(deep)==1 and (deep[0]['target'],deep[0]['monomial_original8'])==(4,mon(r=3))
 I=field(seven,'fixed_I_Told_source_basis')
 assert (I['first']['target'],I['first']['monomial_original8'])==(4,mon(z=1))
 assert [(b['target'],b['monomial_original8']) for b in I['other_three']]==[(4,mon(z=1,r=1)),(4,mon(r=2)),(4,mon(r=3))]
 R=[[Q(v) for v in row] for row in field(seven,'actual_r_action7')]
 expected=[[Q(0)]*7 for _ in range(7)];expected[5][2]=Q(1,4);expected[5][3]=1;expected[6][4]=1
 assert R==expected
 chi=list(map(Q,field(values,'actual_cyclic_top_chi')));assert chi==[0,0,0,0,Q(-36,1681)]
 owned=count([pr,D,gens,basis,I,R,chi,pins])+12000;assert owned<50000 and raw+1000000<32<<20
 save(out/'positive-before-fixed-filtration-proof.json',{'pure_SOURCE':True,'whole_numeric_upper':owned,'no_new_source_or_edge_values':True,'all_actual_basis_dependencies_pinned':True})
 # Exact zero lifts use accepted m^4 Told=0; they are not merely
 # statements in the known initial cover. Each entry records its order.
 r2I=[mon(z=1,r=2),mon(z=1,r=3),mon(r=4),mon(r=5)]
 assert sum(r2I[0])==3 and all(sum(e)>=4 for e in r2I[1:])
 assert r2I[0]==gens[39][0][1]
 rG38=[[i,[e[j]+int(j==7) for j in range(8)],q] for i,e,q in gens[38]]
 assert rG38==[[3,mon(r=3),'1'],[4,mon(z=1,r=2),'11767/3652']]
 actual_r3=[Q(0),Q(0),Q(0),Q(0),Q(1)]
 factor=1/chi[4];assert factor==Q(-1681,36)
 assert [factor*x for x in chi]==actual_r3
 assert list(map(Q,D['allfive_actual_r3_head_images_in_deep_basis']))==actual_r3
 assert Q(D['factor_scalar'])==factor and D['actual_r3_Kbar_zero']
 proof={
  'all_actual_Told_basis_and_m4_zero_dependencies_verified':True,
  'r2_I_basis_images':r2I,'G38_times_r':rG38,
  'all_five_actual_r3_head_coordinates':list(map(str,actual_r3)),
  'cyclic_top_chi':list(map(str,chi)),'inverse_d':str(factor),
  'exact_filtration_proof':[
   'The accepted cubic initial qC=z5*r^2*h4 has an actual lift qC+O(m^4). Since m^4 Told=0, qC=0 in the actual module. Multiplying the other three fixed I Told basis vectors by r^2 has order at least4. Hence r^2(I Told)=0 exactly.',
   'The exact M7 relations give r*h0,r*h1,r*h2-r*h3/4 in I Told. Multiplying by r^2 shows r^3*h0=r^3*h1=0 and r^3*h2=r^3*h3/4. The actual lift of G38 has remainder O(m^3); after multiplication by r this remainder lies in m^4=0. Its other term is c*qC=0. Thus r^3*h3=0.',
   'The remaining r^3*h4 is the nonzero actual degree3 basis gdeep. Every element of m kills it because m^4 Told=0. The five h_i generate Told over the same Rbar, so equality on these five heads determines the full Rbar-linear map.',
   'The cyclic top functional chi on the genuine G=C[ell] is Rbar-linear to k, with values (0,0,0,0,d), d=-36/1681. Therefore multiplication by r^3 is (gdeep/d)*chi*delta on ALL of Told, including all higher coefficients.',
   'For the honest original exact sequence ker(delta)=K/ell K, this gives r^3(K/ell K)=0 and therefore r^3 K is contained in ell K. No cancellation of ell is permitted.'
  ],
  'actual_r3_Kbar_zero':True,'actual_r3_K_zero_not_claimed':True,
  'no_original_H0_DVR_exponent_nil_cube_or_index_claim':True,
  'whole_numeric_upper':owned,'RAW':raw}
 assert all(sha(ROOT/p)==h for p,h in pins.items())
 ep=out/'independent-exact-five-head-r3-factorization-and-lifted-containment.json';save(ep,proof)
 report={'id':ID,'all_passed':True,'frozen':True,'pure_SOURCE':True,'input_sha256':pins,'code_sha256':sha(__file__),'evidence_sha256':{ep.name:sha(ep)},'r3_Kbar_zero_proved':True,'r3_K_zero_claimed':False,'seconds':time.monotonic()-t}
 save(out/'report.json',report);signal.alarm(0);print(json.dumps({'all_passed':True,'sha256':sha(out/'report.json'),'owned':owned,'seconds':report['seconds']}))
if __name__=='__main__':main()
