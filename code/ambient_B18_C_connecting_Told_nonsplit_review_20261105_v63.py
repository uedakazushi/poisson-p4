#!/usr/bin/env python3
"""Independent associated-graded obstruction to the ORIGINAL UCT splitting."""
import hashlib,json,pathlib,resource,signal,time
from fractions import Fraction as Q
from ambient_B18_C3_alpha4_cached_source_review_20261105_v63 import field
from ambient_B18_C_genuine_projection7_source_review_20261105_v63 import count,sha
from ambient_B18_C_connecting_kernel3_split_review_20261105_v63 import c,v,add,sc,mm
ROOT=pathlib.Path(__file__).resolve().parents[1];PRI='AMBIENT_B18_C_CONNECTING_TOLD_NONSPLIT_20261105_v63';ID=PRI+'_REVIEW'
def save(p,d):p.write_text(json.dumps(d,indent=2)+'\n')
def main():
 started=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30))
 out=ROOT/'reports'/ID;out.mkdir(exist_ok=False);pins={};RAW=0
 def pin(p,h=None,parse=True):
  nonlocal RAW
  p=pathlib.Path(p);s=sha(p)
  if h:assert s==h
  k=str(p.relative_to(ROOT))
  if k not in pins:RAW+=p.stat().st_size
  pins[k]=s
  return json.loads(p.read_text()) if parse else None
 pr=pin(ROOT/'reports'/PRI/'report.json')
 for p,h in pr['input_sha256'].items():
  pin(ROOT/p,h,False)
  if p.endswith('_REVIEW/report.json'):assert json.loads((ROOT/p).read_text())['all_passed']
 code=ROOT/'code/ambient_B18_C_connecting_Told_nonsplit_source_20261105_v63.py';pin(code,pr['code_sha256'],False)
 pin(ROOT/'research/notes'/f'{PRI}.md',pr['note_sha256'],False)
 en=next(iter(pr['evidence_sha256']));D=pin(ROOT/'reports'/PRI/en,pr['evidence_sha256'][en])
 def find(s):return ROOT/next(p for p in pr['input_sha256'] if p.endswith(s))
 old=find('known-initial40-and-length-upper17.json');exact=find('exact-Loewy11-and-source-basis-corollary.json')
 seven=find('honest-seven-source-quotient-four-fixed-kernel-classes-and-edge-constraints.json')
 vals=find('phase2-full-coherent-head23-and-low-two-values.json')
 gens=field(old,'homogeneous_module_generators40');gr1=field(old,'standard_degree1');gr2=field(old,'standard_degree2')
 basis=field(exact,'actual_filtered_monomial_basis11')
 assert field(exact,'exact_length')==11 and field(exact,'exact_Loewy_HF')==[5,3,2,1]
 assert gr1==['e3*r','e4*z5','e4*r'] and gr2==['e4*z5*r','e4*r²']
 z5=[0]*8;z5[5]=1;r=[0]*8;r[7]=1;rr=[0]*8;rr[7]=2;z5r=[a+b for a,b in zip(z5,r)]
 assert [(x['target'],x['monomial_original8']) for x in basis if x['degree']==2]==[(4,z5r),(4,rr)]
 relation=gens[38];assert relation==[[3,rr,'1'],[4,z5r,'11767/3652']]
 coeff=Q(relation[1][2]);I=field(seven,'fixed_I_Told_source_basis');assert I['first']['degree']==1 and I['first']['target']==4 and I['first']['monomial_original8']==z5
 assert len(I['other_three'])==3 and all(x['degree']>=2 for x in I['other_three'])
 R=[[Q(q) for q in row] for row in field(seven,'actual_r_action7')]
 expected=[[Q(0)]*7 for _ in range(7)];expected[5][2]=Q(1,4);expected[5][3]=Q(1);expected[6][4]=Q(1)
 assert R==expected
 chi=list(map(Q,field(vals,'actual_cyclic_top_chi')));assert chi==[0,0,0,0,Q(-36,1681)]
 owned=count([pr,D,gens,gr1,gr2,basis,relation,I,R,chi,pins])+10000;assert owned<50000 and RAW+1000000<32<<20
 save(out/'positive-before-formal-gr2-source-check.json',{'pure_SOURCE':True,'three_alpha_unknowns_retained':True,'whole_owned_upper':owned,'no_source_or_edge_values':True})
 # rhi modm² in the ACTUALgr1basis(u,v,w), with everyallowed unknown
 # correction inITold retained. Multiplicationr uses actualG38 ingr2.
 heads=[[{}, {},c(Q(1,4)),c(1),{}],
        [v('alpha0'),v('alpha1'),v('alpha2'),{},{}],
        [{},{},{},{},c(1)]]
 rmap=[[c(-coeff),c(1),{}],[{},{},c(1)]]
 result=mm(rmap,heads)
 assert result[0]==[v('alpha0'),v('alpha1'),add(v('alpha2'),c(-coeff/4)),c(-coeff),{}]
 assert result[1]==[{}, {}, {}, {},c(1)]
 assert [[str(-coeff),'1','0'],['0','0','1']]==D['actual_r_on_gr1_to_gr2']
 assert list(map(Q,D['q2_coordinate_of_r2_allfiveheads']))==[0,0,0,0,1]
 obstruction=1/chi[4];assert obstruction==Q(-1681,36)==Q(D['section_obstruction_coefficient'])
 assert all(sha(ROOT/p)==h for p,h in pins.items())
 proof={'actual_HF5321_and_displayed_gr1_gr2_basis_pinned':True,'initial_G38_gives_actual_gr2_r_u_relation':True,
  'all_five_heads_r2_gr2_q2_row':[0,0,0,0,1],'unknown_alpha0_1_2_never_evaluated':True,
  'exact_nonzero_preimage_obstruction':str(obstruction),
  'full_nonsplitting_proof':[
   'ITold is4-dimensional asQQspace withbasis(v=z5h4,z5rh4,r²h4,r³h4). Onlyvhasdegree1 anditsgr1classisnonzero. Thusitsimageingr1Told isEXACTLYthelinev; arbitraryR0multiplierscannotaddanothergr1direction becauseallIToldelementsalreadyareQQlinearcombinationsofthefourbasisclasses.',
   'The actualquotientrelationsrh0=rh1=0,rh2=rh3/4 inM7 implyrh0,rh1,rh2-rh3/4∈ITold. Hence theirleadingclassesarealpha0v,alpha1v,alpha2v. Multiplicationr sendsu=rh3 to-cq1 byactualinitialG38,vtoq1,andw=rh4toq2. ExactHF5321provesq1,q2independent inactualgr2,notmerelyaknownuppercover.',
   'Everyactualpreimagexofthegenuinecyclictopg hasconstanthead4coefficient1/d becausechiheads0..3=0,chihead4=d andchiannihilatesmG. Thehigherpartofx liesinmTold andr²ofitliesinm³Told. Thusq2coordinateofr²x ingr2 is1/d=-1681/36,nonzero,forEVERYchoiceofpreimage.',
   'AnR0-linearsection wouldsatisfyr²s(g)=s(r²g)=0 sincegenuineGhasr²=0. Thiscontradicts theactualgr2obstruction. ThereforetheORIGINALUCT0→K/ellK→Told11→C[ell]4→0 isNOTR0split.'
  ],
  'quotient_split_is_compatible':'Theobstructionq2belongsITold andvanishesinM7,whereourindependentSPLITproof remainsvalid.',
  'no_Kbar_extension_or_original_H0_niltorsion_claim':True,'owned_upper':owned}
 ep=out/'independent-actual-Told-gr2-preimage-nonsplitting-proof.json';save(ep,proof)
 report={'id':ID,'all_passed':True,'frozen':True,'pure_SOURCE':True,'input_sha256':pins,'code_sha256':sha(__file__),
  'evidence_sha256':{ep.name:sha(ep)},'actual_Told_UCT_R0_nonsplit_proved':True,'M7_split_compatible':True,'seconds':time.monotonic()-started}
 save(out/'report.json',report);signal.alarm(0);print(json.dumps({'all_passed':True,'sha256':sha(out/'report.json'),'owned':owned,'seconds':report['seconds']}))
if __name__=='__main__':main()
