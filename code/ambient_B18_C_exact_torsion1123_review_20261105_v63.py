#!/usr/bin/env python3
"""Independent pure exact-C DVR, filtration and genuine C[ell] corollary."""
import hashlib,json,pathlib,resource,signal,time
ROOT=pathlib.Path(__file__).resolve().parents[1]
PRI='AMBIENT_B18_C_EXACT_TORSION1123_20261105_v63';ID=PRI+'_REVIEW'
PRIMARY='f8719eb1c4d819311800ff25f40dd2e13341926d013deadf5e343f3755138469'
def sha(p):return hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()
def save(p,d):pathlib.Path(p).write_text(json.dumps(d,indent=2)+'\n')
def main():
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30))
 out=ROOT/'reports'/ID;out.mkdir(exist_ok=False);pins={}
 def pin(p,h=None):
  p=pathlib.Path(p);s=sha(p)
  if h:assert s==h,(str(p),s,h)
  pins[str(p.relative_to(ROOT))]=s
  return json.loads(p.read_text()) if p.suffix=='.json' else None
 pr=pin(ROOT/'reports'/PRI/'report.json',PRIMARY);assert pr['all_passed'];docs={}
 for p,h in pr['input_sha256'].items():
  d=pin(ROOT/p,h);docs[pathlib.Path(p).name]=d
  if p.endswith('_REVIEW/report.json'):assert d['all_passed']
 D=pin(ROOT/'reports'/PRI/'exact-comparison-C-torsion-full-HF-and-genuine-target-corollary.json',
   '937053c8184530e5bfe13314e982de1c3fd98db6bdd4da999549b3188058040a')
 pin(ROOT/'code/ambient_B18_C_exact_torsion1123_20261105_v63.py',pr['code_sha256'])
 pin(ROOT/'research/notes'/f'{PRI}.md',pr['note_sha256'])
 assert pr['actual_comparison_C_torsion_exponents']==[1,1,2,3]
 # Explicit abstract Jordan chains provide ranks without using any
 # uncomputed actual C4 action/basis or producer fixture.
 lengths=[4]*4+[1,1,2,3];basis=[(i,j) for i,n in enumerate(lengths) for j in range(n)]
 assert len(basis)==23
 ranks=[sum(max(0,n-k) for n in lengths) for k in range(5)]
 assert ranks==[23,15,9,4,0]
 kernel=[(i,n-1) for i,n in enumerate(lengths)]
 image3=[(i,j) for i,n in enumerate(lengths) for j in range(3,n)]
 assert len(kernel)==8 and len(image3)==4 and set(image3).issubset(kernel)
 genuine=[v for v in kernel if v not in image3];assert len(genuine)==4
 pattern=[1,1,2,3]
 quotients=[[n,4*n+sum(min(n,a) for a in pattern)] for n in range(1,9)]
 assert quotients==D['ell_power_quotient_lengths_first8']
 hf=[3,6,6,6]+[4]*8
 assert [sum(hf[:n]) for n in range(1,9)]==[3,9,15]+[4*n+5 for n in range(4,9)]
 numer=[hf[0]]+[hf[i]-hf[i-1] for i in range(1,len(hf))]
 assert numer==[3,3,0,0,-2]+[0]*7
 proof={'same_actual_comparison_C_not_original_H0_Ktorsion':True,
  'accepted_actual_DVR_type':{'free_rank':4,'torsion_exponents':pattern,'torsion_length':sum(pattern)},
  'full_HF':[3,6,6,6,'4 for every degree>=4'],
  'series':'(3+3t-2t^4)/(1-t)',
  'm2_D_type':'D4 plus (D/ell²)^2','m3_D_type':'D4 plus k²','m4_D_type':'Dfree4',
  'annihilator_proofs':[
    'The exact D-type gives ell³Tor_D C=0.',
    'm4Tor_D C is contained both in m4C and Tor_D C. Since m4C is D-free, their intersection is0. Thus m4Tor_D C=0.',
    'No equality Tor_D(m^nC)=m^nTor_D C, no m3-annihilation, and no original H0 annihilator is assumed.'],
  'positive_n_ell_quotient_lengths':quotients,
  'C4_HF':[3,6,6,6,2,0],'C4_length':23,
  'abstract_Jordan_check_ranks':ranks,'abstract_ker_dimension':8,'abstract_imE3_dimension':4,
  'canonical_genuine_target_proof':[
   'Reduction and quotient induce an R0-linear map C[ell] to ker(E:C4->C4)/im E³.',
   'If ell*x=ell4*y, put t=x-ell³*y. Then ell*t=0. Changing y changes t by ell³C[ell], hence0 since ell³Tor C=0.',
   'Changing the representative x by ell4*z leaves t unchanged with y shifted by ell*z. Changing its ker representative by ell³*w leaves t unchanged with y shifted by w.',
   'This constructs the inverse. Injectivity also follows because C[ell] meets ell³C only in0: ell³C is D-torsionfree after exact exponent3.',
   'Therefore the quotient is canonically genuine C[ell], dimension4. The finite kernel alone is dimension8 and is not the genuine target.'],
  'actual23_actions_basis_and_connecting_values_not_computed':True,
  'no_new_coefficients_physical_gauge_IFT_jet_or_source_replay':True,
  'no_original_H0_nil_cube_index_or_exhaustion_transfer':True}
 ep=out/'independent-exact-C-filtration-annihilators-and-canonical-genuine-target.json';save(ep,proof)
 report={'id':ID,'all_passed':True,'frozen':True,'source_only':True,'input_sha256':pins,
  'code_sha256':sha(__file__),'evidence_sha256':{ep.name:sha(ep)},
  'exact_comparison_corollaries_and_genuine_target_accepted':True,
  'actual_C4_action_interface_not_yet_constructed':True,'elapsed_seconds':time.monotonic()-start}
 save(out/'report.json',report);signal.alarm(0)
 print(json.dumps({'all_passed':True,'report':str(out/'report.json'),'sha256':sha(out/'report.json'),'seconds':report['elapsed_seconds']}))
if __name__=='__main__':main()
