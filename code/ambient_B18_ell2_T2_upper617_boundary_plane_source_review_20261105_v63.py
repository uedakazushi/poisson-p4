#!/usr/bin/env python3
"""Literal617 ancestor/catalog and plane9 source descent; no higher values."""
import json,resource,signal
from fractions import Fraction as Q
from pathlib import Path
from ambient_B18_ell2_T2_mu_review_20261105_v63 import fields
from ambient_B18_ell2_T2_nonactive171_source_review_20261105_v63 import sha,sp
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_ELL2_T2_UPPER617_BOUNDARY_PLANE_SOURCE_20261105_v63_REVIEW'
def main():
 signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024));pins={}
 def read(p,keys=None):pins[p]={'path':p,'sha256':sha(ROOT/p)};return fields(ROOT/p,keys) if keys else json.loads((ROOT/p).read_text())
 b='reports/AMBIENT_B18_ELL2_T2_UPPER617_P8_ANCESTOR_SOURCE_20261105_v63';r=read(b+'/report.json');assert r['all_passed'] and pins[b+'/report.json']['sha256']=='6dfd725a4901e40dac12456060453573f64559cc920ab842c527a60a95e67065'
 for p,h in r['input_sha256'].items():assert sha(ROOT/p)==h
 ep=b+'/literal86-P8-composite-program-and-exact-boundary-catalog.json';d=read(ep);assert pins[ep]['sha256']==r['evidence_sha256'][Path(ep).name]
 nodes=dict(d['selected_original_C0_arithmetic_nodes']);bd={i:n for i,kind,n in d['exact_literal_cache_boundary_nodes']};assert len(nodes)==532 and len(bd)==85 and len(d['root_nodes'])==86 and not d['uncached_or_unrecognized_opaque_boundaries'] and not d['source_inverse_nodes_remain_symbolic']
 needed=dict(nodes);needed.update(bd);cp='reports/AMBIENT_B18_T0_HONEST_C0_20261104_v63/honest-C0-D0-D4-DAG.json';text=(ROOT/cp).read_text();pins[cp]={'path':cp,'sha256':sha(ROOT/cp)};assert len(text.encode())<12*1024*1024
 start=text.index('[',text.index('"DAG_scalar_nodes"'))+1;depth=1;quoted=False;escape=False;idx=-1;node_start=None;found={}
 for pos in range(start,len(text)):
  c=text[pos]
  if quoted:
   if escape:escape=False
   elif c=='\\':escape=True
   elif c=='"':quoted=False
   continue
  if c=='"':quoted=True
  elif c=='[':
   depth+=1
   if depth==2:node_start=pos;idx+=1
  elif c==']':
   if depth==2 and idx in needed:found[idx]=json.loads(text[node_start:pos+1])
   depth-=1
   if depth==1 and idx>=max(needed):break
 assert found==needed;del text
 reached=set()
 def visit(i):
  if i in reached:return
  reached.add(i)
  if i in bd:return
  n=nodes[i];assert n[0] in ['c','a','m','n']
  kids=[] if n[0]=='c' else n[1:]
  for k in kids:assert isinstance(k,int) and k<i;visit(k)
 for i in d['root_nodes']:visit(i)
 assert reached==set(needed)
 f='reports/AMBIENT_B18_ELL2_T2_MU_20261105_v63_REVIEW/phase1-independent-full9-first640-source-proof.json';known79={i for i,v in read(f,['first640'])['first640']['all79_atom_values_full10']}
 f='reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63/actual-source105-ell2-raw33-values.json';known164={i for i,v in read(f,['all_requested164_root_values'])['all_requested164_root_values']}
 physical=sorted(set(bd)&known164);only79=sorted(set(bd)-known164);overlap=sorted(set(bd)&known164&known79);assert len(physical)==45 and len(only79)==40 and len(overlap)==28 and set(only79)<=known79 and not(set(bd)-known79-known164)
 sup='reports/AMBIENT_B18_ELL2_T2_UPPER617_SUPPORT_PREFLIGHT_20261105_v63';s=read(sup+'/factored617-lower-kernel-direct-third-support-and-source-frontier.json',['pruned_typed_source_closure','all617_source_and_direct_third_support_recipes']);atoms={n[1] for i,n in s['pruned_typed_source_closure'] if n[0]=='actual_C0_scalar_function'};assert len(atoms)==59 and len(atoms&known164)==31 and len(atoms-known164)==28 and atoms<=known79|known164
 # Recheck the actual ring map using only fixed conductor-ring constants.
 # No scalar function value, third derivative or proper-gauge coefficient is used.
 raw='reports/AMBIENT_B18_ELL2_RAW33_QUOTIENT_20261105_v63';rr=read(raw+'/report.json');assert rr['all_passed'];assert read(raw+'_REVIEW/report.json')['all_passed']
 a=read(raw+'/actual-ell2-product33-and-raw45-actions.json',['actual33_basis','coordinate_actions9_sparse_columns']);basis=a['actual33_basis'];A=[list(map(sp,cols)) for cols in a['coordinate_actions9_sparse_columns']];assert len(basis)==33 and len(A)==9
 plane=[(deg-e,e) for deg in range(5) for e in range(min(1,deg)+1)];pi={m:i for i,m in enumerate(plane)};P=[{pi[(r,e)]:Q(1)} if sector==0 else {} for sector,r,e in basis];assert len(plane)==9 and P[0]=={0:Q(1)}
 plane_source_indices=[next(i for i,p in enumerate(P) if p=={j:Q(1)}) for j in range(9)];assert len(set(plane_source_indices))==9
 def mapped(v):
  z={}
  for i,x in v.items():
   for j,c in P[i].items():
    y=z.get(j,Q(0))+x*c
    if y:z[j]=y
    else:z.pop(j,None)
  return z
 for coord in range(9):
  for i,p in enumerate(P):
   expected={}
   for j,x in p.items():
    r0,e=plane[j]
    if coord<7:continue
    m=(r0+int(coord==7),e+int(coord==8))
    if m in pi:expected[pi[m]]=x
   assert mapped(A[coord][i])==expected
 # These297 intertwining identities and the cyclic monomial source basis
 # imply all1089 ring products, by applying each basis monomial successively.
 for sector,rpow,epow in basis:assert 0<=sector<=7 and rpow+epow+int(sector>0)<=4 and epow<=1
 scope={'certified':'Onlyliteral164scalarR0functions, orarithmeticcompositeswhoseentireboundaryhascertifieddescent, maybereadtotheplane9quotient. Itpreservespureplaneordinaryr^aell^b witha+b<=4,b<=1, inparticularr^2ell. All297coordinateintertwinersandcyclicmonomialbasisproveafullringmap(all1089products).',
 'uncertified':'The40only79P8boundariesand28only79typedatomsarenotcertifiedhigherfunctionvaluesbycatalogmembershiporabsenceofinversearithmeticnodes. Keephighestthirdsymbolsindependentuntilformalweightsvanishorseparatefunctiondescentisproved. Overlap79-and164hasstronger164certificateandmustnotbemisclassified.',
 'slot_limit':'Thisdoesnotcertifymixedzthirds,pureell^2coefficients,orders>=5,formalregularquotientthirds,nativeN/gaugecoefficients,orfullupperresponse. TheC72/A18m3quotientkillsr^2ellandcannotreplaceplane9forathird-orderquery. No scalarhighercoefficienthasbeensubstituted.'}
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);ve=out/'independent617-literal-catalog-split-and-plane9-ring-source-proof.json';ve.write_text(json.dumps({'all532_arithmetic_and85_boundary_original_nodes_match':True,'all86_root_closure_exact':True,'P8_45_literal164_boundary_ids':physical,'P8_40_only79_boundary_ids':only79,'P8_28_overlap79_and164':overlap,'typed31_literal164_atoms':sorted(atoms&known164),'typed28_only79_atoms':sorted(atoms-known164),'plane9_monomial_order':plane,'raw33_plane9_source_basis_indices':plane_source_indices,'all297_coordinate_intertwiners_checked':True,'all1089_ring_products_follow_from_intertwiners_and_cyclic_basis':True,'no_higher_scalar_or_gauge_coefficients_queried':True,'scope':scope},separators=(',',':'))+'\n');report={'id':STEM,'all_passed':True,'literal_source_and_descent_scope_only':True,'P8_literal164_count':45,'P8_only79_count':40,'typed_literal164_count':31,'typed_only79_count':28,'input_sha256':pins,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ve.name:sha(ve)}};(out/'report.json').write_text(json.dumps(report,separators=(',',':'))+'\n');print(json.dumps({'all_passed':True,'report_sha256':sha(out/'report.json')}))
if __name__=='__main__':main()
