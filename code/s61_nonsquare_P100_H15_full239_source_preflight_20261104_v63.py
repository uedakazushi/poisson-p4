#!/usr/bin/env python3
"""SOURCE/SUPPORT only for all319 Gamma rows on the entire finite H15 algebra."""
import argparse,collections,gzip,hashlib,json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
ID='S61_NONSQUARE_P100_H15_FULL239_SOURCE_PREFLIGHT_20261104_v63'
P=32003
def read(p):return json.loads(p.read_text())
def sha(p):
 h=hashlib.sha256()
 with p.open('rb')as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(parents=True,exist_ok=False);R=ROOT/'reports';paths=[]
 def load(stem,name='certificate.json'):
  p=R/stem/name;paths.append(p);return read(p)
 fac=load('S61_NONSQUARE_P100_H15_MODP_FACTOR_INVENTORY_20261104_v63','native.json');fr=load('S61_NONSQUARE_P100_H15_MODP_FACTOR_INVENTORY_20261104_v63','report.json')
 assert fac['all_passed']and fac['Hbar_squarefree']and max(x['degree']for x in fac['factors'])==4
 typed=load('S61_NONSQUARE_P100_PURE90_TYPED_SPLIT_20261104_v63','operator.json')
 cat=load('S61_NONSQUARE_QUARTIC_319_ROWS_20261103_v63');catalogue=[cat['source_rows'][i]for i in cat['retained_row_indices']]
 assert len(catalogue)==319
 counts=collections.Counter(row[0]for row in catalogue)
 oracle=load('S61_NONSQUARE_LAZY_QUARTIC_ORACLE_20261103_v63','oracle.json');ar=oracle['arithmetic'];nodes=ar['nodes']
 fd=load('S61_NONSQUARE_JOINT_LELL_FREE164_20261104_v63');NF=fd['full_target_projection239_to122'];LN=fd['linear_nodes'];SN=fd['scalar_circuit']['nodes'];SA=fd['scalar_circuit']['atoms_num_den']
 h=load('S61_NONSQUARE_H_SQUARE_H1_PREFLIGHT_20261104_v63')
 ri=load('S61_NONSQUARE_P100_H1_RIGHT_INVERSE_20261104_v63','native.json');L=ri['right_inverse112x23']
 assert ri['all23_source_identities_both_source_frames_verified'] and ri['all_cofactor_denominators_are_whole_curve_units']
 manifest=load('S61_NONSQUARE_P100_PURE98_SOURCE_PREFLIGHT_20261104_v63','manifest.json')
 kernel=load('S61_NONSQUARE_P100_PURE_CUBIC_KERNEL_20261104_v63','source-binding.json')
 labels=load('S61_NONSQUARE_ORACLE_FRAME_ADAPTER_20261103_v63','adapter.json')
 frame=load('S61_NONSQUARE_P100_G18_PURE7_20261104_v63','joint-cubic-frame.json')
 quad=load('S61_BINARY_HYPERPLANE_QUADRATIC_CODIMENSION_20261102', 'unused.json') if False else None
 quadp=R/'S61_BINARY_HYPERPLANE_QUADRATIC_CODIMENSION_20261102.json';paths.append(quadp);quad=read(quadp)
 pre=load('S61_NONSQUARE_P100_PURE70_SOURCE_PREFLIGHT_20261104_v63')
 gamma=typed['nodes'][typed['exports']['full_Gamma319x98']]
 assert gamma['shape']==[319,98]
 assert catalogue[174:286]==h['all112_cubic_source319_definitions']
 assert len(catalogue[:174])==174 and all(row[0]=='quadric'for row in catalogue[:174])
 assert all(row[0]=='cubic_normalized'for row in catalogue[174:286])
 assert len(catalogue[286:])==33 and len(h['all_remaining21_quartic_H2free_labels'])==21
 assert pre['source_quartics21_needed_later']['H2free_source_labels21']==h['all_remaining21_quartic_H2free_labels']
 raw33=oracle['genuine_quartic_original_rows'];roots=oracle['genuine_quartic_roots'];assert len(raw33)==len(roots)==33
 # Complete original source bindings, not an unconstructed70+15+4 kernel.
 sources=dict(original319_catalogue=catalogue,retained_original354_indices=cat['retained_row_indices'],
  counts=dict(counts),actual_source_blocks=dict(quadratic_products=[0,174],cubic_products=[174,286],raw_true_quartics=[286,319]),
  full_Gamma_export=typed['exports']['full_Gamma319x98'],full_Gamma_operator_path='reports/S61_NONSQUARE_P100_PURE90_TYPED_SPLIT_20261104_v63/operator.json',
  physical_to_joint_before_NF=True,pure98_target122_labels=typed['context']['pure98_target122_labels'],
  Gamma_formula='f0-(f2/U2)*chi0-(f1-(f2/U2)*chi1)*(Lactual^T*Ddelta7*PiOrig)',
  H2_pair_actual135=[77,131],U2='(50*a*B*C)^2',delta7='-64000*C^2',
  raw_quartic_actual135=raw33,raw_quartic_polynomial_roots=roots,
  raw_oracle_normalization='Qhat_l=F_l(T34*diag(Tden)*v)/D_H^2=(50/B)^2 F_l(T34*diag(Tden)*v)',
  H2free21_labels=h['all_remaining21_quartic_H2free_labels'],H2free21_cache_indices=list(range(113,134)),
  H2free21_recipe='psi=-b131*Qhat77+b77*Qhat131; others Qhat_l-(b_l/U2)*chi, then actual H1 section with full tails.',
  all33_raw_sources_retained_simpler_for_direct_full_test=True,
  normalized_quadratic6_sparse_path=str(quadp.relative_to(ROOT)),normalized_quadratic6_sparse_key='quadratic_rows_sparse',
  physical_light_inverse_path='reports/S61_NONSQUARE_P100_G18_PURE7_20261104_v63/joint-cubic-frame.json#/physical_light_inverse7x7',
  H2free21_pure210_roots_NOT_prepared=True,
  all_cubic_source_frame_original16_then_Ddelta7=True,no_unconstructed70_plus15_plus4_assumed=True)
 # Exact scalar/linear closure of NF122, with no coefficient arithmetic.
 lneed=set();sroots=set();todo=NF['output122_linear_roots'][:]
 while todo:
  i=todo.pop()
  if i in lneed:continue
  lneed.add(i);n=LN[i]
  if n[0]=='l':
   for c,j in n[1]:sroots.add(c);todo.append(j)
  else:assert n[0]in['z','v']
 sneed=set();atoms=set();todo=list(sroots)
 while todo:
  i=todo.pop()
  if i in sneed:continue
  sneed.add(i);n=SN[i]
  if n[0]=='c':atoms.add(n[1])
  elif n[0]in['a','s','m']:todo+=n[1:]
  else:assert n[0]=='z'
 HS=h['scalar_circuit']['nodes'];HA=h['scalar_circuit']['atoms_num_den'];hneed=set();hatoms=set();todo=[r for row in h['H2_source_pair_change_scalar_roots']for r in row]+[h['H2_known_unit_root']]
 while todo:
  i=todo.pop()
  if i in hneed:continue
  hneed.add(i);n=HS[i]
  if n[0]=='c':hatoms.add(n[1])
  elif n[0]in['a','s','m']:todo+=n[1:]
  else:assert n[0]=='z'
 # Purely combinatorial monomial support of original oracle polynomials.
 # No coefficient node, core value, normal form or field rank is evaluated.
 supp=[];conv_pairs=0;linear_terms=0
 for i,n in enumerate(nodes):
  if n[0]=='z':v=set()
  elif n[0]=='v':v={tuple(int(j==n[1])for j in range(8))}
  elif n[0]=='m':
   a,b=supp[n[1]],supp[n[2]];conv_pairs+=len(a)*len(b);v={tuple(x+y for x,y in zip(m,k))for m in a for k in b}
  else:
   assert n[0]=='l';linear_terms+=sum(len(supp[j])for c,j in n[1]);v=set().union(*(supp[j]for c,j in n[1]))
  assert len(v)<=239 and all(sum(m[j]*ar['variable_weights'][j]for j in range(8))<=4 for m in v);supp.append(v)
 def schedule(rs):
  need=set();use=collections.Counter();todo=rs[:];aa=set()
  while todo:
   i=todo.pop()
   if i in need:continue
   need.add(i);n=nodes[i]
   children=n[1:]if n[0]=='m'else[j for c,j in n[1]]if n[0]=='l'else[]
   if n[0]=='l':aa.update(c for c,j in n[1])
   for j in children:use[j]+=1;todo.append(j)
  use.update(rs);live=peak=0;paircount=addcount=0
  for i in sorted(need):
   live+=len(supp[i]);peak=max(peak,live);n=nodes[i]
   children=n[1:]if n[0]=='m'else[j for c,j in n[1]]if n[0]=='l'else[]
   if n[0]=='m':paircount+=len(supp[n[1]])*len(supp[n[2]])
   if n[0]=='l':addcount+=sum(len(supp[j])for c,j in n[1])
   for j in children:
    use[j]-=1
    if not use[j]:live-=len(supp[j])
  return dict(nodes=len(need),needed_atoms=len(aa),polynomial_coefficient_peak=peak,
   retained_root_coefficients=sum(len(supp[r])for r in rs),multiply_pair_count_upper=paircount,scaled_add_count_upper=addcount)
 combined=schedule(roots);individual=[dict(actual135=i,**schedule([r]))for i,r in zip(raw33,roots)]
 # Stream original57x91 and inspect only original16 actual normalized pure84.
 p=R/'S61_BINARY_NONSQUARE_CUBIC_ROW_BASIS_20261103_v62/remaining.txt.gz';paths.append(p)
 bp=R/'S61_BINARY_NONSQUARE_CUBIC_ROW_BASIS_20261103_v62/report.json';paths.append(bp);br=read(bp)
 sp=R/'S61_BINARY_NONSQUARE_CUBIC_REPLAY_20261103_v62/source-report.json';paths.append(sp);sm=read(sp)['normal_monomials'];pure=[j for j,m in enumerate(sm)if m[5]==0]
 assert len(pure)==84
 hh=hashlib.sha256();terms=0;maxterms=0;qqbad=[];rowterms=[]
 with gzip.open(p,'rb')as zf:
  head=zf.readline();hh.update(head);assert head.strip()==b'57 91 57'
  for row in range(57):
   rt=0
   for col in range(91):
    line=zf.readline();hh.update(line)
    if row>=41 and col in pure:
     a=line.decode().split();n=int(a[0]);assert len(a)==1+4*n;rt+=n;maxterms=max(maxterms,n)
     for k in range(n):
      c=Q(a[1+4*k])
      if c.denominator%P==0:qqbad.append([row,col,k])
   if row>=41:rowterms.append(rt);terms+=rt
  assert not zf.read()
 assert hh.hexdigest()==br['remaining_matrix_sha256'] and terms==139746
 # The finite rank proposal must FIRST check integrality of all literal inputs.
 def badterms(polys):return sum(Q(c).denominator%P==0 for p in polys for m,c in p)
 integer_bad=dict(oracle_atoms=badterms([p for nd in ar['atoms']for p in nd]),
  NF122_atoms=badterms([p for i in atoms for p in SA[i]]),
  H2pair_atoms=badterms([p for i in hatoms for p in HA[i]]),
  original_quadratic6=badterms([p for row in quad['quadratic_rows_sparse']for mon,p in row]),
  H1_rightinverse_RF=badterms([p for row in L for r in row for p in [r['numerator'],r['denominator']]]),
  physical_light_inverse_RF=badterms([p for row in frame['physical_light_inverse7x7']for r in row for p in [r['numerator'],r['denominator']]]),
  cubic_pure1344=len(qqbad))
 # Pole comparison is literal factor-list arithmetic: an irreducible factor
 # can divide a known linear pole only if that factor itself is linear.
 poles=[(1,0),(1,20),(1,60),(203,2700),(1,12)];pole_receipt=[]
 for i,f in enumerate(fac['factors']):
  tests=[]
  for a,b in poles:
   if f['degree']>1:tests.append(dict(linear=[a,b],coprime_by_irreducible_degree=True))
   else:
    c=(b-a*f['coefficients'][0])%P;tests.append(dict(linear=[a,b],remainder_constant=c,unit=c!=0))
  pole_receipt.append(dict(factor_index=i,tests=tests))
 badpole=[i for i,r in enumerate(pole_receipt)if any(t.get('unit')is False for t in r['tests'])]
 dmax=max(f['degree']for f in fac['factors']);section=23*98*dmax;rankbasis=98*98*dmax;minorcopy=rankbasis
 fixed=(len(ar['atoms'])+len(sneed)+len(lneed)+239+122+1344+504+len(atoms)+len(hneed)+len(hatoms)+174)*dmax
 allroots_peak=combined['polynomial_coefficient_peak']*dmax
 one_peak=max(x['polynomial_coefficient_peak']for x in individual)*dmax
 phase1=dict(label='coefficient preparation, no rank basis resident',section23x98=section,
  oracle_all33_polynomial_peak=allroots_peak,oracle_current_one_peak=one_peak,fixed_conservative_field_slots=fixed,
  all33_plus_fixed_upper=section+allroots_peak+fixed,one_root_plus_fixed_upper=section+one_peak+fixed,
  streamed_literal_core_coefficient_reserve=maxterms+512,rank_basis_resident=False)
 phase2=dict(label='source coefficients released; stream Gamma rows into rank basis',basis98x98=rankbasis,
  selected_original98x98_copy=minorcopy,current_two_rows=2*98*dmax,field_inverse_and_temp_reserve=4096,
  full_live_upper=rankbasis+minorcopy+2*98*dmax+4096,
  all_seven_selected_minors_basis_slots=98*98*sum(f['degree']for f in fac['factors']))
 bounds=dict(max_factor_degree=dmax,source_cubic139746_terms_not_retained=True,
  cubic_source_total_terms=terms,cubic_source_max_one_polynomial_terms=maxterms,cubic_source16_row_terms=rowterms,
  NF122_linear_nodes=len(lneed),NF122_scalar_nodes=len(sneed),NF122_atoms=len(atoms),
  H2pair_scalar_nodes=len(hneed),H2pair_atoms=len(hatoms),
  oracle_all33_support_schedule=combined,oracle_individual_support_schedules=individual,
  phase1=phase1,phase2=phase2,
  output_selected_minor_uint16_bytes_upper=2*phase2['all_seven_selected_minors_basis_slots'],
  temporary_all319_Gamma_uint16_bytes_all_factors_upper=2*319*98*15,
  no_all319_field_matrix_resident=True,
  all_coefficients_and_ranks_NOT_evaluated=True)
 cert=dict(source=sources,bounds=bounds,QQ_coefficient_denominator_bad_counts=integer_bad,
  actual_source_denominator_evaluation_NOT_performed=True,
  five_known_curve_pole_factors=poles,pole_factor_receipt=pole_receipt,nonunit_pole_factor_indices=badpole,
  future_gate='Evaluate all required scalar denominator polynomials in each finite factor; certify every inversion. Known curve factors have the factor-list proof above; QQ denominator checks alone do not justify arbitrary parameter inversions.',
  positive_only_theorem='If all seven full Gamma factors have row span Fq^98, CRT gives surjectivity over Fp[z]/Hbar. A monic finite-free Z_(p)[z]/H model and the certified curve denominators permit finite Nakayama, hence QQ full239 at every point of this finite chosen cover. A deficient factor or pole gives no QQ rank failure and cannot be silently discarded.',
  no_unconstructed70_plus15_plus4_kernel=True,no_new_source_coefficient_NF_rank_point_prime_GB=True)
 (out/'certificate.json').write_text(json.dumps(cert,indent=2)+'\n')
 paths += [Path(__file__),ROOT/f'research/notes/{ID}.md',ROOT/'code/s61_nonsquare_oracle_frame_adapter_20261103_v63.py',ROOT/'code/s61_nonsquare_lazy_quartic_oracle_20261103_v63.py']
 report=dict(id=ID,source_only=True,frozen=True,complete_original319_and_H2_H1_binding=True,
  max_factor_degree=dmax,coefficient_or_rank_job_executed=False,all_literal_QQ_denominators_prime_units=not any(integer_bad.values()),
  all_five_curve_poles_avoided_mod_H=not badpole,coefficient_phase_all33_live_bound=phase1['all33_plus_fixed_upper']+maxterms+512,
  coefficient_phase_one_root_live_bound=phase1['one_root_plus_fixed_upper']+maxterms+512,
  rank_phase_live_bound=phase2['full_live_upper'],
  input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},certificate_sha256=sha(out/'certificate.json'),seconds=time.monotonic()-start)
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024;signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256']},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output.resolve())
