#!/usr/bin/env python3
"""Generic scalar-entry linear span and explicit Fitt4-origin certificate."""
import argparse
from fractions import Fraction
import gzip
import itertools
import json
from pathlib import Path
import resource
import subprocess
import time
from s61_affine_normal_global_quintic_20261103_v63 import ROOT,uadd,umul,sha,add,pmul
from s61_affine_fiber_macaulay import load_fiber
P=32003

def main():
 parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);parser.add_argument('--work',type=Path,required=True);args=parser.parse_args();out=args.output.resolve();work=args.work.resolve();out.mkdir(parents=True,exist_ok=False);work.mkdir(parents=True,exist_ok=False);start=time.monotonic()
 dep=ROOT/'reports/S61_AFFINE_NORMAL_GLOBAL_QUINTIC_20261103_v63/report.json';old=json.loads(dep.read_text());assert old['all_passed']
 for rel,want in old['input_sha256'].items():assert sha(ROOT/rel)==want,rel
 for rel,want in old['evidence_sha256'].items():assert sha(dep.parent/rel)==want,rel
 source=dep.parent/'presentation.json.gz';data=json.loads(gzip.decompress(source.read_bytes()));decode=lambda poly:{degree:{tuple(m):Fraction(c)for m,c in coeff}for degree,coeff in poly};G=decode(data['monic_quintic_G']);res=list(map(decode,data['all134_remainders']))
 def reduce(poly):
  quotient={};rem=poly
  while rem and max(rem)>=5:
   k=max(rem);term={k-5:rem[k]};quotient=uadd(quotient,term);rem=uadd(rem,umul(term,G),Fraction(-1))
  assert uadd(umul(quotient,G),rem)==poly;return rem
 rows=[reduce({degree+j:coeff for degree,coeff in r.items()})for r in res for j in range(5)];assert len(rows)==670
 F,_,_,original,pointsource=load_fiber();framepath=ROOT/'reports/S61_AFFINE_RESIDUAL42_20261103_v63/linear129x42-QQ.json.gz';frame=json.loads(gzip.decompress(framepath.read_bytes()));E=[[Fraction(c)for c in row]for row in frame['embedding']];L=[]
 for i,row in enumerate(E):
  poly={}
  if row[1]:poly[42,]=row[1]
  if row[2]:poly[43,]=row[2]
  if i>=3:poly[i-3,]=Fraction(1)
  L.append(poly)
 Q=[]
 for f in F:
  quadratic={}
  for mon,c in f.items():
   if len(mon)==2:
    for ex,d in pmul(L[mon[0]],L[mon[1]]).items():add(quadratic,ex,c*d)
  Q.append(quadratic)
 for i in range(134):
  assert {m:c for m,c in res[i].get(0,{}).items()if len(m)==2 and 44 not in m}==Q[i]
 for index,row in enumerate(rows):
  i,j=divmod(index,5)
  for target in range(5):
   coeff=row.get(target,{});assert not coeff.get((),0)and all(len(m)<=2 for m in coeff)
   assert {m:c for m,c in coeff.items()if len(m)==2 and 44 not in m}==(Q[i]if target==j else {})
 quad=list(itertools.combinations_with_replacement(range(45),2));monomials=quad+[(i,)for i in range(45)];indices={m:i for i,m in enumerate(monomials)};entries=[];observed=set()
 for i,row in enumerate(rows):
  for target,coeff in row.items():
   for mon,c in coeff.items():
    assert c.denominator%P;observed.add(mon);v=c.numerator*pow(c.denominator,-1,P)%P
    if v:entries.append((i*5+target,indices[mon],v))
 assert len(monomials)==1080 and len(entries)==237811
 raw=(f'3350 1035 45 {len(entries)}\n'+'\n'.join(' '.join(map(str,e))for e in entries)+'\n').encode();inputpath=work/'input.txt';inputpath.write_bytes(raw);(out/'input-mod32003.txt.gz').write_bytes(gzip.compress(raw,mtime=0))
 cpp=ROOT/'code/s61_affine_normal_scalar_fitting_20261103_v63.cpp';native=work/'native';flags=subprocess.check_output(['pkg-config','--cflags','--libs','flint'],text=True).split();include=subprocess.check_output(['pkg-config','--variable=includedir','flint'],text=True).strip();subprocess.run(['g++','-O2','-I'+include,str(cpp),'-o',str(native),*flags],check=True)
 result=subprocess.run([str(native),str(inputpath)],capture_output=True,text=True,check=True,timeout=60,preexec_fn=lambda:resource.setrlimit(resource.RLIMIT_AS,(2200*1024*1024,2200*1024*1024)));(out/'native.log').write_text(result.stdout+result.stderr);lines=result.stdout.splitlines();header=lines[0].split();assert header[:5]==['QUADRANK','179','TOTALRANK','224','DET'];assert int(header[5])==21497;selectedrows=list(map(int,lines[1].split()[1:]));selectedcolumns=list(map(int,lines[2].split()[1:]));assert len(selectedrows)==len(set(selectedrows))==224 and len(selectedcolumns)==len(set(selectedcolumns))==224;assert set(range(1035,1080))<=set(selectedcolumns)
 certificate=dict(all3350_exact_QQ_nonH_quadratic_patterns_checked=True,original134_quadratic_forms_L_bound=True,
  universal_structural_quadratic_span_upper=179,reason='134formsQ_i(L)plus45formsH*x_a;nonHquadraticatcomponentkofh^jr_iisdelta_kj*Q_i(L)',
  all1080_monomials_included=True,full_universe_quadratic_monomials=1035,linear_monomials=45,special_QQ_absent_monomials=[list(m)for m in quad if m not in observed],generic_absence_of_six_not_assumed=True,
  minor_rows=selectedrows,minor_columns=selectedcolumns,minor_column_labels=[list(monomials[j])for j in selectedcolumns],minor_det_mod32003=21497,
  exact_generic_scalar_quadratic_rank=179,exact_generic_total_scalar_rank=224,all45_parameter_linears_in_scalar_entry_constant_span=True,
  cofactor_recipe='Invertselected224x224coefficientminoroverK(beta);solveall45linearcoordinatevectorsonselectedcolumns. Structuralrankupper224impliesallremaining1080coefficientsmatch.',
  actual_interpretation='Fitt_4(A)=(n42,S1,S2,H)overformalR45;scalarENTRYspanonly,NOTvectorrowspanorannihilatorA',
  fitting_ideal_formulas={'Fitt_j,j1..4':'m^(5-j)','Fitt0':'m^5 intersect I_B'},
  no_42normal_basechange_constant_span_inference=True,no_parameter_annihilation_or_unit_or_exhaustion=True)
 (out/'certificate.json').write_text(json.dumps(certificate,indent=2)+'\n')
 paths=[Path(__file__),cpp,dep,source,framepath,original,pointsource,ROOT/'research/notes/S61_AFFINE_NORMAL_SCALAR_FITTING_20261103_v63.md'];report=dict(id='S61_AFFINE_NORMAL_SCALAR_FITTING_20261103_v63',all_passed=True,generic_structural_rank_upper179_proved=True,explicit224_punit_minor_saved=True,all45_generic_scalar_entry_linear_span_proved=True,actual_Fitt4_origin_certificate=True,no_A_annihilation_or_unit_claim=True,input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},evidence_sha256={p.name:sha(p)for p in out.iterdir()if p.is_file()},seconds=time.monotonic()-start)
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print('SCALAR FITT4/GENERIC179+45SPAN/224MINOR21497 PASS',report['seconds'],flush=True)
if __name__=='__main__':main()
