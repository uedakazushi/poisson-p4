#!/usr/bin/env python3
"""Certified full quadratic row rank: no parameter-constant linear consequence."""
import argparse
from fractions import Fraction
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import resource
import subprocess
import time
from s61_affine_normal_global_quintic_20261103_v63 import ROOT,uadd,umul,sha
P=32003

def main():
 parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);parser.add_argument('--work',type=Path,required=True);args=parser.parse_args();out=args.output.resolve();work=args.work.resolve();out.mkdir(parents=True,exist_ok=False);work.mkdir(parents=True,exist_ok=False);start=time.monotonic()
 dep=ROOT/'reports/S61_AFFINE_NORMAL_GLOBAL_QUINTIC_20261103_v63/report.json';old=json.loads(dep.read_text());assert old['all_passed']
 for rel,want in old['input_sha256'].items():assert sha(ROOT/rel)==want,rel
 for rel,want in old['evidence_sha256'].items():assert sha(dep.parent/rel)==want,rel
 source=dep.parent/'presentation.json.gz';data=json.loads(gzip.decompress(source.read_bytes()));decode=lambda poly:{degree:{tuple(m):Fraction(c)for m,c in coeff}for degree,coeff in poly};G=decode(data['monic_quintic_G']);res=list(map(decode,data['all134_remainders']))
 def reduce(poly):
  quotient={};remainder=poly
  while remainder and max(remainder)>=5:
   degree=max(remainder);term={degree-5:remainder[degree]};quotient=uadd(quotient,term);remainder=uadd(remainder,umul(term,G),Fraction(-1))
  assert uadd(umul(quotient,G),remainder)==poly;return remainder
 rows=[reduce({degree+j:p for degree,p in row.items()})for row in res for j in range(5)];assert len(rows)==670
 mons=sorted({m for row in rows for coeff in row.values()for m in coeff});assert()not in mons and all(1<=len(m)<=2 for m in mons)
 quad=[(d,m)for d in range(5)for m in mons if len(m)==2];linear=[(d,m)for d in range(5)for m in mons if len(m)==1];columns=quad+linear;index={label:i for i,label in enumerate(columns)};entries=[]
 for i,row in enumerate(rows):
  for degree,coeff in row.items():
   for mon,c in coeff.items():
    assert c.denominator%P;v=c.numerator*pow(c.denominator,-1,P)%P
    if v:entries.append((i,index[degree,mon],v))
 assert len(quad)==5145 and len(linear)==225 and len(entries)==237811
 raw=(f'670 {len(quad)} {len(linear)} {len(entries)}\n'+'\n'.join(' '.join(map(str,row))for row in entries)+'\n').encode();inputpath=work/'input.txt';inputpath.write_bytes(raw);(out/'input-mod32003.txt.gz').write_bytes(gzip.compress(raw,mtime=0))
 cpp=ROOT/'code/s61_affine_normal_global_quintic_rank_20261103_v63.cpp';native=work/'native';flags=subprocess.check_output(['pkg-config','--cflags','--libs','flint'],text=True).split();include=subprocess.check_output(['pkg-config','--variable=includedir','flint'],text=True).strip();subprocess.run(['g++','-O2','-I'+include,str(cpp),'-o',str(native),*flags],check=True)
 result=subprocess.run([str(native),str(inputpath)],capture_output=True,text=True,check=True,timeout=60,preexec_fn=lambda:resource.setrlimit(resource.RLIMIT_AS,(2200*1024*1024,2200*1024*1024)));(out/'native.log').write_text(result.stdout+result.stderr)
 lines=result.stdout.splitlines();header=lines[0].split();assert header[:5]==['RANK_QUADRATIC','670','RANK_TOTAL','670','MINOR_DETERMINANT'];det=int(header[5]);assert det%P;selected=list(map(int,lines[1].split()[1:]));assert len(selected)==len(set(selected))==670 and max(selected)<5145
 certificate=dict(prime=P,module_rows=670,target_basis_h_degrees=list(range(5)),parameter_variables=data['coefficient_symbols'],all_entry_parameter_degree_atmost2=True,no_constant_coefficients=True,quadratic_columns=5145,linear_columns=225,nonzero_modular_entries=237811,quadratic_rank_mod_prime=670,total_rank_mod_prime=670,
  exact_QQ_and_generic_beta_quadratic_rank=670,nonzero_QQ_punit_minor_rows=list(range(670)),minor_columns_zero_based=selected,minor_column_labels=[[d,list(m)]for d,m in(quad[j]for j in selected)],minor_determinant_mod_prime=det,
  exact_scope='No nonzero parameter-constant combination of all670 vectorpolynomials hasdegree<=1;cofactorsh^j withKcoefficientsj<=4included.',
  not_claimed=['parameterdependentcofactorimpossibility','unit','finiteoverbeta','genericfiveexhaustion'])
 (out/'certificate.json').write_text(json.dumps(certificate,indent=2)+'\n')
 paths=[Path(__file__),cpp,dep,source,ROOT/'research/notes/S61_AFFINE_NORMAL_QUINTIC_CONSTANT_BARRIER_20261103_v63.md'];report=dict(id='S61_AFFINE_NORMAL_QUINTIC_CONSTANT_BARRIER_20261103_v63',all_passed=True,explicit670_QQ_punit_minor_from_modular_certificate=True,generic_parameterconstant_purelinear_span_zero=True,original_lowdegree_q_barrier_not_repeated=True,unit_or_exhaustion_claimed=False,input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},evidence_sha256={p.name:sha(p)for p in out.iterdir()if p.is_file()},seconds=time.monotonic()-start)
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print('GLOBAL5x670 QUADRATICRANK670;GENERICCONSTANTLINEARSPAN0 BARRIER PASS',report['seconds'],flush=True)
if __name__=='__main__':main()
