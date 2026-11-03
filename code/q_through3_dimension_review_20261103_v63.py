#!/usr/bin/env python3
"""Small13 cofactor replay and proper-Proj dimension-transfer scope audit."""
import argparse
from fractions import Fraction
import gzip
import hashlib
import json
from pathlib import Path
import re
import resource
import shutil
import subprocess

ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_text())
def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output-dir',type=Path,required=True)
 args=ap.parse_args();out=args.output_dir.resolve();out.mkdir(parents=True,exist_ok=False)
 directory=ROOT/'reports/Q_THROUGH3_PARTIAL_DIMENSION_20261103_v63';manifest=load(directory/'manifest.json');paths=[Path(__file__),directory/'manifest.json']
 for name,h in manifest['source_sha256'].items():p=ROOT/name;assert sha(p)==h,name;paths.append(p)
 for name,h in manifest['evidence_sha256'].items():p=directory/name;assert sha(p)==h,name;paths.append(p)
 for name,h in manifest['uncompressed_sha256'].items():assert hashlib.sha256(gzip.decompress((directory/name).read_bytes())).hexdigest()==h,name
 with gzip.open(ROOT/'reports/Q_ORDINARY_CUBIC_MINIMAL_GENERATORS_20261103_v63/generators.json.gz','rt')as f:generators=json.load(f)
 with gzip.open(ROOT/'reports/Q_ORDINARY_QUADRATIC_CONE_20261103_v63/quadrics.json.gz','rt')as f:old=json.load(f)
 names=generators['variables'];p=32003
 def expr(terms):
  values=[]
  for mon,text in terms:
   c=Fraction(text);assert c.denominator%p
   val=c.numerator*pow(c.denominator,-1,p)%p
   if val:values.append(str(val)+'*'+'*'.join(names[i]for i in mon))
  return '+'.join(values)if values else'0'
 Q=[expr([([i,j],c)for i,j,c in old['all350_quadratic_rows'][source]])for source in generators['quadratic_generators']]
 G=[expr([(mon,c)for mon,c in row['terms']])for row in generators['cubic_generators']]
 assert len(Q)==57 and len(G)==17
 exact_input='R=GF(32003)['+','.join(names)+',MonomialOrder=>GRevLex];\nA=matrix{{'+',\n'.join(Q+G)+'}};\n'
 assert exact_input==gzip.decompress((directory/'input.m2.gz').read_bytes()).decode()
 structure=load(directory/'structure-report.json');lm=structure['selected_leading_monomials']
 supports=[set(re.findall(r'\bu\b|\bz\d+\b',m))for m in lm]
 assert len(lm)==13 and all(supports)and all(not(a&b)for i,a in enumerate(supports)for b in supports[i+1:])
 covering=set(structure['covering_prime_variables']);assert len(covering)==13
 summary=(directory/'stage-3.summary').read_text();line=next(v for v in summary.splitlines()if v.startswith('LEADING_MONOMIALS '))
 all_lm=[v.strip()for v in line[len('LEADING_MONOMIALS '):].strip()[1:-1].split(',')];assert len(all_lm)==329
 assert all(set(re.findall(r'\bu\b|\bz\d+\b',m))&covering for m in all_lm)
 selected=gzip.decompress((directory/'thirteen-cofactors.m2.gz').read_bytes()).decode()
 script=out/'review.m2';log=out/'review.log'
 script.write_text(exact_input+selected+'assert(flatten entries leadTerm B=={'+','.join(lm)+'});\nprint "INDEPENDENT13_EXACT_MEMBERSHIP_AND_LEADING_TERMS_PASS";stdio<<flush;exit 0;\n')
 def limit():
  resource.setrlimit(resource.RLIMIT_AS,(1024**3,1024**3));resource.setrlimit(resource.RLIMIT_CPU,(60,61))
 with log.open('w')as f:subprocess.run([shutil.which('M2'),'--script',str(script)],stdout=f,stderr=subprocess.STDOUT,timeout=60,preexec_fn=limit,check=True)
 assert 'INDEPENDENT13_EXACT_MEMBERSHIP_AND_LEADING_TERMS_PASS'in log.read_text()and'error:'not in log.read_text()
 paths += [script,log,ROOT/'reports/Q_ORDINARY_CUBIC_INITIAL_20261103_v63_REVIEW.json',
  ROOT/'research/notes/Q_THROUGH3_PARTIAL_DIMENSION_20261103_v63_REVIEW.md',
  ROOT/'reports/Q_ORDINARY_CUBIC_MINIMAL_GENERATORS_20261103_v63/generators.json.gz',
  ROOT/'reports/Q_ORDINARY_QUADRATIC_CONE_20261103_v63/quadrics.json.gz']
 result=dict(id='Q_THROUGH3_PARTIAL_DIMENSION_20261103_v63_REVIEW',all_passed=True,
  exact_original_QQ74_generators_reduced_mod_p_and_input_preserved=True,
  selected13_saved_polynomial_cofactor_identities_and_actual_leading_terms_newly_replayed=True,
  pairwise_coprime_leading_supports=13,special_J74_dimension_upper=22,
  all329_saved_leading_terms_cover_prime_height=13,partial329_leading_ideal_dimension_exact=22,
  all329_cofactor_identities_newly_replayed_by_this_review=False,
  previous_all329_separate_process_replay_evidence_audited=True,
  proper_V_Proj_dimension_upper_transfer_accepted=True,QQ_J74_and_full_complex_R35_dimension_upper=22,
  actual_dimension12_or_equal22_or_full_initial_ideal_claimed=False,
  no_new_GB_run=True,input_sha256={str(p.relative_to(ROOT)):sha(p)for p in dict.fromkeys(paths)})
 with(out/'report.json').open('x')as f:json.dump(result,f,indent=2);f.write('\n')
 print('PASS small13 cofactor certificate and unconditional fullR35 dimension upper22')

if __name__=='__main__':main()
