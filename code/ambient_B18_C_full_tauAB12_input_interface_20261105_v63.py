#!/usr/bin/env python3
"""Literal cached39 scalar input reader, grouped by one source file at a time."""
import hashlib,json,re
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def count(x):
 if isinstance(x,(int,float)):return 1
 if isinstance(x,str):return int(bool(re.fullmatch(r'-?\d+(?:/\d+)?',x)))
 if isinstance(x,(list,tuple)):return sum(map(count,x))
 if isinstance(x,dict):return sum(map(count,x.values()))
 return 0
def read39(source,circuit,parse):
 groups={}
 for row in circuit['provided_inputs']:groups.setdefault(row['source_path'],[]).append(row)
 values={};checks=[];peak=raw=0
 for path,rows in sorted(groups.items()):
  p=ROOT/path;assert sha(p)==source['input_sha256'][path]
  raw+=p.stat().st_size;doc=json.loads(p.read_text());peak=max(peak,count(doc))
  for row in rows:
   kind,*args=row['key']
   if kind=='delta':
    physical,index=args;record=next(r for r in doc['physical12_coefficients']if r['physical_row']==physical)
    actual=dict(record['delta']).get(index,'0')
   elif kind=='nu_sparse33':
    root,index=args;record=next(v for i,v in doc['all_requested164_root_values']if i==root)
    actual=dict(record).get(index,'0')
   elif kind=='nu_ell2':actual=doc['actual_nuA_B_pure_ell3'][args[0]][2]
   elif kind=='lambda':actual=next(r['delta_r_ell2']for r in doc['all3_lambda_readbacks']if r['physical_row']==args[0])
   elif kind=='mu03':actual=next(r['delta_ell3']for r in doc['all3_pivot_mu03_readbacks']if r['physical_row']==args[0])
   else:raise AssertionError(kind)
   q=parse(actual);assert q==parse(row['value']);values[row['node']]=q
   checks.append([row['node'],row['key'],str(q),path,row['locator']])
  del doc
 assert len(values)==len(checks)==39
 return values,checks,{'largest_single_parsed_source_numeric_count':peak,'literal_source_RAW_bytes':raw,
  'source_files':list(groups),'provided39_leaves_only_no_circuit_evaluation':True}
