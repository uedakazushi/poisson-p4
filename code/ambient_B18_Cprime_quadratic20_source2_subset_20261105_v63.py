#!/usr/bin/env python3
"""Selective immutable source2 custody; never decode the whole FULL9 file."""
import hashlib,json,re,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_CPRIME_QUADRATIC20_SOURCE2_SUBSET_20261105_v63'
SOURCE='reports/AMBIENT_B18_T3_FULL9_FIRSTJET_20261104_v63/full9-firstjet-certificate.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def numeric(v):
 if isinstance(v,bool)or v is None:return 0
 if isinstance(v,(int,float)):return 1
 if isinstance(v,str):return int(bool(re.fullmatch(r'-?\d+(?:/\d+)?',v)))
 if isinstance(v,list):return sum(map(numeric,v))
 if isinstance(v,dict):return sum(map(numeric,v.values()))
 return 0
def main():
 t=time.monotonic();out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);p=ROOT/SOURCE;raw=p.read_bytes();assert sha(p)=='270db857dfc5188c9dfcfdeb2270ba3e5c02cea7e180542b8dc7f1cfa522f492';text=raw.decode();assert text.isascii();decoder=json.JSONDecoder();offsets={}
 def field(name):
  token=json.dumps(name);assert text.count(token)==1;key=text.index(token);start=text.index(':',key)+1
  while text[start].isspace():start+=1
  value,end=decoder.raw_decode(text,start);offsets[name]={'start_byte':start,'end_byte':end,'source_slice_sha256':hashlib.sha256(raw[start:end]).hexdigest()};return value
 block=field('both57_source2');parse_peak=numeric(block);assert parse_peak<30000
 keep=['unknown57','equations57','F9_nodes','f63','determinant','all57_source2jets','all57_equations_all55_coefficients_verified'];families={f:{k:block[f][k]for k in keep}for f in['CA','CB']};del block
 mons=field('source2_monomials');minimal={'source2_monomials':mons,'both57_source2':families};assert len(mons)==55
 for f in families.values():assert len(f['unknown57'])==len(f['equations57'])==len(f['all57_source2jets'])==57 and all(len(v)==55 for i,v in f['all57_source2jets'])
 retained=numeric(minimal);assert retained<7000;subset=out/'minimal-two57-source2-and-labels.json';subset.write_text(json.dumps(minimal,separators=(',',':'))+'\n');assert len(raw)+subset.stat().st_size<32<<20;del text,raw
 report={'id':STEM,'all_passed':True,'status':'EXACT_selective_source_custody_no_mathematical_values','source_path':SOURCE,'source_sha256':sha(p),'source_raw_bytes':p.stat().st_size,'raw_decode_fields_only':offsets,'largest_parsed_block_numeric_leaves':parse_peak,'retained_numeric_leaves_including_labels':retained,'retained_rational_tensor_entries':6270,'discarded_before_other_working_tables':['fixed_J57','fixed_J57_inverse','pure_function2jets','all unrelated FULL9 top-level coefficient arrays'],'whole_FULL9_JSON_never_decoded':True,'custody_phase_rule':'raw source bytes and one selected under30000-numeric-leaf block only in extraction phase; after pruning, future worker reads only this under7000-leaf subset before Fraction arithmetic. Do not sum discarded parse phase with later118472 arithmetic buffers.','prior_profile_reader_scope':'the earlier SOURCE metadata profile used a general loader; this selective custody is required before any value job, without changing its frozen report or claiming a completed value run.','code_sha256':sha(Path(__file__)),'evidence_sha256':{subset.name:sha(subset)},'subset_raw_bytes':subset.stat().st_size,'seconds':time.monotonic()-t}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
if __name__=='__main__':main()
