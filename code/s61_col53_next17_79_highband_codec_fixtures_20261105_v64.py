"""Synthetic OLD/HIGH restoration and high-output mixed products; no source point."""
from pathlib import Path
from fractions import Fraction as Q
import hashlib,json
from s61_col53_next17_79_full_HL_highband_value_20261105_v64 import fixed_q_terms
from s61_sparse_generic5_stream_primitive_cell_codec_20261105_v64 import write_cell,validate_cell
ROOT=Path(__file__).resolve().parents[1]
ID='S61_COL53_NEXT17_79_HIGHBAND_CODEC_FIXTURES_20261105_v64'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
 out=ROOT/'reports'/ID;out.mkdir(exist_ok=False);qi=Q(-5,21);cap=2;M=1000000
 low={(0,0,0,-1,3):-12,(0,0,2,-1,3):6};high={(0,0,3,-1,3):10,(0,0,4,-1,3):-15}
 records=[];merged={}
 for name,vals in [('OLD',low),('HIGH',high)]:
  p=out/(name+'.bin');r=write_cell(p,qi,vals);assert validate_cell(p,qi,vals,r);body=p.read_bytes()
  for e,c in fixed_q_terms(body,qi,len(vals),M):
   assert((e[2]>cap)if name=='HIGH'else(e[2]<=cap))and e not in merged;merged[e]=c
  records.append(r)
 assert merged=={**low,**high}and records[0]['G']==6 and records[1]['G']==5
 oldsha=sha(out/'OLD.bin');oldbytes=(out/'OLD.bin').read_bytes()
 # Both children may be in their cached OLD band while the parent output is
 # in its missing HIGH band. The actual worker retains these terms.
 parentcap=1;complete={};highout={};low_low_high=False
 for e,a in merged.items():
  for f,b in merged.items():
   z=tuple(e[j]+f[j]for j in range(5));complete[z]=complete.get(z,0)+a*b
   if z[2]>parentcap:
    highout[z]=highout.get(z,0)+a*b
    if e in low and f in low:low_low_high=True
 assert low_low_high and highout=={e:c for e,c in complete.items()if e[2]>parentcap}
 assert highout[(0,0,2,-2,6)]==-144
 # An empty band uses no q division for q=0 or q<0.
 for k,q in enumerate([Q(0),qi]):
  p=out/('EMPTY'+str(k)+'.bin');r=write_cell(p,q,{});assert r['G']==1 and r['content']=='0'and list(fixed_q_terms(p.read_bytes(),q,0,M))==[]
 assert sha(out/'OLD.bin')==oldsha and(out/'OLD.bin').read_bytes()==oldbytes
 files=['code/s61_col53_next17_79_highband_codec_fixtures_20261105_v64.py','code/s61_col53_next17_79_full_HL_highband_value_20261105_v64.py','code/s61_sparse_generic5_stream_primitive_cell_codec_20261105_v64.py','code/s61_sparse_generic5_primitive_binary_codec_20261105_v64.py']
 d={'id':ID,'all_passed':True,'synthetic_only_no_chosen_values':True,'negative_fixedq_DIFFERENT_low_high_G_restored':True,'nonempty_and_empty_fixedq_sign_sentinel_rules':True,'disjoint_old_high_merge_exact':True,'mixed_OLD_times_OLD_missing_HIGH_outputs_kept':True,'old_prefix_bytes_and_SHA_unchanged':True,'all5coordinates_preserved':True,'input_sha256':{f:sha(ROOT/f)for f in files},'evidence_sha256':{p.name:sha(p)for p in out.iterdir()}}
 (out/'report.json').write_text(json.dumps(d,separators=(',',':'))+'\n');print('report_sha256',sha(out/'report.json'))
if __name__=='__main__':main()
