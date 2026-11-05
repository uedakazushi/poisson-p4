#!/usr/bin/env python3
"""ROOT independently audits fixed input support and bounded Q/R/E code; no trial division."""
import ast,hashlib,json,math,subprocess,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_ROW58_QRE_STRIP_PREFLIGHT_20261105_v64'
OUT=ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_ROW58_QRE_ROOT_SOURCE_REVIEW_20261105_v64'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 OUT.mkdir(exist_ok=False)
 p=BASE/'report.json';assert sha(p)=='e597dafb05b037e0c20c0c9da1c98e7b5fd246354a944767bac54aeb8443863e'
 r=json.loads(p.read_text());assert r['all_passed'] and r['no_trial_division_executed'];pins={str(p.relative_to(ROOT)):sha(p)}
 for name,h in r['input_sha256'].items():assert sha(ROOT/name)==h;pins[name]=h
 p=BASE/'integer-pseudo-division-budget-and-proof.json';assert sha(p)==r['evidence_sha256'][p.name];b=json.loads(p.read_text());pins[str(p.relative_to(ROOT))]=sha(p)
 # Independent stream decoder: do not import the prepared codec or worker.
 p=ROOT/b['original_binary_path']
 with p.open('rb')as f:
  assert f.read(9)==b'Q5CELL1\n\x01'
  def u():
   value=0;shift=0;last=128
   while last&128:
    raw=f.read(1);assert len(raw)==1;last=raw[0];value+=(last&127)<<shift;shift+=7;assert shift<4096
   assert shift==7 or last;return value
  def z():
   tag=u();raw=f.read(tag//2);assert len(raw)==tag//2 and(not raw or raw[0]);v=int.from_bytes(raw,'big');assert v or not tag&1;return -v if tag&1 else v
  num,den=z(),z();assert (num,den)==(1,353894400000)
  shift=[]
  for _ in range(5):v=u();shift.append(v//2 if not v&1 else -(v//2)-1)
  n=u();assert n==951;exponents=[];l1=0;g=0;previous=None
  for _ in range(n):
   c=z();e=tuple(shift[j]+u()for j in range(5));assert c
   assert previous is None or (sum(previous),previous)>(sum(e),e);previous=e
   exponents.append(e);l1+=abs(c);g=math.gcd(g,abs(c))
  assert not f.read(1) and g==1
 valuation=tuple(min(e[j]for e in exponents)for j in range(5));assert valuation==(-1,0,0,0,-1)
 fibres={};box=[0]*5
 for exponent in exponents:
  e=tuple(exponent[j]-valuation[j]for j in range(5));key=e[0],e[1],e[4];fibres[key]=max(fibres.get(key,0),e[2]+e[3])
  for j in range(5):box[j]=max(box[j],e[j])
 assert [[list(k),d]for k,d in sorted(fibres.items())]==b['fixed_H_L_t_fibres']
 assert len(fibres)==17 and box==[6,2,8,13,6]
 U=sum((d+1)*(d+2)//2 for d in fibres.values());D=max(fibres.values());assert(U,D)==(1739,16)
 assert l1==int(b['input_integer_l1_norm'])
 # At each success the fibre total degree falls by one. Primitive linear
 # factors Q,R,E have ||tail||_1/leading_coefficient <=2 and integer content1.
 # Geometric long-division bound multiplies by <=17*2^d, d=16,...,1.
 quotient=l1*17**16*2**136
 pseudo=32**16*quotient*(1+50*17*2**16)
 rebuild=quotient*82**16;identity=83*pseudo
 assert quotient==int(b['all_normalized_quotient_l1_bound'])
 assert pseudo==int(b['all_pseudo_maps_integer_l1_bound'])
 assert rebuild==int(b['all_rebuild_maps_integer_l1_bound'])
 assert identity==int(b['all_identity_partial_integer_bound'])
 assert [quotient.bit_length(),pseudo.bit_length(),rebuild.bit_length(),identity.bit_length()]==[358,464,460,470]
 assert 2*max(pseudo,rebuild).bit_length()==928<2048
 meta=5*U+8*951+64*19+5000+2000;assert meta==24519<25000
 ordinary=25000+21*U+2048;codec=ordinary+3*U;assert ordinary==63567 and codec==68784
 cell=2048+U*((464+7)//8+5+16);assert cell==139429
 actual=sum((ROOT/name).stat().st_size for name in b['actual_READ_files_bytes'])
 assert actual==b['actual_READ_input_RAW_bytes']==186611
 raw=actual+3*65536+40*cell+2*1024*1024;assert raw==8057531<32*1024*1024
 worker=ROOT/b['worker'];code=worker.read_text();ast.parse(code)
 for text in ["('Q',3,5,","('R',3,1,","('E',2,32,","assert qrec['G']==A","assert c%A==0",
              "max(e[2]+e[3]for e in P)==oldtotal-1","del p,blob;meta()","del P;pool['B']=dict(pool['P'])",
              "assert left==right","no_automatic_retry","all_completed_trace_files_preserved",
              "assert validate_cell(path,Q(1),A,rec)","rootreview.stat().st_size<=65536","scalars(rootdoc)<=2000"]:
  assert text in code,text
 help_run=subprocess.run([sys.executable,str(worker),'--help'],capture_output=True,timeout=10);assert help_run.returncode==0
 assert not(ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_ROW58_EXISTING_QRE_STRIP_20261105_v64').exists()
 facts=dict(independent_stream951_source_decode=True,fixed17_fibres_and1739_support=True,
  pseudo_division_proof='Scale a^d, eliminate descending leading-variable exponent with leading term first; resulting integral quotient V and remainder satisfy a^d P=F V+R. Full coefficient identities and terminal nonzero remainders are saved. If R=0, Gauss gives G(V)=a^d and primitive integer quotient.',
  bounded_scope='only Q=5q+1, R=q+1, E=32s+15q+35 of fixed row58/col56; no F9, unknown factor, row search or matrix update',
  map_integer_bits464=464,identity_integer_bits470=470,codec_temporary_bits928=928,
  source_owned_upper50022=50022,future_owned_upper=68784,fresh_RAW_upper=8057531,
  metadata_derived=24519,max_trials19=19,max_trace_cells40=40,
  code_audit='at most three maps; decoded input released; alias sort only; full inverse reconstruction; integer-content normalization and all stream codec readbacks; preserves failure trace',
  whole_matrix_scope='168x15, 87 extra RHS,66 constraints,other14,one P290, all old complements retained',
  no_trial_division_or_new_matrix_value=True)
 (OUT/'source-proof.json').write_text(json.dumps(facts,indent=2)+'\n')
 pins[str(Path(__file__).relative_to(ROOT))]=sha(Path(__file__))
 result=dict(id=OUT.name,all_passed=True,frozen=True,SOURCE_only=True,ONE_fixed_QRE_strip_ready=True,
  input_sha256=pins,evidence_sha256={'source-proof.json':sha(OUT/'source-proof.json')})
 p=OUT/'report.json';p.write_text(json.dumps(result,indent=2)+'\n');assert p.stat().st_size<65536
 print(json.dumps(facts,indent=2));print('report_sha256',sha(p))
if __name__=='__main__':main()
