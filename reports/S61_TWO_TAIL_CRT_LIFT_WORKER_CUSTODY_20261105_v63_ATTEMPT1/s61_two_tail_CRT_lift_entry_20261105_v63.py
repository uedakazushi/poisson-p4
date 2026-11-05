#!/usr/bin/env python3
"""ONE fixed40 CRT candidate/exact verification; ONLY after fresh ROOT code/custody audit and GO."""
from pathlib import Path
import argparse,hashlib,json,math,struct,subprocess,time
ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'reports/S61_TWO_TAIL_CRT_LIFT_WORKER_CUSTODY_20261105_v63'
ID='S61_TWO_TAIL_CRT_LIFT_ENTRY_20261105_v63'
WRAD=1<<26
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
class Reader:
 def __init__(self,p):self.f=p.open('rb');self.raw=hashlib.sha256();self.canonical=hashlib.sha256()
 def take(self,n):
  b=self.f.read(n);assert len(b)==n,'truncated binary';self.raw.update(b);return b
 def fixed(self,b):assert self.take(len(b))==b;self.canonical.update(b)
 def num(self,n):
  b=self.take(n);v=int.from_bytes(b,'little');self.canonical.update(v.to_bytes(n,'little'));return v
 def u8(self):return self.num(1)
 def u32(self):return self.num(4)
 def u64(self):return self.num(8)
 def text(self):
  n=self.u32();assert n<=4096;b=self.take(n);s=b.decode();self.canonical.update(s.encode());return s
 def z(self,bits):
  sign=self.u8();assert sign<=1;n=self.u32();assert n<=(bits+7)//8;b=self.take(n);v=int.from_bytes(b,'little');assert v.bit_length()<=bits
  bb=v.to_bytes((v.bit_length()+7)//8,'little');assert b==bb and(sign==0 or v),'noncanonical sign/magnitude';self.canonical.update(bb);return -v if sign else v
 def end(self):
  assert not self.f.read(1);self.f.close();assert self.raw.hexdigest()==self.canonical.hexdigest();return self.raw.hexdigest()
def inspect_prefix(path):
 r=Reader(path);r.fixed(b'TTCRT001');why=r.text();u=[r.u64()for _ in range(64)];z=[r.z(8193)for _ in range(14)];ledger=[[r.u64()for _ in range(8)]for __ in range(40)];paths=[r.text()for _ in range(3)]
 inputs=[]
 for _ in range(2):
  cap=r.u32();n=r.u32();assert n<=cap<=6817;prev=None
  for __ in range(n):
   key=r.u32();c=r.z(301);assert c and(prev is None or prev<key);prev=key
  inputs.append([cap,n])
 fv=[]
 for _ in range(2):
  cap=r.u32();n=r.u32();assert n<=cap<=10545
  for __ in range(n):assert r.u32()<32768
  fv.append([cap,n])
 cc=r.u32();cn=r.u32();assert cn<=cc<=3900
 for _ in range(cn):r.z(600)
 maps=[];max_map_bits=0
 for capbits in [600,8192,8192]:
  n=r.u32();assert n<=12021;prev=None
  for _ in range(n):
   key=r.u64();c=r.z(capbits);assert c and(prev is None or prev<key);prev=key;max_map_bits=max(max_map_bits,abs(c).bit_length())
  maps.append(n)
 owned=16000+2*(inputs[0][0]+inputs[1][0]+sum(maps))+fv[0][0]+fv[1][0]+cc;assert owned<=60000
 digest=r.end();assert digest==sha(path)
 return dict(version='TTCRT001',reason=why,mode=u[0],fixed_prime_index=u[1],all64_cursors=u,
  input_capacity_size=inputs,field_capacity_size=fv,CRT_capacity_size=[cc,cn],H_R_Q=maps,
  current_owned=owned,max_map_integer_bits=max_map_bits,max_scalar_bits=max(abs(x).bit_length()for x in z),
  completed_prime_ledger=ledger,source_paths=paths,binary_roundtrip_SHA256=digest,
  all_state_retained=True,no_auto_resume=True,
  restore='versioned LE typed schema: all64 cursors/14 scalars/40x8 ledger/sourcepaths/inputs/f/g/CRT/H/R/Q; pending scalar may be8193bits before commit, storedcoeffs<=8192')
def map_file(path,magic,which,keep=False):
 r=Reader(path);r.fixed(magic);assert r.u32()==which;n=r.u32();prev=None;stats=dict(terms=n,max_bits=0,max_a=0,max_e=0,max_d=0,max_W=0,total_degree=0);entries=[];content=0
 for _ in range(n):
  key=r.u64();c=r.z(600 if magic==b'TTFAC001'else8192);assert c and(prev is None or prev<key);prev=key
  a,rem=divmod(key%WRAD,1<<20);e,d=divmod(rem,1<<10);w=key//WRAD;assert w==2*a+e
  stats['max_bits']=max(stats['max_bits'],abs(c).bit_length());stats['max_a']=max(stats['max_a'],a);stats['max_e']=max(stats['max_e'],e);stats['max_d']=max(stats['max_d'],d);stats['max_W']=max(stats['max_W'],w);stats['total_degree']=max(stats['total_degree'],a+e+d)
  if magic==b'TTFAC001':assert a<=6 and e<=9 and d<=4 and a+e+d<=11 and(w<13 or(a,e,d)==(6,1,4));content=math.gcd(content,c)
  elif which==0:assert a<=10 and d<=14 and w<=29
  else:assert a<=6 and d<=10 and w<=22
  if keep:entries.append((a,e,d,c))
 digest=r.end();stats['SHA256']=digest
 if keep:assert content==1 and entries[-1][:3]==(6,1,4)and entries[-1][3]>0
 return stats,entries
def main(out):
 start=time.monotonic();out.mkdir(parents=True,exist_ok=False);pre=json.loads((BASE/'report.json').read_text());assert pre['all_passed']and pre['arithmetic_worker_NOT_run']
 pins={str((BASE/'report.json').relative_to(ROOT)):sha(BASE/'report.json')}
 for n,h in pre['input_sha256'].items():
  p=ROOT/n;assert sha(p)==h;pins[n]=h
 for n,h in pre['evidence_sha256'].items():
  p=BASE/n;assert sha(p)==h;pins[str(p.relative_to(ROOT))]=h
 pins[str(Path(__file__).relative_to(ROOT))]=sha(Path(__file__))
 base=sum((ROOT/n).stat().st_size for n in pins)+20000;assert base<=3*2**20
 (out/'phase.json').write_text(json.dumps(dict(phase='BEFORE_FIXED40',source_pins=pins,ONLY_fixed_candidate_lane=True),indent=2)+'\n')
 cmd=[str(BASE/'native-worker'),'--run',*(str(ROOT/x)for x in pre['math_input_paths']),str(out),str(base)];assert len(cmd)==7
 with(out/'native.log').open('w')as log:run=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,timeout=58)
 for n,h in pins.items():assert sha(ROOT/n)==h,n
 np=out/'native-report.json'
 try:n=json.loads(np.read_text())
 except(Exception)as e:n=dict(all_passed=False,common_char0_factor_proved=False,prefix_saved=(out/'TTCRT001-prefix.bin').exists(),native_report_IO_recovery=str(e),returncode=run.returncode)
 if run.returncode or not n.get('all_passed'):
  assert not n.get('common_char0_factor_proved')
  if(out/'TTCRT001-prefix.bin').exists():
   try:n['prefix_inspection']=inspect_prefix(out/'TTCRT001-prefix.bin')
   except(Exception)as e:n['prefix_IO_or_parse_failure']=str(e)
 else:
  assert n['common_char0_factor_proved']and n['candidate_LM_checked']and n['verified_originals']==[1,1]and n['peak_owned']<=60000
  hf,H=map_file(out/'factor.bin',b'TTFAC001',2,True);assert hf['terms']==n['factor_terms']and hf['terms']<=207
  q0,_=map_file(out/'cofactor-N0.bin',b'TTCOF001',0);q2,_=map_file(out/'cofactor-N2.bin',b'TTCOF001',1)
  # No new gcd: verify the reconstructed factor reduces to the EXISTING seed.
  p=32003;lc=H[-1][3]%p;assert lc;il=pow(lc,p-2,p);hm={646*a+19*e+d:(c%p)*il%p for a,e,d,c in H if c%p};seed={}
  with(ROOT/pre['math_input_paths'][2]).open()as f:
   assert int(f.readline())==3899
   for line in f:E,c=map(int,line.split());seed[E]=c
  assert hm==seed
  n['factor_and_cofactor_readback']=dict(factor=hf,N0=q0,N2=q2,primitive=True,seed_reduction_verified=True)
 total=base+sum(p.stat().st_size for p in out.iterdir())+20000;assert total<=32*2**20
 result=dict(id=ID,all_passed=bool(n.get('all_passed')),frozen=True,result=n,
  no_irreducibility_gcd_maximality_Fitting_exhaustion_or_fullrank_claim=True,
  RAW_union_with_inputs=total,input_sha256=pins,evidence_sha256={p.name:sha(p)for p in out.iterdir()},seconds=time.monotonic()-start)
 (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items()if 'sha256'not in k},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path);p.add_argument('--inspect-prefix',type=Path);a=p.parse_args()
 if a.inspect_prefix:print(json.dumps(inspect_prefix(a.inspect_prefix),indent=2))
 elif a.output:main(a.output.resolve())
 else:p.error('--output or --inspect-prefix required')
