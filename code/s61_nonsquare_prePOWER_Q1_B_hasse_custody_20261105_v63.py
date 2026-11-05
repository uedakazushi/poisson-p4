#!/usr/bin/env python3
"""Prepare fixed Q1/B S-cache worker/custody; compile and help ONLY, NO B values."""
from pathlib import Path
import argparse,hashlib,json,resource,signal,subprocess,time
ROOT=Path(__file__).resolve().parents[1]
ID='S61_NONSQUARE_PREPOWER_Q1_B_HASSE_WORKER_CUSTODY_20261105_v63'
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60)
 out.mkdir(parents=True,exist_ok=False);pins={}
 def pin(p):pins[str(p.relative_to(ROOT))]=sha(p);return p
 def load(st,name='report.json'):return json.loads(pin(ROOT/'reports'/st/name).read_text())
 def dump(n,x):(out/n).write_text(json.dumps(x,indent=2)+'\n')
 src=load('S61_NONSQUARE_PREPOWER_Q1_B_HASSE_SOURCE_20261105_v63');assert src['all_passed']and src['source_only']
 pre=load('S61_NONSQUARE_PREPOWER_Q1_B_HASSE_SOURCE_20261105_v63','finite-cache-formula-lifetime-height-and-cost.json')
 assert pre['formula_owned37436']==37436 and pre['ALL_product_upper12768337']==12768337
 lemma=load('S61_NONSQUARE_PREPOWER_Q1_B_HASSE_SOURCE_20261105_v63','ordinary-B-valuation-and-no-unknown-chart.json')
 assert lemma['leading_b1_is_NOT_assumed_FULLU_unit']
 jp=load('S61_NONSQUARE_PREPOWER_Q1_J_ROOT_DHORNER_REVIEW_20261105_v63');assert jp['all_passed']
 peer=load('S61_NONSQUARE_PREPOWER_T1_C4_LONGDIVISION_REVIEW_20261105_v63');assert peer['all_passed']and peer['all4811_quotient_coefficients_match']
 inp=pin(ROOT/'reports/S61_NONSQUARE_PREPOWER_T1_C4_LONGDIVISION_REVIEW_20261105_v63/quotient-original-common-D.txt')
 assert sha(inp)==peer['evidence_sha256'][inp.name]
 L=0;mx=0;n=0;prev=None;hist=[0]*19;body=hashlib.sha256()
 with inp.open()as f,(out/'native-input.txt').open('w')as g:
  assert f.readline().strip()=='T1C4QN001';Dline=f.readline();D=int(Dline);assert D>0 and D.bit_length()==197
  assert int(f.readline())==4811;g.write('Q1BH001 4811 18 0\n');g.write(Dline)
  for line in f:
   key,c=map(int,line.split());assert c and(prev is None or prev<key);prev=key
   a,x=divmod(key,1048576);e,d=divmod(x,1024)
   assert 0<=a<=16 and 0<=e<=33 and 0<=d<=18 and a+e+d<=50
   g.write(line);body.update(line.encode());L+=abs(c);mx=max(mx,abs(c).bit_length());n+=1;hist[d]+=1
 assert n==4811 and L.bit_length()<=301
 H=[sum(n*19*(19-l+k)for l,n in enumerate(hist)if l>=k)for k in range(19)]
 N=[sum(n for l,n in enumerate(hist)if l>=k)for k in range(19)]
 root_build=342+8436+13300
 schedules=[root_build+sum(H[:k+1])+2*H[k]+sum(N[:k+1])+2*N[k]for k in range(19)]
 assert schedules==pre['actual_index_schedule_including_threepasses']and max(schedules)==12768337
 cpp=pin(ROOT/'code/s61_nonsquare_Q1_B_hasse_20261105_v63.cpp')
 drv=pin(ROOT/'code/s61_nonsquare_prePOWER_Q1_B_hasse_entry_20261105_v63.py');pin(Path(__file__))
 note=pin(ROOT/f'research/notes/{ID}.md')
 history=[]
 for folder in sorted((ROOT/'reports').glob(ID+'_ATTEMPT*')):
  for p in sorted(folder.rglob('*')):
   if p.is_file():pin(p);history.append(str(p.relative_to(ROOT)))
 cmd=['bash','-c','source /workspace/cloud-setup/activate.sh && c++ -std=c++17 -O2 "$1" -o "$2" -Wl,-rpath,/workspace/cloud-setup/prefix/usr/lib/x86_64-linux-gnu -lflint','Q1B-build',str(cpp),str(out/'native-worker')]
 with(out/'compile.log').open('w')as f:subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,check=True,timeout=30)
 for name,cmd in [('worker-help',[str(out/'native-worker'),'--help']),('driver-help',[str(ROOT/'.venv/bin/python'),str(drv),'--help'])]:
  r=subprocess.run(cmd,text=True,capture_output=True,check=True,timeout=5);(out/(name+'.log')).write_text(r.stdout+r.stderr)
 proof=dict(prepared_only_NO_Bvalue=True,canonical_Q1='T1/(a3*C4), not the originalS61 quadricq1',
  BODY4811_byteidentical=True,body_sha256=body.hexdigest(),D197=D,input_maxbit=mx,input_L1bit=L.bit_length(),
  B='b1*d-R,b1=50a+18e2-9e-2,R=40-60e',b1_is_NOT_assumed_FULLU_unit=True,
  ordinary_valuation='accepted SOURCE primitive/Gauss proof; denominatorD*b1^18 is SYMBOLIC; no b1=0 FULLU component is discarded and no physicalPGL transformation claimed',
  fixed_cache='S_h=R^h*b1^(18-h),19maps; keys a*128+e, a<=18-h,e<=36-h-2a; total3610 upper',
  exact_Hasse='P_k=sum I_ijl*choose(l,k)*S_(l-k)*a^i*e^j; [z^k]Q1(a,e,R/b1+z)=P_k/(D*b1^18)',
  orders0to18_FIRSTnonzero_ONLY=True,
  output_keybox='a<=34,e<=69,a+e<=86-k,key=a*128+e; e69<128 and allinputoffset+root sums cannot carry',
  input_keybox='a<=16,e<=33,d<=18,total50; literalkey a*1048576+e*1024+d unchanged',
  cache_workspace='R0..18 arrays190; currentP_t=b1^t andnextP; 4literalb1terms. For t0..18 direct updateS_(18-t) withP_t*R_(18-t); no product polynomial temporary.',
  exact_P_lifetime='after completingnextP, clearoldP then swapnextP, counts updated; simultaneouscurrent+next<=324+361=685; ALLR/Pcurrent/Pnext freed BEFORE Hasse accumulator',
  construction_owned34592=34592,main_owned37436=37436,
  ownership='2*(input4811+S3610+oneacc2297)+ALLmetadata16000; rootphase2*(4811+3610+R190+Pold324+Pnext361)+16000',
  metadata16000='choose361,allscalar/constants/counters/exponentdecode/source/cachekeys/cursors,active scaledinputbase,D/L1,temporaryZs,prefix/outputformatbuffers',
  main_source436_hard518=True,root_l1bit_upper120=120,
  source_height='L1I<2^301; correlated79^(18-l)*179^l<=179^18<2^135; all forward or remaining-subset reversepartials<2^436. Rootcache<100^18<2^120.',
  source_index_schedule=schedules,ALL_multiplication_upper12768337=max(schedules),
  scope='includesall19orders throughfirstnonzero plusreverseandforwardrepeat;rootconstruction22078,onechoose*input each eligible source,cacheS products. CPU completion NOT promised.',
  prefix_version='Q1BHP001',prefix_cursor='phase,t,h,Rindex,Rterm,b1term,inputindex,k,Pkey,Skey,base_ready,products,scalings,updates,peak',
  prefix_arrays='input4811partialvector,Rarraysincluding0slots,Pcurrent/Pnext,all19S maps,oneacc,D,L1,activebase; all sizes/owned counts reconstructed byparser',
  pending_semantics='all guards beforemap/arraycommit; pendingcursor recomputes only uncommitted product; activebase stored withbase_ready flag; no source replay; choose is literalmetadata',
  zeroorders='separatecompleted-zero-orders file retained; phase8 clockguard afterall3passes beforeoutput',
  no_IO_doubleprefix='mathSTOP producesONE fullstate; outputdevicefailure producesno secondfullprefix; noatomicdoublebuffer',
  raw_state_upper4642060=4642060,raw_state_reserve5000000=5000000,
  soft45_hard60_RAM2GB=True,history=history,
  no_Bquotient_GCD_otherfactor_tail_or_badlocus=True)
 dump('source-cache-lifetime-height-cursors-and-guards.json',proof)
 files={str(p.relative_to(ROOT)):p for p in out.iterdir()};files.update({p:ROOT/p for p in pins})
 raw=sum(p.stat().st_size for p in files.values())+20000;assert raw<=2800000
 total=raw+5000000+200000;assert total<=8000000
 dump('RAW-closure.json',dict(all_prepared_history_inputs_codes_andmetadata=raw,full_state_OR_output_reserve5000000=5000000,
  actual_allstate_bound4642060=4642060,future_receipts200000=200000,total_upper=total,hard8000000=8000000,
  no_second_fullprefix_or_atomicdoublebuffer=True))
 r=dict(id=ID,all_passed=True,frozen=True,prepared_only=True,value_worker_NOT_run=True,
  canonical4811_BODY_unchanged=True,b1_NOT_assumed_FULLU_unit=True,
  whole_owned37436=37436,source436_hard518=True,ALL_products_upper12768337=12768337,
  RAW_total_upper=total,compiled_help_ONLY=True,fresh_ROOT_source_code_audit_and_ONE_GO_REQUIRED=True,
  input_sha256=pins,evidence_sha256={p.name:sha(p)for p in out.iterdir()},seconds=time.monotonic()-start)
 dump('report.json',r);signal.alarm(0);print(json.dumps({k:v for k,v in r.items()if 'sha256'not in k},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
