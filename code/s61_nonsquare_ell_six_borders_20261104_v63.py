#!/usr/bin/env python3
"""ONE bounded exact six-border test from the same frozen parametric cubic source."""
import argparse,gzip,hashlib,itertools,json,math,os,resource,signal,subprocess,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_ELL_SIX_BORDERS_20261104_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def treebytes(p):return sum(f.stat().st_size for f in p.rglob('*')if f.is_file())
def add(d,k,c):
 v=d.get(k,Q())+c
 if v:d[k]=v
 else:d.pop(k,None)
def writepoly(f,p):f.write(str(len(p))+''.join(' '+str(c)+' '+' '.join(map(str,m))for m,c in sorted(p.items(),reverse=True))+'\n')
def main(out):
 started=time.monotonic();deadline=started+60;resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));out.mkdir(parents=True,exist_ok=False)
 prev=ROOT/'reports/S61_NONSQUARE_ELL_COLON_CUBIC11_20261104_v63';rp=prev/'report.json';pr=json.loads(rp.read_text());assert pr['all_saved_source_and_exact_QQ_rank_checks_passed']
 for rel,h in pr['input_sha256'].items():assert sha(ROOT/rel)==h,rel
 nativep=prev/'native.json';assert sha(nativep)==pr['evidence_sha256']['native.json'];n=json.loads(nativep.read_text())
 r10p=ROOT/'reports/S61_NONSQUARE_ELL_CUBIC_RANK10_20261104_v63/report.json';r10=json.loads(r10p.read_text());assert r10['generic_actual_cubic11_rank_atleast10_proved'];c10p=r10p.parent/'certificate.json';assert sha(c10p)==r10['certificate_sha256'];c10=json.loads(c10p.read_text())
 circ=json.loads((prev/'circuit.json').read_text());sourcegz=ROOT/circ['parametric_source'];data=gzip.decompress(sourcegz.read_bytes());assert hashlib.sha256(data).hexdigest()==circ['source_uncompressed_sha256']
 source=out/'frozen57-source.txt';source.write_bytes(data)
 mons=json.loads((ROOT/circ['normal_monomials_source']).read_text())['normal_monomials']
 pure=[]
 for c in itertools.combinations_with_replacement(range(5),3):pure.append(tuple(c.count(j)for j in range(5)))
 pureidx={m:i for i,m in enumerate(pure)};restriction=[]
 for mon in mons:
  terms=[]
  if not mon[3]:
   for k in range(mon[6]+1):
    ex=(mon[0]+mon[6]-k,mon[1]+k,mon[2],mon[4],mon[7])
    target=pureidx[ex]if not mon[5]else 35+ex.index(1)
    terms.append((target,math.comb(mon[6],k)*(-3)**(mon[6]-k)))
  restriction.append(terms)
 base_rows=c10['parametric_nonzero34minor']['rows'];base_cols=c10['parametric_nonzero34minor']['columns'];leftover=[i for i in range(11)if i not in c10['extracted10_rows']];assert leftover==[4]
 lastrow=20+n['remaining11_genuine_rows'][4];remaining_cols=[n['quotient16_columns'][j]for j in range(16)if j not in c10['fixed_columns_from_completed_RREF']]
 rows=base_rows+[lastrow];cols=base_cols+remaining_cols;assert len(rows)==35 and len(set(rows))==35 and len(cols)==len(set(cols))==40
 selection=out/'selection.txt'
 with selection.open('w')as f:
  f.write('35 40\n'+' '.join(map(str,rows))+'\n'+' '.join(map(str,cols))+'\n91\n')
  for v in restriction:f.write(str(len(v))+''.join(' '+str(j)+' '+str(c)for j,c in v)+'\n')
 qdata=json.loads((ROOT/circ['quad_source']).read_text());qrestricted=[]
 for qid in[0,1,3,5]:
  q={}
  for mon,ts in qdata['quadratic_rows_sparse'][qid]:
   assert mon[5]==0
   if mon[3]:continue
   for k in range(mon[6]+1):
    lm=(mon[0]+mon[6]-k,mon[1]+k,mon[2],mon[4],mon[7]);factor=math.comb(mon[6],k)*(-3)**(mon[6]-k)
    p=q.setdefault(lm,{})
    for core,c in ts:add(p,tuple(core),Q(c)*factor)
  qrestricted.append(q)
 qrows=[]
 for q in qrestricted:
  for j in range(5):
   row=[{}for _ in range(35)]
   for m,p in q.items():ex=list(m);ex[j]+=1;row[pureidx[tuple(ex)]]=p
   qrows.append(row)
 qfile=out/'quad20x35.txt'
 with qfile.open('w')as f:
  f.write('20 35\n')
  for row in qrows:
   for p in row:writepoly(f,p)
 cpp=ROOT/'code/s61_nonsquare_ell_six_borders_20261104_v63.cpp';exe=out/'six-borders';flags=subprocess.check_output(['pkg-config','--cflags','--libs','flint'],text=True).split()
 with(out/'compile.log').open('w')as f:subprocess.run(['c++','-O2','-std=c++17',str(cpp),*flags,'-o',str(exe)],check=True,stdout=f,stderr=subprocess.STDOUT,timeout=max(.1,deadline-time.monotonic()))
 stop=None
 with(out/'process.log').open('w')as f:
  proc=subprocess.Popen([str(exe),str(source),str(qfile),str(selection),str(out/'progress.log'),str(out/'remainders.txt')],stdout=f,stderr=subprocess.STDOUT,start_new_session=True)
  while proc.poll()is None:
   if time.monotonic()>=deadline:stop='shared60s_deadline';break
   if treebytes(out)>32*1024*1024:stop='artifact32MiB_guard';break
   time.sleep(.1)
  if stop:os.killpg(proc.pid,signal.SIGKILL);proc.wait()
 exe.unlink(missing_ok=True)
 output=out/'remainders.txt';status=output.read_text().splitlines()[0]if output.exists()and proc.returncode==0 else'blocked_no_identity_conclusion'
 positive=status in['all_six_zero_exact','at_least_one_border_nonzero_exact']
 evidence={'same_fixed34_polynomial_minor_rows':base_rows,'same_fixed34_polynomial_minor_columns':base_cols,'remaining_original_row':lastrow,'six_remaining_original_columns':remaining_cols,
 'selected35x40_rows':rows,'selected35x40_columns':cols,'base34_nonzero_oldpoint_value':c10['parametric_nonzero34minor']['oldpoint_value'],
 'restriction_rule':circ['restriction_v8'],'actual_genuine16_source_rows':circ['genuine_basis_rows'],'quad_source_rows':[0,1,3,5],
 'border_relation':'full35border_j = DeltaQuad19*DeltaHeavy5*R11border11_j; onlyzero/nonzero tested',
 'generic_row_elimination':'invertible overQQ(a,e,d); exact polynomialcontent divisions checked; unknown factors never declaredFULLUunits',
 'native_status':status,'native_returncode':proc.returncode,'budget_stop_reason':stop,'completed_six_border_test':positive,
 'if_allzero_relation_recipe':'fixed34adjugate/Delta34 expresseslastrow; inR11coordinates row4coefficient1. Lift restrictedcubic relation F in(Z,W), writeF=ZA+WB byfirstfactor, multiplyell andsubtract(ellZ)A+(ellW)B toobtainactualquartic kernel. Noquarticcofactor; allgenuinecubiccofactors divisibleell.',
 'kernel_count_not_increased_without_independent_review':True,'generic_only_no_global_unit_or_allpoint_claim':True,'no_newpoint_prime_rank_GB_or_quartic_expansion':True}
 (out/'certificate.json').write_text(json.dumps(evidence,indent=2)+'\n')
 paths=[Path(__file__),cpp,ROOT/'research/notes'/(ID+'.md'),rp,nativep,r10p,c10p,prev/'circuit.json',sourcegz,ROOT/circ['quad_source'],ROOT/circ['normal_monomials_source']]
 report={'id':ID,'source_binding_and_preparation_passed':True,'six_border_test_completed':positive,'status':status,'generic_exact_rank10_proved':status=='all_six_zero_exact','generic_rank11_proved':status=='at_least_one_border_nonzero_exact','blocked_frontier':not positive,'native_returncode':proc.returncode,'budget_stop_reason':stop,'no_automatic_retry':True,'no_newpoint_prime_GB_or_quartic_expansion':True,'no_allpoint_full239_or_unreviewed_kernel_count_claim':True,'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},'artifact_sha256':{str(p.relative_to(out)):sha(p)for p in out.rglob('*')if p.is_file()},'combined_artifact_bytes':treebytes(out),'seconds':time.monotonic()-started,'frozen':True}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','artifact_sha256']},indent=2),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=ROOT/'reports'/ID);a=p.parse_args();main(a.output.resolve())
