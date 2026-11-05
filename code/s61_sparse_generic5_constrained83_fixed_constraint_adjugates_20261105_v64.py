#!/usr/bin/env python3
"""ONE fixed source-gated two5x5 adjugate coefficient worker."""
from pathlib import Path
from fractions import Fraction as Q
import gc,hashlib,json,resource,signal,subprocess,time
ROOT=Path(__file__).resolve().parents[1]
ID='S61_SPARSE_GENERIC5_CONSTRAINED83_FIXED_CONSTRAINT_ADJUGATES_20261105_v64'
PREF='S61_SPARSE_GENERIC5_CONSTRAINED83_EXACT49_FIXED_CONSTRAINT_MINOR_PREFLIGHT_20261105_v64'
def sha(p):
    h=hashlib.sha256()
    with p.open('rb')as f:
        for b in iter(lambda:f.read(65536),b''):h.update(b)
    return h.hexdigest()
def read(p):return json.loads(p.read_text())
def enc(x):return(json.dumps(x,indent=2)+'\n').encode()
def dump(p,x):p.write_bytes(enc(x))
def tree(p):return sum(x.stat().st_size for x in p.rglob('*')if x.is_file())
def numeric(v):
    if isinstance(v,bool):return 0
    if isinstance(v,(int,float)):return 1
    if isinstance(v,list):return sum(map(numeric,v))
    if isinstance(v,dict):return sum(map(numeric,v.values()))
    return 0
def asdict(p):
    if not p:return {}
    q=Q(p['content']);return{tuple(p['shift'][j]+t[j+1]for j in range(5)):q*int(t[0])for t in p['terms']}
def main():
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000))
    out=ROOT/'reports'/ID;work=ROOT/'work'/ID;out.mkdir(parents=True,exist_ok=False);work.mkdir(parents=True,exist_ok=False)
    pd=ROOT/'reports'/PREF;pr=pd/'report.json';assert sha(pr)=='15ec0fa2870a5985257a5a9e62500141ade5619b2ba1afe7b50a9bddcafea06e';r=read(pr);assert r['all_passed']and r['source_only']
    paths={ROOT/p:h for p,h in r['input_sha256'].items()}
    paths[pr]=sha(pr)
    for p,h in r['evidence_sha256'].items():paths[pd/p]=h
    for p,h in paths.items():assert sha(p)==h
    inp=pd/'two-actual-fixed-constraint-minors.json';pp=pd/'support-and-worker-plan.json';bp=pd/'whole-budget-and-scope.json'
    data,plan,budget=map(read,[inp,pp,bp]);assert budget['future_owned_including_metadata_upper']==85813 and budget['future_RAW_including_actual_source_support_receipts_all52_outputs_and3MiB_other_upper']==10622833
    assert plan['actual_metadata_slots']<=25000
    cpp=Path(__file__).with_suffix('.cpp');extra={Path(__file__),cpp,ROOT/f'research/notes/{ID}.md',ROOT/'code/s61_sparse_generic5_row_stream_kernel_20261104_v63.hpp'}
    extra={p for p in extra if p not in paths};external=sum(p.stat().st_size for p in extra);base=sum(p.stat().st_size for p in paths)
    control=work/'control.txt'
    with control.open('w')as f:
        f.write(f'{base} {external}\n')
        for block in data['blocks']:
            for i,scale in enumerate(block['row_clearings']):
                rho=int(scale['positive_integer_denominator']);shift=scale['parameter_monomial_exponent']
                for p in block['literal_H5x5'][i]+block['literal_B5x2'][i]:
                    values=asdict(p);f.write(str(len(values))+'\n')
                    for e,c in values.items():
                        z=c*rho;assert z.denominator==1;ee=[e[j]+shift[j]for j in range(5)];assert min(ee)>=0
                        f.write(str(z.numerator)+' '+' '.join(map(str,ee))+'\n')
        for b,block in enumerate(plan['plans']):
            for j,p in enumerate(block['outputs']):
                f.write(f"{b} {j} {p['literal_output_RAW_upper']}\n")
                rr,cc=p['rows'],p['columns'];f.write(str(len(rr))+' '+' '.join(map(str,rr+cc))+'\n');f.write(str(len(p['DP_nodes']))+'\n')
                for mask,n,terms in p['DP_nodes']:
                    f.write(f'{mask} {n} {len(terms)}\n')
                    for x in terms:f.write(' '.join(map(str,x))+'\n')
        for b,block in enumerate(plan['plans']):
            for order,i,j,terms in block['fixed_adjugate_identity_support_plans']:
                f.write(f"{b} {int(order=='AdjH')} {i} {j}\n")
                for x in terms:f.write(' '.join(map(str,x))+'\n')
    outputcap=budget['literal52_coefficient_outputs_RAW_upper'];othercap=10622833-base-outputcap
    # Exact helper/source versions and controller are included in the ledger.
    # Reserved future bytes cover70 source-clearing files (20k each), compact
    # DP/identity logs, native log/summary, final checks/report and receipts.
    future_other=70*20000+200000+100000+20000+100000+300000
    pre={'id':ID,'all_passed':True,'one_fixed_coefficient_worker':True,'source_preflight_sha256':sha(pr),
         'actual_READ_input_RAW_bytes':base,'worker_code_header_note_charge':external,'fixed_metadata_upper':25000,
         'whole_owned_upper':85813,'whole_RAW_upper':10622833,'all52_output_RAW_upper':outputcap,'other_RAW_cap':othercap,
         'future_other_artifact_upper':future_other,'hard_seconds':60,'hard_address_space_bytes':2000000000,
         'input_sha256':{**{str(p.relative_to(ROOT)):h for p,h in paths.items()},**{str(p.relative_to(ROOT)):sha(p)for p in extra}},
         'control_sha256':sha(control),'all_original_source_versions_preserved':True}
    dump(out/'preflight.json',pre)
    del data,plan,budget,r,paths;gc.collect()
    binary=work/'native'
    try:
        with(out/'compile.log').open('wb')as log:
            subprocess.run(['bash','-c','source /workspace/cloud-setup/activate.sh; exec c++ -O2 -std=c++17 "$1" -o "$2" $(pkg-config --cflags --libs flint) -lgmp','compile',str(cpp),str(binary)],stdout=log,stderr=subprocess.STDOUT,check=True,timeout=max(1,59-(time.monotonic()-start)))
        assert tree(out)+tree(work)+external+future_other<=othercap
        with(out/'native.log').open('wb')as log:
            subprocess.run([str(binary),str(control),str(out),str(work)],stdout=log,stderr=subprocess.STDOUT,check=True,timeout=max(1,59-(time.monotonic()-start)))
        n=read(out/'native-summary.json');assert n['all_passed']and n['actual_determinants_nonzero']and n['all52_outputs']==52 and n['all100_left_right_identities']==100
        # Independent Python source clearing and coefficient-support/height
        # readback, one polynomial at a time; no determinant recomputation.
        data=read(inp);clears=0
        for block in data['blocks']:
            for i,scale in enumerate(block['row_clearings']):
                for j,p in enumerate(block['literal_H5x5'][i]+block['literal_B5x2'][i]):
                    expected={tuple(e[k]+scale['parameter_monomial_exponent'][k]for k in range(5)):v*int(scale['positive_integer_denominator'])for e,v in asdict(p).items()}
                    actual=asdict(read(out/f"{block['name']}-cleared-{i}-{j}.json"));assert expected==actual and all(v.denominator==1 for v in actual.values());clears+=1
        del data;gc.collect();full=read(pp)
        thin=[(b['name'],p['output'],p['support_upper'],p['support_path'],p['support_sha256'],p['integer_coefficient_abs_upper'],p['literal_output_RAW_upper'])for b in full['plans']for p in b['outputs']]
        del full;gc.collect();results=[];post_peak=0
        for name,label,upper,path,h,height,cap in thin:
            filename=f'{name}-det.json'if label=='det'else f"{name}-adj-{label.split('-')[2]}-{label.split('-')[1]}.json"
            p=out/filename;q=read(p);value=asdict(q);sp=ROOT/path;assert sha(sp)==h;supp={tuple(x)for x in read(sp)}
            post_peak=max(post_peak,numeric(q)+5*len(supp)+len(value)+numeric(pre)+numeric(thin)+4096);assert post_peak<85813
            assert set(value)<=supp and len(value)<=upper and all(v.denominator==1 and abs(v.numerator)<=int(height)for v in value.values());assert p.stat().st_size<=cap
            if label=='det':assert value
            results.append({'file':filename,'sha256':sha(p),'terms':len(value),'source_support_and_integer_height_checked':True,'bytes':p.stat().st_size})
        ids=[json.loads(s)for s in(out/'all100-adjugate-identities.jsonl').read_text().splitlines()];assert len(ids)==100 and all(x[-1]is True for x in ids)
        for p,h in pre['input_sha256'].items():assert sha(ROOT/p)==h
        assert sha(control)==pre['control_sha256']and not list(out.rglob('*.partial'))
        check={'all_passed':True,'literal_source_clearing_checks':clears,'output_files':results,'post_readback_owned_numeric_bound':post_peak,
               'adjugate_convention':'adj-i-j is signed cofactor deleting source row j and source column i',
               'inverse_formula':'x_p=Adj(Htilde)*(row_scaled_rhs-Btilde*x_f)/det(Htilde)',
               'all_parameters_only':['H','L','s0','q0','t0'],'normal_unknown_inverted':False,'top111x24_not_substituted':True,
               'specialization_scope':'Generic K5; require both cleared determinant numerators and old monomial denominators nonzero. Complementary divisors remain open.'}
        dump(out/'source-clearings-and-all52-output-checks.json',check)
        report={'id':ID,'all_passed':True,'actual_fixed_determinants_nonzero':True,'all52_outputs':52,'all100_left_right_adjugate_identities':100,
                'all70_literal_source_clearings':clears,'fixed_rows_oldL':[6,14,16,41,42],'determinants_base_only':True,
                'owned_whole_peak':max(n['whole_owned_upper_checked'],n['prospective_owned_upper'],post_peak),'fixed_metadata_charge':25000,
                'actual_polynomial_output_terms':n['actual_output_terms'],'future_top_remains_factored':True,
                'no_rank_search_GCD_point_prime_or_complete_inverse':True,'input_sha256':pre['input_sha256'],
                'evidence_sha256':{str(p.relative_to(out)):sha(p)for p in out.rglob('*')if p.is_file()},
                'ignored_executable_sha256':sha(binary),'native_seconds':n['seconds'],'total_seconds':time.monotonic()-start,
                'whole_RAW_upper':10622833,'actual_fresh_RAW':0}
        rp=out/'report.json';report['actual_fresh_RAW']=base+external+tree(out)+tree(work)+len(enc(report))
        for _ in range(5):
            z=base+external+tree(out)+tree(work)+len(enc(report))
            if report['actual_fresh_RAW']==z:break
            report['actual_fresh_RAW']=z
        assert report['actual_fresh_RAW']<=10622833 and report['owned_whole_peak']<=85813 and report['total_seconds']<60
        dump(rp,report);assert base+external+tree(out)+tree(work)==report['actual_fresh_RAW']
        print(json.dumps({k:v for k,v in report.items()if k not in('input_sha256','evidence_sha256')}));print('report_sha256',sha(rp))
    except Exception as e:
        dump(out/'failure.json',{'all_passed':False,'failure':repr(e),'seconds':time.monotonic()-start,'all_complete_outputs_preserved':True,'no_automatic_retry':True});raise
    finally:signal.alarm(0)
if __name__=='__main__':main()
