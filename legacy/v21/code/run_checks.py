"""Reproduce v21 checks in disposable copies; archived editions stay read-only."""
from pathlib import Path
import argparse
import json
import shutil
import subprocess
import sys
import tempfile
import time

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--output-directory',type=Path,required=True)
    ap.add_argument('--regenerate-linear-models',action='store_true')
    args=ap.parse_args()
    bundle=Path(__file__).resolve().parents[1]
    out=args.output_directory.resolve();out.mkdir(exist_ok=True,parents=True)
    (out/'logs').mkdir(exist_ok=True)
    reports=[]
    with tempfile.TemporaryDirectory(prefix='p4_v21_checks_') as td:
        work=Path(td)
        shutil.copytree(bundle/'proofs',work/'proofs')
        shutil.copytree(bundle/'data/tree_linear_models',work/'proofs/reducible/tree_linear_models')
        prior=bundle/'previous_v20/data/quartic_models'
        p=work/'proofs'
        jobs=[
          ('nodal_nilpotency',[sys.executable,str(p/'nodal/verify_nodal_nilpotency.py'),'--model',str(prior/'nodal_model.json')]),
          ('cuspidal_nilpotency',[sys.executable,str(p/'cuspidal/verify_cusp_nilpotency.py'),'--model',str(prior/'cuspidal_model.json')]),
          ('ci_31',[sys.executable,str(p/'nodal/verify_ci_31.py')]),
          ('ci_22',[sys.executable,str(p/'cuspidal/verify_two_conics.py')]),
          ('ci_211_1111',[sys.executable,str(p/'nodal/verify_211_1111.py')]),
          ('independent_star',[sys.executable,str(p/'cuspidal/verify_star_integrability.py')]),
          ('ci_representatives',[sys.executable,str(p/'reducible/verify_ci_representatives.py')]),
          ('linear_tree_audit',[sys.executable,str(p/'reducible/audit_tree_models.py')]),
          ('cusp_unit_certificate',[sys.executable,str(p/'cuspidal/verify_unit_certificate.py'),'--model',str(prior/'cuspidal_model.json')]),
          ('stdlib_full_jacobi',[sys.executable,str(bundle/'code/audit_tree_jacobi_sparse.py'),str(bundle/'data/tree_jacobi_models'),'--output-directory',str(out/'full_jacobi')]),
        ]
        if args.regenerate_linear_models:
            jobs.append(('regenerate_linear_models',[sys.executable,str(bundle/'code/build_p4_tree_linear.py'),'--output-directory',str(work/'regenerated')]))
        for name,command in jobs:
            start=time.monotonic()
            r=subprocess.run(command,cwd=work,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
            (out/'logs'/f'{name}.txt').write_text(r.stdout)
            report={'name':name,'exit_code':r.returncode,'seconds':round(time.monotonic()-start,3),'passed':r.returncode==0}
            reports.append(report);print(json.dumps(report),flush=True)
            if r.returncode:raise RuntimeError('Verification failed: '+name)
            try:
                value=json.loads(r.stdout)
                if isinstance(value,dict):(out/(name+'.json')).write_text(json.dumps(value,indent=2)+'\n')
            except json.JSONDecodeError:pass
        if args.regenerate_linear_models:
            count=0
            for path in (work/'regenerated').glob('*.json'):
                if path.name=='summary.json':continue
                assert json.loads(path.read_text())==json.loads((bundle/'data/tree_linear_models'/path.name).read_text()),path.name
                count+=1
            assert count==8
        shutil.copy2(p/'reducible/tree_linear_models/direct_audit.json',out/'linear_tree_audit.json')
    summary={'all_passed':all(r['passed'] for r in reports),'groups':reports,
      'regenerated_linear_models':args.regenerate_linear_models,
      'scope':'Exact algebra and stored-model verification; existence and component decomposition for the eight residual supports remain open.'}
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps({'all_passed':summary['all_passed'],'groups':len(reports)}))

if __name__=='__main__':main()
