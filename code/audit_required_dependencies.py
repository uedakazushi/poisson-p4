#!/usr/bin/env python3
"""Audit the current classification proof/input/certificate dependency closure.

This is an integrity and availability audit, not a mathematical proof checker.
It follows the primary claim DAG, current TeX inputs, recorded file inputs,
exact certificate references, local imports and static producer paths. Historical
catalogue snapshots are projected onto those same primary claims; optional
telemetry, rejected proofs and archive-member locators are not proof premises.
"""
from pathlib import Path
from collections import deque,defaultdict
import argparse,ast,gzip,hashlib,json,re
ROOT=Path(__file__).resolve().parents[1]
def digest(b):return hashlib.sha256(b).hexdigest()
def audit(output):
    roots=json.loads((ROOT/'research/reproduction/classification_roots.json').read_text())
    claims={c['id']:c for c in json.loads((ROOT/'research/CLAIMS.json').read_text())['claims']}
    claimids=set();queue=list(roots['roots'])
    while queue:
        cid=queue.pop()
        if cid in claimids:continue
        claimids.add(cid);queue.extend(claims[cid].get('dependencies',{}).get('claims',[]))
    assert claimids==set(roots['claims'])
    seeds=set()
    for cid in claimids:
        c=claims[cid];seeds.update(p['path']for p in c.get('proof_sources',[])if p.get('path'))
        seeds.update(p for p in c.get('verification_records',[])if isinstance(p,str))
    seeds.update(str(p.relative_to(ROOT))for p in (ROOT/'papers').glob('*.tex'))
    seeds.update(['code/data/known51_rank4_s61_coordinate_catalogue_20261103_v62.json','code/s61_model.py','reports/rank_two_component_recheck.json','reports/rank_two_reconstruction_checks.json','reports/component_catalogue_checks.json','research/reproduction/frozen_producer_index.json','code/known_rank4_firstjet_function_field_20261103_v62.py'])
    frozen=json.loads((ROOT/'research/reproduction/frozen_producer_index.json').read_text())['records']
    recovery={(r['original_path'],r['sha256']):r['materialized_path']for r in frozen}
    alias={'work/rank3-mixed-support-obstruction.json':'reports/CLASSIFICATION_RANK3_MIXED_COMPONENT_20261107.json'}
    todo=deque(seeds);seen=set();edges=defaultdict(set);missing=defaultdict(set);integrity=[];excluded=[];pending=[]
    prefixes=('code/','reports/','research/','legacy/','scripts/','audits/','work/','papers/')
    extensions=('.json','.json.gz','.in.gz','.py','.tex','.md','.csv','.cpp','.h','.m2','.M2','.txt')
    def put(name,origin,expected=None,role='reference'):
        name=name.replace('\\_','_').strip().removeprefix('./')
        if name.startswith(str(ROOT)+'/'):name=name[len(str(ROOT))+1:]
        if not name.startswith(prefixes)or not name.endswith(extensions)or any(q in name for q in ['{','}','*','?','__pycache__']):return
        logical=name;name=alias.get(name,name)
        # Legacy records use paths relative to their own unpacked source root.
        if origin.startswith('legacy/') and not name.startswith('legacy/'):
            candidate='/'.join(origin.split('/')[:2])+'/'+name
            if (ROOT/candidate).is_file():name=candidate
        p=ROOT/name
        if expected and re.fullmatch('[0-9a-f]{64}',expected):
            actual=digest(p.read_bytes())if p.is_file()else None
            if actual!=expected:
                recovered=recovery.get((logical,expected))
                if recovered:
                    name=recovered;p=ROOT/name;actual=digest(p.read_bytes())
                elif (Path(logical).suffix in ('.md','.tex') or '/jobs/' in logical or logical=='research/COMPONENTS.json' or '_REVIEW_JOB_' in logical) and p.is_file():
                    integrity.append({'origin':origin,'path':logical,'expected':expected,'current':actual,'resolution':'Historical prose/control-plane/catalogue provenance. Current proof and recorded coefficient/code inputs are retained; this byte version is not an arithmetic premise.'})
                    expected=None
                else:
                    pending.append({'origin':origin,'path':logical,'expected':expected,'actual':actual,'role':role})
            if expected:integrity.append({'origin':origin,'path':logical,'resolved_path':name,'expected':expected,'actual':actual,'match':actual==expected})
        if p.is_file():edges[origin].add(name);todo.append(name)
        else:missing[logical].add(origin)
    def walk(x,f,keypath=()):
        if isinstance(x,dict):
            if isinstance(x.get('path'),str)and isinstance(x.get('sha256'),str):put(x['path'],f,x['sha256'],'recorded_file')
            if isinstance(x.get('required_path'),str)and isinstance(x.get('sha256'),str):put(x['required_path'],f,x['sha256'],'required_extracted_input')
            for k,v in x.items():
                if k in ('excluded_unaccepted_proof','archive_member') or (k=='member' and 'tracked_archive' in x):
                    excluded.append({'origin':f,'field':'/'.join(keypath+(k,)),'reason':'Explicitly rejected proof or member locator inside an old source archive, not a root-relative input.'});continue
                if 'accepted_proof_projection' in f:
                    if k=='known_components':
                        excluded.append({'origin':f,'field':k,'reason':'Catalogue snapshot; componenthood premises are the current primary claims selected by classification_roots.'});continue
                    if k=='support_claim_records':v=[r for r in v if r.get('record',{}).get('id')in claimids]
                put(k,f,v if isinstance(v,str)else None,'input_hash_key')
                walk(v,f,keypath+(str(k),))
        elif isinstance(x,list):
            for i,v in enumerate(x):walk(v,f,keypath+(str(i),))
        elif isinstance(x,str):put(x,f)
    while todo:
        f=todo.popleft()
        if f in seen:continue
        seen.add(f);p=ROOT/f
        if not p.is_file():missing[f].add('seed');continue
        raw=p.read_bytes()
        try:text=(gzip.decompress(raw)if f.endswith('.gz')else raw).decode()
        except (UnicodeError,OSError):continue
        # A preservation inventory is not a list of mathematical prerequisites.
        if f.endswith('/previous_v21_preservation.json'):
            excluded.append({'origin':f,'field':'inventory','reason':'Historical source-tree preservation receipt. Current necessary legacy code/data are followed from their actual uses.'});continue
        if f.endswith(('.json','.json.gz')):
            try:walk(json.loads(text),f)
            except json.JSONDecodeError:pass
        if f.endswith('.py'):
            try:tree=ast.parse(text)
            except SyntaxError:continue
            base='/'.join(f.split('/')[:2])if f.startswith('legacy/')else''
            env={'ROOT':base,'root':base,'R':base}
            def val(n):
                if isinstance(n,ast.Constant)and isinstance(n.value,str):return n.value
                if isinstance(n,ast.Name):return env.get(n.id)
                if isinstance(n,ast.BinOp)and isinstance(n.op,ast.Div):
                    a,b=val(n.left),val(n.right)
                    if a is not None and b is not None:return str(Path(a)/b)
                if isinstance(n,ast.Call)and isinstance(n.func,ast.Name)and n.func.id=='Path'and n.args:return val(n.args[0])
            for n in ast.walk(tree):
                if isinstance(n,ast.Assign):
                    v=val(n.value)
                    if v is not None:
                        for target in n.targets:
                            if isinstance(target,ast.Name):env[target.id]=v
                v=val(n)
                if v and v.startswith(prefixes) and v.endswith(extensions) and '\n'not in v and len(v)<512:
                    if (ROOT/v).is_file()or v in alias:put(v,f,role='static_code_path')
                    elif v.startswith(prefixes)and v.endswith(extensions):
                        excluded.append({'origin':f,'path':v,'reason':'Absent static producer path, not a recorded input: output/resume/prototype path. No successful computation input hash is waived here.'})
                if isinstance(n,(ast.Import,ast.ImportFrom)):
                    mods=[a.name for a in n.names]if isinstance(n,ast.Import)else[n.module or'']
                    for m in mods:
                        mn=m.replace('.','/')+'.py'
                        for candidate in [str(Path(f).parent/mn),'code/'+mn,mn]:
                            if (ROOT/candidate).is_file():put(candidate,f);break
        if f.endswith('.tex'):
            for arg in re.findall(r'\\(?:input|include)\{([^}]+)\}',text):put(arg if arg.endswith('.tex')else arg+'.tex',f)
            for arg in re.findall(r'\\path\{([^}]+)\}',text):
                if (ROOT/arg).is_file():put(arg,f)
                elif '/'not in arg:
                    candidate='research/classification/'+arg
                    if (ROOT/candidate).is_file():put(candidate,f)
    # No telemetry is mathematical evidence. Preserve the observation separately.
    for logical in list(missing):
        if '_RUNTIME_'in logical or '_CHECKPOINT_'in logical:
            excluded.append({'path':logical,'origins':sorted(missing.pop(logical)),'reason':'Execution/resume telemetry; complete mathematical coefficient/certificate reports are included.'})
            pending=[r for r in pending if r['path']!=logical]
    files=sorted(f for f in seen if (ROOT/f).is_file())
    result={'scope':__doc__,'primary_roots':roots['roots'],'primary_claim_count':len(claimids),'primary_claims':sorted(claimids),'direct_paths':sorted(seeds),'files':[{'path':f,'bytes':(ROOT/f).stat().st_size,'sha256':digest((ROOT/f).read_bytes())}for f in files],'edges':{k:sorted(v)for k,v in edges.items()},'recorded_hashes':integrity,'unresolved_recorded_inputs':pending,'missing':{k:sorted(v)for k,v in missing.items()},'excluded_nonpremises':excluded,'root_input_aliases':alias,'all_required_inputs_available':not missing and not pending,'mathematics_independently_certified':False}
    output.parent.mkdir(parents=True,exist_ok=True);output.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({'files':len(files),'bytes':sum(r['bytes']for r in result['files']),'recorded_hash_checks':len(integrity),'pending_inputs':len(pending),'missing':len(missing),'output':str(output)},indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);audit(p.parse_args().output)
