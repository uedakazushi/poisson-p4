#!/usr/bin/env python3
"""Package the current necessary proof/input/certificate closure in Library volumes.

Content-addressed chunks deduplicate exact bytes. JSON.gz and native .in.gz are
individual exact data inputs, not historical ZIP snapshots or nested bundles.
"""
from pathlib import Path
import argparse,hashlib,json,subprocess,zipfile,zlib
ROOT=Path(__file__).resolve().parents[1]
CHUNK=4*1024*1024;LIMIT=9*1024*1024

def sha(b):return hashlib.sha256(b).hexdigest()
def main(output):
    output.mkdir(parents=True,exist_ok=False)
    audit=json.loads((ROOT/'reports/required_proof_dependency_closure.json').read_text());assert audit['all_required_inputs_available']
    paths={r['path']for r in audit['files']}
    paths.update(['README.md','AGENTS.md','RESEARCH_POLICY.md','requirements.txt','research/CLAIMS.json','research/COMPONENTS.json','research/REFERENCES.md','research/STATE.md','research/conceptual/FINAL_TARGET.md','research/conceptual/GLOBAL_PLAN.md','research/classification/CURRENT.md','research/conceptual/component_entrances.csv','code/data/component_manuscript_catalogue.json','code/build_component_catalogue.py','code/map_component_entrances.py','code/check_component_catalogue.py','code/audit_required_dependencies.py','scripts/build_structural_manuscripts.sh','scripts/make_structural_bundle.py','scripts/materialize_frozen_proof_inputs.py','scripts/unpack_current_bundle.py','research/reproduction/classification_roots.json','research/reproduction/frozen_producer_index.json','research/reproduction/README.md','reports/required_proof_dependency_closure.json','reports/dependency_reconstruction_checks.json','reports/final_manuscript_validation.json','reports/component_catalogue_integration_validation.json','reports/rank_two_integration_validation.json','reports/normalization_residue_validation.json','reports/conductor_conormal_validation.json','reports/quadric_residue_validation.json','code/check_global_structure.py','code/check_marked_reconstruction.py','code/check_quadric_residue_reconstruction.py','code/check_conductor_conormal_reconstruction.py','code/check_normalization_residue_reconstruction.py','code/check_rank_two_reconstruction.py','reports/marked_reconstruction_checks.json','reports/quadric_residue_checks.json','reports/conductor_conormal_checks.json','reports/normalization_residue_checks.json'])
    paths.update(str(p.relative_to(ROOT))for p in (ROOT/'papers').glob('*.tex'))
    paths.update(['papers/poisson_p4_classification_ja.pdf','papers/poisson_divisor_symmetries_ja.pdf','papers/normal_actual_divisor_recognition_ja.pdf','papers/s61_first_jet_reduction_ja.pdf'])
    files=[];chunks={};symlinks={};logical_paths={}
    expected={r['path']:r['sha256']for r in audit['files']}
    physical=set()
    for name in sorted(paths):
        p=ROOT/name
        if any(v in Path(name).parts for v in ['work','__pycache__','.git'])or name.endswith(('.zip','.tar','.tar.gz','.log','.aux','.pyc')):raise ValueError(('forbidden',name))
        actual=p.resolve();assert actual.is_relative_to(ROOT)
        canonical=str(actual.relative_to(ROOT));logical_paths[name]=canonical;physical.add(canonical)
        if name in expected:assert sha(p.read_bytes())==expected[name],('Input changed after audit',name)
        for i in range(1,len(Path(name).parts)+1):
            prefix=ROOT.joinpath(*Path(name).parts[:i])
            if prefix.is_symlink():
                key=str((prefix.parent.resolve()/prefix.name).relative_to(ROOT))
                symlinks[key]=__import__('os').readlink(prefix)
    for name in sorted(physical):
        p=ROOT/name;data=p.read_bytes();digest=sha(data);ids=[]
        for offset in range(0,len(data),CHUNK):
            b=data[offset:offset+CHUNK];h=sha(b);ids.append(h)
            if h not in chunks:chunks[h]={'source':name,'offset':offset,'bytes':len(b),'compressed_estimate':len(zlib.compress(b,6))+180}
        files.append({'path':name,'bytes':len(data),'sha256':digest,'chunks':ids})
    parts=[[]];size=0
    for h,r in sorted(chunks.items()):
        if size+r['compressed_estimate']>LIMIT:parts.append([]);size=0
        parts[-1].append(h);size+=r['compressed_estimate'];r['part']=len(parts)
    commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    manifest={'version':'v72-current-dependency-closure','source_commit':commit,'repository':'https://github.com/uedakazushi/poisson-p4','scope':'Current main/structural manuscripts and necessary primary proof, exact input, code and certificate closure.','mathematical_verification':'Integrity checks do not independently certify the classification theorem. See the scoped proof/read/reexecution reports.','part_count':len(parts),'files':files,'symlinks':symlinks,'logical_paths':logical_paths,'chunks':{h:{k:r[k]for k in ['bytes','part']}for h,r in chunks.items()},'reproduction':'Read research/reproduction/README.md after reconstruction. Run UNPACK.py with all ZIP volumes and a new destination.'}
    manifest['bundle_id']=sha(json.dumps(manifest,sort_keys=True,separators=(',',':')).encode())
    encoded=(json.dumps(manifest,ensure_ascii=False,indent=2)+'\n').encode();receipts=[]
    readme='All volumes are disjoint parts of ONE current research bundle. They contain no nested ZIP.\nExtract UNPACK.py from any volume, then run:\npython UNPACK.py --parts poisson-p4-v72-part*.zip --destination poisson-p4-current\nAll chunks and reconstructed files are checked by SHA256.\n'
    for number,ids in enumerate(parts,1):
        path=output/f'poisson-p4-v72-part{number:02d}-of{len(parts):02d}.zip'
        with zipfile.ZipFile(path,'x',compression=zipfile.ZIP_DEFLATED,compresslevel=6)as z:
            z.writestr('BUNDLE_MANIFEST.json',encoded);z.writestr('README.txt',readme);z.writestr('UNPACK.py',(ROOT/'scripts/unpack_current_bundle.py').read_bytes())
            for h in ids:
                r=chunks[h]
                with (ROOT/r['source']).open('rb')as stream:stream.seek(r['offset']);b=stream.read(r['bytes'])
                assert sha(b)==h;z.writestr('payload/'+h,b)
        assert path.stat().st_size<10*1024*1024
        with zipfile.ZipFile(path)as z:assert z.testzip()is None
        receipts.append({'path':str(path),'bytes':path.stat().st_size,'sha256':sha(path.read_bytes())});print(json.dumps(receipts[-1]),flush=True)
    result={'source_commit':commit,'bundle_id':manifest['bundle_id'],'files':len(files),'symlinks':len(symlinks),'logical_paths':len(logical_paths),'unique_chunks':len(chunks),'uncompressed_file_bytes':sum(r['bytes']for r in files),'parts':receipts,'zip_crc':'PASS','nested_archives':False,'library_save':'Separate persistence receipt required.'}
    (output/'VALIDATION.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
