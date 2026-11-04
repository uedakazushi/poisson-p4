#!/usr/bin/env python3
"""Source-bound corrected K2 cached84 API, seven products and exact NF values."""
import hashlib,json,resource,signal,subprocess,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_P100_K2_PURE7_20261104_v63_REVIEW'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main(out):
    started=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));out.mkdir(exist_ok=False,parents=True)
    base=ROOT/'reports';pb=base/'S61_NONSQUARE_P100_K2_PURE7_20261104_v63';r=read(pb/'report.json');assert r['complete']
    paths={Path(__file__),Path(__file__).with_suffix('.cpp'),ROOT/'research/notes'/(ID+'.md'),pb/'report.json'}
    for rel,h in r['input_sha256'].items():assert sha(ROOT/rel)==h,rel;paths.add(ROOT/rel)
    for rel,h in r['evidence_sha256'].items():assert sha(pb/rel)==h,rel;paths.add(pb/rel)
    kp=base/'S61_NONSQUARE_P100_PURE_CUBIC_KERNEL_20261104_v63/pure-cubic84-2.json';parent=read(kp)
    kernelpeer=base/'S61_NONSQUARE_P100_PURE_CUBIC_KERNEL_20261104_v63_REVIEW/report.json';assert read(kernelpeer)['all_passed'];paths.add(kernelpeer)
    gp=base/'S61_NONSQUARE_P100_G18_PURE8_20261104_v63_REVIEW/report.json';gr=read(gp);assert gr['all_passed'];paths.add(gp)
    binding=read(pb/'source-binding.json');assert binding['kernel_label']==2 and binding['actual_nonpivot135']==23 and parent['nonpivot_genuine_index']==2
    assert binding['old16_source_cofactors']==parent['old16_source_cofactors'] and binding['actual_delta7_source_cofactors']==parent['actual_delta7_source_cofactors']
    assert binding['physical_lights']==[0,1,2,3,4,6,7]
    kb=base/'S61_NONSQUARE_P100_PURE_CUBIC_KERNEL_20261104_v63/source-binding.json';kbinding=read(kb);assert binding['pure_monomials84']==kbinding['pure_monomials84']
    mons=kbinding['pure_monomials84'];physical=[0,1,2,3,4,6,7]
    frame=read(base/'S61_NONSQUARE_P100_G18_PURE7_20261104_v63/joint-cubic-frame.json');jc=read(pb/'joint-cubic-frame.json');pc=read(pb/'physical-cubic84.json')
    assert pc['values84']==parent['pure84_values']
    fdpath=base/'S61_NONSQUARE_JOINT_LELL_FREE164_20261104_v63/certificate.json';hp=base/'S61_NONSQUARE_H_SQUARE_H1_PREFLIGHT_20261104_v63/certificate.json';mp=base/'S61_NONSQUARE_P100_PURE98_SOURCE_PREFLIGHT_20261104_v63/manifest.json'
    for p in [fdpath,hp,mp]:assert sha(p)==gr['input_sha256'][str(p.relative_to(ROOT))]
    fd,h,m=map(read,[fdpath,hp,mp]);nf=fd['full_target_projection239_to122'];sn=fd['scalar_circuit']['nodes'];atoms=fd['scalar_circuit']['atoms_num_den'];ln=fd['linear_nodes']
    for i in range(16):
        for j,p in enumerate(physical):assert h['all112_cubic_source319_definitions'][7*i+j][0]=='cubic_normalized' and h['all112_cubic_source319_definitions'][7*i+j][3]==p
    rows=[read(pb/f'pure98-row{i}.json')for i in range(7)];purecols=h['target122_Hdegree_classes']['0']
    assert all(row['target122_columns']==purecols for row in rows)
    work=ROOT/'work/P100-K2-pure7-independent-review';work.mkdir(exist_ok=True,parents=True);inp=work/'fixed-input.txt'
    with inp.open('w')as f:
        def pw(p):f.write(str(len(p))+''.join(f' {mm[0]} {Q(c)}'for mm,c in p)+'\n')
        def rw(v):pw(v['numerator']);pw(v['denominator'])
        def cw(p):f.write(str(len(p))+''.join(f' {mm[0]} {mm[1]} {mm[2]} {Q(c)}'for mm,c in p)+'\n')
        def nw(ids,array):
            f.write(str(len(ids))+'\n')
            for i in ids:
                n=array[i]
                f.write(f'{i} l {len(n[1])}'+''.join(f' {c} {j}'for c,j in n[1])+'\n'if n[0]=='l'else str(i)+' '+' '.join(map(str,n))+'\n')
        for c in [Q(0),Q(20),Q(60),Q(2700,203),Q(12)]:pw([[[1],'1']]+([[[0],str(c)]]if c else[]))
        for v in parent['pure84_values']:rw(v)
        for mon in mons:f.write(' '.join(str(mon[i])for i in physical)+'\n')
        for v in pc['values84']:rw(v)
        for row in frame['physical_light_inverse7x7']:
            for v in row:rw(v)
        for row in jc['physical_light_inverse7x7']:
            for v in row:rw(v)
        for mon in jc['joint_cubic_monomials84']:f.write(' '.join(map(str,mon))+'\n')
        for v in jc['joint_cubic_values84']:rw(v)
        ids=m['NF_needed_atom_ids'];f.write(str(len(ids))+'\n')
        for i in ids:f.write(str(i)+'\n');cw(atoms[i][0]);cw(atoms[i][1])
        nw(m['NF_needed_scalar_ids'],sn);nw(m['NF_needed_linear_ids'],ln)
        for values in [m['NF_pure98_linear_roots'],purecols,nf['input_variables']]:f.write(' '.join(map(str,values))+'\n')
        for mon in h['parent_joint_target239_monomials']:f.write(' '.join(map(str,mon))+'\n')
        for row in rows:
            for v in row['values98']:rw(v)
    cpp=Path(__file__).with_suffix('.cpp');arithmetic=ROOT/'code/s61_nonsquare_P100_pure98_review_20261104_v63.cpp';paths.add(arithmetic);exe=work/'verify-cached-parent'
    with(out/'compile.log').open('w')as log:subprocess.run(['bash','-lc','source /workspace/cloud-setup/activate.sh\nc++ -O2 -std=c++17 "$1" -o "$2" $(pkg-config --cflags --libs flint)','compile',str(cpp),str(exe)],check=True,stdout=log,stderr=subprocess.STDOUT,timeout=max(1,55-(time.monotonic()-started)))
    with inp.open()as f:run=subprocess.run(['bash','-lc','source /workspace/cloud-setup/activate.sh\nexec "$@"','cached-parent-peer',str(exe)],stdin=f,capture_output=True,text=True,timeout=max(1,55-(time.monotonic()-started)))
    (out/'native.log').write_text(run.stdout);(out/'native.stderr').write_text(run.stderr);assert run.returncode==0,run.stderr
    n=read(out/'native.log');assert n['all_passed']and n['peak_polynomial_slots']<=150000 and n['nonzero_counts']==[55,89,89,36,45,89,89]
    proof={'all_passed':True,'parent84':'Exact cached K2 values are taken from the independently accepted ten-kernel packet. K2 is corrected C23 minus its pivot combination, not the literal C23.',
      'source_cofactors_both_original16_and_actualdelta7_parent_match':True,
      'actual_quartic_source_recipe':'For physical lightpositionj, multiply the entire parent source combination. Original16 cofactor s goes to original319 row174+7*s+j; no Hlift/pure tail is dropped. Actualdelta7 cofactors belong to the separately preconditioned H1 cache only.',
      'all84_physical_parent_joint_and49_cached_inverse_source_checks':True,'all1470_products_all686_NF_outputs':True,
      'Hdegree_and_NF_source_operator':'The identical FD/H/manifest hashes pin the already independently reconstructed all122 H-degree closure and exact244scalar/547linear NF operator. Accepted K2 has H0 and noH2; only pure98 outputs remain.',
      'scope':'Only seven source/value rows. No rank,selectedminor,gcd,70-product block,quartic lower-bound change or full239 inference. Whole142 and open149 bounds unchanged.',
      'new_source_oracle_rank_point_prime_GB_jobs':0,'native':n,'local_input':{'path':str(inp.relative_to(ROOT)),'sha256':sha(inp),'bytes':inp.stat().st_size}}
    (out/'independent-proof.json').write_text(json.dumps(proof,indent=2)+'\n');exe.unlink(missing_ok=True)
    report={'id':ID,'all_passed':True,'corrected_K2_parent84_and_seven_actual_NF_source_rows':True,'no_rank_or_bounds_update':True,
      'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in sorted(paths)},'evidence_sha256':{p.name:sha(p)for p in out.iterdir()},
      'seconds':time.monotonic()-started,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'frozen':True}
    assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024;signal.alarm(0)
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2),flush=True)
if __name__=='__main__':
    import argparse
    pa=argparse.ArgumentParser();pa.add_argument('--output',type=Path,default=ROOT/'reports'/ID);main(pa.parse_args().output.resolve())
