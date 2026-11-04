#!/usr/bin/env python3
"""Independent full-G18 source, seven NF rows and fixed8-minor audit; no search."""
import gzip,hashlib,json,resource,signal,subprocess,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_P100_G18_PURE8_20261104_v63_REVIEW'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def core(p):return {tuple(m):Q(v)for m,v in p if Q(v)}
def native_matrix(p):
    lines=p.read_text().splitlines();n,m=map(int,lines[0].split()[:2]);assert len(lines)==1+n*m;ans=[]
    for line in lines[1:]:
        t=line.split();k=int(t[0]);assert len(t)==1+4*k
        ans.append([[list(map(int,t[2+4*i:5+4*i])),str(Q(t[1+4*i]))]for i in range(k)])
    return n,m,ans
def main(out):
    started=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));out.mkdir(exist_ok=False,parents=True)
    base=ROOT/'reports';paths={Path(__file__),Path(__file__).with_suffix('.cpp'),ROOT/'research/notes'/(ID+'.md')}
    def audit(p,key=None):
        r=read(p);paths.add(p)
        if key:assert r[key]
        for rel,h in r.get('input_sha256',{}).items():assert sha(ROOT/rel)==h,rel;paths.add(ROOT/rel)
        for rel,h in r.get('evidence_sha256',{}).items():assert sha(p.parent/rel)==h,rel;paths.add(p.parent/rel)
        return r
    gb=base/'S61_NONSQUARE_P100_G18_PURE7_20261104_v63';audit(gb/'report.json','complete')
    rb=base/'S61_NONSQUARE_P100_PURE8_GENERIC_RANK_20261104_v63';audit(rb/'report.json','complete')
    pf=base/'S61_NONSQUARE_P100_PURE21_SOURCE_PREFLIGHT_20261104_v63';pfr=audit(pf/'report.json','full_source_input_bound')
    for stem in ['S61_NONSQUARE_P100_Q5_PURE98_ROW_20261104_v63_REVIEW','S61_NONSQUARE_ORACLE_FRAME_ADAPTER_20261103_v63_REVIEW',
      'S61_NONSQUARE_JOINT_LELL_FREE164_20261104_v63_REVIEW','S61_NONSQUARE_P100_H1_20261104_v63_REVIEW',
      'S61_NONSQUARE_H_SQUARE_H1_PREFLIGHT_20261104_v63_REVIEW','S61_BINARY_NONSQUARE_CUBIC_ROW_BASIS_20261103_v62_REVIEW']:
        audit(base/stem/'report.json','all_passed')
    pre8=base/'S61_NONSQUARE_P100_PURE8_RANK_SOURCE_PREFLIGHT_20261104_v63';pre8r=audit(pre8/'report.json');pre8c=read(pre8/'certificate.json')
    assert sha(pre8/'certificate.json')==pre8r['certificate_sha256']
    src=read(pf/'source-input.json');assert sha(pf/'source-input.json')==pfr['source_input_sha256'];source=src['three_full_sources'][0]
    assert source['basis57_position']==41 and source['genuine_index']==0 and source['original135_equation']==18 and source['original63_normalized_row']==42
    # Bind the actual FULL remaining57x91 row, including every H coefficient.
    rowbase=base/'S61_BINARY_NONSQUARE_CUBIC_ROW_BASIS_20261103_v62';basis=read(rowbase/'report.json');matrix=rowbase/'remaining.txt.gz'
    raw=gzip.decompress(matrix.read_bytes());assert hashlib.sha256(raw).hexdigest()==basis['remaining_matrix_sha256'];lines=raw.decode().splitlines();assert lines[0]=='57 91 57'
    full=[]
    for line in lines[1+41*91:1+42*91]:
        tt=line.split();n=int(tt[0]);assert len(tt)==1+4*n
        full.append([[list(map(int,tt[2+4*i:5+4*i])),str(Q(tt[1+4*i]))]for i in range(n)])
    sr=read(base/'S61_BINARY_NONSQUARE_CUBIC_REPLAY_20261103_v62/source-report.json');mons=sr['normal_monomials']
    pure=[i for i,m in enumerate(mons)if m[5]==0];heavy=[i for i,m in enumerate(mons)if m[5]]
    assert len(full)==91 and len(pure)==84 and len(heavy)==7 and all(not full[i]for i in heavy)
    for i,p in enumerate(source['full_pure84_coefficients']):assert core(p)==core(full[pure[i]])
    assert src['physical_pure_cubic_monomials84']==[mons[i]for i in pure] and sum(map(len,source['full_pure84_coefficients']))==1858
    assert basis['remaining_original_row_indices'][41]==42 and sr['labels'][42][0]!='product' and sr['labels'][42][1]==18
    adapter=read(base/'S61_NONSQUARE_ORACLE_FRAME_ADAPTER_20261103_v63/adapter.json')
    assert adapter['basis57_source_definitions'][41]=={'kind':'genuine_cubic','original_equation':18}
    assert source['oracle_to_saved_row_unit']==adapter['basis57_encoded_source_units'][41]
    h=read(base/'S61_NONSQUARE_H_SQUARE_H1_PREFLIGHT_20261104_v63/certificate.json');fd=read(base/'S61_NONSQUARE_JOINT_LELL_FREE164_20261104_v63/certificate.json')
    physical=[0,1,2,3,4,6,7]
    for i,j in enumerate(physical):
        assert h['all112_cubic_source319_definitions'][i]==['cubic_normalized',42,18,j]
        assert src['twentyone_actual_source319_products'][i]['source319_index']==174+i
    # Derive the complete pure98 operator closure from the accepted full NF117.
    NF=fd['full_target_projection239_to122'];ln=fd['linear_nodes'];sn=fd['scalar_circuit']['nodes'];atoms=fd['scalar_circuit']['atoms_num_den']
    purecols=h['target122_Hdegree_classes']['0'];roots=[NF['output122_linear_roots'][i]for i in purecols]
    Lneed=set();Sroots=set();todo=roots[:]
    while todo:
        i=todo.pop()
        if i in Lneed:continue
        Lneed.add(i);n=ln[i]
        if n[0]=='l':
            for c,j in n[1]:Sroots.add(c);todo.append(j)
        else:assert n[0]in ['z','v']
    Sneed=set();SA=set();todo=list(Sroots)
    while todo:
        i=todo.pop()
        if i in Sneed:continue
        Sneed.add(i);n=sn[i]
        if n[0]=='c':SA.add(n[1])
        elif n[0]in ['a','s','m']:todo.extend(n[1:])
        else:assert n[0]=='z'
    manifest=read(base/'S61_NONSQUARE_P100_PURE98_SOURCE_PREFLIGHT_20261104_v63/manifest.json')
    assert sorted(Lneed)==manifest['NF_needed_linear_ids'] and sorted(Sneed)==manifest['NF_needed_scalar_ids'] and sorted(SA)==manifest['NF_needed_atom_ids']
    assert [len(Lneed),len(Sneed),len(SA)]==[547,244,122]
    # Audit H-degree for ALL122 output roots, not merely the selected pure98.
    index={x:i for i,x in enumerate(NF['input_variables'])};allneed=set();todo=NF['output122_linear_roots'][:]
    while todo:
        i=todo.pop()
        if i in allneed:continue
        allneed.add(i)
        if ln[i][0]=='l':todo.extend(j for c,j in ln[i][1])
    degrees={}
    for i in sorted(allneed):
        n=ln[i]
        if n[0]=='v':degrees[i]={h['parent_joint_target239_monomials'][index[n[1]]][7]}
        elif n[0]=='z':degrees[i]=set()
        else:
            assert n[0]=='l';degrees[i]=set().union(*(degrees[j]for c,j in n[1]))if n[1]else set()
    for deg,cols in h['target122_Hdegree_classes'].items():
        for i in cols:assert not degrees[NF['output122_linear_roots'][i]]or degrees[NF['output122_linear_roots'][i]]=={int(deg)}
    native=read(rb/'native.json');labels=read(rb/'source-labels.json');selrows=native['selected_rows'];selcols=native['selected_columns']
    assert selrows==[0,2,3,5,7,6,1,4] and selcols==[0,1,2,3,4,5,6,9]
    assert [purecols[i]for i in selcols]==[0,1,2,4,5,6,7,10] and labels['target122_labels98']==purecols
    q5p=base/'S61_NONSQUARE_P100_Q5_PURE98_ROW_20261104_v63/pure98.json';q5=read(q5p)
    rows=[q5]+[read(gb/f'pure98-row{i}.json')for i in range(7)]
    for i,row in enumerate(rows):
        assert row['target122_columns']==purecols
        if i:assert row['actual_source319_index']==173+i and row['physical_light']==physical[i-1]
    assert sum(len(v['numerator'])+len(v['denominator'])for row in rows for v in row['values98'])==pre8c['RF_input_terms']==18193
    jc=read(gb/'joint-cubic-frame.json');pc=read(gb/'physical-cubic84.json');assert pc['actual135']==18
    inv=base/'S61_NONSQUARE_JOINT_LELL_PREFLIGHT_20261104_v63/coordinate-inverse-numerators7x7.txt';fwd=inv.with_name('coordinate-forward7x7.txt')
    n,m,inverse=native_matrix(inv);assert(n,m)==(7,7);n,m,forward=native_matrix(fwd);assert(n,m)==(7,7);paths.update([inv,fwd,matrix,rowbase/'report.json',base/'S61_BINARY_NONSQUARE_CUBIC_REPLAY_20261103_v62/source-report.json'])
    rawframe=read(base/'S61_NONSQUARE_RAW_H_LINEAR_UNITS_20261104_v63/preflight.json')
    work=ROOT/'work/P100-G18-pure8-independent-review';work.mkdir(exist_ok=True,parents=True);ip=work/'fixed-input.txt'
    with ip.open('w')as stream:
        def pw(p):
            p={int(m[0]):Q(v)for m,v in p}if isinstance(p,list)else p
            stream.write(str(len(p))+''.join(f' {i} {v}'for i,v in sorted(p.items()))+'\n')
        def cw(p):stream.write(str(len(p))+''.join(f' {m[0]} {m[1]} {m[2]} {v}'for m,v in p)+'\n')
        def rw(v):pw(v['numerator']);pw(v['denominator'])
        def aw(ids):
            stream.write(str(len(ids))+'\n')
            for i in ids:stream.write(str(i)+'\n');cw(atoms[i][0]);cw(atoms[i][1])
        def nw(ids,array):
            stream.write(str(len(ids))+'\n')
            for i in ids:
                n=array[i]
                if n[0]=='l':stream.write(f'{i} l {len(n[1])}'+''.join(f' {c} {j}'for c,j in n[1])+'\n')
                else:stream.write(str(i)+' '+' '.join(map(str,n))+'\n')
        for c in [Q(0),Q(20),Q(60),Q(2700,203),Q(12)]:pw({0:c,1:Q(1)}if c else{1:Q(1)})
        stream.write(' '.join(map(str,selrows))+'\n'+' '.join(map(str,selcols))+'\n')
        for powers in labels['row_LCD_powers']:stream.write(' '.join(map(str,powers))+'\n')
        for p in source['full_pure84_coefficients']:cw(p)
        for mon in src['physical_pure_cubic_monomials84']:stream.write(' '.join(str(mon[i])for i in physical)+'\n')
        for v in pc['values84']:rw(v)
        for p in inverse:cw(p)
        for p in forward:cw(p)
        for row in rawframe['constant_inverse7x7']:stream.write(' '.join(map(str,row))+'\n')
        for row in jc['physical_light_inverse7x7']:
            for v in row:rw(v)
        for mon in jc['joint_cubic_monomials84']:stream.write(' '.join(map(str,mon))+'\n')
        for v in jc['joint_cubic_values84']:rw(v)
        aw(sorted(SA));nw(sorted(Sneed),sn);nw(sorted(Lneed),ln)
        for a in [roots,purecols,NF['input_variables']]:stream.write(' '.join(map(str,a))+'\n')
        for mon in h['parent_joint_target239_monomials']:stream.write(' '.join(map(str,mon))+'\n')
        for v in rows[0]['values98']:rw(v)
        for row in rows[1:]:
            for v in row['values98']:rw(v)
        for row in native['selected_square']:
            for p in row:pw(p)
        pw(native['determinant_nonzero'])
        for row in native['final_adjugate']:
            for p in row:pw(p)
    assert ip.stat().st_size<16*1024*1024
    cpp=Path(__file__).with_suffix('.cpp');peer_cpp=ROOT/'code/s61_nonsquare_P100_pure98_review_20261104_v63.cpp';paths.add(peer_cpp);exe=work/'verify-native'
    with(out/'compile.log').open('w')as log:subprocess.run(['bash','-lc','source /workspace/cloud-setup/activate.sh\nc++ -O2 -std=c++17 "$1" -o "$2" $(pkg-config --cflags --libs flint)','compile',str(cpp),str(exe)],check=True,stdout=log,stderr=subprocess.STDOUT,timeout=max(1,55-(time.monotonic()-started)))
    with ip.open()as f:run=subprocess.run(['bash','-lc','source /workspace/cloud-setup/activate.sh\nexec "$@"','fixed-peer',str(exe)],stdin=f,capture_output=True,text=True,timeout=max(1,55-(time.monotonic()-started)))
    (out/'native.log').write_text(run.stdout);(out/'native.stderr').write_text(run.stderr)
    assert run.returncode==0,run.stderr;n=read(out/'native.log');assert n['all_passed'] and n['peak_polynomial_slots']<=150000 and n['det_degree']==200
    proof={'all_passed':True,'literal_FULL57x91_G18_row41_all91coefficients':True,'seven_Hlinear_coefficients_identically_zero':True,
      'actual135_row':18,'old_normalized63_row':42,'Tn_normalization':'Saved row is already divided by189a/F; no second row-unit operation. Full source adapter/basis receipts independently accepted.',
      'whole_curve_source84_joint84_inverse49_and_all7_actual_products_NF98':True,'all122_NF_outputs_Hdegree_source_audit':True,
      'knownunit784_rowLCD_and_fixed64_selected_entries':True,'independent_fixed8_determinant_and_both128_adj_products':True,
      'selected_source_rows':selrows,'selected_target98_columns':selcols,'selected_original122_targets':[purecols[i]for i in selcols],
      'exact_generic_rank_original8':8,'residual7x97_rank':7,'residual_rank_reason':'The independently accepted q5 row is unimodular. The accepted explicit GL2 target completion and source clear give a block1 direct_sum residual7x97. Thus rank(original8)=1+rank(residual). No22MB residual matrix replay is needed for this formal identity.',
      'original_quartic_lower_on_chart_D_delta':149,'whole_P100_lower':142,
      'scope':'The polynomial delta is nonzero but its unknown zeros are not removed. LCDs are whole-curve units, so the8 actual pure sources are independent only on D(delta). Together with the accepted117quad+23H1+1H2 this gives149 there. No uniform149/full239/allFULLU claim.',
      'new_source_oracle_point_prime_ranksearch_GB_jobs':0,'native':n,'input_local':{'path':str(ip.relative_to(ROOT)),'sha256':sha(ip),'bytes':ip.stat().st_size}}
    (out/'independent-proof.json').write_text(json.dumps(proof,indent=2)+'\n');exe.unlink(missing_ok=True)
    report={'id':ID,'all_passed':True,'source_G18_full57_seven_rows_and_fixed8_certificate':True,'generic_rank8_and_open_quartic149_accepted':True,'whole142_unchanged_not_uniform149':True,
      'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in sorted(paths)},'evidence_sha256':{p.name:sha(p)for p in out.iterdir()},
      'seconds':time.monotonic()-started,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'frozen':True}
    assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024;signal.alarm(0)
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2),flush=True)
if __name__=='__main__':
    import argparse
    pa=argparse.ArgumentParser();pa.add_argument('--output',type=Path,default=ROOT/'reports'/ID);main(pa.parse_args().output.resolve())
