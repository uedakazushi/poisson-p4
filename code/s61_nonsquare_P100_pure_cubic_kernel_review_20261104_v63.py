#!/usr/bin/env python3
"""Stream full57 source rows and audit ten pure-H kernels, no new source query."""
import gzip,hashlib,json,resource,signal,subprocess,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_P100_PURE_CUBIC_KERNEL_20261104_v63_REVIEW'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def canonical(p):return {tuple(m):Q(c)for m,c in p if Q(c)}
def main(out):
    started=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));out.mkdir(exist_ok=False,parents=True)
    base=ROOT/'reports';pb=base/'S61_NONSQUARE_P100_PURE_CUBIC_KERNEL_20261104_v63';r=read(pb/'report.json');assert r['complete']
    paths={Path(__file__),Path(__file__).with_suffix('.cpp'),ROOT/'research/notes'/(ID+'.md'),pb/'report.json'}
    for rel,h in r['input_sha256'].items():assert sha(ROOT/rel)==h,rel;paths.add(ROOT/rel)
    for rel,h in r['evidence_sha256'].items():assert sha(pb/rel)==h,rel;paths.add(pb/rel)
    for stem in ['S61_NONSQUARE_P100_Q5_PURE98_ROW_20261104_v63_REVIEW','S61_NONSQUARE_P100_G18_PURE8_20261104_v63_REVIEW',
      'S61_NONSQUARE_P100_H1_20261104_v63_REVIEW','S61_NONSQUARE_ORACLE_FRAME_ADAPTER_20261103_v63_REVIEW',
      'S61_BINARY_NONSQUARE_CUBIC_ROW_BASIS_20261103_v62_REVIEW']:
        p=base/stem/'report.json';assert read(p)['all_passed'];paths.add(p)
    Pids=[1,6,8,9,11,14];ids=[0,2,3,4,5,7,10,12,13,15]
    raw=read(base/'S61_NONSQUARE_RAW_H_LINEAR_UNITS_20261104_v63/preflight.json')['raw16x7']
    sigmafile=base/'S61_NONSQUARE_P100_H1_RIGHT_INVERSE_20261104_v63/native.json';sigma=read(sigmafile)['linear_lift_coefficients7x16']
    eta=read(base/'S61_NONSQUARE_P100_H1_DUAL_20261104_v63/native.json')['eta_physical7']
    ci=read(base/'S61_NONSQUARE_P100_PURE98_SOURCE_PREFLIGHT_20261104_v63/cubic-coefficient-input.json');mons=ci['old_physical_monomials84']
    pivotdata=read(base/'S61_NONSQUARE_P100_Q5_PURE98_ROW_20261104_v63/cubic504.json');pivot=pivotdata['values6x84']
    assert pivotdata['normalized_Tn_FrameAdapter_ROWUNIT_DIVISION'] and [s['genuine_index']for s in ci['sources6']]==Pids
    basis=read(base/'S61_BINARY_NONSQUARE_CUBIC_ROW_BASIS_20261103_v62/report.json');sr=read(base/'S61_BINARY_NONSQUARE_CUBIC_REPLAY_20261103_v62/source-report.json')
    adapter=read(base/'S61_NONSQUARE_ORACLE_FRAME_ADAPTER_20261103_v63/adapter.json');binding=read(pb/'source-binding.json')
    fullmons=sr['normal_monomials'];positions={tuple(m):i for i,m in enumerate(fullmons)};purecols=[positions[tuple(m)]for m in mons]
    hcols=[]
    for j in [0,1,2,3,4,6,7]:
        m=[0]*8;m[5]=1;m[j]+=1;hcols.append(positions[tuple(m)])
    assert len(purecols)==84 and len(hcols)==7 and binding['pure_monomials84']==mons
    for i,p in enumerate(Pids):
        source=ci['sources6'][i];assert source['basis57_position']==41+p and source['actual135_row']==adapter['basis57_source_definitions'][41+p]['original_equation']
    g18=read(base/'S61_NONSQUARE_P100_G18_PURE7_20261104_v63/physical-cubic84.json')['values84']
    litmons=read(base/'S61_NONSQUARE_P100_PURE21_SOURCE_PREFLIGHT_20261104_v63/source-input.json')['physical_pure_cubic_monomials84']
    g18=[g18[litmons.index(m)]for m in mons]
    work=ROOT/'work/P100-pure-cubic-kernel-independent-review';work.mkdir(exist_ok=True,parents=True);inp=work/'streamed-input.txt'
    matrix=base/'S61_BINARY_NONSQUARE_CUBIC_ROW_BASIS_20261103_v62/remaining.txt.gz';hmatched=0;source_terms=0
    with inp.open('w')as f:
        def pw(p):f.write(str(len(p))+''.join(f' {m[0]} {Q(c)}'for m,c in p)+'\n')
        def cw(p):f.write(str(len(p))+''.join(f' {m[0]} {m[1]} {m[2]} {Q(c)}'for m,c in p)+'\n')
        def rw(v):pw(v['numerator']);pw(v['denominator'])
        for c in [Q(0),Q(20),Q(60),Q(2700,203),Q(12)]:
            pp=[[[1],'1']]+([[[0],str(c)]]if c else[]);pw(pp)
        for row in raw:
            for p in row:cw(p)
        for row in sigma:
            for v in row:rw(v)
        for v in eta:rw(v)
        for row in pivot:
            for v in row:rw(v)
        digest=hashlib.sha256()
        with gzip.open(matrix,'rb')as stream:
            header=stream.readline();digest.update(header);assert header.strip()==b'57 91 57'
            for pos in range(57):
                genuine=pos-41;store=genuine in ids and genuine!=0;polys={}
                for col in range(91):
                    line=stream.readline();assert line;digest.update(line)
                    if 0<=genuine<16 and(col in hcols or(store and col in purecols)):
                        tt=line.decode().split();n=int(tt[0]);assert len(tt)==1+4*n
                        p=[[list(map(int,tt[2+4*i:5+4*i])),str(Q(tt[1+4*i]))]for i in range(n)]
                        if col in hcols:assert canonical(p)==canonical(raw[genuine][hcols.index(col)]);hmatched+=1
                        if store and col in purecols:polys[col]=p
                if genuine in ids:
                    original=basis['remaining_original_row_indices'][pos];equation=adapter['basis57_source_definitions'][pos]['original_equation']
                    assert sr['labels'][original][0]!='product'and sr['labels'][original][1]==equation
                    saved=next(q for q in binding['nonpivot_sources10']if q['genuine_index']==genuine)
                    assert(saved['basis57_position'],saved['actual135'],saved['original63'])==(pos,equation,original)
                    f.write(str(genuine)+' '+str(int(genuine==0))+'\n')
                    if genuine==0:
                        for v in g18:rw(v)
                    else:
                        assert len(polys)==84
                        for col in purecols:cw(polys[col]);source_terms+=len(polys[col])
                    savedrow=read(pb/f'pure-cubic84-{genuine}.json');assert savedrow['nonpivot_genuine_index']==genuine
                    for key,n in [('old16_source_cofactors',16),('actual_delta7_source_cofactors',16),('pure84_values',84)]:
                        assert len(savedrow[key])==n
                        for v in savedrow[key]:rw(v)
                polys.clear()
            assert not stream.read()
        assert digest.hexdigest()==basis['remaining_matrix_sha256'] and hmatched==112 and source_terms==79559
    order=Pids+ids;inversions=sum(order[i]>order[j]for i in range(16)for j in range(i+1,16));assert inversions%2==0
    cpp=Path(__file__).with_suffix('.cpp');helper=ROOT/'code/s61_nonsquare_P100_pure98_review_20261104_v63.cpp';paths.add(helper);exe=work/'verify-stream'
    with(out/'compile.log').open('w')as log:subprocess.run(['bash','-lc','source /workspace/cloud-setup/activate.sh\nc++ -O2 -std=c++17 "$1" -o "$2" $(pkg-config --cflags --libs flint)','compile',str(cpp),str(exe)],check=True,stdout=log,stderr=subprocess.STDOUT,timeout=max(1,55-(time.monotonic()-started)))
    with inp.open()as stream:run=subprocess.run(['bash','-lc','source /workspace/cloud-setup/activate.sh\nexec "$@"','independent-stream',str(exe)],stdin=stream,capture_output=True,text=True,timeout=max(1,55-(time.monotonic()-started)))
    (out/'native.log').write_text(run.stdout);(out/'native.stderr').write_text(run.stderr);assert run.returncode==0,run.stderr
    native=read(out/'native.log');assert native['all_passed'] and native['Hkernel_identity_checks']==140 and native['pure_source_frame_identity_checks']==1680 and native['peak_polynomial_slots']<=150000
    proof={'all_passed':True,'full_saved57x91_source_H112_direct_matches':True,'new_nine_nonpivot_core_terms_streamed':source_terms,
      'cached504_pivot_values_and84_G18_reused_from_independent_source_receipts':True,'no139746_source_JSON_retained':True,
      'all16_Heta_49sigma_36pivot_and_sigma_row3_zero':True,'all140_Hkernel_100_identityblock_1680_pure_frame_identities':True,
      'formal_kernel_basis':'Since H_i*eta=0 and sigma*H=I-eta*e3^T, alpha_i=H_i*sigma gives H_i=sum_p alpha_ip H_p. The sixpivot alpha block isI6, hence their H images are independent. Every H-kernel vector subtracts its ten nonpivot coordinates times K_i, leaving a pivot-supported vector, which is zero by that independence. Thus the10 displayed vectors are an actual free basis of the formal source H-kernel on the whole P100 unit ring, not only at a point.',
      'full_source_basis_determinants':'Pivot-first permutation has even sign; adjoining the original six pivot vectors gives block triangular source change determinant1. The actualdelta7 frame change only replaces C6 bydelta*(C6-C1), so coefficient conversion has determinant1/delta; delta=-64000*C² is a whole-curve unit.',
      'normalization':'All complete source rows are the saved normalized Tn57 basis, already divided by the accepted source row unit. No new row-unit multiplication/division.',
      'scope':'Ten formal source H-kernel basis vectors and their pure84 values. Their actual cubic images may have relations modulo the quadratic ideal. No image rank10,70quartic products,NF117,quartic lower-bound update or full239 statement.',
      'native':native,'input_local':{'path':str(inp.relative_to(ROOT)),'sha256':sha(inp),'bytes':inp.stat().st_size}}
    (out/'independent-proof.json').write_text(json.dumps(proof,indent=2)+'\n');exe.unlink(missing_ok=True)
    report={'id':ID,'all_passed':True,'ten_full_source_Hkernel_basis_and_pure84_values_accepted':True,'no_image_rank_or_quartic_inference':True,
      'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in sorted(paths)},'evidence_sha256':{p.name:sha(p)for p in out.iterdir()},
      'seconds':time.monotonic()-started,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'frozen':True}
    assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024;signal.alarm(0)
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2),flush=True)
if __name__=='__main__':
    import argparse
    pa=argparse.ArgumentParser();pa.add_argument('--output',type=Path,default=ROOT/'reports'/ID);main(pa.parse_args().output.resolve())
