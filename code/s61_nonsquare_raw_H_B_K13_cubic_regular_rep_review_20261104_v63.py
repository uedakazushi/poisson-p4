#!/usr/bin/env python3
"""Combined exact408 source receipt and fixed finite coefficient-minor proof."""
import hashlib,json,resource,signal,subprocess,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_RAW_H_B_K13_CUBIC_REGULAR_REP_20261104_v63_REVIEW';P=32003
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main(out):
    started=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));out.mkdir(exist_ok=False,parents=True)
    pb=ROOT/'reports/S61_NONSQUARE_RAW_H_B_K13_CUBIC_REGULAR_REP_20261104_v63';pr=read(pb/'report.json');assert pr['all_passed'] and pr['rank_modp']==221
    paths={Path(__file__),Path(__file__).with_suffix('.cpp'),ROOT/'research/notes'/(ID+'.md'),pb/'report.json'}
    for rel,h in pr['input_sha256'].items():assert sha(ROOT/rel)==h,rel;paths.add(ROOT/rel)
    for rel,h in pr['evidence_sha256'].items():assert sha(pb/rel)==h,rel;paths.add(pb/rel)
    sourcepeer=ROOT/'reports/S61_NONSQUARE_RAW_H_B_K13_CUBIC17_COEFFICIENTS_20261104_v63_REVIEW/report.json';assert read(sourcepeer)['all_passed'];paths.add(sourcepeer)
    cb=ROOT/'reports/S61_NONSQUARE_RAW_H_B_K13_CUBIC17_COEFFICIENTS_20261104_v63';values=[json.loads(s)for s in(cb/'matrix17x24.jsonl').read_text().splitlines()]
    dg=ROOT/'reports/S61_NONSQUARE_RAW_H_B_K13_D_EUCLID_20261104_v63/checkpoint.json';g={int(i):int(v)for i,v in read(dg)['g13']};assert len(values)==408
    gate=read(pb/'denominator-punit-gate.json');assert gate['prime']==P and gate['passed']
    lc=g[13]%P;assert lc==gate['g13_leading_coefficient_residue']==16763
    denoms=[int(v['value']['denominator'])%P for v in values];assert all(denoms)and denoms==gate['all408_QQ_common_denominator_residues']
    ilc=pow(lc,P-2,P);assert lc*ilc%P==1;monic=[g.get(i,0)%P*ilc%P for i in range(13)];coeffs=[]
    for k,r in enumerate(values):
        assert(r['index'],r['row17'],r['column24'])==(k,k//24,k%24)
        d=denoms[k];di=pow(d,P-2,P);assert d*di%P==1
        num={int(i):int(v)for i,v in r['value']['numerator']};assert all(0<=i<13 for i in num)
        coeffs.extend(num.get(i,0)%P*di%P for i in range(13))
    original_input=list(map(int,(pb/'regular-input.txt').read_text().split()));assert original_input==[P,17,24,13]+monic+coeffs
    n=read(pb/'native-certificate.json');assert n['prime']==P and n['regular_shape']==[221,312]and n['rank_modp']==221
    cols=n['selected_columns'];perm=n['row_permutation_PA_equals_LU'];diag=n['diagonal_U'];lu=n['combined_LU221x221']
    assert len(cols)==len(set(cols))==221 and all(0<=i<312 for i in cols)
    assert sorted(perm)==list(range(221)) and len(diag)==len(lu)==221 and all(len(row)==221 for row in lu)
    parity=sum(perm[i]>perm[j]for i in range(221)for j in range(i+1,221))%2;assert parity==n['row_swaps']%2
    work=ROOT/'work/K13-cubic-regular-independent-review';work.mkdir(exist_ok=True,parents=True);inp=work/'fixed-input.txt'
    with inp.open('w')as f:
        f.write('32003 17 24 13\n'+' '.join(map(str,monic))+'\n'+' '.join(map(str,coeffs))+'\n')
        for a in [cols,perm,diag]:f.write(' '.join(map(str,a))+'\n')
        for row in lu:f.write(' '.join(map(str,row))+'\n')
        f.write(str(n['actual_selected_minor_determinant'])+'\n')
    cpp=Path(__file__).with_suffix('.cpp');exe=work/'verify-fixed-minor'
    with(out/'compile.log').open('w')as log:subprocess.run(['bash','-lc','source /workspace/cloud-setup/activate.sh\nc++ -O2 -std=c++17 "$1" -o "$2" $(pkg-config --cflags --libs flint)','compile',str(cpp),str(exe)],check=True,stdout=log,stderr=subprocess.STDOUT,timeout=max(1,55-(time.monotonic()-started)))
    with inp.open()as f:run=subprocess.run(['bash','-lc','source /workspace/cloud-setup/activate.sh\nexec "$@"','fixed-peer',str(exe)],stdin=f,capture_output=True,text=True,timeout=max(1,55-(time.monotonic()-started)))
    (out/'native.log').write_text(run.stdout);(out/'native.stderr').write_text(run.stderr);assert run.returncode==0,run.stderr
    native=read(out/'native.log');assert native['all_passed']and native['all_PLU_entries']==48841 and native['independent_FLINT_fixed_determinant']==30542 and native['explicit_scalar_peak_bound']<=150000
    # Bind the genuine full-source criterion, not only the numerical shape.
    gates=['S61_NONSQUARE_CUBIC_ZW25_UNITS_20261104_v63_REVIEW',
      'S61_BINARY_NONSQUARE_CUBIC_ROW_BASIS_20261103_v62_REVIEW',
      'S61_NONSQUARE_CUBIC_PROJECTED15_UNIT_PAIR_20261104_v63_REVIEW',
      'S61_NONSQUARE_RAW_H_B_K13_ACTUAL_SOURCE_SCHEME_20261104_v63_REVIEW']
    for stem in gates:
        p=ROOT/'reports'/stem/'report.json';assert read(p)['all_passed'];paths.add(p)
    for p in [ROOT/'reports/S61_NONSQUARE_THETA_FULL_CUBIC_SOURCE_PREFLIGHT_20261104_v63/report.json',
      ROOT/'reports/S61_NONSQUARE_CUBIC_ZW25_AUGMENTED_INPUT_20261104_v63/input.json',
      ROOT/'reports/S61_NONSQUARE_PRODUCT41_FULL_OPEN_20261103_v63/report.json',dg,cb/'matrix17x24.jsonl']:
        paths.add(p)
    proof={'all_passed':True,'combined_exact_source408_receipt':str(sourcepeer.relative_to(ROOT)),
      'all408_denominators_and_g13_LC_punits_independently_checked':True,'fixed_prime':32003,
      'regular_basis':'1,a,...,a^12. Own monomial reduction table forms each selected column from C_ij*a^k, rowblock i*13+l and fullcolumn j*13+k. No a/root/beta point evaluation.',
      'all48841_PLU_identities_and_pivot_product':True,'independent_FLINT_fixed221_determinant':30542,
      'characteristic_zero_lift':'V=Z_(32003), A_V=V[a]/monic(g13) is free13 because the g13 leading coefficient is a V-unit. All408 C_ij lie in A_V by the full denominator check. The selected regular221minor reduces to30542!=0, so the same rational minor is nonzero. Thus C:Ka13^24→Ka13^17 is surjective as a QQ-linear map and therefore as a Ka13-linear map. Right-exact base change yields rank17 at EVERY residue field; no field/irreducibility/reducedness assumption or socle shortcut is needed.',
      'actual_full_source_transfer':'The independent exact408 source audit binds the universal pre-curve flat matrix byte-for-byte to the actual25x32 source checkpoint, including its genuine cubic rows. C=cLCM*aB*M with a,B units on the whole actual Ka13 bad algebra. The old8prefix determinant−a⁴FB²C²/4608 is a QQ/FULLU unit. Together with the accepted product41 and Cdelta primitive-source quotient criterion (actual full57 source, not an arbitrary leading section), rank17 of M gives actual full cubic57 at every residue field of this whole finite bad algebra. GF unitness of the oldprefix/original cubic is not assumed: the inference uses characteristic0 source units.',
      'whole_actual_Bbad_full_cubic57':True,'scope':'Only the actual Ka13 bad locus. No wholeFULLU cubic57, quartic239, P100149 or genericPoisson exhaustion claim.',
      'new_rank_search_QQ_Gaussian_prime_point_or_GB_jobs':0,'native':native,
      'local_input':{'path':str(inp.relative_to(ROOT)),'sha256':sha(inp),'bytes':inp.stat().st_size}}
    (out/'independent-proof.json').write_text(json.dumps(proof,indent=2)+'\n');exe.unlink(missing_ok=True)
    report={'id':ID,'all_passed':True,'exact408_source_and_fixed_PLU_minor_combined':True,'whole_actual_Ka13_Bbad_full_cubic57':True,'no_allFULLU_or_quartic_inference':True,
      'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in sorted(paths)},'evidence_sha256':{p.name:sha(p)for p in out.iterdir()},
      'seconds':time.monotonic()-started,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'frozen':True}
    assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024;signal.alarm(0)
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2),flush=True)
if __name__=='__main__':
    import argparse
    pa=argparse.ArgumentParser();pa.add_argument('--output',type=Path,default=ROOT/'reports'/ID);main(pa.parse_args().output.resolve())
