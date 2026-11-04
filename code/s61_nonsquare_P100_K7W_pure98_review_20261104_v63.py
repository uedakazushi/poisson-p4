#!/usr/bin/env python3
"""Independent one-row C68*W source/frame/NF98 replay from pinned cached84."""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import resource
import signal
import subprocess
import time

ROOT=Path(__file__).resolve().parents[1]
ID='S61_NONSQUARE_P100_K7W_PURE98_20261104_v63_REVIEW'


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main(out):
    start=time.monotonic();signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
    resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024))
    out=out.resolve();out.mkdir(parents=True,exist_ok=False);hashes={}
    def read(p):hashes[str(p.relative_to(ROOT))]=sha(p);return json.loads(p.read_text())
    def gate(stem,passed=True):
        p=ROOT/'reports'/stem/'report.json';r=read(p)
        if passed:assert r['all_passed']
        for name,h in r.get('input_sha256',{}).items():assert sha(ROOT/name)==h,name;hashes[name]=h
        for name,h in r.get('evidence_sha256',{}).items():
            q=p.parent/name;assert sha(q)==h,name;hashes[str(q.relative_to(ROOT))]=h
        return r
    b=ROOT/'reports';base=b/'S61_NONSQUARE_P100_K7W_PURE98_20261104_v63'
    r=gate(base.name,False);assert r['complete']
    assert sha(base/'report.json')=='42a14c31455101dced25ef66d38e315a846407522144979575676af6163525ee'
    gate('S61_NONSQUARE_P100_PURE_CUBIC_KERNEL_20261104_v63_REVIEW')
    gate('S61_NONSQUARE_P100_G18_PURE8_20261104_v63_REVIEW')
    gate('S61_NONSQUARE_JOINT_LELL_FREE164_20261104_v63_REVIEW')
    gate('S61_NONSQUARE_H_SQUARE_H1_PREFLIGHT_20261104_v63_REVIEW')
    kb=b/'S61_NONSQUARE_P100_PURE_CUBIC_KERNEL_20261104_v63'
    parent=read(kb/'pure-cubic84-7.json');pbind=read(kb/'source-binding.json')
    physical=read(base/'physical-cubic84.json');joint=read(base/'joint-cubic-frame.json')
    binding=read(base/'source-binding.json');row=read(base/'pure98-row3.json')
    frame=read(b/'S61_NONSQUARE_P100_G18_PURE7_20261104_v63/joint-cubic-frame.json')
    h=read(b/'S61_NONSQUARE_H_SQUARE_H1_PREFLIGHT_20261104_v63/certificate.json')
    fd=read(b/'S61_NONSQUARE_JOINT_LELL_FREE164_20261104_v63/certificate.json')
    manifest=read(b/'S61_NONSQUARE_P100_PURE98_SOURCE_PREFLIGHT_20261104_v63/manifest.json')
    one=dict(numerator=[[[0],'1']],denominator=[[[0],'1']]);zero=dict(numerator=[],denominator=[[[0],'1']])
    literal=[one if i==7 else zero for i in range(16)]
    assert parent['nonpivot_genuine_index']==7 and parent['old16_source_cofactors']==parent['actual_delta7_source_cofactors']==literal
    assert parent['all7_Hlinear_coefficients_zero'] and parent['both_full_source_frame_identities_verified']
    assert binding['old16_source_cofactors']==binding['actual_delta7_source_cofactors']==literal
    assert binding['original319_source_index']==226 and binding['source_label']==h['all112_cubic_source319_definitions'][7*7+3]==['cubic_normalized',50,68,3]
    assert binding['normalization_ROWUNIT_DIVISION_already_applied'] and binding['no_H1_H2_cleanup_or_truequartic_Qhat_rescale']
    assert physical['values84']==parent['pure84_values'] and binding['pure_monomials84']==pbind['pure_monomials84']
    assert frame['physical_light_inverse7x7'][3]==[one if j==6 else zero for j in range(7)]
    assert joint['physical_light_inverse7x7']==frame['physical_light_inverse7x7']
    assert row['physical_light']==3 and row['target122_columns']==h['target122_Hdegree_classes']['0']
    assert row['nonzero_count']==40 and row['individual_known_unit_count']==0
    # Rebuild the accepted exact pure NF dependency closures, not just trust a
    # manifest list of coefficients or substitute H=0 into the original ideal.
    nf=fd['full_target_projection239_to122'];ln=fd['linear_nodes'];sn=fd['scalar_circuit']['nodes'];atoms=fd['scalar_circuit']['atoms_num_den']
    roots=[nf['output122_linear_roots'][i] for i in row['target122_columns']]
    Lneed=set();Sroots=set();todo=roots[:]
    while todo:
        i=todo.pop()
        if i in Lneed:continue
        Lneed.add(i);op,*a=ln[i]
        if op=='l':
            for c,j in a[0]:Sroots.add(c);todo.append(j)
        else:assert op in ['v','z']
    Sneed=set();SA=set();todo=list(Sroots)
    while todo:
        i=todo.pop()
        if i in Sneed:continue
        Sneed.add(i);op,*a=sn[i]
        if op=='c':SA.add(a[0])
        elif op in ['a','s','m']:todo.extend(a)
        else:assert op=='z'
    assert [len(SA),len(Sneed),len(Lneed)]==[122,244,547]
    assert sorted(SA)==manifest['NF_needed_atom_ids'] and sorted(Sneed)==manifest['NF_needed_scalar_ids'] and sorted(Lneed)==manifest['NF_needed_linear_ids']
    assert roots==manifest['NF_pure98_linear_roots']
    (out/'preflight-and-source-binding.json').write_text(json.dumps(dict(
        only_source319_row226=True,actual_source='normalized C68, literal old16/actualdelta7 unitvector7',
        only_cached84_coefficients=True,physical_light3_exact_jointW=True,
        NF_closure=[122,244,547],one_product210_and_NF98=True,
        no_other_row_source_query_rank_or_factor_inventory=True,
        limits=dict(seconds=60,address_bytes=2_000_000_000,output_bytes=32*1024*1024,live_polynomial_slots=150000)),indent=2)+'\n')
    work=ROOT/'work/P100-K7W-independent-review';work.mkdir(parents=True,exist_ok=True);inp=work/'one-pinned-row-input.txt'
    with inp.open('w') as f:
        def pw(p):f.write(str(len(p))+''.join(f' {m[0]} {c}' for m,c in p)+'\n')
        def rw(v):pw(v['numerator']);pw(v['denominator'])
        def cw(p):f.write(str(len(p))+''.join(f' {m[0]} {m[1]} {m[2]} {c}' for m,c in p)+'\n')
        def nw(ids,array):
            f.write(str(len(ids))+'\n')
            for i in ids:
                n=array[i]
                if n[0]=='l':f.write(f'{i} l {len(n[1])}'+''.join(f' {c} {j}' for c,j in n[1])+'\n')
                else:f.write(str(i)+' '+' '.join(map(str,n))+'\n')
        for c in [F(0),F(20),F(60),F(2700,203),F(12)]:pw([[[1],'1']]+([[[0],str(c)]] if c else []))
        for v in parent['pure84_values']:rw(v)
        for m in binding['pure_monomials84']:f.write(' '.join(str(m[i]) for i in [0,1,2,3,4,6,7])+'\n')
        for v in physical['values84']:rw(v)
        for matrix in [frame['physical_light_inverse7x7'],joint['physical_light_inverse7x7']]:
            for rr in matrix:
                for v in rr:rw(v)
        for m in joint['joint_cubic_monomials84']:f.write(' '.join(map(str,m))+'\n')
        for v in joint['joint_cubic_values84']:rw(v)
        f.write(str(len(SA))+'\n')
        for i in sorted(SA):f.write(str(i)+'\n');cw(atoms[i][0]);cw(atoms[i][1])
        nw(sorted(Sneed),sn);nw(sorted(Lneed),ln)
        for values in [roots,row['target122_columns'],nf['input_variables']]:f.write(' '.join(map(str,values))+'\n')
        for m in h['parent_joint_target239_monomials']:f.write(' '.join(map(str,m))+'\n')
        for v in row['values98']:rw(v)
    assert inp.stat().st_size<2*1024*1024
    cpp=Path(__file__).with_suffix('.cpp');helper=ROOT/'code/s61_nonsquare_P100_pure98_review_20261104_v63.cpp';exe=work/'verify-one-source'
    for p in [Path(__file__),cpp,helper,ROOT/f'research/notes/{ID}.md']:hashes[str(p.relative_to(ROOT))]=sha(p)
    flags=subprocess.check_output(['pkg-config','--cflags','--libs','flint'],text=True).split()
    with (out/'compile.log').open('w') as log:
        subprocess.run(['c++','-O2','-std=c++17',str(cpp),'-o',str(exe),*flags],stdout=log,stderr=subprocess.STDOUT,check=True,timeout=max(1,55-(time.monotonic()-start)))
    with inp.open() as f:proc=subprocess.run([str(exe)],stdin=f,text=True,capture_output=True,timeout=max(1,55-(time.monotonic()-start)))
    (out/'native.log').write_text(proc.stdout);(out/'native.stderr').write_text(proc.stderr)
    if proc.returncode:
        (out/'failure.json').write_text(json.dumps(dict(all_passed=False,reason=proc.stderr,returncode=proc.returncode,no_automatic_larger_retry=True),indent=2)+'\n');raise SystemExit(proc.returncode)
    result=json.loads(proc.stdout);assert result['all_passed'] and result['peak_polynomial_slots']<=150000
    proof=dict(all_passed=True,accepted_cached_literal_C68_source84_only=True,
        old16_and_actualdelta7_source_vectors_exact_unitvector7=True,
        actual319_source_row226_with_no_new_rowunit_or_Qhat_scaling=True,
        physical_v3_equals_joint_W_coefficient1_no_tail=True,
        all84_physical_and_joint_coefficients_and49_frame_entries_checked=True,
        all210_product_coefficients_and98_NF_outputs_checked=True,
        no_H1_H2_subtraction_needed='The accepted literal cubic has all7 H-linear coefficients zero; multiplying a pure light keeps Hdegree0. NF preserves Hdegree. Its difference from the full source is only accepted quadratic rows.',
        all_knownunit_denominators=True,nonzero_count=40,individual_named_unit_count=0,
        no_rank_or_newbound_other_source_NF_query_factor_inventory_point_prime_or_GB=True,checker=result)
    (out/'independent-proof.json').write_text(json.dumps(proof,indent=2)+'\n')
    report=dict(id=ID,all_passed=True,actual_C68_times_W_row226_source_and_NF98_accepted=True,
        one_row_only_no_rank_or_bound_update=True,input_sha256=hashes,
        evidence_sha256={p.name:sha(p) for p in out.iterdir()},ignored_native_input_sha256=sha(inp),ignored_executable_sha256=sha(exe),
        seconds=time.monotonic()-start,max_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
    print(json.dumps(dict(all_passed=True,seconds=report['seconds'],report_sha256=sha(out/'report.json'))))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    main(p.parse_args().output)
