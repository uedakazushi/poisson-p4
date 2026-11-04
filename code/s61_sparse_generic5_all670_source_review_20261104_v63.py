#!/usr/bin/env python3
"""Stream every stored plane-selector polynomial against its original source."""
import argparse,gzip,json,resource,signal,subprocess,time
from fractions import Fraction as F
from pathlib import Path
from s61_sparse_generic5_plane_rank3_review_20261104_v63 import read,sha,dec
ROOT=Path(__file__).resolve().parents[1]
BASE='S61_SPARSE_GENERIC5_PLANE_FIVE_SELECTOR_20261104_v63';ID='S61_SPARSE_GENERIC5_ALL670_SOURCE_BINDING_20261104_v63_REVIEW'
def main(out,work):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False);work.mkdir(parents=True,exist_ok=False)
    folder=ROOT/'reports'/BASE;pr=read(folder/'report.json');mp=folder/'five-column-matrix.json.gz';matrix=read(mp);assert matrix['all_passed'];assert not pr['all_passed'] and pr['complete_matrix_saved']
    assert matrix['column_order']==['P*F(C)','J(C)*W0','J(C)*W1','polarQ(W0,d)','polarQ(W1,d)'];assert len(matrix['matrix134x5'])==134 and all(len(r)==5 for r in matrix['matrix134x5'])
    counts=[sum(len(p['terms']) for p in row) for row in matrix['matrix134x5']];maxentry=max(len(p['terms']) for row in matrix['matrix134x5'] for p in row);assert sum(counts)==matrix['total_terms']==692953 and max(counts)<=10415 and maxentry==1685<=2083
    hashes={}
    for name,h in pr['input_sha256'].items():assert sha(ROOT/name)==h;hashes[name]=h
    assert pr['evidence_sha256']['five-column-matrix.json.gz']==sha(mp)
    source=ROOT/'reports/S61_SPARSE_GENERIC5_INPUT_PREFLIGHT_20261104_v63';wp=ROOT/'reports/S61_SPARSE_GENERIC5_PLANE_RANK3_CIRCUIT_20261104_v63/plane-W-original.json.gz';sp=ROOT/'reports/S61_SPARSE_GENERIC5_SMALL_TAIL_FIELD_20261104_v63/shared-denominator-kernel-plane.json.gz';fp=source/'pulled134-input.json.gz';kp=source/'known-five-input.json.gz'
    data=read(fp);known=read(kp);original=read(wp);P=dec(read(sp)['common_denominator_P']);C=list(map(dec,original['center45']));d=list(map(dec,known['Q52_minus_center']));W=[[dec(row[j]) for j in range(2)] for row in original['W45x3']]
    assert P==dec(original['common_P'])
    gates=['S61_SPARSE_GENERIC5_PLANE_RANK3_CIRCUIT_20261104_v63_REVIEW/report.json','S61_SPARSE_GENERIC5_FULL134_NORMAL_CIRCUIT_20261104_v63_REVIEW/report.json','S61_SPARSE_GENERIC5_PLANE_FIVE_UNIVARIATE_20261104_v63_REVIEW/report.json']
    for name in gates:
        p=ROOT/'reports'/name;assert read(p)['all_passed'];hashes[str(p.relative_to(ROOT))]=sha(p)
    inp=work/'original-source-only.txt';groups=0
    with inp.open('w') as f:
        def put(p):
            f.write(str(len(p))+'\n')
            for e,c in sorted(p.items()):f.write(str(c)+' '+' '.join(map(str,e))+'\n')
        f.write('134 45\n');put(P)
        for p in C+d:put(p)
        for row in W:
            for p in row:put(p)
        for encoded in data['equations']:
            row={}
            for m,e,c in encoded:
                p=row.setdefault(tuple(m),{});e=tuple(e);p[e]=p.get(e,F(0))+F(c)
            row={m:{e:c for e,c in p.items() if c} for m,p in row.items()};groups+=len(row);f.write(str(len(row))+'\n')
            for m,p in sorted(row.items()):assert len(m)<=2;f.write(' '.join(map(str,(len(m),*m)))+'\n');put(p)
    assert inp.stat().st_size<32*1024*1024
    preflight=dict(actual_source_qmonomial_groups=groups,expected_total_terms=692953,max_expected_entry_terms=maxentry,max_expected_row_terms=max(counts),producer_intermediate_bound2083_is_not_final_entry_max=True,matrix_compressed_bytes=mp.stat().st_size,source_only_work_input_bytes=inp.stat().st_size,source_only_work_input_sha256=sha(inp),expected_matrix_not_duplicated_as_artifact=True,one_native_checker_no_field_arithmetic=True)
    (out/'preflight.json').write_text(json.dumps(preflight,indent=2)+'\n')
    cpp=ROOT/'code/s61_sparse_generic5_all670_source_review_20261104_v63.cpp';binary=work/'verify-all670'
    subprocess.run(['bash','-c','source /workspace/cloud-setup/activate.sh; c++ -O2 -std=c++17 "$1" -o "$2" $(pkg-config --cflags --libs flint) -lgmp','compile',str(cpp),str(binary)],check=True,timeout=15)
    native_start=time.monotonic()
    with (out/'native.log').open('wb') as log:
        proc=subprocess.Popen(['bash','-c','source /workspace/cloud-setup/activate.sh; exec "$1" "$2"','check',str(binary),str(inp)],stdin=subprocess.PIPE,stdout=log,stderr=subprocess.STDOUT)
        try:
            for row in matrix['matrix134x5']:
                lines=[]
                for p in row:
                    lines.append(f"{p['content']} {len(p['terms'])}\n")
                    lines.extend(' '.join(map(str,t))+'\n' for t in p['terms'])
                proc.stdin.write(''.join(lines).encode())
            proc.stdin.close();returncode=proc.wait(timeout=max(1,45-(time.monotonic()-native_start)));assert returncode==0
        finally:
            if proc.poll() is None:proc.kill();proc.wait()
    result=json.loads((out/'native.log').read_text());assert result['all_passed'] and result['literal_QQ_polynomial_equalities']==670 and result['all_expected_terms_compared']==692953
    proof=dict(all670_actual_five_variable_QQ_coefficients_source_bound=True,source_functionals=['P F(C)','J_F(C) W0','J_F(C) W1','polar_QF(W0,d)','polar_QF(W1,d)'],actual_original134_rows_preserved=True,all_q_squares_both_polarization_terms_retained=True,actual_25_selected_gate_extended_to_every_other645_entry=True,streamed_expected_rows_discarded_without_new_matrix_artifact=True,previous_692953_guard_failure_not_retried_or_promoted=True,native=result,no_field_Gaussian_point_prime_GB_or_univariate_only_inference=True,scope='Source binding of the full five-functional matrix; residual129 plane-zero still also needs its separate algebraic circuit identities')
    (out/'independent-proof.json').write_text(json.dumps(proof,indent=2)+'\n')
    for p in [Path(__file__),cpp,ROOT/'code/s61_sparse_generic5_plane_rank3_review_20261104_v63.py',ROOT/f'research/notes/{ID}.md',folder/'report.json',mp,fp,kp,wp,sp]:hashes[str(p.relative_to(ROOT))]=sha(p)
    report=dict(id=ID,all_passed=True,all670_five_variable_QQ_source_polynomial_equalities_independently_verified=True,no_new_full_matrix_field_elimination_or_source_point=True,input_sha256=hashes,evidence_sha256={p.name:sha(p) for p in out.iterdir()},seconds=time.monotonic()-start,max_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(dict(all_passed=True,seconds=report['seconds'],rss=report['max_rss_kib'],report_sha256=sha(out/'report.json'))))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--work',type=Path,required=True);a=p.parse_args();main(a.output,a.work)
