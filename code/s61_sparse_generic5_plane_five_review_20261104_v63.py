#!/usr/bin/env python3
"""Independent selected-five source rows and formal QQ[T] determinant audit."""
import argparse,gzip,hashlib,itertools,json,resource,signal,subprocess,time
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
BASE='S61_SPARSE_GENERIC5_PLANE_FIVE_UNIVARIATE_20261104_v63';ID=BASE+'_REVIEW'
WEIGHTS=[1,2,3,5,7]
def read(p):
    with gzip.open(p,'rt') if str(p).endswith('.gz') else open(p) as f:return json.load(f)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def collapse(terms,content=F(1)):
    out={}
    for coefficient,exponents in terms:
        e=sum(a*b for a,b in zip(WEIGHTS,exponents));out[e]=out.get(e,F(0))+content*F(coefficient)
    return {e:c for e,c in out.items() if c}
def decode(p):
    if isinstance(p,dict):return collapse([(t[0],t[1:]) for t in p['terms']],F(p['content']))
    if p and isinstance(p[0][0],list):return collapse([(c,e) for e,c in p])
    return collapse([(t[0],t[1:]) for t in p])
def rank(A):
    A=[row[:] for row in A];i=0
    for j in range(len(A[0])):
        k=next((k for k in range(i,len(A)) if A[k][j]),None)
        if k is None:continue
        A[k],A[i]=A[i],A[k];p=A[i][j];A[i]=[v/p for v in A[i]]
        for k in range(len(A)):
            if k!=i and A[k][j]:
                q=A[k][j];A[k]=[x-q*y for x,y in zip(A[k],A[i])]
        i+=1
        if i==len(A):break
    return i
def functionals():
    # A universal normal form for a centered noncollinear triangle and the two
    # fixed vertical points, not new numerical beta/source samples.
    points=[(F(1),F(0),F(0)),(F(0),F(1),F(0)),(F(-1),F(-1),F(0)),(F(0),F(0),F(-1,2)),(F(0),F(0),F(1))]
    evaluation=[[1,a,b,c,a*a,a*b,a*c,b*b,b*c,c*c] for a,b,c in points];assert rank(evaluation)==5
    ideal=[[1,0,0,1,-1,1,0,-1,0,-2],[0,1,0,0,-1,2,0,0,0,0],[0,0,1,0,0,2,0,-1,0,0],[0,0,0,0,0,0,1,0,0,0],[0,0,0,0,0,0,0,0,1,0]]
    assert all(sum(F(a)*b for a,b in zip(g,p))==0 for g in ideal for p in evaluation)
    selected=[0,1,2,6,8];assert [[g[j] for j in selected] for g in ideal]==[[int(i==j) for j in range(5)] for i in range(5)]
    return dict(evaluation_rank=5,dimension_quadrics_vanishing_on_five=5,independent_universal_ideal_basis=ideal,functional_coefficient_indices=selected,functionals_restrict_to_identity=True,Tr_c0_and_GL2_triangle_coordinate_change_required=True)
def main(out,work):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False);work.mkdir(parents=True,exist_ok=False)
    folder=ROOT/'reports'/BASE;pr=read(folder/'report.json');assert pr['all_passed'];hashes={}
    for p,h in pr['input_sha256'].items():assert sha(ROOT/p)==h;hashes[p]=h
    for p,h in pr['evidence_sha256'].items():assert sha(folder/p)==h;hashes[str((folder/p).relative_to(ROOT))]=h
    cert=read(folder/'native-certificate.json.gz');rows=cert['selected_original134_rows'];assert rows==[50,10,52,51,15] and len(set(rows))==5
    assert cert['rank']==5 and cert['pivot_column_order']==list(range(5))
    gates=['S61_SPARSE_GENERIC5_PLANE_RANK3_CIRCUIT_20261104_v63_REVIEW/report.json','S61_SPARSE_GENERIC5_FULL134_NORMAL_CIRCUIT_20261104_v63_REVIEW/report.json']
    for name in gates:
        p=ROOT/'reports'/name;g=read(p);assert g['all_passed'];hashes[str(p.relative_to(ROOT))]=sha(p)
    source=ROOT/'reports/S61_SPARSE_GENERIC5_INPUT_PREFLIGHT_20261104_v63';plane=ROOT/'reports/S61_SPARSE_GENERIC5_PLANE_RANK3_CIRCUIT_20261104_v63';shared=ROOT/'reports/S61_SPARSE_GENERIC5_SMALL_TAIL_FIELD_20261104_v63/shared-denominator-kernel-plane.json.gz'
    sf=source/'pulled134-input.json.gz';kf=source/'known-five-input.json.gz';wf=plane/'plane-W-original.json.gz';data=read(sf);known=read(kf);original=read(wf);den=read(shared)['common_denominator_P']
    C=list(map(decode,original['center45']));d=list(map(decode,known['Q52_minus_center']));W=[[decode(row[j]) for j in range(2)] for row in original['W45x3']];P=decode(den)
    assert P and decode(original['common_P'])==P
    equations=[]
    for r in rows:
        equation={}
        for qm,exp,c in data['equations'][r]:
            assert len(qm)<=2;p=equation.setdefault(tuple(qm),{});e=sum(a*b for a,b in zip(WEIGHTS,exp));p[e]=p.get(e,F(0))+F(c)
        equations.append({m:{e:c for e,c in p.items() if c} for m,p in equation.items()})
    # The original failed five-variable field run is not retried or promoted.
    # Its 692953-term matrix is not reconstructed or used as the arithmetic input.
    matrix=cert['selected_original_matrix5x5'];assert len(matrix)==5 and all(len(r)==5 for r in matrix)
    inp=out/'independent-selected-source-input.txt'
    with inp.open('w') as f:
        def put(p):
            f.write(str(len(p))+'\n')
            for e,c in sorted(p.items()):f.write(f'{e} {c}\n')
        put(P)
        for v in C+d:put(v)
        for r in W:
            for p in r:put(p)
        for eq in equations:
            f.write(str(len(eq))+'\n')
            for m,p in sorted(eq.items()):f.write(' '.join(map(str,(len(m),*m)))+'\n');put(p)
        for row in matrix:
            for p in row:put({int(e):F(c) for e,c in p})
        put({int(e):F(c) for e,c in cert['determinant_polynomial']})
    assert inp.stat().st_size<32*1024*1024
    pf=functionals();tracepath=ROOT/'reports/TORIC_E2_SCALAR_INVERSE_20261101.json';tr=read(tracepath);assert tr['all_passed'] and str(tr['formal_trace_X_over_derivative'])=='0' and tr['universal_toric_barycenter']=='Tr(c)=0; Tr(q_toric)=3*q_center=2*Q51+Q52'
    fivepath=ROOT/'reports/S61_KNOWN_FIVE_POINT_SECTION_20261101.json';fr=read(fivepath);assert fr['all_passed'] and fr['affine_span_dimension']==3 and fr['intersection_length']==5 and fr['intersection_reduced']
    cpp=ROOT/'code/s61_sparse_generic5_plane_five_review_20261104_v63.cpp';native=work/'verify-selected-source'
    subprocess.run(['bash','-c','source /workspace/cloud-setup/activate.sh; c++ -std=c++17 -O2 "$1" -o "$2" $(pkg-config --cflags --libs flint) -lgmp','compile',str(cpp),str(native)],check=True,timeout=15)
    result=subprocess.run(['bash','-c','source /workspace/cloud-setup/activate.sh; exec "$1" "$2"','check',str(native),str(inp)],capture_output=True,check=True,timeout=40)
    (out/'native.log').write_bytes(result.stdout+result.stderr);arithmetic=json.loads(result.stdout);assert arithmetic['all_passed'] and arithmetic['determinant_degree']==cert['determinant_degree']==448
    proof=dict(formal_parameter_exponents=WEIGHTS,selected_original134_rows=rows,only_selected_five_source_rows_reconstructed=True,all25_actual_entries_bound=True,nonzero_determinant_from_independent_integer_polynomial_matrix=True,arithmetic=arithmetic,universal_functional_injectivity=pf,scalar_P_multiplies_all_five_true_functionals=True,Tr_c0_trace_gate=True,invariant_genericity_via_Sp4_and_Gm_dominance=True,actual_plane_restriction_quadrics_rank=5,not_numeric_sourcepoint_or_prime=True,no_full134x5_reconstruction_field_Gaussian_cross126_or_exhaustion=True)
    (out/'independent-proof.json').write_text(json.dumps(proof,indent=2)+'\n')
    for p in [Path(__file__),cpp,ROOT/f'research/notes/{ID}.md',folder/'report.json',sf,kf,wf,shared,tracepath,fivepath]:hashes[str(p.relative_to(ROOT))]=sha(p)
    report=dict(id=ID,all_passed=True,selected5_original_source_rows_and_QQ_T_determinant_independently_verified=True,plane_quadrics_rank5_not_full_fiber_exhaustion=True,no_author_solver_newpoint_prime_GB_or_full_matrix_retry=True,input_sha256=hashes,evidence_sha256={p.name:sha(p) for p in out.iterdir()},seconds=time.monotonic()-start,max_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(dict(all_passed=True,seconds=report['seconds'],rss=report['max_rss_kib'],report_sha256=sha(out/'report.json'))))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--work',type=Path,required=True);a=p.parse_args();main(a.output,a.work)
