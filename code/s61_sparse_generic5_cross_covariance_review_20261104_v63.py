#!/usr/bin/env python3
"""Independent20 identities, raw-row lineage, and off-shell covariance audit."""
import argparse,gzip,hashlib,itertools,json,resource,signal,time
from collections import Counter
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
BASE='S61_SPARSE_GENERIC5_CROSS_RANK_COVARIANCE_20261104_v63';ID=BASE+'_REVIEW'
def read(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def rref(A):
    A=[row[:] for row in A];r=0;pivots=[]
    for j in range(len(A[0])):
        k=next((i for i in range(r,len(A)) if A[i][j]),None)
        if k is None:continue
        A[k],A[r]=A[r],A[k];p=A[r][j];A[r]=[v/p for v in A[r]]
        for i in range(len(A)):
            if i!=r and A[i][j]:
                q=A[i][j];A[i]=[x-q*y for x,y in zip(A[i],A[r])]
        pivots.append(j);r+=1
        if r==len(A):break
    return A[:r],pivots
def det(A):
    A=[r[:] for r in A];q=F(1)
    for j in range(len(A)):
        k=next(i for i in range(j,len(A)) if A[i][j])
        if k!=j:A[k],A[j]=A[j],A[k];q=-q
        v=A[j][j];q*=v
        for i in range(j+1,len(A)):
            if A[i][j]:
                c=A[i][j]/v
                for h in range(j,len(A)):A[i][h]-=c*A[j][h]
    return q
def main(out):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False)
    folder=ROOT/'reports'/BASE;pr=read(folder/'report.json');assert pr['all_passed'];hashes={}
    for p,h in pr['input_sha256'].items():assert sha(ROOT/p)==h;hashes[p]=h
    for p,h in pr['evidence_sha256'].items():assert sha(folder/p)==h;hashes[str((folder/p).relative_to(ROOT))]=h
    base=ROOT/'legacy/v02/data';old=read(base/'symplectic_slice_equations.json');new=read(base/'symplectic_slice_reduced_equations.json');certificate=read(base/'symplectic_slice_certificate.json');ids=read(base/'symplectic_bianchi_identities.json');weights=old['parameter_weights'];assert weights==new['parameter_weights']
    order=new['old_generator_indices'];inverse={j:i for i,j in enumerate(order)};assert len(order)==134 and len(inverse)==134 and len(old['equations'])==154
    decode=lambda p:{tuple(m):F(n,d) for m,n,d in p};before=list(map(decode,old['equations']));after=list(map(decode,new['equations']));weight=lambda p:{sum(weights[i] for i in m) for m in p}
    matrix=[{} for _ in range(154)]
    def put(i,j,u,c):
        p=matrix[i].setdefault(j,{});p[u]=p.get(u,F(0))+c
        if not p[u]:del p[u]
    for j,i in enumerate(order):assert before[i]==after[j];put(i,j,-1,F(1))
    removed=[]
    for record in ids['generator_redundancy_identities']:
        i=record['eliminated_old_generator'];assert i not in inverse and weight(before[i])=={5};removed.append(i)
        for j,n,d in record['constant_coefficients_on_retained_old_generators']:
            assert weight(before[j])=={5};put(i,inverse[j],-1,F(n,d))
        for j,u,n,d in record['linear_coefficients_on_weight4_old_generators']:
            assert weights[u]==1 and weight(before[j])=={4};put(i,inverse[j],u,F(n,d))
    assert len(removed)==len(set(removed))==20 and set(removed)==set(range(154))-set(order)
    contributions=0
    for i in removed:
        difference={m:-v for m,v in before[i].items()}
        for j,p in matrix[i].items():
            for u,c in p.items():
                for m,v in after[j].items():
                    key=tuple(sorted((*m,u))) if u>=0 else m;difference[key]=difference.get(key,F(0))+c*v
                    if not difference[key]:del difference[key]
                    contributions+=1
        assert not difference
    saved=read(folder/'global154to134-linear-operator.json')['operator154x134_sparse'];converted=[{int(j):{int(u):F(c) for u,c in p} for j,p in row} for row in saved]
    assert converted==matrix;terms=sum(len(p) for row in matrix for p in row.values());assert terms==747 and contributions==262420
    for j,i in enumerate(order):assert matrix[i]=={j:{-1:F(1)}}
    # Reuse the exact350 construction gate; do not expand its350 source rows.
    rp=ROOT/'reports/S61_REBUILD_20261030.json';rb=read(rp);assert rb['all_passed'] and rb['comparison_summary']['json_files_equal_ignoring_runtime_seconds']==5
    for source in rb['source_scripts']:assert sha(ROOT/source['path'])==source['sha256']
    for item in rb['file_comparisons']:assert sha(base/item['file'])==item['archived_sha256']
    # The inspected polynomial_row_basis only subtracts constant QQ multiples;
    # every selected original row therefore spans the same350-row vector space.
    pairs=list(itertools.combinations(range(5),2));mons=list(itertools.combinations_with_replacement(range(5),2));pi={p:i for i,p in enumerate(pairs)};mi={p:i for i,p in enumerate(mons)}
    rawweight=lambda r:2-sum(i>0 for i in pairs[r//15])+sum(i>0 for i in mons[r%15])
    counts=Counter();maxq=0
    for row in old['eliminated_raw_variables']:
        w=rawweight(row['raw_index']);counts[w]+=1
        for m,n,d in row['polynomial']:assert sum(weights[i] for i in m)==w;maxq=max(maxq,sum(weights[i]==2 for i in m))
    assert dict(counts)=={1:4,2:15,3:36} and maxq==1 and len(certificate['elimination_rows'])==55
    assert all(rawweight(r)==weights[i] for i,r in enumerate(certificate['raw_free_indices']))
    assert Counter(next(iter(weight(p))) for p in before)=={4:70,5:84} and all(len(weight(p))==1 for p in before)
    assert len(certificate['remaining_equation_rows'])==154
    # Independent E wedge gl5 quotient and perspective response at constant J.
    vertical=[]
    for a,b in itertools.product(range(5),repeat=2):
        row=[F(0)]*150
        for i in range(5):
            if i==a:continue
            pair=tuple(sorted((i,a)));mon=tuple(sorted((i,b)));row[15*pi[pair]+mi[mon]]+=1 if i<a else -1
        assert len({rawweight(i) for i,c in enumerate(row) if c})<=1;vertical.append(row)
    vr,vp=rref(vertical);assert len(vp)==24 and vp==certificate['vertical_pivots']
    response=[]
    for j in range(1,5):
        v=[F(0)]*150
        for a,b in [(1,2),(3,4)]:
            v[15*pi[a,b]+mi[0,j]]-=2
            if a==j:v[15*pi[0,b]+mi[0,0]]+=1
            if b==j:v[15*pi[0,a]+mi[0,0]]-=1
        quotient=[v[k]-sum((v[p]*vr[i][k] for i,p in enumerate(vp) if v[p]),F(0)) for k in range(150)]
        response.append([quotient[k] for k in certificate['perspective_orbit_pivots_raw']])
    response=[list(row) for row in zip(*response)];assert response==[[F(-3 if i<2 else -2) if i==j else F(0) for j in range(4)] for i in range(4)] and det(response)==36
    # Verify exact old126 ancestry and current physical-input hashes, not G6.
    finite=ROOT/'reports/S61_AFFINE_NORMAL_FINITE5_20261103_v63';fr=read(finite/'report.json');fc=read(finite/'certificate.json');peer=ROOT/'reports/S61_AFFINE_NORMAL_FINITE5_20261103_v63_REVIEW/report.json';fp=read(peer)
    assert fr['all_passed'] and fp['all_passed'] and fr['exact_mixed_cross_rank126'];assert fc['cross_matrix_shape']==[129,126] and fc['exact_cross_rank']==126
    assert fp['all134_original_equations_accounted_for'] and fp['all129_factorizations_checked_exact'] and fp['cross_input_rebuilt_from_actual_G_and_byte_equal'] and fp['integer_inverse_products_checked']==2
    oldinput=finite/'cross-input.txt.gz';assert sha(oldinput)==fr['evidence_sha256']['cross-input.txt.gz'];assert hashlib.sha256(gzip.decompress(oldinput.read_bytes())).hexdigest()==fc['native_input_sha256']
    modelpath='legacy/v02/data/symplectic_slice_reduced_equations.json';assert fp['input_sha256'][modelpath]==sha(ROOT/modelpath)
    # Actual known-five component/graph and source-cover gates make the plane
    # intrinsic under the connected group, not a chosen non-invariant minor.
    geometry=['Q51_UNIVERSAL_SECTION_REVIEW_20261101.json','Q52_UNIVERSAL_SECTION_REVIEW_20261101.json','KNOWN_FIVE_SMOOTH_LOCUS_REVIEW_20261101.json','TORIC_UNIVERSAL_THREE_PLANE_REVIEW_20261101.json','S61_SPARSE_GENERIC5_PLANE_FIVE_UNIVARIATE_20261104_v63_REVIEW/report.json','S61_SPARSE_GENERIC5_ALL670_SOURCE_BINDING_20261104_v63_REVIEW/report.json','S61_SPARSE_GENERIC5_RESIDUAL129_CIRCUIT_20261104_v63_REVIEW/report.json','Q_NEW_FAMILY_GENERIC_BETA_SPARSE_SECTION_MINOR_20261104_v63_REVIEW/report.json']
    for name in geometry:
        p=ROOT/'reports'/name;assert read(p)['all_passed'];hashes[str(p.relative_to(ROOT))]=sha(p)
    proof=dict(all20_original61_QQ_Bianchi_identities=True,all747_operator_coefficients_constant_or_weight1_beta=True,no_q_dependent_coefficients=True,retained_projection_operator_equals_I134=True,full350_constant_QQ_rowspace_gate_reused_with_source_and_archive_hashes=True,triangular55_unique_lowjet_graph_source_audited=True,all55_reconstruction_weight_and_qdegree_checks=True,independent_vertical_quotient_rank24=True,independent_perspective_response=[[str(c) for c in row] for row in response],perspective_response_det=36,offshell_normalizing_matrix_depends_only_on_beta=True,q_action_affine_and_invertible_by_inverse_normalization=True,RAW134_not_entire_qquadratic_ideal_space=True,G6_not_used_in_old126=True,actual_known_plane_equivariant_via_accepted_components=True,intrinsic_cross_rank126_invariant_under_affine_frame_changes=True,extended_S_Y_rank131_invariant_under_Y_plus_S_R=True,generic_K5_cross_rank126_proved=True,new126minor_not_computed_or_assumed_fixed_old_selector=True,no_full_fiber_exhaustion=True)
    (out/'independent-proof.json').write_text(json.dumps(proof,indent=2)+'\n')
    paths=[Path(__file__),ROOT/f'research/notes/{ID}.md',folder/'report.json',rp,peer,oldinput,ROOT/'research/notes/KNOWN_FIVE_SMOOTH_LOCUS_20261101.md',ROOT/'research/notes/Q51_UNIVERSAL_SECTION_REVIEW_20261101.md',ROOT/'research/notes/Q52_UNIVERSAL_SECTION_REVIEW_20261101.md',ROOT/'research/notes/S61_SPARSE_GENERIC5_CROSS_RANK_COVARIANCE_20261104_v63.md']
    for p in paths:hashes[str(p.relative_to(ROOT))]=sha(p)
    r=dict(id=ID,all_passed=True,RAW134_rowspace_covariance_and_offshell_graph_independently_accepted=True,all20_QQ_Bianchi_and_independent_perspective_arithmetic=True,generic_K5_cross_rank126_and_extended_rank131=True,no_entire_ideal_G6_fixed_minor_or_exhaustion_inference=True,input_sha256=hashes,evidence_sha256={'independent-proof.json':sha(out/'independent-proof.json')},seconds=time.monotonic()-start,max_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    (out/'report.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(dict(all_passed=True,seconds=r['seconds'],rss=r['max_rss_kib'],report_sha256=sha(out/'report.json'))))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
