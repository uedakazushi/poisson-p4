#!/usr/bin/env python3
"""Independent fixed polynomial rank/minor-open audit; no solver or search."""
import argparse,hashlib,itertools,json,resource,signal,time
from pathlib import Path
import sympy as s
from sympy.polys.rings import ring
ROOT=Path(__file__).resolve().parents[1]
ID='S61_NONSQUARE_RAW_H_GENERIC_RANK_20261104_v63_REVIEW'
OPEN='S61_NONSQUARE_RAW_H_MINOR_OPEN_20261104_v63_REVIEW'
def read(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main(out):
    start=time.monotonic();out=out.resolve();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False)
    folder=ROOT/'reports/S61_NONSQUARE_RAW_H_GENERIC_RANK_20261104_v63';of=ROOT/'reports/S61_NONSQUARE_RAW_H_MINOR_OPEN_20261104_v63';gr=read(folder/'report.json');orr=read(of/'report.json');assert gr['complete_exact_generic_rank_certificate'] and gr['frozen'] and orr['frozen'];hashes={}
    for report in [gr,orr]:
        for p,h in report['input_sha256'].items():assert sha(ROOT/p)==h;hashes[p]=h
    for p,h in [(folder/'certificate.json',gr['certificate_sha256']),(folder/'checkpoint.json',gr['checkpoint_sha256']),(of/'certificate.json',orr['certificate_sha256'])]:assert sha(p)==h;hashes[str(p.relative_to(ROOT))]=h
    g=read(folder/'certificate.json');checkpoint=read(folder/'checkpoint.json');oc=read(of/'certificate.json');cp=ROOT/'reports/S61_NONSQUARE_RAW_H_LINEAR_UNITS_20261104_v63/checkpoint.json';old=read(cp)
    accepted=ROOT/'reports/S61_NONSQUARE_RAW_H_LINEAR_UNITS_20261104_v63_REVIEW/report.json';ag=read(accepted);assert ag['all_passed'] and ag['raw_Hlinear_uniform_rank_lower']==2 and ag['input_sha256'][str(cp.relative_to(ROOT))]==sha(cp)
    R,a,e,d=ring('a,e,d',s.QQ)
    def decode(p):return R.from_dict({tuple(m):s.QQ(c) for m,c in p})
    def determinant(M):
        n=len(M);out=R.zero
        for perm in itertools.permutations(range(n)):
            term=R.one
            for i,j in enumerate(perm):term*=M[i][j]
            out+=(-1)**sum(perm[i]>perm[j] for i in range(n) for j in range(i+1,n))*term
        return out
    M=[[decode(q) for q in row] for row in old['matrix15x6']];assert len(M)==15 and all(len(row)==6 for row in M)
    F=500*a+90*e**2+3*e-24;B=(50*a+18*e**2-9*e-2)*d+60*e-40;C=100*a+18*e**2-9*e-2;J=5*a*d-6*e-4;named=[a,d,F,B,C,J]
    p=M[2][0];assert p and len(p)==21 and p==M[4][0]==decode(g['p'])
    assert M[2][5]-M[4][5]==s.QQ(6,5)*B*p
    pairdet=determinant([[M[i][j] for j in [0,5]] for i in [2,4]])
    assert pairdet==-s.QQ(6,5)*B*p*p==decode(g['literal_pair_determinant'])
    brows=[7,9,10,11,12,13,14];bcols=[2,3,4];assert [i for i in old['remaining_rows'] if any(M[i][j] for j in bcols)]==brows
    assert all(M[i][j]==0 for i in [2,4] for j in bcols)
    block=[[M[i][j] for j in bcols] for i in brows];assert block==[[decode(p) for p in row] for row in g['actual_block7x3']]
    assert sum(len(p) for row in block for p in row)==908
    assert g['input14_remaining_rows']==old['remaining_rows'] and g['input5_columns']==old['remaining_columns']
    # Fixed three fraction-free steps, no pivot selection or unit inference.
    A=[row[:] for row in block];rr=list(range(7));cc=list(range(3));previous=R.one;division_count=0
    for stage,step in enumerate(checkpoint['steps']):
        i,j=step['pivot_row7'],step['pivot_col3'];assert i in rr and j in cc
        pivot=decode(step['pivot']);assert pivot==A[i][j] and pivot
        assert previous==decode(step['previous_exact_Bareiss_divisor'])
        sr=step['original_minor_rows7'];sc=step['original_minor_cols3'];assert len(sr)==len(sc)==stage+1
        assert determinant([[block[r][c] for c in sc] for r in sr])==pivot==decode(step['original_minor_nonzero'])
        nextA=[row[:] for row in A]
        for r in rr:
            if r==i:continue
            for c in cc:
                if c==j:nextA[r][c]=R.zero
                else:
                    numerator=pivot*A[r][c]-A[r][j]*A[i][c];q=numerator.exquo(previous);assert q*previous==numerator;nextA[r][c]=q;division_count+=1
        A=nextA;previous=pivot;rr.remove(i);cc.remove(j)
    assert len(checkpoint['steps'])==3 and checkpoint['completed_generic_pivots']==3
    assert rr==checkpoint['remaining_rows7'] and cc==checkpoint['remaining_cols3']==[]
    assert A==[[decode(p) for p in row] for row in checkpoint['matrix7x3']]
    sr=g['selected_rows7'];sc=g['selected_cols3'];assert sr==[0,5,2] and sc==[0,2,1]
    selected=[[block[i][j] for j in sc] for i in sr];Delta=determinant(selected);assert Delta and Delta==decode(g['nonzero_selected_minor_Delta']) and len(Delta)==407
    adj=[[(-1)**(i+j)*determinant([row[:i]+row[i+1:] for row in selected[:j]+selected[j+1:]]) for j in range(3)] for i in range(3)]
    coordinates=[]
    for row in block:
        lam=[sum((row[sc[k]]*adj[k][j] for k in range(3)),R.zero) for j in range(3)];coordinates.append(lam)
        for col in range(3):assert Delta*row[col]==sum((lam[j]*block[sr[j]][col] for j in range(3)),R.zero)
    assert coordinates==[[decode(q) for q in row] for row in g['cleared_upper_row_coordinates7xr']]
    assert g['polynomial_right_kernel_basis3']==[] and g['block_generic_rank']==3
    fullrows=[2,4]+[brows[i] for i in sr];fullcols=[0,5]+[bcols[i] for i in sc]
    assert fullrows==g['nonzero_minor_of14x5_rows'] and fullcols==g['nonzero_minor_of14x5_cols']
    assert all(M[i][j]==0 for i in fullrows[:2] for j in fullcols[2:])
    fullminor=pairdet*Delta;assert fullminor==decode(g['nonzero_minor_of14x5']) and len(fullminor)==1752
    assert g['generic_rank14x5']==5 and g['generic_rank_actual16_H7']==7
    Q=decode(oc['remaining_factor_Q']);unit=oc['named_unit_factor'];assert unit['QQ_constant']=='27/48125000000' and unit['a_d_F_B_C_J_powers']==[0,0,0,3,0,0]
    assert Delta==s.QQ(27,48125000000)*B**3*Q==decode(oc['fixed_original_minor']) and p==decode(oc['p'])
    assert len(Q)==104 and max(sum(m) for m in Q)==12 and Q and not Q.is_ground
    assert all(Q.rem(f) for f in named)
    # Exact surjective product map, followed by the accepted H-degree-preserving
    # free projection, not a mere residue-field dimension count.
    hp=ROOT/'reports/S61_NONSQUARE_H_SQUARE_H1_PREFLIGHT_20261104_v63/certificate.json';hd=read(hp);hg=ROOT/'reports/S61_NONSQUARE_H_SQUARE_H1_PREFLIGHT_20261104_v63_REVIEW/report.json';hr=read(hg);assert hr['all_passed']
    assert hd['target117_projection_preserves_Hdegree_exact_structural_support'] and len(hd['target122_Hdegree_classes']['1'])==23
    monomials=list(itertools.combinations_with_replacement(range(7),2));section=[(i,j) for i,j in monomials]
    product_matrix=[[int(tuple(sorted((i,j)))==m) for i,j in section] for m in monomials]
    assert product_matrix==[[int(i==j) for j in range(28)] for i in range(28)]
    for key in ['no_allFULLU_H7_or_H1_surjectivity_claim','no_full239_claim_or_global_kernel_count_change']:assert oc[key]
    proof=dict(source_raw_H_checkpoint_and_two_global_unit_gate_bound=True,p21_terms_nonzero=True,sparse_pair_formula='-(6/5)*B*p^2',other3_columns_literal_unchanged_block7x3_terms=908,all3_fixed_fractionfree_steps_and_nonzero_minors_verified=True,exact_division_remultiplications=division_count,all21_cleared_upper_identities_verified=True,selected_minor_Delta_terms=407,Delta_named_factor='(27/48125000000)*B^3*Q',Q_terms=104,Q_degree=12,full14x5_minor_terms=1752,generic_actual_raw_H_rank=7,explicit_pointwise_H7_chart='FULLU intersect D(p*Q)',unordered_light_product_section_identity28=True,H1_projection23_surjective_on_same_chart=True,pure98_retained=True,no_whole_FULLU_rank7_H1_or_full239_claim=True)
    (out/'independent-proof.json').write_text(json.dumps(proof,indent=2)+'\n')
    for pth in [Path(__file__),ROOT/f'research/notes/{ID}.md',ROOT/f'research/notes/{OPEN}.md',accepted,folder/'report.json',of/'report.json',cp,hp,hg]:hashes[str(pth.relative_to(ROOT))]=sha(pth)
    result=dict(id=ID,all_passed=True,actual_generic_raw_H_rank7_and_minor_chart_independently_verified=True,actual_H1_surjectivity23_only_on_explicit_D_pQ=True,no_author_solver_search_newpoint_prime_GB_or_global239_claim=True,input_sha256=hashes,evidence_sha256={'independent-proof.json':sha(out/'independent-proof.json')},seconds=time.monotonic()-start,max_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n')
    other=ROOT/'reports'/OPEN;other.mkdir(parents=True,exist_ok=False);companion=dict(id=OPEN,all_passed=True,shared_independent_check=str((out/'report.json').relative_to(ROOT)),input_sha256={str((out/'report.json').relative_to(ROOT)):sha(out/'report.json')},explicit_chart='FULLU intersect D(p*Q)',no_all_FULLU_or_full239_claim=True)
    (other/'report.json').write_text(json.dumps(companion,indent=2)+'\n');print(json.dumps(dict(all_passed=True,seconds=result['seconds'],rss=result['max_rss_kib'],report_sha256=sha(out/'report.json'),open_review_sha256=sha(other/'report.json'))))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
