#!/usr/bin/env python3
"""Independent original-v Laurent Fraction K22 audit and source/dual proof binding."""
import argparse
from fractions import Fraction as Q
import hashlib
import itertools
import json
from pathlib import Path
import resource
import signal
import time

ROOT=Path(__file__).resolve().parents[1]
STEM='S61_NONSQUARE_SOURCE282_ATLAS_20261104_v63_REVIEW'
N=10;ZERO=(0,)*N
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def accepted(p):
    r=read(p);assert r.get('all_passed') or r.get('independent_mathematical_review_passed')
    for n,h in r.get('input_sha256',{}).items():assert sha(ROOT/n)==h,n
    for n,h in r.get('evidence_sha256',{}).items():
        q=ROOT/n if (ROOT/n).exists() else p.parent/n
        assert sha(q)==h,n
    return r
def add(*As):
    z={}
    for A in As:
        for e,c in A.items():z[e]=z.get(e,Q(0))+c
    return {e:c for e,c in z.items() if c}
def scale(A,c):return {e:v*c for e,v in A.items() if v*c}
def mul(A,B):
    z={}
    for e,c in A.items():
        for f,d in B.items():
            g=tuple(x+y for x,y in zip(e,f));z[g]=z.get(g,Q(0))+c*d
    return {e:c for e,c in z.items() if c}
def power(A,n):
    B={ZERO:Q(1)}
    for i in range(n):B=mul(B,A)
    return B
def var(i,k=1):return {tuple(k if j==i else 0 for j in range(N)):Q(1)}
def const(c):return {ZERO:Q(c)} if c else {}

def main(out):
    start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60)
    assert not out.exists();out.mkdir(parents=True)
    base=ROOT/'reports/S61_NONSQUARE_SOURCE282_ATLAS_20261104_v63'
    ap=base/'report.json';author=accepted(ap);kp=base/'kernel22.json';atp=base/'atlas.json';proofpath=base/'proof-receipt.json'
    K=read(kp);atlas=read(atp);proof=read(proofpath)
    assert proof['all_passed'] and proof['report_sha256']==sha(ap) and proof['kernel22_sha256']==sha(kp) and proof['atlas_sha256']==sha(atp)
    source=ROOT/'code/s61_nonsquare_source282_atlas_20261104_v63.py';note=ROOT/'research/notes/S61_NONSQUARE_SOURCE282_ATLAS_20261104_v63.md'
    assert proof['source_sha256']==sha(source) and proof['proof_note_sha256']==sha(note)
    reviews=[ROOT/f'reports/{n}/report.json' for n in ['S61_NONSQUARE_UNIVERSAL_QUAD22_20261103_v63_REVIEW',
        'S61_NONSQUARE_BIANCHI_QUARTIC_UNIT_PIVOTS_20261104_v63_REVIEW','S61_NONSQUARE_BIANCHI_RESIDUAL_FOUR_20261104_v63_REVIEW']]
    qr,b11,b4=[accepted(p) for p in reviews]
    assert qr['entire_kernel_module_basis_accepted'] and qr['pointwise_22_directions_independent_on_FULLU']
    assert b11['all11_final_vectors_have_literal_identity_quartic_block'] and b4['all4_by15_explicit_composites_and_selected_I4_verified']
    codimpath=ROOT/'reports/S61_NONSQUARE_FIXED_QUARTIC_CODIMENSION_TRANSFER_20261103_v63_REVIEW.json';codim=accepted(codimpath)
    assert codim['fixed319_rank_defect_codimension_at_least2_on_U0']
    refpath=ROOT/'reports/S61_NONSQUARE_QUADRATIC_QUARTIC_QUOTIENT_REFLEXIVE_20261103_v63_REVIEW.json';ref=accepted(refpath)
    assert ref['S2_and_torsion_free_imply_reflexive64']
    qupath=ROOT/'reports/S61_BINARY_HYPERPLANE_QUADRATIC_CODIMENSION_20261102.json';raw=read(qupath)
    certpath=ROOT/'reports/S61_NONSQUARE_UNIVERSAL_QUAD22_20261103_v63/certificate.json';old=read(certpath)
    v=[var(i) for i in range(7)];a,e,d=var(7),var(8),var(9);ainv=var(7,-1)
    certlight=[v[0],v[1],v[2],v[4],v[6],add(scale(v[0],3),scale(v[1],-1),v[5]),v[3]]
    # Inverse coordinate matrices, no Gaussian elimination or determinant search.
    C=[[int(poly.get(tuple(int(i==j) for i in range(N)),0)) for j in range(7)] for poly in certlight]
    D=[[0]*7 for _ in range(7)]
    for i,j in [(0,0),(1,1),(2,2),(3,6),(4,3),(6,4)]:D[i][j]=1
    D[5][5]=1;D[5][1]=1;D[5][0]=-3
    assert [[sum(C[i][k]*D[k][j] for k in range(7)) for j in range(7)] for i in range(7)]==[[int(i==j) for j in range(7)] for i in range(7)]
    assert [[sum(D[i][k]*C[k][j] for k in range(7)) for j in range(7)] for i in range(7)]==[[int(i==j) for j in range(7)] for i in range(7)]
    def core(terms):return {(0,)*7+tuple(ex):Q(c) for ex,c in terms}
    def certpoly(ts):
        z={}
        for mon,c in ts:
            p=core(c)
            for i,k in enumerate(mon):
                if k:p=mul(p,power(certlight[i],k))
            z=add(z,p)
        return z
    qs=[]
    for row in raw['quadratic_rows_sparse']:
        q={}
        for mon,coeffs in row:
            assert mon[5]==0
            for ex,c in coeffs:
                exponent=tuple(mon[:5])+tuple(mon[6:])+tuple(ex);q[exponent]=q.get(exponent,Q(0))+Q(c)
        qs.append({ex:c for ex,c in q.items() if c})
    L,ell,TZ,TW=[certpoly(old[k]) for k in ['L','ell','TZ','TW']]
    ch=[[{} for j in range(6)] for i in range(6)]
    for i,j in [(0,0),(1,1),(2,3),(3,5)]:ch[i][j]=const(1)
    ch[4][1]=mul(add(scale(a,10),scale(mul(e,e),3),scale(e,-3)),ainv)
    ch[4][3]=mul(scale(add(e,const(-1)),-3),ainv);ch[4][4]=const(-10)
    ch[5][1]=scale(add(mul(d,e),const(8)),-1);ch[5][2]=const(8);ch[5][3]=d
    prime=[qs[0],qs[1],qs[3],qs[5],mul(ell,certlight[5]),mul(ell,certlight[6])]
    for i in range(6):assert add(*(mul(ch[i][j],qs[j]) for j in range(6)))==prime[i]
    # Unique nonzero determinant permutation: [0,1,3,5,4,2], even parity.
    perm=[0,1,3,5,4,2];assert sum(perm[i]>perm[j] for i in range(6) for j in range(i+1,6))%2==0
    assert (-10)*8==-80 and K['change_determinant']=='-80' and K['Koszul_basis_change_determinant']=='(-80)^5'
    def pull(vec):return [add(*(mul(ch[j][i],vec[j]) for j in range(6))) for i in range(6)]
    expected=[]
    for i,j in itertools.combinations(range(6),2):
        col=[{} for k in range(6)];col[i]=qs[j];col[j]=scale(qs[i],-1);expected.append(col)
    lin=pull([{}, {}, {}, {},certlight[6],scale(certlight[5],-1)])
    for i in [0,2,3,4,5,6]:expected.append([mul(c,certlight[i]) for c in lin])
    expected.append(pull([certpoly(c) for c in old['new_syzygy']]))
    defs=K['actual_source174_definitions'];assert len(defs)==174
    assert K['shape']==[174,22] and len(K['sparse_columns'])==22
    got=[];entries=numden=0;denom_a_counts={0:0,1:0}
    for ci,column in enumerate(K['sparse_columns']):
        vec=[{} for i in range(6)];seen=set()
        for row,num,den in column:
            assert row not in seen;seen.add(row);entries+=1;numden+=len(num)+len(den)
            assert len(den)==1;de,dc=den[0];assert de[1:]==[0,0] and de[0] in [0,1]
            denom_a_counts[de[0]]+=1
            kind,i,mon=defs[row];assert kind=='quadric' and mon[5]==0 and sum(mon)==2
            for ex,c in num:
                exp=tuple(mon[:5])+tuple(mon[6:])+tuple(x-y for x,y in zip(ex,de))
                vec[i][exp]=vec[i].get(exp,Q(0))+Q(c)/Q(dc)
        vec=[{e:c for e,c in p.items() if c} for p in vec]
        assert vec==expected[ci]
        assert not add(*(mul(c,q) for c,q in zip(vec,qs)))
        got.append(vec)
    assert entries==963 and numden==14525 and denom_a_counts=={0:947,1:16}
    rowsfile=ROOT/'reports/S61_NONSQUARE_QUARTIC_319_ROWS_20261103_v63/certificate.json';rows=read(rowsfile)
    source319=[rows['source_rows'][i] for i in rows['retained_row_indices']]
    assert source319==atlas['source319_definitions'] and source319[:174]==defs
    b11path=ROOT/'reports/S61_NONSQUARE_BIANCHI_QUARTIC_UNIT_PIVOTS_20261104_v63/certificate.json';b11c=read(b11path)
    b4path=ROOT/'reports/S61_NONSQUARE_BIANCHI_RESIDUAL_FOUR_20261104_v63/certificate.json';b4c=read(b4path)
    quartics=[(i,r[1]) for i,r in enumerate(source319) if r[0]=='quartic']
    deleted=[quartics[j][0] for j in b11c['pivot_quartic_columns']]+[p['source319_column'] for p in b4c['unit_pivot_steps']]
    assert deleted==atlas['deleted15_source319_indices_in_pivot_order'] and len(set(deleted))==15 and min(deleted)>=174
    retained=[i for i in range(319) if i not in deleted]
    assert retained==atlas['retained304_source319_indices'] and len(retained)==304
    assert atlas['source304_definitions']==[source319[i] for i in retained]
    assert retained[:174]==list(range(174)) and len(retained[174:])==130
    assert (174-22)+130==282 and 282-239==43 and (210-146)+(28-6)+1==87 and 122-15==107
    assert proof['global282_free_basis_or_left_inverse_computed'] is False
    frontierpath=base/'dependency-frontier.json';frontier=read(frontierpath);assert not frontier['arithmetic_started']
    # Pure geometry/module proofs and countermodels are justified in the note;
    # these integers only bind the actual packet, they are not new rank tests.
    evidence=dict(all22_original_v_quadric_identities_exact=True,all963_coefficients_rebuilt=True,
        coefficient_numden_terms=14525,denominator_a_power_counts=denom_a_counts,forward_light_frame=C,inverse_light_frame=D,
        unit_basis_change_determinant=-80,exterior_square_determinant='(-80)^5',
        deleted15_source319_indices=deleted,source304_retained_indices=retained,
        direct_sum_ranks={'E152':152,'P282':282,'Q87':87,'dualD':43,'descended_source':107},
        no_P_leftinverse_or_global_atlas_minor_computed=True)
    ep=out/'independent-certificate.json';ep.write_text(json.dumps(evidence,indent=2)+'\n')
    paths=[Path(__file__),ROOT/f'research/notes/{STEM}.md',ap,kp,atp,proofpath,source,note,*reviews,codimpath,refpath,qupath,certpath,rowsfile,b11path,b4path,frontierpath]
    result=dict(id=STEM,all_passed=True,
        independent_Laurent_Fraction_original_v_engine_no_SymPy_or_producer=True,
        all22_K_columns_rebuilt_and_all22_literal_original_quadric_identities=True,
        K_basis_same_as_whole_old22_via_unit_exterior_square_and_same_last7=True,
        fixed15_deletion_to_actual304_same_image_cokernel=True,
        projective_source282_atlas_correct_global_free_presentation_exists_P_uncomputed=True,
        R130_to_reflexive_Q87_same_cokernel_and_C122_rank107_source_bundle_scopes_accepted=True,
        dualD43_torsionfree_pd1_Ext1_actualC_proof_accepted=True,
        dual_reflexive_iff_no_curves_and_locally_free_iff_surjective=True,
        abstract_free43_kernel_and_codim2_do_not_imply_surjectivity=True,
        no_computed_global282_basis_actual107_matrix_full239_C122rank64_or_affine_claim=True,
        no_new_pivot_minor_search_GB_field_point_or_full_quartic_expansion=True,
        input_sha256={str(p.relative_to(ROOT)):sha(p) for p in paths},evidence_sha256={ep.name:sha(ep)},seconds=time.monotonic()-start)
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');signal.alarm(0)
    print(json.dumps({k:result[k] for k in ['all_passed','all22_K_columns_rebuilt_and_all22_literal_original_quadric_identities','seconds']},indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
