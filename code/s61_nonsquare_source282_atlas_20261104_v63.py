#!/usr/bin/env python3
"""One sparse source preprocessing: actual174x22 and fixed15 source deletions.

No maximal minor, pivot search, actual quartic matrix or rank is evaluated.
"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import resource
import signal
import time
import sympy as s
from sympy.polys.rings import ring

ROOT=Path(__file__).resolve().parents[1]
ID='S61_NONSQUARE_SOURCE282_ATLAS_20261104_v63'
def read(p):return json.loads(p.read_text())
def sha(p):
    h=hashlib.sha256()
    with p.open('rb')as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()

def main(out):
    start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60)
    out.mkdir(parents=True,exist_ok=False)
    paths={
      'quads':ROOT/'reports/S61_BINARY_HYPERPLANE_QUADRATIC_CODIMENSION_20261102.json',
      'quad22':ROOT/'reports/S61_NONSQUARE_UNIVERSAL_QUAD22_20261103_v63/certificate.json',
      'quad22_report':ROOT/'reports/S61_NONSQUARE_UNIVERSAL_QUAD22_20261103_v63/report.json',
      'quad22_review':ROOT/'reports/S61_NONSQUARE_UNIVERSAL_QUAD22_20261103_v63_REVIEW/report.json',
      'source319':ROOT/'reports/S61_NONSQUARE_QUARTIC_319_ROWS_20261103_v63/certificate.json',
      'quartic11':ROOT/'reports/S61_NONSQUARE_BIANCHI_QUARTIC_UNIT_PIVOTS_20261104_v63/certificate.json',
      'quartic11_review':ROOT/'reports/S61_NONSQUARE_BIANCHI_QUARTIC_UNIT_PIVOTS_20261104_v63_REVIEW/report.json',
      'cubic4':ROOT/'reports/S61_NONSQUARE_BIANCHI_RESIDUAL_FOUR_20261104_v63/certificate.json',
      'cubic4_review':ROOT/'reports/S61_NONSQUARE_BIANCHI_RESIDUAL_FOUR_20261104_v63_REVIEW/report.json',
      'containment':ROOT/'reports/S61_NONSQUARE_BIANCHI_WEIGHT4_CONTAINMENT_20261104_v63/report.json',
      'codim2_review':ROOT/'reports/S61_NONSQUARE_FIXED_QUARTIC_CODIMENSION_TRANSFER_20261103_v63_REVIEW.json',
      'reflexive64':ROOT/'reports/S61_NONSQUARE_QUADRATIC_QUARTIC_QUOTIENT_REFLEXIVE_20261103_v63.json',
      'reflexive64_review':ROOT/'reports/S61_NONSQUARE_QUADRATIC_QUARTIC_QUOTIENT_REFLEXIVE_20261103_v63_REVIEW.json',
      'circuit122':ROOT/'reports/S61_NONSQUARE_QUARTIC_122_CIRCUIT_20261103_v63/circuit.json'}
    D={k:read(p)for k,p in paths.items()}
    for key in ['quad22_review','quartic11_review','cubic4_review','codim2_review']:assert D[key]['all_passed']
    R,a,e,d=ring('a,e,d',s.QQ);K=R.to_field()
    T,x0,x1,x2,x3,x4,Z,W=ring('x0,x1,x2,x3,x4,Z,W',K);xs=[x0,x1,x2,x3,x4,Z,W]
    V,*vs=ring('v0,v1,v2,v3,v4,v6,v7',K)
    loadpoly=lambda ts:K(R.from_dict({tuple(ex):s.QQ(c)for ex,c in ts}))
    def dec(p):return T.from_dict({tuple(mon):loadpoly(c)for mon,c in p})
    Q=D['quad22'];L=dec(Q['L']);ell=dec(Q['ell']);TZ=dec(Q['TZ']);TW=dec(Q['TW'])
    images=[x0,x1,x2,W,x3,T.zero,Z+x1-3*x0,x4]
    qs=[]
    for row in D['quads']['quadratic_rows_sparse']:
        q=T.zero
        for mon,c in row:
            z=T(loadpoly(c))
            for im,k in zip(images,mon):
                if k:z*=im**k
            q+=z
        qs.append(q)
    assert len(qs)==6
    ch=[[K.zero for _ in range(6)]for _ in range(6)]
    ch[0][0]=ch[1][1]=ch[2][3]=ch[3][5]=K.one
    ch[4][1]=K(10*a+3*e**2-3*e)/K(a);ch[4][3]=-K(3*(e-1))/K(a);ch[4][4]=K(-10)
    ch[5][1]=K(-(d*e+8));ch[5][2]=K(8);ch[5][3]=K(d)
    prime=[qs[0],qs[1],qs[3],qs[5],ell*Z,ell*W]
    assert all(sum((T(c)*q for c,q in zip(row,qs)),T.zero)==p for row,p in zip(ch,prime))
    def pullback(vec):return[sum((T(ch[j][i])*vec[j]for j in range(6)),T.zero)for i in range(6)]
    columns=[];labels=[]
    for i,j in itertools.combinations(range(6),2):
        v=[T.zero]*6;v[i]=qs[j];v[j]=-qs[i];columns.append(v);labels.append(['original_quadric_Koszul',i,j])
    lin=pullback([T.zero,T.zero,T.zero,T.zero,W,-Z])
    for j in[0,2,3,4,5,6]:columns.append([z*xs[j]for z in lin]);labels.append(['linear_relation_times_cert_light',j,str(xs[j])])
    new=pullback([T.zero,ell*(T(K(e))*L-T(K(400*a))*x0),-ell*L,T.zero,-TZ,-TW])
    columns.append(new);labels.append(['twenty_second_actual_relation'])
    assert len(columns)==22
    # Literal identities in the ACTUAL six original quadrics, including each
    # of the15 Koszul and the full newly lifted22nd relation.
    for col in columns:assert sum((c*q for c,q in zip(col,qs)),T.zero)==0
    # Explicit inverse light-coordinate map into the frozen original Tn frame.
    toV=[vs[0],vs[1],vs[2],vs[4],vs[6],3*vs[0]-vs[1]+vs[5],vs[3]]
    def transport(p):
        q=V.zero
        for mon,c in p.items():
            z=V(c)
            for im,k in zip(toV,mon):
                if k:z*=im**k
            q+=z
        return q
    C=D['source319'];source319=[C['source_rows'][i]for i in C['retained_row_indices']]
    assert len(source319)==319 and all(row[0]=='quadric'for row in source319[:174])
    look={(row[1],tuple(row[2])):j for j,row in enumerate(source319[:174])}
    encode=lambda p:[[list(ex),str(c)]for ex,c in p.items()]
    encfrac=lambda c:[encode(c.numer),encode(c.denom)]
    known=[a,500*a+90*e**2+3*e-24,(50*a+18*e**2-9*e-2)*d+60*e-40,100*a+18*e**2-9*e-2]
    sparse=[];terms=0
    for col in columns:
        entries=[]
        for qi,p in enumerate(col):
            for mon,c in transport(p).items():
                assert sum(mon)==2
                mon8=tuple([*mon[:5],0,*mon[5:]]);j=look[qi,mon8]
                den=c.denom
                for f in known:
                    while not den.rem(f):den=den.exquo(f)
                assert den.is_ground and den
                entries.append([j,*encfrac(c)]);terms+=len(c.numer)+len(c.denom)
        assert len({x[0]for x in entries})==len(entries)
        assert all(source319[x[0]][2][5]==0 for x in entries)
        sparse.append(sorted(entries,key=lambda x:x[0]))
    # Actual fixed15 SOURCE-coordinate deletion: its15x15 coefficient minor
    # is triangular I11/I4 by the two accepted arithmetic packets.
    parent=D['quartic11'];qlookup={row[1]:i for i,row in enumerate(source319)if row[0]=='quartic'}
    q_original_rows=[source319[i][1]for i in range(286,319)]
    delete11=[qlookup[q_original_rows[j]]for j in parent['pivot_quartic_columns']]
    delete4=[p['source319_column']for p in D['cubic4']['unit_pivot_steps']]
    assert len(set(delete11+delete4))==15
    assert all(i>=174 for i in delete11+delete4)
    retained=[i for i in range(319)if i not in delete11+delete4]
    assert len(retained)==304 and retained[:174]==list(range(174))
    kernel=dict(shape=[174,22],coefficient_core_order=['a','e','d'],
        coefficient_encoding='column list of[source174row,numerator polynomial,denominator polynomial]',
        sparse_columns=sparse,column_labels=labels,
        exact_all22_actual_original_quadric_identities_passed=True,
        six_quadric_change_to_old_prime_basis=[[str(c)for c in row]for row in ch],
        change_determinant='-80',Koszul_basis_change_determinant='(-80)^5',
        light_frame='cert(x0,x1,x2,x3,x4,Z,W)=(v0,v1,v2,v4,v7,3v0-v1+v6,v3)',
        actual_source174_definitions=source319[:174])
    (out/'kernel22.json').write_text(json.dumps(kernel,separators=(',',':'))+'\n')
    atlas=dict(source319_definitions=source319,deleted15_source319_indices_in_pivot_order=delete11+delete4,
        deleted15_coefficient_minor='[[I11,*],[0,I4]], determinant1',
        retained304_source319_indices=retained,source304_definitions=[source319[i]for i in retained],
        source304_first174_quadratic_indices=list(range(174)),source304_remaining130_indices=list(range(174,304)),
        source_bundle='E152=coker(K:R22->R174); P282=E152 direct_sum R130',target_bundle='R239',
        source304_coefficients='Precisely the unchanged319 rows with these15 source rows deleted; same actual image/cokernel',
        minor_atlas='For every22-subset sigma of[0,174), delta_sigma=detK[sigma,:]. OnD(delta_sigma), delete these22 source304 rows, obtaining actual282x239.',
        local_source282_index_recipe='retained304_source319_indices with the22 quadratic indices sigma removed',
        no_selected_global22minor_or_global_basis_saved=True,
        global_split_recipe='Any P with P*K=I22 exists. Split R130=R22+R108; if M304=[A B C], a global282 presentation is[A+B*P,C]. P not computed.',
        all22_minor_opens_cover_FULLU_by_accepted_pointwise_kernel22=True,
        C122_recipe='On accepted DeltaH DeltaQ charts, quotient sourceR122 by descended15 subbundle gives rank107->R64; no explicit107 matrix evaluated.',
        C122_kernel_projection='After accepted target/source Schur map(p,z)->(p+A*z,C122*z), a kernel vector maps toz, with inversez->(-A*z,z).',
        no_actual_fullquartic_or_C122_coefficients_evaluated=True)
    (out/'atlas.json').write_text(json.dumps(atlas,indent=2)+'\n')
    report=dict(id=ID,all_passed=True,all22_literal_original_quadric_identities=True,
        sparse_K_shape=[174,22],sparse_nonzero_entries=sum(map(len,sparse)),core_numden_terms=terms,
        fixed_source_deletions=15,actual_fixed_reduced_source_shape=[304,239],
        projective_source_bundle_rank=282,explicit_single_global282_basis_not_claimed=True,
        no_minor_or_pivot_search_no_quartic_matrix_expansion=True,
        no_new_rank_field_point_GB_or_full239_claim=True,
        input_sha256={str(p.relative_to(ROOT)):sha(p)for p in [Path(__file__),*paths.values()]},
        evidence_sha256={n:sha(out/n)for n in ['kernel22.json','atlas.json']},
        seconds=time.monotonic()-start,frozen=True,independent_review_pending=True)
    for n in ['kernel22.json','atlas.json']:assert(out/n).stat().st_size<32*1024*1024
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');signal.alarm(0)
    print(json.dumps(report,indent=2),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
