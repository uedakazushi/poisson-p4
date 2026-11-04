#!/usr/bin/env python3
"""Independent full135 ell restriction and fixed QQ Bareiss minor audit."""
import argparse
from fractions import Fraction as Q
import hashlib
import json
from math import gcd
from pathlib import Path
import resource
import signal
import time

import sympy as s

ROOT=Path(__file__).resolve().parents[1]
ID='S61_NONSQUARE_BIANCHI_ELL_WEIGHT5_20261104_v63_REVIEW'


def sha(p):
    h=hashlib.sha256()
    with p.open('rb')as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()


def read(p):return json.loads(p.read_text())


def bareiss_rational(A):
    """Clear each rational row, then use fraction-free integer elimination."""
    B=[];denprod=1
    for row in A:
        d=1
        for c in row:d=d*c.denominator//gcd(d,c.denominator)
        denprod*=d;B.append([int(c*d)for c in row])
    n=len(B);previous=1;sign=1
    for k in range(n-1):
        pivot=next((i for i in range(k,n)if B[i][k]),None)
        if pivot is None:return Q(),B,denprod
        if pivot!=k:B[k],B[pivot]=B[pivot],B[k];sign=-sign
        p=B[k][k]
        for i in range(k+1,n):
            for j in range(k+1,n):
                numerator=p*B[i][j]-B[i][k]*B[k][j]
                assert numerator%previous==0
                B[i][j]=numerator//previous
        for i in range(k+1,n):B[i][k]=0
        previous=p
    return Q(sign*B[-1][-1],denprod),B,denprod


def main(out):
    start=time.monotonic();signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
    if out.exists():
        assert not any(out.iterdir()), 'refuse to overwrite evidence'
    else:out.mkdir(parents=True)
    base=ROOT/'reports/S61_NONSQUARE_BIANCHI_ELL_WEIGHT5_20261104_v63'
    author=read(base/'report.json');cert=read(base/'certificate.json');pre=read(base/'preflight.json')
    assert author['all_exact_source_frame_and_oldpoint_checks_passed']
    assert author['no_generic_core_combination_divisible_by_this_ell_proved']
    for rel,h in author['input_sha256'].items():assert sha(ROOT/rel)==h,rel
    for rel,h in author['evidence_sha256'].items():assert sha(base/rel)==h,rel
    opath=ROOT/'reports/S61_NONSQUARE_FULL_BIANCHI_OPERATOR_20261104_v63/operator.json'
    op=read(opath);chartpath=ROOT/'reports/S61_BINARY_HYPERPLANE_WEIGHTED_CHART_20261102.json'
    chart=read(chartpath)
    adapterpath=ROOT/'reports/S61_NONSQUARE_ORACLE_FRAME_ADAPTER_20261103_v63/adapter.json'
    adapter=read(adapterpath)
    anchorpath=ROOT/'reports/S61_NONSQUARE_UNIVERSAL_QUAD22_20261103_v63/certificate.json'
    anchor=read(anchorpath)
    gp=ROOT/'reports/S61_BINARY_HYPERPLANE_GENERIC_RANKS_20261102.json'
    assert read(gp)['all_passed']and read(gp)['point']==[1,2,3]
    opreview=ROOT/'reports/S61_NONSQUARE_FULL_BIANCHI_OPERATOR_20261104_v63_REVIEW/report.json'
    adreview=ROOT/'reports/S61_NONSQUARE_ORACLE_FRAME_ADAPTER_20261103_v63_REVIEW/report.json'
    assert read(opreview)['all_passed']and read(adreview)['all_passed']
    w=op['local_coordinate_weights'];ew=op['equation_weights']
    assert w==chart['local_coordinate_weights']and ew==chart['equation_weights']
    assert len(w)==45 and len(ew)==135
    assert [i for i,z in enumerate(w)if z==0]==[19,23]
    T,U,normal=chart['T34'],chart['U34'],chart['normal34']
    assert sorted(normal+[1,2,5])==list(range(11))
    for i in(19,23):
        assert all(T[i][j]=='0'for j in normal+[1])and all(z=='0'for z in U[i])
    assert [[T[i][j]for j in(2,5)]for i in(19,23)]==[['1','0'],['0','1']]
    pieces=[(i,q)for i,q in enumerate(op['local_weight_pieces'])if q['weight']==5]
    assert len(pieces)==20 and [i for i,_ in pieces]==pre['source_piece_indices']
    assert op['equation_weights'][134]==3
    assert [q['weight']for q in op['local_weight_pieces']if q['spatial_degree']=='extra_tG6']==[6]
    # Reconstruct the coordinate permutation/shear from its coordinate
    # meaning, rather than copying the primary's fixed8 matrix.
    zz=s.symbols('z0:8')
    oldv=[zz[0],zz[1],zz[2],zz[6],zz[3],zz[7],zz[5]+zz[1]-3*zz[0],zz[4]]
    M=[[int(s.diff(q,z))for z in zz]for q in oldv]
    assert M==pre['M_oldv8_from_newnormal8']
    a,e,d=s.symbols('a e d');point={a:s.Integer(1),e:s.Integer(2),d:s.Integer(3)}
    cache={}
    def ev(expr):
        if expr not in cache:
            q=s.cancel(s.sympify(expr).subs(point,simultaneous=True))
            assert q.is_Rational
            cache[expr]=Q(int(q.p),int(q.q))
        return cache[expr]
    def core(terms):
        # Exact sparse integer-power arithmetic at the SAME saved old point.
        answer=Q()
        for ex,co in terms:
            assert len(ex)==3 and all(i>=0 for i in ex)
            answer+=Q(co)*(2**ex[1])*(3**ex[2])
        return answer
    diagonal=list(map(ev,adapter['Tden']))
    frame=[[Q()for _ in range(43)]for _ in range(45)]
    for i in range(45):
        for j in range(8):frame[i][j]=sum((ev(T[i][normal[k]])*diagonal[k]*M[k][j]for k in range(8)),Q())
        frame[i][8]=ev(T[i][1])
        for j in range(34):frame[i][9+j]=ev(U[i][j])
    assert all(not q for i in(19,23)for q in frame[i])
    newweights=[1]*7+[2,1]
    for j in range(34):
        support={w[i]for i in range(45)if U[i][j]!='0'}
        assert len(support)==1;newweights.append(support.pop())
    assert newweights==pre['new_variable_weights']
    block_evidence=[]
    for weight in sorted(set(newweights)):
        rr=[i for i in range(45)if w[i]==weight];cc=[j for j in range(43)if newweights[j]==weight]
        assert len(rr)==len(cc)
        determinant,_,_=bareiss_rational([[frame[i][j]for j in cc]for i in rr]);assert determinant
        saved=next(q for q in cert['frame_block_proofs']if q['weight']==weight)
        assert determinant==Q(saved['exact_determinant'])
        block_evidence.append(dict(weight=weight,size=len(rr),determinant=str(determinant)))
    ell=[Q()for _ in range(8)]
    for mon,terms in anchor['ell']:
        assert sum(mon)==1;ell[mon.index(1)]+=core(terms)
    assert list(map(str,ell))==cert['oldpoint_ell_newnormal8']
    C=100+18*2**2-9*2-2;J=5*3-6*2-4
    assert ell[1]==-Q(16,243)*C*J and ell[1]and ell[7]==0
    # Reverse the primary's order: FIRST restrict every raw-coordinate
    # functional to ell=0, THEN apply each full135 cofactor vector.
    restricted=[row[:]for row in frame]
    for i in range(45):
        for j in range(7):
            if j!=1:restricted[i][j]-=frame[i][1]*ell[j]/ell[1]
        restricted[i][1]=Q()
    nonzero_frame=[[(j,q)for j,q in enumerate(row)if q]for row in restricted]
    constraints={}
    def put(key,column,value):
        if not value:return
        row=constraints.setdefault(key,{})
        row[column]=row.get(column,Q())+value
        if not row[column]:del row[column]
        if not row:del constraints[key]
    raw_records=0
    for column,(_,piece)in enumerate(pieces):
        raw={}
        for equation,var,terms in piece['coefficients']:
            assert 0<=equation<135 and -1<=var<45
            assert ew[equation]+(0 if var==-1 else w[var])==5
            raw[equation,var]=raw.get((equation,var),Q())+core(terms);raw_records+=1
        for(equation,var),co in raw.items():
            if var==-1:put((equation,-1),column,co)
            else:
                for j,q in nonzero_frame[var]:put((equation,j),column,co*q)
    stored={tuple(label):{j:Q(q)for j,q in row}for label,row in cert['all_nonzero_oldpoint_constraints']}
    assert constraints==stored and len(constraints)==1166
    assert all(var!=1 and 0<=equation<135 for equation,var in constraints)
    assert all(var==-1 and ew[equation]==5 or var>=0 and ew[equation]+newweights[var]==5 for equation,var in constraints)
    assert sum(ew[g]==5 and var==-1 for g,var in constraints)==27
    assert all(g!=134 for g,var in constraints)
    keys=[tuple(q)for q in cert['selected_constraint_keys']]
    minor=[[constraints[key].get(j,Q())for j in range(20)]for key in keys]
    assert minor==[[Q(q)for q in row]for row in cert['selected20_minor']]
    determinant,integer_final,row_denominators=bareiss_rational(minor)
    assert determinant==Q(cert['selected20_minor_determinant'])and determinant
    evidence=dict(
        arithmetic='own exact rational point substitution + ell-restricted raw frame FIRST + integer-cleared Bareiss determinant',
        same_accepted_core_point=[1,2,3],full135_equations_retained=True,
        independent_variables=43,only_eliminated_nuisance_rows=[19,23],
        frame_weight_blocks=block_evidence,ell=list(map(str,ell)),ell_unit='-16*C*J/243',
        raw_cofactor_records=raw_records,dense_constraint_shape=[135*43,20],
        equation_weight5_total_count=ew.count(5),nonzero_weight5_constant_rows_at_saved_point=27,
        all_nonzero_constraint_rows_matched=1166,selected_constraint_keys=[list(q)for q in keys],
        exact20_minor_determinant=str(determinant),integer_row_denominator_product=str(row_denominators),
        numerator_bits=abs(determinant.numerator).bit_length(),denominator_bits=determinant.denominator.bit_length(),
        source_column134_G6h_retained='local equation weight3; the20 weight5 pieces have this cofactor identically zero, checked directly',
        extra_operator_row175_tG6='true local identity weight6, not eligible for the20 weight5 pieces; distinct from source column134',
        negative_proof='A nonzero exact QQ specialization of a regular20-column minor implies a nonzero QQ(a,e,d) minor. Restriction to ell=0 is injective on the20-dimensional core-constant candidate space.',
        scope='ONLY the20 frozen true-weight5 pieces; coefficients constant in all43 local variables and rational in core; divisibility by this actual ell. No higher-weight or general polynomial-cofactor conclusion.')
    ep=out/'independent-certificate.json';ep.write_text(json.dumps(evidence,indent=2)+'\n')
    inputs=[Path(__file__),ROOT/'research/notes'/(ID+'.md'),base/'report.json',base/'certificate.json',base/'preflight.json',
            ROOT/'code/s61_nonsquare_bianchi_ell_weight5_20261104_v63.py',
            ROOT/'research/notes/S61_NONSQUARE_BIANCHI_ELL_WEIGHT5_20261104_v63.md',opath,chartpath,adapterpath,anchorpath,gp,opreview,adreview]
    report=dict(id=ID,all_passed=True,
        accepted_full_Bianchi_operator_and_frame_source_hashes_bound=True,
        independent_exact43_frame_and_actual_ell_restriction_before_cofactor_action=True,
        all135_F5_F6_extra_tG6_equation_slots_retained=True,
        all1166_nonzero_rows_independently_reconstructed_and_compared=True,
        independent_integer_Bareiss_fixed20_minor=True,
        exact_QQ_minor=str(determinant),no_nonzero_generic_core_constant_weight5_combination_divisible_by_this_ell=True,
        no_general_polynomial_cofactor_higher_weight_kernel_or_uniform_rank_claim=True,
        no_new_point_prime_GB_author_rerun_source_regeneration_or319_transport=True,
        input_sha256={str(q.relative_to(ROOT)):sha(q)for q in inputs},
        evidence_sha256={ep.name:sha(ep)},seconds=time.monotonic()-start,frozen=True)
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');signal.alarm(0)
    print(json.dumps({k:report[k]for k in('all_passed','exact_QQ_minor','seconds')},indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,default=ROOT/'reports'/ID)
    main(p.parse_args().output.resolve())
