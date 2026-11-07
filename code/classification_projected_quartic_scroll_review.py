#!/usr/bin/env python3
"""Independent all150/Euler/conormal/fixed-factor review of four scroll orbits."""
import argparse
from fractions import Fraction as F
import hashlib
import itertools as it
import json
from pathlib import Path
import re
import resource
import time

import classification_degenerate_quadric_q4_nil31_I_review as q

ROOT=Path(__file__).resolve().parents[1]

def parse(text):
    text=re.sub(r'\bx0\b','h',text)
    text=re.sub(r'\bz([0-2])\b',lambda m:'tau'+m[1],text)
    return q.parse(text)

def decode(vector):
    out={}
    for c,(ij,e)in zip(vector,q.COLS):
        if c:out[ij]=q.add(out.get(ij,{}),{e+(0,0,0):F(c)})
    return q.antisymmetric(out)

def div(P):return [q.add(*(q.der(P.get((i,j),{}),i)for i in range(5)))for j in range(5)]

def canonical(P):
    d=div(P)
    out=q.antisymmetric({(i,j):q.add(P.get((i,j),{}),q.scale(q.add(q.mul(q.atom(1,i),d[j]),q.scale(q.mul(q.atom(1,j),d[i]),-1)),-F(1,5)))for i,j in q.PAIRS})
    assert all(not p for p in div(out))
    return out

def quintic(P):
    out={}
    for missing in range(5):
        i,j,k,l=[a for a in range(5)if a!=missing]
        pf=q.add(q.mul(P.get((i,j),{}),P.get((k,l),{})),q.scale(q.mul(P.get((i,k),{}),P.get((j,l),{})),-1),q.mul(P.get((i,l),{}),P.get((j,k),{})))
        out=q.add(out,q.scale(q.mul(q.atom(1,missing),pf),(-1)**missing))
    return out

def matrix(data):return [[F(c)for c in row]for row in data]

def kills(A,B):
    cols=q.transpose(B)
    return all(not sum((v*col[j]for j,v in enumerate(row)if v and col[j]),F())for row in A for col in cols)

def minor(A,certificate):
    rows,cols=certificate['rows'],certificate['columns']
    value=q.determinant([[A[i][j]for j in cols]for i in rows])
    assert value==F(certificate['determinant'])and value
    return dict(rank=len(rows),rows=rows,columns=cols,determinant=str(value))

def bracket(P,left,right):return q.add(*(q.mul(q.mul(q.der(left,i),P.get((i,j),{})),q.der(right,j))for i,j in it.product(range(5),repeat=2)))

def poly_division(dividend,divisor):
    leading=max(divisor);co=divisor[leading];remainder=dict(dividend);quotient={}
    while remainder:
        mon=max(remainder)
        if not all(a>=b for a,b in zip(mon,leading)):break
        term={tuple(a-b for a,b in zip(mon,leading)):remainder[mon]/co}
        quotient=q.add(quotient,term);remainder=q.add(remainder,q.scale(q.mul(term,divisor),-1))
    assert not remainder
    assert q.mul(divisor,quotient)==dividend
    return quotient

def coeff(poly,coordinate_exponents):
    return {(0,)*5+e[5:]:c for e,c in poly.items()if e[:5]==coordinate_exponents}

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args();started=time.monotonic()
    resource.setrlimit(resource.RLIMIT_CPU,(170,175));resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,)*2)
    factorpath=ROOT/'reports/CLASSIFICATION_PROJECTED_QUARTIC_SCROLL_FACTORS_20261107.json'
    factors=json.loads(factorpath.read_text());paths=[factorpath,
        ROOT/'research/classification/PROJECTED_QUARTIC_SCROLL_ISOLATED_MODELS.md',
        ROOT/'research/classification/RANK_ZERO_SURFACE_GENERIC_FIRST_JET.md',
        ROOT/'reports/CLASSIFICATION_REGULAR_PENCIL_CI22_REVIEW_V2_20261107.json',
        ROOT/'reports/CLASSIFICATION_ALL_DEGENERATE_CI22_REVIEW_20261107.json']
    records=[]
    for factor in factors['cases']:
        path=ROOT/factor['input_report'];paths.append(path)
        assert hashlib.sha256(path.read_bytes()).hexdigest()==factor['input_sha256']
        saved=json.loads(path.read_text());model=saved['model'];params=saved['parametrization_variables']
        aliases={name:['h','x1','x2','x3'][i]for i,name in enumerate(params)}
        def paramparse(text):return q.parse(re.sub(r'\b(ss|tt|u|v|rr)\b',lambda m:aliases[m[1]],text))
        Qs=list(map(paramparse,saved['surface_parametrization']))
        pull=lambda p:q.substitute(p,Qs+[q.atom(1,i)for i in range(5,8)])
        columnlabels=[(tuple(c['pair']),tuple(c['exponents']))for c in saved['full150_coefficient_columns']]
        assert columnlabels==q.COLS
        labels=[(tuple(r['triple']),tuple(r['exponents']))for r in saved['all280_zero_on_surface_rows']]
        lookup={(tr,e+(0,)*(8-len(e))):i for i,(tr,e)in enumerate(labels)}
        M=[[F()for _ in range(150)]for _ in labels]
        for j,(ij,ex)in enumerate(q.COLS):
            image=pull({ex+(0,0,0):F(1)})
            for tri in it.combinations(range(5),3):
                if set(ij)<=set(tri):
                    k=next(i for i in tri if i not in ij)
                    for e,c in q.scale(q.mul(Qs[k],image),(-1)**tri.index(k)).items():M[lookup[tri,e]][j]=c
        claimM=saved.get('complete_zero_matrix',saved.get('complete280_by150_zero_matrix'))
        assert M==[[F(claimM[str(i)].get(str(j),'0'))for j in range(150)]for i in range(len(M))]
        K=matrix(saved['full_zero_kernel_basis']);B=matrix(saved['complete150_by_projective_source_matrix'])
        assert len(K)==len(B)==150 and len(K[0])==35 and len(B[0])==11
        assert kills(M,K)and kills(M,B)and q.rank(K)[0]==35 and q.rank(B)[0]==11
        mcert=minor(M,saved['zero_rank_minor']);assert mcert['rank']==115
        # Full24 Euler quotient and the extra eleventh class are retained.
        Euler=[]
        for a,b in it.product(range(5),repeat=2):
            Euler.append(q.vec({(i,j):q.add(q.atom(1,i,b)if j==a else{},q.atom(-1,j,b)if i==a else{})for i,j in q.PAIRS}))
        E=q.transpose(Euler);assert q.rank(E)[0]==24 and kills(M,E)
        assert q.rank([e+b for e,b in zip(E,B)])[0]==35
        polys=[decode(v)for v in q.transpose(B)]
        assert all(all(not p for p in div(P))for P in polys)
        for vector in q.transpose(K):
            can=q.vec(canonical(decode(vector)))
            assert q.rank([row+[v]for row,v in zip(B,can)])[0]==11
        F2,F3=list(map(parse,saved['generic_ideal_generators']))
        assert not pull(F2)and not pull(F3)
        normal_minor=None
        for i,j in it.combinations(range(5),2):
            poly=pull(q.add(q.mul(q.der(F2,i),q.der(F3,j)),q.scale(q.mul(q.der(F2,j),q.der(F3,i)),-1)))
            if poly:normal_minor=[i,j];break
        assert normal_minor is not None and Qs[0]
        # Verify the unique quadric and its10 constant-skew classes do not exhaust11.
        restrictions=[pull({ex+(0,0,0):F(1)})for ex in q.monomials(2)]
        monomials=sorted({e for p in restrictions for e in p},reverse=True)
        quadric_matrix=[[p.get(e,F())for p in restrictions]for e in monomials]
        assert q.rank(quadric_matrix)[0]==14
        quadric_columns=[]
        for ij in q.PAIRS:quadric_columns.append(q.vec(canonical(q.antisymmetric({ij:F2}))))
        assert q.rank(q.transpose(quadric_columns))[0]==10
        assert q.rank([row+extra for row,extra in zip(B,q.transpose(quadric_columns))])[0]==11
        clabels=[(r['ambient_partial'],tuple(r['exponents']))for r in saved['all330_conormal_firstpartial_labels']]
        cindex={(i,e+(0,)*(8-len(e))):at for at,(i,e)in enumerate(clabels)}
        C=[[F()for _ in range(11)]for _ in clabels]
        for j,P in enumerate(polys):
            numerator=q.add(q.mul(q.atom(1,0),bracket(P,F2,F3)),q.scale(q.mul(F3,bracket(P,F2,q.atom(1,0))),-3),q.scale(q.mul(F2,bracket(P,q.atom(1,0),F3)),-2))
            assert numerator==parse(saved['projective_conormal_bracket_numerators'][j])and not pull(numerator)
            for i in range(5):
                for ex,c in pull(q.der(numerator,i)).items():C[cindex[i,ex]][j]=c
        savedC=saved.get('complete_abelian_matrix',saved.get('complete330_by_projective_source_abelian_matrix'))
        assert C==matrix(savedC)
        needed=matrix(saved['complete_projective_by_abelian_source_matrix'])
        final=matrix(saved['complete150_by_abelian_source_matrix']);d=len(needed[0])
        assert kills(C,needed)and q.rank(needed)[0]==d
        ccert=minor(C,saved['abelian_rank_minor']);assert ccert['rank']+d==11
        assert q.matmul(B,needed)==final and q.rank(final)[0]==d
        P={}
        # Literal full parameter tensor; tau0,tau1,tau2 are z0,z1,z2.
        for j,(ij,ex)in enumerate(q.COLS):
            for i in range(d):
                if final[j][i]:P[ij]=q.add(P.get(ij,{}),{ex+tuple(int(k==i)for k in range(3)):final[j][i]})
        P=q.antisymmetric(P)
        assert all(P.get(ij,{})==parse(saved['full_reduced_source'].get(','.join(map(str,ij)),'0'))for ij in q.PAIRS)
        actualJ=q.jacobi(P);J=[coeff(actualJ[quad],ex)for quad,ex in q.ROWS]
        assert J==list(map(parse,saved['all350_Jacobi']))
        parameter_mons=[]
        for i,j in it.combinations_with_replacement(range(d),2):
            e=[0]*8;e[5+i]+=1;e[5+j]+=1;parameter_mons.append(tuple(e))
        JJ=[[p.get(e,F())for e in parameter_mons]for p in J]
        assert JJ==[[F(r.get(str(j),'0'))for j in range(len(parameter_mons))]for r in saved['all350_parameter_quadric_rows']]
        assert q.rank(JJ)[0]==saved['Jacobi_row_rank']
        D=quintic(P);assert D==parse(saved['actual_D'])
        DD=[[coeff(D,ex).get(e,F())for e in parameter_mons]for ex in q.monomials(5)]
        assert DD==[[F(r.get(str(j),'0'))for j in range(len(parameter_mons))]for r in saved['all126_actual_D_parameter_quadric_rows']]
        quotient=poly_division(D,F2)
        assert F2==parse(factor['unique_quadric_factor'])and quotient==parse(factor['exact_cubic_quotient_of_actual_D'])
        square_records=[]
        for cert in factor['S13_square_parameter_functionals']:
            identity=q.add(*(q.scale(J[int(k)],F(c))for k,c in cert['original350_Jacobi_functional'].items()))
            assert identity==parse(cert['squared_parameter']+'**2');square_records.append(cert['squared_parameter'])
        if square_records:assert not {ex:c for ex,c in D.items()if not any(ex[5+i]for i in range(1,d))}
        records.append(dict(model=model,full_zero_matrix_shape=[len(M),150],zero_minor=mcert,
            complete_raw_kernel_dimension=35,full_Euler_dimension=24,complete_canonical_projective_source_dimension=11,
            unique_quadric_space_dimension=1,quadric_constant_skew_source_dimension=10,extra_Euler_quotient_class_retained=True,
            generic_normal_basis_jacobian_minor_indices=normal_minor,
            full_projective_dehomogenized_conormal_numerator_and_all_firstpartials_verified=True,
            complete_conormal_matrix_shape=[len(C),11],conormal_minor=ccert,whole_necessary_source_dimension=d,
            all350_Jacobi_and126_actual_D_coefficients_reconstructed=True,
            exact_whole_source_D_quadric_factor=True,independent_cubic_quotient=qformat(quotient),
            S13_squared_parameter_functionals_verified=square_records,irreducible_quintic_branch_empty=True))
    result=dict(task_id='G',all_passed=True,no_CAS_import=True,no_producer_import=True,
        all_four_model_full_source_reviews=records,all_four_nonACM_projection_orbits_excluded=True,
        center_orbit_coverage='Separate mathematical proof audits both S22 rank2 pencils and S13 nonpure cubic centers, all actual normalization automorphisms, with excluded centers on the surface removed and the other centers routed by primary Brodmann-Park6.2 to ACM/CI22.',
        ACM_centers_routed_to_accepted_all_CI22_G_bound=True,S04_cone_parameter_family_explicitly_retained=True,
        input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in paths},
        mathematical_note_snapshot=paths[1].read_text(),
        reviewer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),reviewer_bytes=Path(__file__).read_text(),
        seconds=time.monotonic()-started,maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k]for k in ['all_passed','seconds','maximum_RSS_KiB']}))

def qformat(poly):return [[list(e),str(c)]for e,c in sorted(poly.items())]
if __name__=='__main__':main()
