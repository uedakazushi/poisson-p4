#!/usr/bin/env python3
"""Separate-author full F1 strict-Casimir source, common vertex factor.
Different sparse Laurent polynomial implementation from root polarizations.
"""
from fractions import Fraction as F
from itertools import combinations,combinations_with_replacement
from pathlib import Path
import json,hashlib,time,resource,signal
from classification_minimal_quartic_review import LaurentRing
import classification_degenerate_quadric_q4_nil31_I_review as q
ROOT=Path(__file__).resolve().parents[1]
JOB='research/classification/jobs/plane_triple_root_F1_separate_review.json'
OUT='reports/CLASSIFICATION_PLANE_TRIPLE_ROOT_F1_SEPARATE_REVIEW_20261107.json'

def main():
    signal.alarm(60);t=time.monotonic();job=json.loads((ROOT/JOB).read_text());assert all(hashlib.sha256((ROOT/r['path']).read_bytes()).hexdigest()==r['sha256']for r in job['exact_inputs_and_hashes'])
    raw=json.loads((ROOT/'reports/CLASSIFICATION_PLANE_M4_TRIPLE_ROOT_QUINTICS_20261107.json').read_text());row=raw['cases'][1];assert raw['coordinates']==['A','B','u','v','w']and row['epsilon']==1
    labels=raw['full150_column_labels'];assert labels==[dict(pair=list(ij),monomial=list(e))for ij,e in q.COLS]
    R=LaurentRing(('c0','c1','c2','c3'),('A','B','u','v','w'));p=R.parse;n=R.n
    Ff=p('u**4*A+u**3*v*B+u**2*v**3+v**5');assert Ff=={(0,)*n+tuple(r['exponent']):F(r['coefficient'])for r in row['F']}
    K=[[F(c)for c in r]for r in row['full150_kernel']];assert len(K)==4 and q.rank(K)[0]==4
    # Reconstruct the entire original strict-Casimir/divergence linear map.
    columns=[]
    for label in labels:
        i,j=label['pair'];mon={(0,)*n+tuple(label['monomial']):F(1)};cas=[{}for k in range(5)];cas[i]=R.mul(mon,R.derivative(Ff,j));cas[j]=R.scale(R.mul(mon,R.derivative(Ff,i)),-1)
        div=[{}for k in range(5)];div[i]=R.scale(R.derivative(mon,j),-1);div[j]=R.derivative(mon,i)
        columns.append([P.get((0,)*n+e,F(0))for P in cas for e in q.monomials(6)]+[P.get((0,)*n+e,F(0))for P in div for e in q.monomials(1)])
    M=q.transpose(columns);saved=[[F(v)for v in r]for r in row['full1075x150_strict_Casimir_matrix']];assert M==saved
    C=row['rank_certificate'];rows,cols=C['rows'],C['columns'];assert C['rank']==146
    d=q.determinant([[M[i][j]for j in cols]for i in rows]);assert d and d==F(C['minor_determinant'])
    assert all(not sum(a*b for a,b in zip(r,v))for r in M for v in K)
    # Complete source, all four independent parameters in one polynomial tensor.
    P={(i,j):{}for i in range(5)for j in range(5)}
    for k,v in enumerate(K):
        for label,c in zip(labels,v):
            if c:
                i,j=label['pair'];e=(0,)*n+tuple(label['monomial']);term=R.mul(p('c'+str(k)),{e:c});P[i,j]=R.add(P[i,j],term);P[j,i]=R.scale(P[i,j],-1)
    expected=[['-u*v','u**2','0','0'],['-v**2','u*v','0','0'],['-4*A*u-5*B*v','2*B*u-2*v**2','u**2','0'],['-5*A*v-B*u-3*v**2','15*A*u+15*B*v/2','-5*u*v/2','u**2']]
    assert all(P[i,4]==R.add(*(R.mul(p('c'+str(k)),p(expected[k][i]))for k in range(4)))for i in range(4))
    assert all(not P[i,j]for i,j in combinations(range(4),2))and all(not R.derivative(P[i,4],4)for i in range(4))
    strict={tr:R.add(*(R.mul(P[a,l],R.derivative(P[b,c],l))for a,b,c in((tr[0],tr[1],tr[2]),(tr[1],tr[2],tr[0]),(tr[2],tr[0],tr[1]))for l in range(5)))for tr in combinations(range(5),3)}
    assert not any(strict.values());x=[R.variable(s)for s in('A','B','u','v','w')]
    projective={quad:R.add(*(R.scale(R.mul(x[quad[i]],strict[quad[:i]+quad[i+1:]]),(-1)**i)for i in range(4)))for quad in combinations(range(5),4)}
    pf=[]
    for missing in range(5):
        i,j,k,l=[a for a in range(5)if a!=missing];pf.append(R.add(R.mul(P[i,j],P[k,l]),R.scale(R.mul(P[i,k],P[j,l]),-1),R.mul(P[i,l],P[j,k])))
    assert not any(pf)
    pairs=list(combinations_with_replacement(range(4),2));parmon=[tuple(int(i==k)+int(j==k)for k in range(4))for i,j in pairs]
    Jrows=[[poly.get(par+e,F(0))for par in parmon]for poly in projective.values()for e in q.monomials(4)]
    Pfrows=[[poly.get(par+e,F(0))for par in parmon]for poly in pf for e in q.monomials(4)]
    root=json.loads((ROOT/'reports/CLASSIFICATION_PLANE_TRIPLE_ROOT_F1_SUPPORT_REVIEW_20261107.json').read_text());assert Jrows==[[F(v)for v in r]for r in root['full350_Jacobi_quadratic_matrix']]
    assert Pfrows==[[F(v)for v in r]for r in root['full5_by70_Pfaffian_quadratic_matrix']]
    assert root['all_raw_Pfaffians_identically_zero']and root['all_actualD_identically_zero']
    out=dict(task_id='E_reduced_nonnormal',all_passed=True,whole1075x150_original_map_reconstructed=True,rank=146,independent146minor=str(d),complete_kernel_dimension=4,all_four_parameters_and_all_zero_complements_retained=True,exact_common_vertex_source_basis=expected,entire_source_form='Pi=Z(A,B,u,v) wedge partial_w; partial_w Z=0',all10_strict_Jacobi_zero=True,all350_by10_projective_polarizations_zero=True,all5_by70_by10_raw_Pfaffian_coefficients_zero=True,root_full_polarization_matrices_independently_equal=True,all_actualD_zero=True,full_F1_rank4_source_empty=True,whole31_rank4_empty_bridge='All actualD normal forms F0/F1 are obtained by invertible linear shears with kappa!=0 from irreducibility, residualalpha0/1 both retained. Canonical divergencezero strict lift is GL5 equivariant and has Casimir actualD. F0 complete source rank<=2 accepted; F1 fullL4 is decomposable here. No R reduced/log condition.',no_new_anchor_jet_or_GB=True,inputs=job['exact_inputs_and_hashes'],reviewer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),reviewer_source_bytes=Path(__file__).read_text(),seconds=time.monotonic()-t,maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    path=ROOT/OUT;assert not path.exists();path.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:out[k]for k in('all_passed','seconds','maximum_RSS_KiB','full_F1_rank4_source_empty')}))
if __name__=='__main__':main()
