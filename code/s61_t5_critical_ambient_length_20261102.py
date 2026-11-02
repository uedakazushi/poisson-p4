"""Full 18-parameter incidence transverse and actual sixth Macaulay space."""
import itertools
import json
from functools import lru_cache
from pathlib import Path

import sympy as sym
from sympy import QQ
from sympy.polys.rings import ring

from s61_flag_t5_core_geometry_20261102 import fixed_representative
from s61_model import S61,X,J
from weight1_structure import sp_basis

ROOT=Path(__file__).resolve().parents[1]
EXTRA="reports/S61_UNIVERSAL_WEIGHT6_20261031.json"


def ambient_certificate():
    model=S61()
    tensor,H,beta=fixed_representative(model)
    tensor=tensor.subs(sym.Symbol("e"),1)
    ambient=X[0]**2*J
    ambient[1:,1:]+=tensor
    one=[i for i,w in enumerate(model.reduced["parameter_weights"]) if w==1]
    two=[i for i,w in enumerate(model.reduced["parameter_weights"]) if w==2]
    q=sym.Matrix(model.normalize(ambient)["parameters"])[two,:]
    assert q[21]==2 and q[25]==-1
    spatial=sym.Matrix(X[1:])
    orbit_cols=[]
    for g in sp_basis():
        velocity=g*spatial
        variation=g*tensor+tensor*g.T-sum((velocity[i]*tensor.diff(spatial[i]) for i in range(4)),sym.zeros(4))
        changed=X[0]**2*J
        changed[1:,1:]+=variation
        col=sym.Matrix(model.normalize(changed)["parameters"])[two,:]
        orbit_cols.append(col+col[25]*q)
    orbit=sym.Matrix.hstack(*orbit_cols)
    beta_cut=list(H.T.to_DM().rref()[1])
    q_cut=list(orbit.T.to_DM().rref()[1])
    gcols=list(orbit[q_cut,:].to_DM().rref()[1])
    assert len(beta_cut)==12 and len(q_cut)==len(gcols)==6
    incidence_det=H[beta_cut,:].det()*orbit[q_cut,gcols].det()
    assert incidence_det
    cuts={two[25],*[one[i] for i in beta_cut],*[two[i] for i in q_cut]}
    variables=[i for i in range(62) if i not in cuts]
    loc={j:i for i,j in enumerate(variables)}
    point=dict(zip(one,beta))|dict(zip(two,q))|{61:sym.S.Zero}
    assert len(variables)==43
    source=model.reduced["equations"]+[json.loads((ROOT/EXTRA).read_text())["universal_new_equation"]]
    linear=sym.zeros(135,43)
    higher=[]
    for row,polynomial in enumerate(source):
        expanded={}
        for mon,numerator,denominator in polynomial:
            qdegree=sum(i in two for i in mon)
            terms={():sym.Rational(numerator,denominator)}
            for j in [*mon,*([61]*(2-qdegree))]:
                next_terms={}
                for ids,c in terms.items():
                    if point[j]:
                        next_terms[ids]=next_terms.get(ids,0)+point[j]*c
                    if j in loc:
                        key=tuple(sorted((*ids,loc[j])))
                        next_terms[key]=next_terms.get(key,0)+c
                terms=next_terms
            for key,c in terms.items():
                expanded[key]=expanded.get(key,0)+c
        expanded={key:c for key,c in expanded.items() if c}
        assert () not in expanded
        for key,c in expanded.items():
            if len(key)==1:
                linear[row,key[0]]=c
        higher.append({key:c for key,c in expanded.items() if len(key)>=2})
    kernel=sym.Matrix.hstack(*linear.nullspace())
    piv=list(linear.to_DM().rref()[1])
    rows=list(linear[:,piv].T.to_DM().rref()[1])
    assert len(piv)==40 and kernel.cols==3
    block=linear[rows,piv]
    inverse=block.inv()
    right=sym.zeros(43,40)
    for i,k in enumerate(piv):
        right[k,:]=inverse[i,:]
    assert linear[rows,:]*right==sym.eye(40)
    algebra,*v=ring("v0,v1,v2",QQ)
    scalar=lambda z:QQ.from_sympy(sym.Rational(z))
    def apply(mat,vec):
        return [sum((scalar(mat[i,j])*vec[j] for j in range(mat.cols) if mat[i,j] and vec[j]),algebra.zero) for i in range(mat.rows)]
    lifts={1:[sum((scalar(kernel[i,j])*v[j] for j in range(3)),algebra.zero) for i in range(43)]}
    residuals=[algebra.zero for _ in range(135)]
    for degree in range(2,7):
        @lru_cache(None)
        def product(ids,n):
            if not ids:
                return algebra.one if n==0 else algebra.zero
            return sum((lifts[k][ids[0]]*product(ids[1:],n-k) for k in range(1,min(n,degree-1)+1) if lifts[k][ids[0]]),algebra.zero)
        raw=[sum((scalar(c)*product(ids,degree) for ids,c in row.items() if len(ids)<=degree),algebra.zero) for row in higher]
        lifts[degree]=[-z for z in apply(right,[raw[i] for i in rows])]
        terms=[z+y for z,y in zip(raw,apply(linear,lifts[degree]))]
        assert not any(terms[i] for i in rows)
        residuals=[z+y for z,y in zip(residuals,terms)]
    # Start at lowest total degree and retain every later coefficient of each row.
    monomials=[(i,j,d-i-j) for d in range(7) for i in range(d+1) for j in range(d-i+1)]
    index={exp:i for i,exp in enumerate(monomials)}
    basis={}
    selected=[]
    for row_index,polynomial in enumerate(residuals):
        if not polynomial:
            continue
        for shift in monomials:
            vec={index[tuple(x+y for x,y in zip(exp,shift))]:c for exp,c in polynomial.items() if sum(exp)+sum(shift)<=6}
            while vec:
                j=min(vec)
                pivot=vec[j]
                if j not in basis:
                    basis[j]={k:c/pivot for k,c in vec.items()}
                    selected.append([row_index,list(shift)])
                    break
                for k,c in basis[j].items():
                    val=vec.get(k,QQ.zero)-pivot*c
                    if val:
                        vec[k]=val
                    else:
                        vec.pop(k,None)
    hilbert=[sum(1 for j,a in enumerate(monomials) if sum(a)==d and j not in basis) for d in range(7)]
    assert hilbert==[1,3,5,6,5,3,0] and len(basis)==61 and len(monomials)==84
    return dict(H12_beta=list(map(str,beta)),ambient_q_support={str(i):str(z) for i,z in enumerate(q) if z},
                actual_incidence_dimension=18,cut_beta_indices=beta_cut,cut_q_indices=q_cut,orbit_columns=gcols,
                incidence_minor=str(incidence_det),transverse_raw_variables=variables,rank=40,
                jacobian_minor_rows=rows,jacobian_minor_columns=piv,jacobian_minor=str(block.det(method="domain-ge")),
                normal_basis=[[str(z) for z in a] for a in kernel.tolist()],truncated_monomials=monomials,
                independent_Macaulay_rows=selected,pivot_columns=list(basis),Macaulay_rank=61,
                actual_transverse_Hilbert=hilbert,actual_transverse_length=sum(hilbert),
                ambient_free_rank_over_18_parameters=23,homogenizer_power_vanishing=23)
