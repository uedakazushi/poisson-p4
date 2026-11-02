"""Mixed-characteristic transverse bounds at all nonzero rays on the marked a=0 test."""
import itertools
from functools import lru_cache
import sympy as sym
from sympy.polys.rings import ring
from s61_binary_endpoint_full_transverse_20261102 import PRIME, ROOT, load_report
from s61_model import S61, X
from toric_closed_jet_coordinates_20261102 import coords, pairs, m2
from weight1_structure import sp_basis

def unit(z):
    z=sym.Rational(z)
    return int(z.p)%PRIME != 0 and int(z.q)%PRIME != 0

def red(z):
    z=sym.Rational(z)
    return int(z.p)*pow(int(z.q),-1,PRIME)%PRIME

def transverse(one,two,beta,q,variables,source):
    field=sym.GF(PRIME)
    loc={j:i for i,j in enumerate(variables)}
    point=dict(zip(one,beta))|dict(zip(two,q))|{61:sym.S.Zero}
    N=len(variables)

    linear=sym.zeros(135,N)
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
    n = len(variables)-len(piv)
    block=linear[rows,piv]
    jdet=block.det(method="domain-ge")
    assert unit(jdet)
    if not n:return [1], str(jdet), rows, piv
    assert kernel.cols==n
    inverse=block.inv()
    right=sym.zeros(N,len(piv))
    for i,k in enumerate(piv):
        right[k,:]=inverse[i,:]
    assert linear[rows,:]*right==sym.eye(len(piv))
    algebra,*v=ring(",".join(f"v{i}" for i in range(n)),field)
    scalar=lambda z:field(red(z))
    def apply(mat,vec):
        return [sum((scalar(mat[i,j])*vec[j] for j in range(mat.cols) if mat[i,j] and vec[j]),algebra.zero) for i in range(mat.rows)]
    lifts={1:[sum((scalar(kernel[i,j])*v[j] for j in range(n)),algebra.zero) for i in range(N)]}
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
    monomials=[a for d in range(7) for a in itertools.product(range(d+1),repeat=n) if sum(a)==d]
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
                    val=vec.get(k,field.zero)-pivot*c
                    if val:
                        vec[k]=val
                    else:
                        vec.pop(k,None)
    hilbert=[sum(1 for j,a in enumerate(monomials) if sum(a)==d and j not in basis) for d in range(7)]
    return hilbert,str(jdet),rows,piv

def incidence_and_transverse(model, data, source, entry):
    one=[i for i,w in enumerate(model.reduced['parameter_weights']) if w==1]
    two=[i for i,w in enumerate(model.reduced['parameter_weights']) if w==2]
    q=sym.zeros(45,1); beta=sym.zeros(16,1)
    for i,c in entry['q'].items():q[int(i)]=sym.Rational(c)
    for i,c in entry['beta'].items():beta[int(i)]=sym.Rational(c)
    fixed=21
    assert q[fixed]==1
    pos=dict(zip(one,beta))|dict(zip(two,q))
    highest=sym.zeros(135,16); quadratic=sym.zeros(16)
    for k,row in enumerate(source):
        for mon,n,d in row:
            if sum(i in two for i in mon)!=2:continue
            bmon=[one.index(i) for i in mon if i in one]
            coeff=sym.Rational(n,d)*sym.prod(pos[i] for i in mon if i in two)
            if len(bmon)==1:highest[k,bmon[0]]+=coeff
            if len(bmon)==2:quadratic[bmon[0],bmon[1]]+=coeff
    H=sym.Matrix.hstack(*highest.nullspace())
    assert highest*beta==sym.zeros(135,1)
    assert H.T*(quadratic+quadratic.T)*H==sym.zeros(H.cols)
    assert all(int(sym.denom(c))%PRIME for c in H)
    y=sym.Matrix(X[1:]);physical=data['E']*q
    Q=sym.zeros(4)
    for k,(i,j) in enumerate(pairs):
        Q[i,j]=sum(physical[k*10+h]*y[v]*y[w] for h,(v,w) in enumerate(m2))
    Q-=Q.T
    columns=[]
    for g in sp_basis():
        dQ=g*Q+Q*g.T-sum(((g*y)[i]*Q.diff(y[i]) for i in range(4)),sym.zeros(4))
        coefficients=sym.Matrix([sym.Poly(dQ[i,j],*y).coeff_monomial(y[k]*y[l])
                                 for i,j in pairs for k,l in m2])
        v=data['Pi']*coefficients
        columns.append(v-q*v[fixed])
    O=sym.Matrix.hstack(*columns)
    bc=list(H.T.to_DM().rref()[1]);qc=list(O.T.to_DM().rref()[1])
    gc=list(O[qc,:].to_DM().rref()[1])
    det=H[bc,:].det()*O[qc,gc].det()
    assert unit(det)
    frozen={two[fixed],*[two[j] for j in qc],*[one[j] for j in bc]}
    variables=[i for i in range(62) if i not in frozen]
    hf,jdet,rows,piv=transverse(one,two,beta,q,variables,source)
    length={'xzw':2,'T2':2,'double root':3,'smooth':1,'T5':23,'T6':2}[entry['type']]
    assert sum(hf)==length and (hf==[1] or hf[-1]==0)
    return dict(**entry,tensor=[str(sym.factor(Q[i,j])) for i,j in pairs],incidence_dimension=len(qc)+len(bc),
                orbit_dimension=len(qc),H_dimension=len(bc),incidence_minor=str(det),beta_cuts=bc,
                q_cuts=qc,orbit_columns=gc,Jacobian_minor=jdet,Jacobian_rows=rows,Jacobian_columns=piv,
                transverse_hilbert_mod_32003=hf,transverse_length=length)

def all_rays(points):
    model=S61()
    source=model.reduced['equations']+[load_report('S61_UNIVERSAL_WEIGHT6_20261031')['universal_new_equation']]
    data=coords()
    return [incidence_and_transverse(model,data,source,pt) for pt in points]
