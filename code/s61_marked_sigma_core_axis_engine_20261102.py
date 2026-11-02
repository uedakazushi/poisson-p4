"""Build the exact positive-weight axis equations for a specified homogeneous core."""
from collections import defaultdict

import sympy as s
from sympy.polys.rings import ring

from s61_model import S61
from s61_flag_transverse_obstruction_20261101 import specialize


def build_axis(core, orbit_cuts, weights):
    m=S61()
    one=[i for i,w in enumerate(m.reduced["parameter_weights"]) if w==1]
    beta={j:value for j,value in zip(one,[0,*core])}
    r=s.Symbol("rho")
    beta[4]=r
    _,forms=specialize(m,beta,True)
    chart=[i for i in range(45) if i!=25]
    cut=[chart.index(i) for i in [*orbit_cuts,32]]
    columns=[i for i in range(45) if i not in cut]
    linear=s.Matrix([[-f.get(tuple(sorted((25,j))),0) for j in chart]+[-f.get((25,),0)] for f in forms])
    tangent=linear.subs(r,0)
    rr,pivot=tangent[:,columns].to_DM().to_field().rref()
    rr=rr.to_Matrix()
    pc=[columns[j] for j in pivot]
    free=[j for j in range(len(columns)) if j not in pivot]
    rows=list(tangent[:,pc].T.to_DM().rref()[1])
    inverse=tangent[rows,pc].inv(method="DM")
    names=[f"v{i}" for i in range(len(free))]+["rho"]
    ws=[weights[columns[j]] for j in free]+[1]
    assert min(ws)>0
    R,*v=ring(",".join(names),s.QQ)
    rho=v[-1]
    pols=[R.zero for i in range(45)]
    for k,j in enumerate(free):
        pols[columns[j]]=v[k]
    for k,j in enumerate(pc):
        pols[j]=-sum((s.QQ.convert(rr[k,z])*v[h] for h,z in enumerate(free) if rr[k,z]),R.zero)

    # Separate the constant Jacobian from every positive-weight product,
    # including uniformizer powers multiplying a single spatial coordinate.
    higher=[]
    rw=[]
    for i,f in enumerate(forms):
        terms=defaultdict(lambda:R.zero)
        monomials=set()
        for mon,c in f.items():
            c*=(-1)**mon.count(25)
            indices=tuple([chart.index(j) for j in mon if j!=25]+[44]*(2-len(mon)))
            if not indices:
                assert c==0
                continue
            for (k,),a in s.Poly(c,r).terms():
                if not a: continue
                monomials.add(k+sum(weights[j] for j in indices))
                if len(indices)>1 or k:
                    terms[indices]+=s.QQ.convert(a)*rho**k
        assert len(monomials)<=1
        rw.append(next(iter(monomials),None))
        higher.append(terms)

    def residual(i):
        return sum((c*pols[ids[0]] if len(ids)==1 else c*(pols[ids[0]]*pols[ids[1]])
                    for ids,c in higher[i].items()),R.zero)

    for d in sorted({weights[i] for i in pc}):
        chosen=[i for i in rows if rw[i]==d]
        raw={i:residual(i) for i in chosen}
        for k,j in enumerate(pc):
            if weights[j]!=d: continue
            assert all(rw[i]==d for h,i in enumerate(rows) if inverse[k,h])
            pols[j]-=sum((s.QQ.convert(inverse[k,h])*raw[i] for h,i in enumerate(rows)
                          if inverse[k,h]),R.zero)
    for j,p in enumerate(pols):
        assert all(sum(a*b for a,b in zip(mon,ws))==weights[j] for mon in p)
    output=[sum((s.QQ.convert(tangent[i,j])*pols[j] for j in columns if tangent[i,j]),R.zero)
            +residual(i) for i in range(135)]
    assert not any(output[i] for i in rows)
    assert all(sum(a*b for a,b in zip(mon,ws))==rw[i]
               for i,p in enumerate(output) for mon in p)
    assert pols[44]==v[-2]
    return names,ws,output,dict(Jacobian_rank=len(pc),selected_rows=rows,selected_columns=pc,
                                free_chart_columns=[columns[i] for i in free],uniformizer=names[-2],
                                normal_weights_with_rho=ws,full135_residual_terms=sum(map(len,output)))
