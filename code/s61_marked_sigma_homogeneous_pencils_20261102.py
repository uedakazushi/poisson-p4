"""Classify all torus-invariant affine pencils in the marked flag direction."""
import argparse
import hashlib
import json
from pathlib import Path
import time

import sympy as s

from s61_model import S61
from s61_flag_transverse_obstruction_20261101 import specialize

ROOT=Path(__file__).resolve().parents[1]
STEM="S61_MARKED_SIGMA_HOMOGENEOUS_PENCILS_20261102"
NOTE=f"research/notes/{STEM}.md"


def verify():
    start=time.monotonic()
    prior_path="reports/S61_MARKED_SIGMA_PROJECTIVE_CORE_20261102.json"
    weighted_path="reports/S61_MARKED_SIGMA_RESIDUAL_FINITE_20261102.json"
    deps=[NOTE,str(Path(__file__).relative_to(ROOT)),prior_path,weighted_path]
    def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
    for path in [prior_path,weighted_path]:
        old=json.loads((ROOT/path).read_text());assert old["all_passed"]
        for p,digest in old["input_sha256"].items():assert sha(p)==digest,p
        deps.extend(old["input_sha256"])
    weights=json.loads((ROOT/weighted_path).read_text())["frame"]["q_weights"]
    chart=[j for j in range(45)if j!=25]
    cols=[chart[j]for j in range(44) if weights[j]==weights[44]]
    m=S61();one=[j for j,w in enumerate(m.reduced["parameter_weights"])if w==1]
    c,h,inv=s.symbols("c h inv")
    records=[]
    def case(a,cc,hh,expected_minor,expected,parameters=()):
        beta=dict.fromkeys(one,s.S.Zero)
        beta.update({11:1,25:a,26:hh,50:cc})
        _,F=specialize(m,beta,True)
        # Coordinate q25 equals the free affine line parameter.
        assert all(not row.get((25,25),0) for row in F)
        A=s.Matrix([[*[row.get(tuple(sorted((25,j))),0)for j in cols],-row.get((25,),0)] for row in F])
        rr,pc=A.to_DM().to_field().rref();rr=rr.to_Matrix()
        rows=list(A[:,list(pc)].T.to_DM().to_field().rref()[1])
        minor=s.factor(A[rows,list(pc)].det(method="domain-ge"))
        assert s.expand(minor-expected_minor)==0,(minor,expected_minor)
        record=dict(core=list(map(str,[a,cc,hh])),rows=rows,pivots=list(pc),minor=str(minor))
        records.append(record)
        if 12 in pc:
            assert expected=="inconsistent"
            record["augmented_pivot"]=True
            return
        fr=[j for j in range(12)if j not in pc]
        variables=s.symbols(f"z0:{len(fr)}")
        q=s.zeros(45,1)
        for j,z in zip(fr,variables):q[cols[j]]=z
        for i,j in enumerate(pc):q[cols[j]]=s.factor(rr[i,-1]-sum(rr[i,k]*q[cols[k]] for k in fr))
        assert (A[:,:12]*q[cols,:]-A[:,-1]).applyfunc(s.cancel)==s.zeros(135,1)
        residual=[s.factor(sum(v*s.prod(q[j]for j in mon)for mon,v in row.items()))for row in F]
        generators=[s.fraction(s.cancel(z))[0] for z in residual if z]
        field=s.QQ.algebraic_field(s.sqrt(-2)) if s.sympify(hh).has(s.I) else s.QQ
        G=s.groebner([inv*minor-1,*generators],inv,*variables,*parameters,
                     domain=field,order="lex")
        gb=[f.as_expr()for f in G.polys]
        record["constant_equation_Groebner_basis"]=list(map(str,gb))
        if expected=="empty":
            assert gb==[1]
        elif expected=="c3":
            assert gb==[inv+s.Rational(1,4608),variables[0],c-3]
            assert [s.factor(x.subs({variables[0]:0,c:3})) for x in q]==[
                -s.Rational(5,2) if j==39 else 0 for j in range(45)]
        elif expected=="plane":
            assert gb==[inv-s.Rational(1,288)]
            assert [s.factor(x)for x in q]==[
                {2:s.Rational(3,4),13:variables[0],15:s.Rational(1,2),44:s.Rational(1,4)}.get(j,0)
                for j in range(45)]
        else:raise AssertionError(expected)
    case(1,c,h,-576*h*(c-5),"empty",(c,h))
    case(1,c,0,1152*(c-5)*(c-1),"c3",(c,))
    case(1,1,0,18432,"inconsistent")
    case(1,5,h,144*h*(h*h+32),"inconsistent",(h,))
    case(1,5,0,-36864,"inconsistent")
    case(1,5,4*s.sqrt(-2),576*s.sqrt(-2),"empty")
    case(0,1,h,288*h,"empty",(h,))
    case(0,1,0,288,"plane")
    case(0,0,1,144,"inconsistent")
    return dict(id=STEM,all_passed=True,field="QQ, one specialization in QQ(sqrt(-2))",
                zero_weight_affine_coordinates=cols,exact_cases=records,
                projective_loci_with_homogeneous_affine_pencil=["[1:3:0]","[0:1:0]"],
                generic_beta_counterexample_claimed=False,other_affine_arcs_classified=False,
                proof_note=NOTE,input_sha256={d:sha(d) for d in dict.fromkeys(deps)},seconds=time.monotonic()-start)

if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument("--output",type=Path,required=True)
    out=parser.parse_args().output;result=verify()
    with out.open("x")as stream:json.dump(result,stream,indent=2);stream.write("\n")
    print(f"PASS {result['seconds']:.3f}s / {len(result['input_sha256'])} hashes")
