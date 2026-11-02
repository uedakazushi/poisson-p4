"""A projective two-parameter criterion for all remaining marked flags."""
import argparse
import hashlib
import json
from pathlib import Path
import time
import sympy as s
from s61_model import S61
from s61_flag_transverse_obstruction_20261101 import specialize

ROOT=Path(__file__).resolve().parents[1]
STEM="S61_MARKED_SIGMA_PROJECTIVE_CORE_20261102"
NOTE=f"research/notes/{STEM}.md"


def verify():
    start=time.monotonic()
    deps=[NOTE,str(Path(__file__).relative_to(ROOT)),"code/s61_model.py","code/s61_flag_transverse_obstruction_20261101.py",
          "reports/S61_UNIVERSAL_WEIGHT6_20261031.json","reports/Q51_UNIVERSAL_SECTION_20261101.json",
          "reports/Q52_UNIVERSAL_SECTION_20261101.json","reports/S61_FLAG_UNIFORM_OPEN_20261101.json",
          "reports/S61_FLAG_BOUNDARY_NORMAL_FORMS_20261101.json","reports/S61_MARKED_SIGMA_RESIDUAL_FINITE_20261102.json",
          "reports/S61_MARKED_SIGMA_INTERSECTION_REMOVAL_20261102.json"]
    for p in deps[4:]:
        d=json.loads((ROOT/p).read_text())
        assert d["all_passed"]
    flag=json.loads((ROOT/"reports/S61_FLAG_UNIFORM_OPEN_20261101.json").read_text())
    b=s.symbols("b0:15")
    D=s.Matrix(flag["frame_change"]["beta_derivative"]).applyfunc(s.sympify)
    E=D[:,3].jacobian(b);assert E**4==s.zeros(15)
    transform=lambda t:s.eye(15)+t*E+t*t/2*E**2+t**3/6*E**3
    marked=s.Matrix(b);marked[0]=0;marked[3]=4*b[12]/5
    t=3*b[12]/(5*b[1]);moved=(transform(t)*marked).applyfunc(s.factor)
    assert moved[0]==moved[3]==moved[12]==0 and moved[1]==b[1]
    assert (transform(-t)*moved-marked).applyfunc(s.factor)==s.zeros(15,1)
    core=[s.factor(moved[i]/b[1])for i in [5,14,6]]
    assert s.factor(core[0]-(50*b[1]*b[5]-9*b[12]**2)/(50*b[1]**2))==0
    assert s.factor(core[1]-(10*b[1]*b[14]-9*b[12]**2)/(10*b[1]**2))==0
    assert s.factor(core[2]-(125*b[1]**2*b[6]-75*b[1]*b[12]*(b[14]+b[5])+54*b[12]**3)/(125*b[1]**3))==0
    d1=[s.cancel(D[i,1]/b[i]) for i in range(15)]
    d2=[s.cancel(D[i,5]/b[i]) for i in range(15)]
    contraction=[int(1+d1[i]+d2[i]) for i in range(15)]
    nuisance=[2,4,7,8,9,10,11,13]
    assert [contraction[i]for i in nuisance]==[2,2,4,2,4,2,2,2]
    assert all(contraction[i]==0 for i in [1,5,14,6])
    assert [2+d1[i] for i in [1,5,14,6]]==[0,2,2,3]
    # Test whether translating either universal elliptic lift in the flag
    # direction gives an affine line, on the exact zero-nuisance core.
    m=S61();one=[i for i,v in enumerate(m.reduced["parameter_weights"]) if v==1]
    a,c,h=s.symbols("a c h")
    values=dict.fromkeys(one,s.S.Zero);values.update({11:1,25:a,50:c,26:h})
    _,rows=specialize(m,values,True)
    assert all(not row.get((25,25),0) for row in rows)
    subs={s.Symbol(f"u{i+1}"):val for i,val in values.items()}
    tests={}
    for n in ["Q51","Q52"]:
        data=json.loads((ROOT/f"reports/{n}_UNIVERSAL_SECTION_20261101.json").read_text())
        q=[s.sympify(f).subs(subs) for f in data["quadratic_section_q"]]
        linear=[]
        for row in rows:
            constant=s.S.Zero
            direction=s.sympify(row.get((25,),0))
            for mon,z in row.items():
                constant+=z*s.prod(q[i] for i in mon)
                if len(mon)==2:
                    i,j=mon
                    if i==25:direction+=z*q[j]
                    if j==25:direction+=z*q[i]
            assert s.expand(constant)==0
            linear.append(s.factor(direction))
        gb=s.groebner([f for f in linear if f],a,c,h,domain=s.QQ)
        assert [f.as_expr()for f in gb.polys]==[a,c,h]
        tests[n]=dict(linear_coefficients={str(j):str(f) for j,f in enumerate(linear) if f},ideal=["a","c","h"])
    return dict(id=STEM,all_passed=True,core_coordinates=list(map(str,core)),projective_weights=[2,2,3],
                nuisance_indices=nuisance,nuisance_weights=[contraction[i] for i in nuisance],
                zero_nuisance_elliptic_line_tests=tests,finite_core_subset_beta_image_dimension_bound=15,
                zero_core_beta_image_dimension_bound=14,complete_ample_curve_elimination_claimed=False,
                proof_note=NOTE,input_sha256={d:hashlib.sha256((ROOT/d).read_bytes()).hexdigest() for d in deps},
                seconds=time.monotonic()-start)

if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("--output",required=True,type=Path)
    out=p.parse_args().output;d=verify()
    with out.open("x")as f: json.dump(d,f,indent=2);f.write("\n")
    print(f"PASS {d['seconds']:.3f}s / {len(d['input_sha256'])} hashes")
