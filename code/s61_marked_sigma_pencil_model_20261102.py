"""A denominator-free relative model for all finite h on the pencil c=2a."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import time
import sympy as s
from sympy.polys.matrices.normalforms import smith_normal_decomp
from s61_model import S61
from s61_flag_transverse_obstruction_20261101 import specialize

ROOT=Path(__file__).resolve().parents[1]
STEM="S61_MARKED_SIGMA_PENCIL_MODEL_20261102"
NOTE=f"research/notes/{STEM}.md"
NATIVE="code/s61_marked_sigma_pencil_model_20261102.cpp"


def prepare(folder):
    m=S61();one=[i for i,w in enumerate(m.reduced["parameter_weights"])if w==1]
    h,r=s.symbols("h rho");values=dict.fromkeys(one,s.S.Zero);values.update({4:r,11:1,25:1,26:h,50:2})
    _,forms=specialize(m,values,True)
    chart=[i for i in range(45)if i!=25]
    get_linear=lambda f:s.Matrix([[-eq.get(tuple(sorted((25,j))),0)for j in chart]+[-eq.get((25,),0)] for eq in f]).subs(r,0)
    J=get_linear(forms)
    ref=json.loads((ROOT/"reports/S61_MARKED_SIGMA_RESIDUAL_FINITE_20261102.json").read_text())
    w=ref["frame"]["q_weights"];cols=[i for i in range(45)if i not in [chart.index(j)for j in [0,19,27,32]]]
    domain=s.QQ.poly_ring(h);new=forms[:];blocks=[];nw=[];selected=[]
    T=s.zeros(45,13);U=s.zeros(45,28)
    for weight in range(1,5):
        cc=[j for j in cols if w[j]==weight];rr=[i for i in range(135)if any(J[i,j] for j in cc)]
        mat=J[rr,cc].to_DM().convert_to(domain)
        diagonal,left,right=smith_normal_decomp(mat);assert left*mat*right==diagonal
        assert domain.is_unit(left.det()) and domain.is_unit(right.det())
        dm,sm,tm=[z.to_Matrix() for z in [diagonal,left,right]]
        diag=[dm[i,i]for i in range(min(dm.shape))if dm[i,i]]
        assert all(z.is_Rational for z in diag)
        k=len(diag)
        for i,row in enumerate(rr):
            eq={}
            for j,oldrow in enumerate(rr):
                if sm[i,j]:
                    for mon,c in forms[oldrow].items():eq[mon]=eq.get(mon,0)+sm[i,j]*c
            new[row]={mon:coeff for mon,v in eq.items()if(coeff:=s.expand(v))}
        for j in range(len(cc)-k):
            for i,z in enumerate(cc):T[z,len(nw)+j]=tm[i,k+j]
        for j in range(k):
            for i,z in enumerate(cc):U[z,len(selected)+j]=tm[i,j]/diag[j]
        nw.extend([weight]*(len(cc)-k));selected.extend(rr[:k])
        blocks.append(dict(weight=weight,rows=rr,columns=cc,invariants=list(map(str,diag)),
                           row_transform=[[str(z)for z in r0]for r0 in sm.tolist()],
                           column_transform=[[str(z)for z in r0]for r0 in tm.tolist()]))
    J=get_linear(new)
    assert (J*T).applyfunc(s.expand)==s.zeros(135,13) and (J[selected,:]*U).applyfunc(s.expand)==s.eye(28)
    assert [len(b["invariants"])for b in blocks]==[3,11,11,3]
    assert nw==[1]*7+[2]*4+[3]*2
    def encode(v):
        coeff=[(m[0],a)for m,a in s.Poly(v,h,domain=s.QQ).terms()if a]
        return " ".join([str(len(coeff)),*[f"{a} {k}"for k,a in coeff]])
    lines=["14 28"," ".join(map(str,w))," ".join(map(str,[*nw,1]))," ".join(map(str,selected))]
    for matrix in [J,T,U]:lines.extend(encode(v)for v in matrix)
    rw=[]
    for row in new:
        terms=[];ws=set()
        for mon,c in row.items():
            ids=[chart.index(j)for j in mon if j!=25]+[44]*(2-len(mon));c*=(-1)**mon.count(25)
            for (power,),coeff in s.Poly(c,r).terms():
                if not coeff:continue
                ws.add(power+sum(w[j]for j in ids));assert ids
                if len(ids)==1 and not power:continue
                terms.append(f"{ids[0]} {ids[1]if len(ids)==2 else -1} {power} {encode(coeff)}")
        assert len(ws)<=1
        rw.append(next(iter(ws),-1));lines.extend([str(len(terms)),*terms])
    lines.append(" ".join(map(str,rw)))
    (folder/"input").write_text("\n".join(lines)+"\n")
    return blocks,nw,selected


def verify(export=None):
    start=time.monotonic();paths=[NOTE,NATIVE,str(Path(__file__).relative_to(ROOT))]
    for stem in ["S61_MARKED_SIGMA_RESIDUAL_FINITE_20261102","S61_MARKED_SIGMA_PROJECTIVE_CORE_20261102",
                 "S61_MARKED_SIGMA_PURE_ENDPOINT_REMOVAL_20261102","S61_MARKED_SIGMA_ZERO_H_REMOVAL_20261102"]:
        p=f"reports/{stem}.json";d=json.loads((ROOT/p).read_text());assert d["all_passed"]
        for f,v in d["input_sha256"].items():assert hashlib.sha256((ROOT/f).read_bytes()).hexdigest()==v,f
        paths.extend([p,d["proof_note"],*d["input_sha256"]])
    hashes={f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in dict.fromkeys(paths)}
    with tempfile.TemporaryDirectory(prefix="marked-square-pencil-")as tmp:
        wd=Path(tmp);blocks,nw,sel=prepare(wd);exe=wd/"recursion"
        subprocess.run(["c++","-O3","-std=c++17","-I"+str(ROOT/"work/flint-headers/flint/3.5.0/include"),
                        "-I"+str(ROOT/"work/flint-headers/mpfr/4.2.2/include"),"-I/opt/homebrew/include",
                        str(ROOT/NATIVE),str(ROOT/"work/m2-runtime/lib/libflint.24.0.dylib"),
                        "-L/opt/homebrew/lib","-lgmp","-o",str(exe)],check=True)
        env=os.environ|{"DYLD_LIBRARY_PATH":str(ROOT/"work/m2-runtime/lib")+":/opt/homebrew/lib"}
        run=subprocess.run([str(exe),str(wd/"input"),str(wd/"rows")],check=True,
                           capture_output=True,text=True,env=env,timeout=100)
        total,observable=map(int,run.stdout.split())
        rows_hash=hashlib.sha256((wd/"rows").read_bytes()).hexdigest()
        if export:
            assert not export.exists()
            shutil.copyfile(wd/"rows",export)
    return dict(id=STEM,all_passed=True,field="QQ[h]",finite_pencil="a=1,c=2,h arbitrary",
                normal_weights=nw,positive_parameter="rho=u5, weight one",Smith_blocks=blocks,
                selected_rows=sel,all135_polynomial_terms=total,homogenizer_polynomial_terms=observable,
                row_output_sha256=rows_hash,no_parameter_factor_inverted=True,full_pencil_removal_claimed=False,
                proof_note=NOTE,input_sha256=hashes,seconds=time.monotonic()-start)

if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("--output",type=Path,required=True)
    p.add_argument("--rows",type=Path);a=p.parse_args();r=verify(a.rows)
    with a.output.open("x")as f:json.dump(r,f,indent=2);f.write("\n")
    print(f"PASS {r['seconds']:.3f}s, {r['all135_polynomial_terms']} terms, {len(r['input_sha256'])} hashes")
