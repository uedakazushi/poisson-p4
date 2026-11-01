#!/usr/bin/env python3
"""Bounded initial-ideal and proof audit of four double-root torsion divisors."""
import argparse,hashlib,json,time
from pathlib import Path
from itertools import combinations_with_replacement
import sympy as s
from s61_minimal_junction_review_20261101 import standards
ROOT=Path(__file__).resolve().parents[1]
M=lambda a:s.Matrix(a).applyfunc(s.sympify)
def run():
    start=time.time();path=ROOT/'reports/S61_NONBINARY_DOUBLE_ROOT_TORSION_DIVISORS_20261101.json';src=json.loads(path.read_text());assert src['all_passed'];out=[]
    for name,rec in src['cases'].items():
        A2=M(rec['A2']);A3=M(rec['A3']);n=len(rec['local']['tangent_basis'][0]);mons={d:list(combinations_with_replacement(range(n),d)) for d in [1,2,3]};allmons=mons[2]+mons[3];ix={mm:i for i,mm in enumerate(allmons)}
        rr,pv=A2.row_join(A3).to_DM().rref();base=rr.to_Matrix()[:len(pv),:];multiples=[]
        for row in base.tolist():
            for j in range(n):
                val=[s.S.Zero]*len(allmons)
                for k,mm in enumerate(mons[2]):
                    if row[k]:val[ix[tuple(sorted(mm+(j,)))]]=row[k]
                if any(val):multiples.append(val)
        rr,pv=base.col_join(s.Matrix(multiples)).to_DM().rref();rr=rr.to_Matrix();v=s.symbols('v0:'+str(n));forms=[];ranks={}
        for row,pc in enumerate(pv):
            d=len(allmons[pc]);ranks[d]=ranks.get(d,0)+1;forms.append(sum(rr[row,j]*s.prod(v[i] for i in mm) for j,mm in enumerate(allmons) if len(mm)==d))
        gb0=s.groebner([f.subs(v[0],0) for f in forms],*reversed(v[1:]),order='grevlex');st0,hf0=standards(gb0);gb1=s.groebner([f.subs(v[0],1) for f in forms],*reversed(v[1:]),order='grevlex');st1,hf1=standards(gb1);assert len(st0)==4 and len(st1)==3 and len(hf0)<=4
        sec=[mm for d in [1,2,3] for mm in mons[d] if 0 not in mm];ix={mm:i for i,mm in enumerate(sec)};rows=[]
        for row in range(135):
            val=[s.S.Zero]*len(sec)
            for d,A in [(2,A2),(3,A3)]:
                for k,mm in enumerate(mons[d]):
                    if 0 not in mm:val[ix[mm]]=A[row,k]
            rows.append(val)
            for j in range(1,n):
                val=[s.S.Zero]*len(sec)
                for k,mm in enumerate(mons[2]):
                    if 0 not in mm:val[ix[tuple(sorted(mm+(j,)))]]=A2[row,k]
                rows.append(val)
        rank=s.Matrix(rows).to_DM().rank();assert 1+len(sec)-rank==4
        tangent=M(src['geometry'][name]['main_tangent']);assert tangent[0]!=0
        assert all(s.expand(f.subs(dict(zip(v,tangent))))==0 for f in forms)
        RB=M(src['geometry'][name]['restricted_beta_differential']);assert RB.rank()==16
        out.append(dict(label=name,truncated_initial_ranks=ranks,zero_section_standard_degrees=hf0,localized_standard_degrees=hf1,transverse_length=4,generic_rank_upper=3,actual_curve_tangent_nonzero_on_parameter=True,restricted_beta_rank=16))
    paths=[Path(__file__),path,ROOT/'code/s61_nonbinary_double_root_torsion_divisors_20261101.py',ROOT/'research/notes/S61_NONBINARY_DOUBLE_ROOT_TORSION_DIVISORS_20261101.md',ROOT/'code/s61_local_cubic_certificate_20261101.py',ROOT/'reports/S61_NONBINARY_DOUBLE_ROOT_REVIEW_20261101.json',ROOT/'reports/S61_NONBINARY_DOUBLE_ROOT_ENDPOINT_REVIEW_20261101.json',ROOT/'code/s61_minimal_junction_review_20261101.py']
    return dict(id='S61_NONBINARY_DOUBLE_ROOT_TORSION_DIVISORS_REVIEW_20261101',all_passed=True,cases=out,accepted_generic_rank3_transverse4=True,accepted_t_nilpotent=True,accepted_whole_generic_beta_strata_removed=True,scope='Independent initial/mixed relation audit from saved A2/A3 plus conceptual review; not a separate raw135 derivative or Sp4-action reconstruction.',seconds=time.time()-start,input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);a=p.parse_args();r=run()
    with a.output.open('x') as f:json.dump(r,f,indent=2);f.write('\n')
    print(json.dumps({'all_passed':True,'seconds':r['seconds']}))
