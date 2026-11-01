#!/usr/bin/env python3
"""Finite initial-ideal multiplicity bounds remove eight minimal junctions."""
import sys,json,time,argparse,hashlib
from pathlib import Path
from itertools import combinations_with_replacement
import sympy as s
from s61_model import S61
from s61_flag_transverse_obstruction_20261101 import specialize
ROOT=Path(__file__).resolve().parents[1]
def run():
    start_total=time.time();m=S61();source_path=ROOT/'reports/S61_MINIMAL_JUNCTION_BRANCHES_20261101.json';source=json.loads(source_path.read_text());out=[]
    for record in source['cases']:
        r=record['local'];order_bound=4 if record['kind']=='M1' and record['k']=='-1/3' else 3
        start=time.time();p=list(map(s.Rational,r['point61']));beta={i:p[i] for i,w in enumerate(m.reduced['parameter_weights']) if w==1};_,F=specialize(m,beta,True);cols=r['chart_columns'];fixed=next(i for i in range(45) if i not in cols);J=s.Matrix(r['jacobian']).applyfunc(s.Rational);T=s.Matrix(r['tangent_basis']).applyfunc(s.Rational);n=T.cols;rows=list(J.T.to_DM().rref()[1]);piv=list(J.to_DM().rref()[1]);U=s.zeros(45,len(piv));iv=J.extract(rows,piv).inv(method='DM')
        for i,c in enumerate(piv):U[c,:]=iv[i,:]
        jac_piv=list(piv);JU=J*U;mon2=list(combinations_with_replacement(range(n),2));mon3=list(combinations_with_replacement(range(n),3));i3={x:i for i,x in enumerate(mon3)};Q=s.zeros(135,len(mon2));terms=[]
        for row,f in enumerate(F):
         ts=[]
         for mon,c in f.items():
          if fixed in mon:continue
          i,j=[cols.index(k) for k in mon]+[44]*(2-len(mon));ts.append((i,j,c))
          for h,(a,b) in enumerate(mon2):Q[row,h]+=c*(T[i,a]*T[j,b]+(T[i,b]*T[j,a] if a!=b else 0))
         terms.append(ts)
        w2=-U*Q[rows,:];A2=Q-JU*Q[rows,:];Q3=s.zeros(135,len(mon3))
        for row,ts in enumerate(terms):
         for i,j,c in ts:
          for a in range(n):
           for h,mon in enumerate(mon2):Q3[row,i3[tuple(sorted((a,)+mon))]]+=c*(T[i,a]*w2[j,h]+T[j,a]*w2[i,h])
        if order_bound==4:
            A3=Q3-JU*Q3[rows,:];w3=-U*Q3[rows,:];mon4=list(combinations_with_replacement(range(n),4));i4={x:i for i,x in enumerate(mon4)};Q4=s.zeros(135,len(mon4))
            for row,ts in enumerate(terms):
             for i,j,c in ts:
              for h,ma in enumerate(mon2):
               for l,mb in enumerate(mon2):Q4[row,i4[tuple(sorted(ma+mb))]]+=c*w2[i,h]*w2[j,l]
              for a in range(n):
               for h,mon in enumerate(mon3):Q4[row,i4[tuple(sorted((a,)+mon))]]+=c*(T[i,a]*w3[j,h]+T[j,a]*w3[i,h])
            A4=Q4-JU*Q4[rows,:];mons=mon2+mon3+mon4;ix={mon:i for i,mon in enumerate(mons)};rr,rp=s.Matrix.hstack(A2,A3,A4).to_DM().rref();R=rr.to_Matrix()[:len(rp),:];mult=[]
            for row in R.tolist():
             for fac in [(i,) for i in range(n)]+mon2:
              z=[s.S.Zero]*len(mons)
              for j,mon in enumerate(mons):
               if len(mon)+len(fac)<=4 and row[j]:z[ix[tuple(sorted(mon+fac))]]+=row[j]
              if any(z):mult.append(z)
            full=R.col_join(s.Matrix(mult));print('macaulay',full.shape,flush=True);rr,piv=full.to_DM().rref();rr=rr.to_Matrix();v=s.symbols(f'v0:{n}');forms=[];ranks={2:0,3:0,4:0}
            for row,pc in enumerate(piv):
             degree=len(mons[pc]);ranks[degree]+=1;forms.append(s.expand(sum(rr[row,j]*s.prod(v[i] for i in mon) for j,mon in enumerate(mons) if len(mon)==degree)))
            ar=list(range(ranks[2]));br=list(range(ranks[3]));print('initial ranks',ranks,flush=True)
        else:
            A3=Q3-JU*Q3[rows,:];ar=list(A2.T.to_DM().rref()[1]);ac=list(A2.to_DM().rref()[1]);C=A2[:,ac]*A2.extract(ar,ac).inv(method='DM') if ar else s.zeros(135,0);B=A3-C*A3[ar,:];br=list(B.T.to_DM().rref()[1]);v=s.symbols(f'v0:{n}');forms=[s.expand(sum(c*s.prod(v[i] for i in mon) for c,mon in zip(row,mons))) for M,mons in [(A2[ar,:],mon2),(B[br,:],mon3)] for row in M.tolist()]
        choices=[[1]+[0]*(n-1),list(range(1,n+1)),[1]*(n)]
        result=None
        for coeff in choices:
         substitute=-sum(coeff[i]*v[i] for i in range(1,n));section=s.groebner([f.subs(v[0],substitute) for f in forms],*v[1:],order='grevlex')
         if not section.is_zero_dimensional:continue
         gb=s.groebner([f.subs(v[0],1+substitute) for f in forms],*v[1:],order='grevlex');assert gb.is_zero_dimensional;leads=[p.LM(order=gb.order).exponents for p in gb.polys];counts=[]
         for degree in range(20):
          count=0
          for inds in combinations_with_replacement(range(n-1),degree):
           e=tuple(inds.count(i) for i in range(n-1))
           if not any(all(a<=b for a,b in zip(l,e)) for l in leads):count+=1
          counts.append(count)
          if count==0:break
         result={'parameter_coefficients':coeff,'multiplicity_upper_bound':sum(counts),'localized_basis':list(map(str,gb)),'standard_degrees':counts,'initial_forms':list(map(str,forms)),'section_groebner':list(map(str,section))};break
        assert result is not None and result['multiplicity_upper_bound']==record['length_lower']
        coeff=result['parameter_coefficients'];values=[]
        for branch in record['branches']:
            value=sum(s.Rational(a)*b for a,b in zip(branch['tangent'],coeff));assert value!=0;values.append(str(value))
        result.update(kind=record['kind'],k=record['k'],initial_degree_bound=order_bound,generic_rank=record['length_lower'],parameter_on_branch_tangents=values,t_nilpotent=True,seconds=time.time()-start,jacobian_pivot_rows=rows,jacobian_pivot_columns=jac_piv)
        out.append(result);print(json.dumps({k:result[k] for k in ['kind','k','generic_rank','initial_degree_bound','seconds']}),flush=True)
    paths=[Path(__file__),source_path,ROOT/'code/s61_model.py',ROOT/'code/s61_flag_transverse_obstruction_20261101.py',ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json',m.source/'symplectic_slice_reduced_equations.json']
    return dict(id='S61_MINIMAL_JUNCTION_SATURATION_20261101',all_passed=True,field='QQ',cases=out,conclusion='At all eight specified fixed-beta junctions, the homogenizing coordinate is nilpotent; these points are removed by t-saturation.',limitations=['No flatness or completed-ring isomorphism is claimed at these junctions.','The pointwise proofs extend on their justified opens; no entire generic-beta stratum removal is claimed here.','No full affine-fiber finiteness or exhaustion conclusion.'],proof_note='research/notes/S61_MINIMAL_JUNCTION_SATURATION_20261101.md',seconds=time.time()-start_total,input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);a=p.parse_args();r=run()
    with a.output.open('x') as f:json.dump(r,f,indent=2);f.write('\n')
    print(json.dumps({'all_passed':True,'cases':len(r['cases']),'seconds':r['seconds']}))
if __name__=='__main__':main()
