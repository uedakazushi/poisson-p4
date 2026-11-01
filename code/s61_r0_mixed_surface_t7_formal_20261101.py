#!/usr/bin/env python3
"""Specified rank4 thickening of the mixed R0 surface T7 and an ambient open."""
import argparse,json,time,hashlib
from pathlib import Path
from itertools import combinations_with_replacement
import sympy as s
from s61_model import S61
from s61_flag_transverse_obstruction_20261101 import specialize
from s61_minimal_infinity_local_data_20261101 import calculate
ROOT=Path(__file__).resolve().parents[1]
def enc(M):return [[str(v) for v in row] for row in M.tolist()]
def universal():
 r=json.loads((ROOT/'reports/S61_R0_MIXED_ORBIT_REDUCTION_20261101.json').read_text())['orbits']['T7'];m=S61();one=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==1];two=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==2];H=s.Matrix(r['beta_kernel']).applyfunc(s.Rational);bb=s.symbols('b:'+str(H.cols));q0=s.Matrix([s.Rational(r['point61'][i]) for i in two]);_,F=specialize(m,dict(zip(one,H*s.Matrix(bb))),True);fixed=next(i for i,v in enumerate(q0) if v);cols=[i for i in range(45) if i!=fixed];J=s.zeros(135,45);start=time.time()
 for row,f in enumerate(F):
  for mon,c in f.items():
   if len(mon)==2:
    i,j=mon
    if i!=fixed:J[row,cols.index(i)]+=c*q0[j]
    if j!=fixed:J[row,cols.index(j)]+=c*q0[i]
   elif len(mon)==1:J[row,44]+=c*q0[mon[0]]
 rr,piv=J.to_DM().rref();rr=rr.to_Matrix();free=[i for i in range(45) if i not in piv];T=s.zeros(45,len(free))
 for j,f in enumerate(free):
  T[f,j]=1
  for k,p in enumerate(piv):T[p,j]=-rr[k,f]
 print('rank',len(piv),'tangent',T.cols,'tzero',T[-1,:]==s.zeros(1,T.cols),time.time()-start,flush=True)
 enc=lambda M:[[str(v) for v in row] for row in M.tolist()]
 out=dict(beta_kernel=enc(H),jacobian=enc(J),tangent=enc(T),pivot_columns=list(piv),free_columns=free)
 pass
 n=T.cols;ms={d:list(combinations_with_replacement(range(n),d)) for d in [2,3]};ix={mon:i for i,mon in enumerate(ms[3])};Q=s.zeros(135,len(ms[2]));terms=[]
 for row,f in enumerate(F):
  ts=[]
  for mon,c in f.items():
   if fixed in mon:continue
   i,j=[cols.index(k) for k in mon]+[44]*(2-len(mon));ts.append((i,j,c))
   for k,(a,b) in enumerate(ms[2]):Q[row,k]+=c*(T[i,a]*T[j,b]+(T[i,b]*T[j,a] if a!=b else 0))
  terms.append(ts)
 aug,ap=J.row_join(-Q).to_DM().rref();assert len(ap)==len(piv);aug=aug.to_Matrix();W2=s.zeros(45,len(ms[2]))
 for k,p in enumerate(piv):W2[p,:]=aug[k,45:]
 assert (J.to_DM()*W2.to_DM()+Q.to_DM()).is_zero_matrix
 print('quad allzero','t2zero',W2[-1,:]==s.zeros(1,len(ms[2])),time.time()-start,flush=True)
 C3=s.zeros(135,len(ms[3]))
 for row,ts in enumerate(terms):
  for i,j,c in ts:
   for a in range(n):
    for k,mon in enumerate(ms[2]):C3[row,ix[tuple(sorted((a,)+mon))]]+=c*(T[i,a]*W2[j,k]+T[j,a]*W2[i,k])
 aug,ap=J.row_join(-C3).to_DM().rref();assert len(ap)==len(piv);aug=aug.to_Matrix();W3=s.zeros(45,len(ms[3]))
 for k,p in enumerate(piv):W3[p,:]=aug[k,45:]
 assert (J.to_DM()*W3.to_DM()+C3.to_DM()).is_zero_matrix
 print('cubic allzero','t3zero',W3[-1,:]==s.zeros(1,len(ms[3])),time.time()-start,flush=True)
 out.update(quadratic_correction=enc(W2),cubic_correction=enc(W3),quadratic_vanishes=True,cubic_vanishes=True,t_quadratic_zero=W2[-1,:]==s.zeros(1,len(ms[2])),t_cubic_zero=W3[-1,:]==s.zeros(1,len(ms[3])),seconds=time.time()-start);pass

 assert len(piv)==42 and T[-1,:]==s.zeros(1,3) and W2[-1,:]==s.zeros(1,6) and W3[-1,:]==s.zeros(1,10)
 sub=dict(zip(bb,range(2,13)));denoms=set(s.fraction(v)[1] for M in [T,W2,W3] for v in M);assert all(d.subs(sub)!=0 for d in denoms);out['all_frame_denominators_nonzero_at_point']=True
 return out

def point_certificate():
 r=json.loads((ROOT/'reports/S61_R0_MIXED_ORBIT_REDUCTION_20261101.json').read_text())['orbits']['T7'];m=S61();p=s.Matrix(r['point61']).applyfunc(s.Rational);one=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==1];local=calculate(m,True,p);assert local['homogenized_jacobian_rank']==42 and local['quadratic_span_rank']==0
 J=s.Matrix(local['jacobian']).applyfunc(s.Rational);T=s.Matrix(local['tangent_basis']).applyfunc(s.Rational);H=s.Matrix(r['beta_kernel']).applyfunc(s.Rational);BD=s.Matrix(r['source_beta_differential']).applyfunc(s.Rational);QD=s.Matrix(r['q_orbit_differential']).applyfunc(s.Rational);curve=s.zeros(44,H.cols).row_join(QD)*s.Matrix.hstack(*BD.nullspace());cp=list(curve.to_DM().rref()[1]);assert len(cp)==2;curve=curve[:,cp].col_join(s.zeros(1,2));tr=list(T.T.to_DM().rref()[1]);ct=T[tr,:].inv()*curve[tr,:];assert T*ct==curve;params=list(ct.T.to_DM().rref()[1]);assert params==[0,1]
 rows=list(J.T.to_DM().rref()[1]);cols=list(J.to_DM().rref()[1]);U=s.zeros(45,len(rows));iv=J[rows,cols].inv(method='DM')
 for k,c in enumerate(cols):U[c,:]=iv[k,:]
 _,F=specialize(m,{i:p[i] for i in one},True);chart=local['chart_columns'];fixed=next(i for i in range(45) if i not in chart);terms=[]
 for f in F:terms.append([([chart.index(i) for i in mon]+[44]*(2-len(mon)),co) for mon,co in f.items() if fixed not in mon])
 W={1:T[:,2]};A={}
 for degree in [2,3,4]:
  Q=s.zeros(135,1)
  for row,ts in enumerate(terms):
   for (i,j),co in ts:
    for k in range(1,degree):Q[row]+=co*W[k][i]*W[degree-k][j]
  W[degree]=-U*Q[rows,:];A[degree]=Q+J*W[degree]
 assert A[2]==A[3]==s.zeros(135,1) and A[4]!=s.zeros(135,1)
 assert all(W[d][-1]==0 for d in [1,2,3,4]);row=next(i for i,v in enumerate(A[4]) if v);columns=list(BD.to_DM().rref()[1]);assert len(columns)==16;dm=BD[:,columns].to_DM();minor=dm.domain.to_sympy(dm.det());assert minor
 return dict(local=local,actual_surface_tangent=enc(ct),parameter_rows=params,transverse_normal_column=2,transverse_corrections={str(d):list(map(str,W[d])) for d in [2,3,4]},transverse_residuals={str(d):list(map(str,A[d])) for d in [2,3,4]},quartic_nonzero_row=row,quartic_coefficient=str(A[4][row]),transverse_length_upper=4,source_beta_rank=16,beta_minor_columns=columns,beta_minor=str(minor),relative_known_surface_dimension=13)
def verify():
 start=time.time();u=universal();p=point_certificate();paths=[Path(__file__),ROOT/'reports/S61_R0_MIXED_ORBIT_REDUCTION_20261101.json',ROOT/'code/s61_model.py',ROOT/'code/s61_flag_transverse_obstruction_20261101.py',ROOT/'code/s61_minimal_infinity_local_data_20261101.py',ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json',ROOT/'research/notes/S61_HOMOGENIZED_EQUIVARIANCE_20261101.md']
 return dict(id='S61_R0_MIXED_SURFACE_T7_FORMAL_20261101',all_passed=True,field='QQ and QQ(b0,...,b10)',universal=u,point=p,completed_ring='C[[u0,u1,e]]/(e^4)',free_rank=4,homogenizer_zero=True,relative_free_rank=4,ambient_open_removed=True,whole_generic_beta_surface_removed=False,limitations=['Specified point and a nonempty ambient open of the18D T7 stratum.','Its17D complement can still dominate beta; the entire orbit is not removed.','T5 and the four mixed curve orbits need their own analysis.'],input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},seconds=time.time()-start)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args();r=verify()
 with a.output.open('x') as f:json.dump(r,f,indent=2);f.write('\n')
 print(json.dumps({k:r[k] for k in ['all_passed','free_rank','homogenizer_zero','seconds']}))
