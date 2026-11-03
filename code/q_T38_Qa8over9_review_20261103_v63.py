#!/usr/bin/env python3
"""Independent quartic-field T38 frame, actual150 arc, edge and seed audit.

All field and tensor arithmetic is from an independent rational4-vector engine.
SymPy is used only to parse the old rational T38 source's ten symbolic terms.
"""
import argparse,gzip,hashlib,itertools,json,time
from fractions import Fraction as F
from pathlib import Path
import sympy as sp
from q_T38_quartic_arithmetic_review_20261103_v63 import Q4,ZERO,ONE,THETA,inverse,norm,determinant,pc,pa,pn,ps,pm,mzero,mid,ma,ms,mm,wedge,push,physical_checks,PAIRS,MONS,BASIS,POS
ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load_report(folder):
 d=json.loads((folder/'report.json').read_text());assert d['all_passed']
 for name,h in d['input_sha256'].items():assert sha(ROOT/name)==h,name
 for name,h in d['evidence_sha256'].items():assert sha(folder/name)==h,name
 return d
def decodepoly(p):return {k:Q4(c)for k,c in p if Q4(c)}
def constmatrix(A):return [[pc(c)for c in row]for row in A]
def fieldmatrix_product(A,B):return [[sum((A[i][k]*B[k][j]for k in range(len(B))),ZERO)for j in range(len(B[0]))]for i in range(len(A))]
def raw_tensor(pair,mon,value=ONE):
 d=[{}for _ in BASIS];d[POS[pair,tuple(sorted(mon))]]=pc(value);return d
def scalar_values(p):
 assert all(not q or set(q)=={0}for q in p);return [q.get(0,ZERO)for q in p]
def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',required=True,type=Path);args=ap.parse_args();out=args.output.resolve();out.mkdir(parents=True,exist_ok=False);start=time.monotonic()
 frames=ROOT/'reports/Q_T38_QA8OVER9_COMPACT_FRAME_20261103_v63';edges=ROOT/'reports/Q_T38_QA8OVER9_EDGE_WEIGHT_TEST_20261103_v63';seeds=ROOT/'reports/Q_T38_QA8OVER9_SEED_WEIGHT_TEST_20261103_v63'
 for folder in [frames,edges,seeds]:load_report(folder)
 with gzip.open(frames/'certificate.json.gz','rt')as f:d=json.load(f)
 Edata=json.loads((edges/'certificate.json').read_text());Sdata=json.loads((seeds/'certificate.json').read_text())
 assert d['field_polynomial_theta']==[305,0,-50,0,1];assert THETA**4-50*THETA**2+305==ZERO
 assert 50%5==305%5==0 and 305%25!=0 #5-Eisenstein, genuine degree4 field.
 rho=(33-THETA**2)/16;z=(THETA-7)*(rho-1)/4;x=-(3*rho+1)*z-(2*rho+1);y=-(2*rho+1)*z-(3*rho+1)
 assert rho*rho-rho-1==ZERO and 2*rho*(z*z+1)+7*z==ZERO
 ds=[ZERO,ONE,z,x,y];alpha=x*x+rho*z*(x+y)+rho*x*y;beta=-rho*(z+x+y)-x;s=3*(rho-1)/2;mu=-4/(9*s*s*(s+3));a=8*rho/(9*(rho*rho-1));assert a==Q4(F(8,9))
 for name,value in [('rho',rho),('z',z),('x',x),('y',y),('alpha',alpha),('beta',beta),('s',s),('mu',mu),('a',a)]:assert Q4(d[name])==value
 assert list(map(Q4,d['diagonal_eigenvalues']))==ds
 gam={}
 for i,j in itertools.combinations(range(1,5),2):
  u,v=ds[i],ds[j];gam[i,j]=u*v*(v-u)*(alpha+beta*(u+v)+rho*(u*u+u*v+v*v)+u*v);assert gam[i,j]==Q4(d['gamma'][f'{i},{j}'])
 p=gam[1,2];pf=gam[1,2]*gam[3,4]-gam[1,3]*gam[2,4]+gam[1,4]*gam[2,3];assert p and pf==p*p/4
 period=[gam[3,4],-gam[2,4],gam[2,3],gam[1,4],-gam[1,3],gam[1,2]];assert period==list(map(Q4,d['period_vector_Pfaffian_dual']))==[ZERO,p/2,ZERO,ZERO,-p/2,p];assert period==[p*Q4(v)for v in d['T38_period_basis']]
 vand=ONE
 for value in ds[1:]:vand=vand*value
 for i,j in itertools.combinations(range(1,5),2):vand=vand*(ds[j]-ds[i])
 assert pf==rho*vand
 unit_values={f'd{i}-d{j}':ds[i]-ds[j]for i in range(5)for j in range(i)};unit_values.update(gamma12=p,gamma_Pfaffian=pf)
 for record in d['unit_certificates']:
  value=unit_values[record['name']];assert value==Q4(record['value'])and value*Q4(record['inverse'])==ONE;assert norm(value)==F(record['norm'])!=0
 assert len(unit_values)==len(d['unit_certificates'])==12
 for value in [s,s+3,mu]:assert value*inverse(value)==ONE and norm(value)!=0
 links=[ONE,s,(s+3)/2,-9/(s+3)];assert links==list(map(Q4,d['N_subdiagonal_links']))
 I=mid();Nt=mzero()
 for i in range(5):Nt[i][i]=pc(ds[i],1)
 for i,c in enumerate(links):Nt[i+1][i]=pc(c)
 powers=[I]
 for k in range(1,5):powers.append(mm(powers[-1],Nt))
 W=[]
 for k in range(1,5):
  trace={}
  for i in range(5):trace=pa(trace,powers[k][i][i])
  correction=[[ps(trace,-F(1,5))if i==j else{}for j in range(5)]for i in range(5)];W.append(ma(powers[k],correction))
 arc=[{}for _ in BASIS]
 for i,j,c in [(0,1,pc(mu*alpha,2)),(0,2,pc(mu*beta,1)),(0,3,pc(mu*rho)),(1,2,pc(mu))]:arc=[pa(a,pm(c,b))for a,b in zip(arc,wedge(W[i],W[j]))]
 assert arc==list(map(decodepoly,d['canonical150_arc']))and all(k>=0 for a in arc for k in a)
 target=[ZERO]*150
 for key,value in {((1,4),(0,0)):Q4(F(8,9)),((2,3),(0,0)):Q4(F(-2,9)),((2,4),(0,1)):Q4(F(4,3)),((3,4),(0,2)):Q4(F(2,3)),((3,4),(1,1)):ONE}.items():target[POS[key]]=value
 assert [a.get(0,ZERO)for a in arc]==target==list(map(Q4,d['canonical150_limit']))
 div,jac=physical_checks(arc);assert not any(div)and not any(jac)
 # Rebuild the unit lower-triangular eigenframe and both Laurent inverse products.
 S0=[[Q4(i==j)for j in range(5)]for i in range(5)]
 for i in range(5):
  for j in range(i):S0[i][j]=-links[i-1]*S0[i-1][j]/(ds[i]-ds[j])
 assert S0==[[Q4(v)for v in row]for row in d['eigenframe_coefficient_S0']]
 S=[[pc(S0[i][j],j-i)for j in range(5)]for i in range(5)];U=ma(S,ms(I,-1));Sinv=I;power=I
 for k in range(1,5):power=mm(power,U);Sinv=ma(Sinv,ms(power,(-1)**k))
 assert mm(S,Sinv)==I and mm(Sinv,S)==I
 diag=mzero()
 for i in range(5):diag[i][i]=pc(ds[i],1)
 assert mm(Nt,S)==mm(S,diag)
 for k in range(1,5):
  expected=mzero();mean=sum((q**k for q in ds),ZERO)/5
  for i in range(5):expected[i][i]=pc(ds[i]**k-mean,k)
  assert mm(mm(Sinv,W[k-1]),S)==expected
 # Actual old T38 source, not a T01 relabeling: compare ALL its symbolic terms.
 oldpath=ROOT/'reports/MPS_SMOOTHING_COMPLETION_20261103_v62/T38.json';old=json.loads(oldpath.read_text());assert old['all_passed']and old['period_pfaffian_open']=='p0**2/4'
 pp=sp.Symbol('p0');tau=sp.symbols('t0:4');xx=sp.symbols('x0:5');local={str(q):q for q in (pp,*tau,*xx)}
 expected={(0,1):tau[0]*xx[2]*xx[4]-2*tau[0]*tau[3]*xx[3]**2/pp,(0,2):tau[1]*xx[1]*xx[3]+2*tau[1]*tau[2]*xx[4]**2/pp,(1,2):pp*xx[1]*xx[2]+4*tau[2]*tau[3]*xx[3]*xx[4]/pp,(1,3):pp*xx[1]*xx[3]/2+tau[2]*xx[4]**2,(2,4):-pp*xx[2]*xx[4]/2+tau[3]*xx[3]**2}
 assert set(old['homogeneous_bivector'])=={','.join(map(str,k))for k in expected}
 for pair,formula in expected.items():assert sp.expand(pp*(sp.sympify(old['homogeneous_bivector'][','.join(map(str,pair))],locals=local)-formula))==0
 log=[{}for _ in BASIS]
 for pair,value in [((1,2),p),((1,3),p/2),((2,4),-p/2)]:log[POS[pair,pair]]=pc(value)
 div,_=physical_checks(log,strict=False);Dmat=mzero()
 for j,row in enumerate(div):
  for mon,c in row.items():assert len(mon)==1;Dmat[j][mon[0]]=c
 vertical=wedge(I,Dmat);canonical=[pa(c,ps(v,-F(1,5)))for c,v in zip(log,vertical)];assert scalar_values(canonical)==list(map(Q4,d['canonical_T38_log150']))
 actual_push=push(canonical,S,Sinv);assert [pm(pc(mu,5),v)for v in actual_push]==arc
 # Push the complete four linear edges and all three quadratic corrections.
 Sconst=constmatrix(S0);U0=ma(Sconst,ms(I,-1));inv0=I;power=I
 for k in range(1,5):power=mm(power,U0);inv0=ma(inv0,ms(power,(-1)**k))
 assert mm(Sconst,inv0)==I and mm(inv0,Sconst)==I
 assert Sinv==[[pc(inv0[i][j].get(0,ZERO),j-i)for j in range(5)]for i in range(5)]
 records=[(0,1,2,4),(0,2,1,3),(1,3,4,4),(2,4,3,3)];weights=[i+j-u-v for(i,j),(u,v)in BASIS];source_weights=[a+b-u-v for a,b,u,v in records];assert source_weights==Edata['source_weights']==[-5,-2,-4,0]and weights==Edata['target_weights']
 pushed=[]
 for a0,b0,u,v in records:
  tensor=push(raw_tensor((a0,b0),(u,v)),Sconst,inv0);dv,_=physical_checks(tensor,strict=False);assert not any(dv);pushed.append(scalar_values(tensor))
 E=[[pushed[j][i]for j in range(4)]for i in range(150)];assert E==[[Q4(v)for v in row]for row in Edata['canonical150_by4_edge_pushforward']]
 selected=Edata['selected4_rows'];assert len(selected)==4 and all(weights[i]>=0 for i in selected);minor=determinant([E[i]for i in selected]);assert minor==Q4(Edata['exact4minor_theta_basis'])and minor*Q4(Edata['exact4minor_inverse'])==ONE;assert norm(minor)==F(Edata['exact4minor_norm'])!=0
 assert [max(W for W,c in zip(weights,col)if c)for col in pushed]==[7]*4
 assert all(not q or set(q)=={5-W}for W,q in zip(weights,arc))
 cross_records=[(0,1,3,3),(0,2,4,4),(1,2,3,4)];cross=[]
 for a0,b0,u,v in cross_records:
  tensor=push(raw_tensor((a0,b0),(u,v)),Sconst,inv0);dv,_=physical_checks(tensor,strict=False);assert not any(dv);cross.append(scalar_values(tensor))
 for pair,record in zip([(0,3),(1,2),(2,3)],cross_records):assert record[0]+record[1]-record[2]-record[3]==sum(source_weights[k]for k in pair)
 assert Edata['distinct_levels_in_full_T38_source']=={'toric':5,'linear':-1,'cross':-7}
 # Cheap seed test: ALL8x7 actual equations, both inverse products, norm/unit.
 seed_rows=[i for i,W in enumerate(weights)if W>=5];assert seed_rows==Sdata['actual_raw_rows']==[90,105,120,121,135,136,137,140];scale=mu*p
 A=[E[i]+[Q4(c)*cross[j][i]/scale for j,c in enumerate([-2,2,4])]for i in seed_rows];assert A==[[Q4(c)for c in row]for row in Sdata['exact8_by7_constant_coefficient_matrix']]
 assert Sdata['selected7_raw_rows']==seed_rows[:7]
 inverse7=[[Q4(c)for c in row]for row in Sdata['first7_equations_inverse']];I7=[[Q4(i==j)for j in range(7)]for i in range(7)];assert fieldmatrix_product(inverse7,A[:7])==I7 and fieldmatrix_product(A[:7],inverse7)==I7
 minor7=determinant(A[:7]);assert minor7==Q4(Sdata['exact7minor'])and minor7*Q4(Sdata['exact7minor_inverse'])==ONE;assert norm(minor7)==F(Sdata['exact7minor_norm'])!=0
 assert Sdata['monomials']==['c0','c1','c2','c3','c0*c3','c1*c2','c2*c3']
 small={'independent_field_arithmetic':'Four rational coefficients; reduce theta^4=50theta^2-305; inverse/norm from rational4x4 multiplication matrix','quartic_irreducible_5_Eisenstein':True,'twelve_unit_norms':{r['name']:r['norm']for r in d['unit_certificates']},'all150_actual_arc_and_eigenframe_pushforward_equal':True,'all25_divergence_and350_Jacobi_coefficients_zero':True,'period_Pfaffian_dual_in_identity_labels':[q.encode()for q in period],'edge_selected_rows':selected,'edge_minor_norm':str(norm(minor)),'seed_actual_rows':seed_rows,'seed_minor_norm':str(norm(minor7)),'levels_fixed_row':{'toric':5,'linear':-1,'cross':-7},'seed_levels':{'toric':5,'linear':5,'cross':5}}
 (out/'certificate.json').write_text(json.dumps(small,indent=2)+'\n')
 files=[Path(__file__),ROOT/'code/q_T38_quartic_arithmetic_review_20261103_v63.py',oldpath]
 for folder in [frames,edges,seeds]:files.extend([folder/'report.json',folder/('certificate.json.gz'if folder==frames else'certificate.json')])
 report={'id':'Q_T38_QA8OVER9_20261103_V63_REVIEW','all_passed':True,'independent_QQ_quartic_vector_arithmetic':True,'all_unit_norms_inverses_and_admissibility_checked':True,'all150_actual_polynomial_arc_coefficients_rebuilt':True,'all150_t0_physical_Qa8over9_coefficients_equal':True,'physical_25_divergence_350_strict_Jacobi_checked':True,'old_actual_T38_full_source_and_Pfaffian_dual_period_bound':True,'all150_Laurent_eigenframe_pushforward_and_inverse_checked':True,'all150x4_edge_coefficients_independently_rebuilt':True,'full_Laurent_frame_inverse_exponent_pattern_checked':True,'all3_cross_tensors_canonical_divergence_zero':True,'edge_nonnegative_weight_minor_rank4_checked':True,'aligned_level_minus1_negative_only_accepted':True,'all8x7_seed_equations_and_both7inverse_products_checked':True,'fixedframe_level5_seed_ideal_exactly_c0_c1_c2_c3':True,'Pplus_P2_T38_closure_membership_or_exclusion_claim':False,'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in files},'seconds':time.monotonic()-start}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:report[k]for k in ['all_passed','all150_actual_polynomial_arc_coefficients_rebuilt','physical_25_divergence_350_strict_Jacobi_checked','all8x7_seed_equations_and_both7inverse_products_checked','seconds']},indent=2))
if __name__=='__main__':main()
