#!/usr/bin/env python3
"""Exact Fraction structural identities; no producer or SymPy is imported."""
import hashlib,itertools,json,resource,time
from pathlib import Path
import classification_degenerate_quadric_q4_nil_support as p
p.ZERO=(0,)*32
p.PARAMETERS=['n'+str(i) for i in range(27)]
p.variable=lambda j:{tuple(int(i==j) for i in range(32)):p.F(1)}
ROOT=Path(__file__).resolve().parents[1]
CASES=['AZERO_CNIL2_ISOTROPIC','AZERO_CNIL2_NONISOTROPIC','AZERO_CNIL2_ZERO','AZERO_CZERO_RANK1_NONISOTROPIC','AZERO_CZERO_RANK1_ISOTROPIC','AZERO_CZERO_RANK2_DEGENERATE','AZERO_CZERO_RANK2_NONDEGENERATE']
def src(case):
 path='reports/CLASSIFICATION_DEGENERATE_QUADRIC_Q3_'+case+'_20261107.json';r=json.loads((ROOT/path).read_text());return path,r,p.decode(r['brackets'])
def main():
 resource.setrlimit(resource.RLIMIT_CPU,(55,56));resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2));start=time.monotonic();records=[];paths=[];q=p.parse('2*x0*x2+x1**2');a,b,c,r,z=[p.variable(i) for i in range(5)]
 for case in CASES:
  path,data,P=src(case);paths.append(path);D=p.potential(P);assert D==p.parse(data['euler_contracted_pfaffian_not_yet_actual_divisor']);rec={'case':case,'full_source_dimension':data['source_dimension'],'Euler_Pfaffian_reconstructed':True}
  if case in CASES[:2]:
   linear=a if case==CASES[0] else b;quot=p.divide_exact(D,p.multiply(q,linear));rec['universal_extra_linear_factor']='x0' if linear==a else 'x1';rec['quotient_spatial_degree']=max(sum(e[:5]) for e in quot)
  if case==CASES[2]:
   A=p.parse('-2*n6*x0+2*n9*x1-2*n1*x2');B=p.parse('n19');ell=p.add(A,p.multiply(B,r));J={(0,1):a,(0,2):p.scale(b,-1),(1,2):c};W=[p.scale(p.add(*(p.multiply(p.entry(J,i,j),p.derivative(A,j)) for j in range(3))),p.F(-1,2)) for i in range(3)]
   assert all(p.entry(P,i,j)==p.multiply(ell,p.entry(J,i,j)) for i,j in itertools.combinations(range(3),2));assert all(p.entry(P,i,3)==p.multiply(r,W[i]) for i in range(3))
   rec.update(base_block_identity=True,A='-2*n6*x0+2*n9*x1-2*n1*x2',B='n19',block='R=(A+B*r)J+r W wedge d/dr; W=-J#dA/2')
  if case==CASES[3]:
   assert all(not p.substitute(p.entry(P,1,j),{1:{}}) for j in [0,2,4]);S={5+i:{} for i in [11,17,7]};RR={ij:p.substitute(p.entry(P,*ij),S) for ij in itertools.combinations([0,1,2,4],2)}
   expected={(0,4):p.parse('-n6*x0*x1'),(2,4):p.parse('n6*x1*x2')};assert all(v==expected.get(ij,{}) for ij,v in RR.items());rec.update(nonisotropic_linear_factor_Poisson_normal=True,after_necessary_H_restriction='n11=n17=n7=0',remaining_base_block='n6*x1*(-x0*d/dx0+x2*d/dx2) wedge d/dx4')
  if case==CASES[4]:
   alpha=p.parse('-n7/2');C=p.parse('-n7*x0**2*x4/2-n11*x0**2*x1-n6*x0**3/3+n15*x0**2*x2-n18*(x0*x1*x2+x1**3/3)');J={(0,1):a,(0,2):p.scale(b,-1),(1,2):c}
   assert all(p.entry(P,i,j)==p.multiply(p.multiply(alpha,a),p.entry(J,i,j)) for i,j in itertools.combinations(range(3),2));assert all(p.multiply(a,p.entry(P,i,4))==p.scale(p.add(*(p.multiply(p.entry(J,i,j),p.derivative(C,j)) for j in range(3))),-1) for i in range(3));rec.update(full_base_Jacobian_identity=True,base_cubic='-n7*x0**2*x4/2-n11*x0**2*x1-n6*x0**3/3+n15*x0**2*x2-n18*(x0*x1*x2+x1**3/3)')
   fixed={5+i:p.scalar(v) for i,v in {7:-2,18:-3,11:0,6:0,15:0}.items()};PS={ij:p.substitute(v,fixed) for ij,v in P.items()};C=p.parse('x0**2*x4+x1**3+3*x0*x1*x2');indices=[5+i for i in range(22) if 5+i not in fixed];slots=[]
   for j in range(5):
    f=p.add(*(p.multiply(p.derivative(C,i),p.entry(PS,i,j)) for i in range(5)),p.multiply(C,a) if j==3 else {})
    for mon in sorted(set(e[:5] for e in f)):
     coeff={(0,)*5+e[5:]:v for e,v in f.items() if e[:5]==mon};assert all(sum(e)<=1 for e in coeff);row=[coeff.get(tuple(int(k==i) for k in range(32)),p.F(0)) for i in indices]+[coeff.get(p.ZERO,p.F(0))];slots.append((j,mon,row))
   matrix=[v[2][:] for v in slots];trans=[[p.F(i==j) for j in range(len(slots))] for i in range(len(slots))];at=0
   for col in range(len(indices)):
    pivot=next((i for i in range(at,len(matrix)) if matrix[i][col]),None)
    if pivot is None:continue
    matrix[at],matrix[pivot]=matrix[pivot],matrix[at];trans[at],trans[pivot]=trans[pivot],trans[at];v=matrix[at][col];matrix[at]=[x/v for x in matrix[at]];trans[at]=[x/v for x in trans[at]]
    for i in range(len(matrix)):
     if i!=at and matrix[i][col]:
      v=matrix[i][col];matrix[i]=[x-v*y for x,y in zip(matrix[i],matrix[at])];trans[i]=[x-v*y for x,y in zip(trans[i],trans[at])]
    at+=1
   inconsistent=next(i for i,row in enumerate(matrix) if not any(row[:-1]) and row[-1]);lam=[x/matrix[inconsistent][-1] for x in trans[inconsistent]]
   assert all(sum((v*slots[j][2][i] for j,v in enumerate(lam)),p.F(0))==p.F(i==len(indices)) for i in range(len(indices)+1))
   rec['fiber_cubic_inconsistency']={'remaining_parameter_indices':[i-5 for i in indices],'linear_rows_QQ':at,'augmented_rank_QQ':at+1,'all_equations':[{'component':j,'spatial_monomial':list(mon),'coefficients_then_constant':list(map(str,row))} for j,mon,row in slots],'left_relation_equal_one':{str(i):str(v) for i,v in enumerate(lam) if v},'exact_left_times_variable_columns_zero_and_constant_one':True}
  if case in CASES[5:]:
   assert all(not p.entry(P,i,j) for i,j in itertools.combinations(range(3),2));assert all(not p.derivative(D,i) for i in [3,4]);assert all(not p.derivative(p.entry(P,i,j),k) for i in range(3) for j in [3,4] for k in [3,4]);rec.update(all_base_base_brackets_zero=True,actual_divisor_depends_only_on_base=True,mixed_columns_depend_only_on_base=True)
  records.append(rec)
 result={'task_id':'E_reduced_nonnormal','all_passed':True,'cases':records,'input_sha256':{f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in paths},'independent_arithmetic_module':'code/classification_degenerate_quadric_q4_nil_support.py','independent_arithmetic_module_sha256':hashlib.sha256(Path(p.__file__).read_bytes()).hexdigest(),'producer_source':Path(__file__).read_text(),'seconds':time.monotonic()-start,'scope':'Full coefficient structural identities and an exact affine inconsistency certificate; geometric Casimir-ring and free-curve arguments are separate proof obligations.'}
 out=ROOT/'reports/CLASSIFICATION_DEGENERATE_QUADRIC_Q3_AZERO_REDUCTION_REVIEW_20261107.json'
 with out.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
 print(json.dumps({'all_passed':True,'seconds':result['seconds']},indent=2))
if __name__=='__main__':main()
