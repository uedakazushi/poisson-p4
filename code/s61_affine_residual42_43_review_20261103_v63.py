#!/usr/bin/env python3
"""Independent QQ coefficient and finite-field audit of scalar residual42/43."""
import argparse,gzip,hashlib,itertools,json,subprocess,sys,time
from fractions import Fraction as Q
from pathlib import Path
import numpy as np
import sympy as s
from s61_affine_fiber_macaulay import load_fiber
from s61_affine_polarization_splitter_review_20261103_v63 import detmod,cv
ROOT=Path(__file__).resolve().parents[1];P=32003

def main(out):
 out.mkdir(parents=True,exist_ok=False);started=time.monotonic();inputs=[Path(__file__)];reports={}
 for degree in [43,42]:
  producer=ROOT/f'code/s61_affine_residual{degree}_20261103_v63.py';rp=ROOT/f'reports/S61_AFFINE_RESIDUAL{degree}_20261103_v63/report.json';r=json.loads(rp.read_text());assert r['all_passed']
  for f,h in r['input_sha256'].items():assert hashlib.sha256((ROOT/f).read_bytes()).hexdigest()==h,f
  for f,h in r['evidence_sha256'].items():assert hashlib.sha256((rp.parent/f).read_bytes()).hexdigest()==h,f
  subprocess.run([sys.executable,str(producer),'--output',str(out/f'producer{degree}-replay')],cwd=ROOT,check=True,capture_output=True,text=True)
  reports[degree]=rp;inputs += [producer,rp,ROOT/f'research/notes/S61_AFFINE_RESIDUAL{degree}_20261103_v63.md']
 F,qids,center,*_=load_fiber();fp=ROOT/'reports/S61_KNOWN_FIVE_POINT_SECTION_20261101.json';five=json.loads(fp.read_text());vs=s.symbols('v0:3');E=[[Q(s.sympify(f).coeff(v))for v in vs]for f in five['affine_plane_substitution']]
 qpairs=list(itertools.combinations_with_replacement(range(45),2));index={():0}|{(i,):1+i for i in range(45)}|{m:46+k for k,m in enumerate(qpairs)};width=46+len(qpairs)
 images=[{j:c for j,c in enumerate(row)if c}for row in E]
 for i in range(3,45):images[i][i]=Q(1)
 cache={():[(0,Q(1))]}
 def transformed_monomial(mon):
  if mon in cache:return cache[mon]
  result={}
  if len(mon)==1:
   for j,c in images[mon[0]].items():result[index[(j,)]]=c
  else:
   for i,a in images[mon[0]].items():
    for j,b in images[mon[1]].items():
     k=index[tuple(sorted((i,j)))];result[k]=result.get(k,Q(0))+a*b
  cache[mon]=list(result.items());return cache[mon]
 T=[[Q(0)]*width for _ in range(134)]
 for i,row in enumerate(F):
  for mon,c in row.items():
   for j,v in transformed_monomial(mon):T[i][j]+=c*v
 change=[[Q(v)for v in row]for row in five['all134_restricted_span_coefficients']];G=[]
 for i in range(5,134):G.append([T[i][k]-sum(change[i][j]*T[j][k]for j in range(5))for k in range(width)])
 archive=reports[42].parent/'linear129x42-QQ.json.gz';stored=json.loads(gzip.decompress(archive.read_bytes()));assert stored['embedding']==[[str(c)for c in row]for row in E]
 reconstructed=[]
 for row in stored['all129_linear_factor_rows']:
  coefficients=[Q(0)]*width
  for j,column in enumerate(row):
   for mon,c in column:
    product=tuple(sorted([j+3,*mon]));coefficients[index[product]]+=Q(c)
  reconstructed.append(coefficients)
 assert G==reconstructed and len(G)==129
 knownpath=ROOT/'reports/S61_AFFINE_ENLARGED_SECTION_20261103_v62/report.json';points=json.loads(knownpath.read_text())['known_five_mod_32003'];c42path=reports[42].parent/'residual42-certificate.json';c42=json.loads(c42path.read_text());selected=c42['selected_residual_rows_zero_based'];factor=stored['all129_linear_factor_rows']
 def factor_matrix(z,n,constant=True):
  variables=[*z,*n];M=[]
  for i in selected:
   row=[]
   for col in factor[i]:
    v=0
    for mon,c in col:
     if not mon:v+=cv(c)if constant else 0
     else:v+=cv(c)*variables[mon[0]]
    row.append(v%P)
   M.append(row)
  return M
 determinants=[detmod(factor_matrix(point[1:],[0]*42))for point in points];assert determinants==c42['known_five_determinants_mod_prime'];norm42=1
 for d in determinants:assert d;norm42=norm42*d%P
 assert norm42==c42['normal_determinant_norm_mod_prime'];top42=detmod(factor_matrix(c42['leading_degree42_direction_z'][1:],c42['leading_degree42_direction_n'],False));assert top42==31743
 c43path=reports[43].parent/'residual43-certificate.json';c43=json.loads(c43path.read_text());sep=json.loads((ROOT/'reports/S61_AFFINE_POLARIZATION_SPLITTER_20261103_v63/known-five-separator.json').read_text());hs=sep['separating_h_values'];V=np.array([[pow(h,j,P)for j in range(5)]for h in hs],dtype=np.int64);VI=np.array(s.Matrix(V.tolist()).inv_mod(P)).astype(np.int64);ell=np.array(c43['functional_on_h_coefficient_basis'],dtype=np.int64);weights=ell@VI%P;assert weights.tolist()==c43['split_functional_weights']and np.all(weights)
 slopes=np.array([[cv(c)for c in row]for row in E],dtype=np.int64);ps=[slopes@np.array(point[1:],dtype=np.int64)%P for point in points];absolute=np.array([(p+np.array(list(map(cv,center)),dtype=np.int64))%P for p in ps]);coords=VI@absolute%P;constraints=np.vstack([np.array([1,0,0,0,0]),coords.T]);assert not np.any(constraints@ell%P)
 rr=c43['selected_function_rows_zero_based'];assert detmod(constraints[rr,:4])==c43['selected_constraint4minor_mod_prime']!=0
 rows=sep['selected_original_rows_zero_based'];L=np.zeros((45,45),dtype=np.int64);H=np.zeros((45,45,45),dtype=np.int64)
 for i,ri in enumerate(rows):
  for mon,co in F[ri].items():
   if len(mon)==1:L[i,mon[0]]=(L[i,mon[0]]+cv(co))%P
   else:
    u,v=mon;H[i,u,v]=(H[i,u,v]+cv(co))%P;H[i,v,u]=(H[i,v,u]+cv(co))%P
 direction=np.array(c43['q_direction'],dtype=np.int64);leading=np.einsum('ijk,k->ij',H,direction)*pow(2,-1,P)%P;lead=detmod(leading);assert lead==c43['leading_H_over2_determinant']!=0;inverse=np.array(s.Matrix(leading.tolist()).inv_mod(P)).astype(np.int64);top=[]
 for p in ps:
  constant=(L+np.einsum('ijk,k->ij',H,p)*pow(2,-1,P))%P;M=inverse@constant%P;tr=int(np.trace(M))%P;tr2=int(np.trace(M@M%P))%P;top.append([lead,lead*tr%P,lead*(tr*tr-tr2)*pow(2,-1,P)%P])
 assert top==c43['individual_D_coefficients_at_q_degrees45_44_43'];elltop=[int(weights@np.array([v[j]for v in top],dtype=np.int64))%P for j in range(3)];assert elltop==[0,0,15015];rvalues=[int(w)*int(d)%P for w,d in zip(weights,sep['known_five_jacobian_determinants'])];assert rvalues==c43['degree43_separator_values_at_known_five'];norm43=1
 for value in rvalues:assert value;norm43=norm43*value%P
 assert norm43==17591 and lead*sum(pow(int(v),-1,P)for v in sep['known_five_jacobian_determinants'])%P==10048
 inputs += [fp,knownpath,archive,c42path,c43path,ROOT/'code/s61_affine_polarization_splitter_review_20261103_v63.py',ROOT/'research/notes/TORIC_UNIVERSAL_THREE_PLANE_20261101.md',ROOT/'research/notes/S61_FIVE_COVER_PRODUCT_20261101.md'];result=dict(id='S61_AFFINE_RESIDUAL42_43_20261103_v63_REVIEW',all_passed=True,both_fresh_producer_replays=True,all129_QQ_factorizations_independently_reconstructed=True,all134_input_rows_accounted_for=True,rank42_determinants=determinants,rank42_norm=norm42,rank42_leading_coefficient=top42,functional43_h_basis=ell.tolist(),functional43_split_weights=weights.tolist(),degree43_leading_three_coefficients=elltop,rank43_norm=norm43,raw_trace_degree45_coefficient=10048,scheme_ideal_proofs_accepted=True,accepted_equalities='I+(R42)=I+(R43)=I+(D0,...,D4)=I+(e); beta-only units; all unknown scheme structure retained',residual_empty_or_finite_or_reduced_claimed=False,input_sha256={str(f.relative_to(ROOT)):hashlib.sha256(f.read_bytes()).hexdigest()for f in inputs},seconds=time.monotonic()-started)
 (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');print('RESIDUAL42+43 INDEPENDENT EXACT REVIEW PASS',result['seconds'],flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
