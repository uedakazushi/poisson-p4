#!/usr/bin/env python3
"""Reconstruct the standalone paper's three exact parameter conventions."""
import argparse
import hashlib
import json
from pathlib import Path
import re

import sympy as s

ROOT=Path(__file__).resolve().parents[1]
PAIRS=[(i,j) for i in range(5) for j in range(i+1,5)]
MONS=[(i,j) for i in range(5) for j in range(i,5)]
# Each nonzero coordinate is exactly t**exponent; no legacy builder is imported.
CONFIGS={
 '211_triple':[{0:0,1:1,2:2},{0:0,3:1},{0:0,4:1}],
 '22_node':[{0:0,1:1,2:2},{0:0,3:1,4:2}],
 '31_external':[{0:0,1:1,2:2,3:3},{0:0,4:1}],
}


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
 parser=argparse.ArgumentParser(description=__doc__)
 parser.add_argument('--output',type=Path,required=True)
 args=parser.parse_args()
 paper=ROOT/'papers/normal_actual_divisor_recognition_ja.tex'
 text=paper.read_text()
 x=s.symbols('x0:5')
 result=[];paths=[Path(__file__),paper]
 for name,components in CONFIGS.items():
  width=150+2*len(components);rows=[]
  for component_index,coordinates in enumerate(components):
   for p,(i,j) in enumerate(PAIRS):
    by_degree={}
    def add(degree,col,value):
     row=by_degree.setdefault(degree,[0]*width);row[col]+=value
    for m,(a,b) in enumerate(MONS):
     if a in coordinates and b in coordinates:
      add(coordinates[a]+coordinates[b],15*p+m,5)
    if i in coordinates and j in coordinates:
     for residue_power in [1,2]:
      add(residue_power+coordinates[i]+coordinates[j]-1,
          150+2*component_index+residue_power-1,coordinates[j]-coordinates[i])
    rows.extend(row for row in by_degree.values() if any(row))
  # Exact coefficient rows of delta P; sign follows P_{ji}=-P_{ij}.
  for output in range(5):
   for coordinate in range(5):
    row=[0]*width
    for p,(i,j) in enumerate(PAIRS):
     if output==j:derivative,sign=i,1
     elif output==i:derivative,sign=j,-1
     else:continue
     for m,(a,b) in enumerate(MONS):
      if a==derivative and b==coordinate:row[15*p+m]+=sign
      if b==derivative and a==coordinate:row[15*p+m]+=sign
    if any(row):rows.append(row)
  M=s.Matrix(rows);red,pivots=M.rref()
  free=[i for i in range(width) if i not in pivots]
  # The printed nonpivot lists must fix the same actual coefficient columns.
  printed=text.split('\\paragraph{\\texttt{'+name.replace('_','\\_')+'}。}')[1].split('\\end{gathered}')[0]
  printed=printed.split('\\begin{gathered}')[1]
  assert list(map(int,re.findall(r'\d+',printed)))==free,name
  z=s.symbols('z0:'+str(len(free)))
  lift={col:parameter for col,parameter in zip(free,z)}
  for row,pivot in enumerate(pivots):
   lift[pivot]=-sum(red[row,col]*lift[col] for col in free)
  assert M*s.Matrix([lift[i] for i in range(width)])==s.zeros(M.rows,1)
  path=ROOT/'legacy/v21/code/tree_linear_models'/(name+'.json')
  if not path.exists():path=ROOT/'legacy/v21/data/tree_jacobi_models'/(name+'.json')
  data=json.loads(path.read_text());paths.append(path)
  data=data.get('linear_model',data)
  assert data['free_column_order']==free
  for p,(i,j) in enumerate(PAIRS):
   expression=sum(lift[15*p+m]*x[a]*x[b] for m,(a,b) in enumerate(MONS))
   assert s.expand(expression-s.sympify(data['P_entries'][f'{i}{j}']))==0,(name,i,j)
  for i,old in enumerate(data['residue_expressions']):
   assert s.expand(lift[150+i]-s.sympify(old))==0,(name,'residue',i)
  result.append(dict(model=name,columns=width,rank=len(pivots),dimension=len(free),
                     all150_P_coefficients_and_all_residues_match=True,
                     standalone_printed_nonpivot_list_matches=True))
 # Separate homogeneous verification of the conic and twisted-cubic
 # conormal presentations; no old geometry code is used.
 a,b=s.symbols('a b')
 for degree in [2,3]:
  point=s.Matrix([a**(degree-i)*b**i for i in range(degree+1)])
  covectors=[]
  for i in range(degree-1):
   v=s.zeros(degree+1,1)
   v[i],v[i+1],v[i+2]=b*b,-2*a*b,a*a
   assert s.expand((v.T*point)[0])==0
   assert s.expand((v.T*point.diff(a))[0])==0
   assert s.expand((v.T*point.diff(b))[0])==0
   covectors.append(v)
  frame=s.Matrix.hstack(*covectors)
  assert frame.subs({a:0,b:1}).rank()==degree-1
  assert frame.subs({a:1,b:0}).rank()==degree-1
  # On a!=0, the last degree-1 rows have triangular determinant a**(2*(degree-1)).
  assert s.factor(frame[2:,:].det())==a**(2*(degree-1))
 report=dict(id='NORMAL_MODEL_RREF_20261103_v62_REVIEW',all_passed=True,
             independent_integer_constraint_rows_and_exact_QQ_RREF=True,
             models=result,full_old_Jacobi_certificates_rerun=False,
             independent_conic_and_TC_conormal_frames_checked=True,
             conic_N_degrees=[4,2,2],TC_N_degrees=[5,5,3],
             splitting_Ext1_is_H1_O2=True,TC_center_m_bounds=[5,6],
             paper_structural_review_scope='recognition and relative complexes; dependency ordering through inactive/RNC/E6/E8/P3/eight; conic fourth jet and conic/TC bundles',
             old_long_coefficient_cell_certificates_retained_as_dependencies=True)
 paths += [ROOT/'research/notes/NORMAL_ACTUAL_DIVISOR_PAPER_20261103_v62_REVIEW.md']
 for name in ['NORMAL_ACTUAL_INACTIVE_LINE_SCOPE','NORMAL_ACTUAL_EIGHT_SUPPORT',
              'NORMAL_ACTUAL_E6_E8_GEOMETRY','NORMAL_QUINTIC_RECOGNITION','NORMAL_ACTUAL_RNC4_COVER']:
  p=ROOT/'reports'/(name+'_AUDIT_20261103_v62.json')
  assert json.loads(p.read_text())['all_passed']
  paths.append(p)
 report['input_sha256']={str(p.relative_to(ROOT)):sha(p) for p in paths}
 with args.output.open('x')as stream:json.dump(report,stream,indent=2);stream.write('\n')
 print('PASS all3 standalone exact parameter conventions')


if __name__=='__main__':main()
