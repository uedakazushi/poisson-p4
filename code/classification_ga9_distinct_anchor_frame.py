"""Load the complete validated QQ frame of the Ga9 distinct-root source."""
from pathlib import Path
import json,time
import sympy as s
from classification_rank3_mixed_component_review import COLUMNS,ROWS
ROOT=Path(__file__).resolve().parents[1]
def build_frame():
 started=time.monotonic();path=ROOT/'reports/CLASSIFICATION_GA9_DISTINCT_WEIGHT_BLOCKS_20261107.json';packet=json.loads(path.read_text());assert packet['all_passed']and(packet['D_rank'],packet['actual_allfive_F_rank'],packet['normal_count'])==(97,47,6)
 assert packet['all150_directsum_verified']and packet['true_normal_relative_weights']==[1]*6
 xx=s.symbols('h x y z w');env=dict(zip(map(str,xx),xx));b=s.Symbol('b');K=s.QQ.frac_field(b);P=s.zeros(5)
 for ij,f in packet['canonical_point'].items():
  i,j=map(int,ij.split(','));P[i,j]=s.sympify(f,locals=env);P[j,i]=-P[i,j]
 normals=[]
 for encoded in packet['true_normal_vectors']:
  v=s.zeros(150,1)
  for i,a in encoded.items():v[int(i)]=s.Rational(a)
  normals.append(v)
 solvers={}
 for block in packet['complete_blocks']:
  if block['implicit_columns']:solvers[block['relative_weight']]=(block['implicit_columns'],block['implicit_rows'],[[K.from_sympy(s.Rational(a))for a in row]for row in block['implicit_inverse']])
 weights=packet['weights'];c=packet['canonical_tensor_weight'];rw=[sum(a*w for a,w in zip(ex,weights))-sum(weights[i]for i in quad)-2*c for quad,ex in ROWS]
 return dict(started=started,prior_path=path,b=b,K=K,P=P,xx=xx,normal_weights=[1]*6,normals=normals,solvers=solvers,selected_rows=packet['implicit_rows'],row_relative_weights=rw,implicit=packet['implicit_columns'])
