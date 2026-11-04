#!/usr/bin/env python3
"""One bounded weight5 FULL135 cofactor divisibility test for the actual ell."""
import argparse,ast,collections,hashlib,json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_BIANCHI_ELL_WEIGHT5_20261104_v63'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def add(p,k,c):
 z=p.get(k,Q())+c
 if z:p[k]=z
 else:p.pop(k,None)
def determinant(A):
 A=[r[:]for r in A];ans=Q(1)
 for j in range(len(A)):
  k=next((i for i in range(j,len(A))if A[i][j]),None)
  if k is None:return Q()
  if k!=j:A[k],A[j]=A[j],A[k];ans=-ans
  v=A[j][j];ans*=v
  for i in range(j+1,len(A)):
   c=A[i][j]/v
   for k in range(j+1,len(A)):A[i][k]-=c*A[j][k]
   A[i][j]=Q()
 return ans
def kernel(A,n):
 B={};selected=[]
 for label,row in sorted(A.items()):
  row=row.copy()
  while row:
   j=min(row)
   if j not in B:
    c=row[j];B[j]={k:v/c for k,v in row.items()};selected.append(label);break
   c=row[j]
   for k,v in B[j].items():add(row,k,-c*v)
  if len(B)==n:break
 free=[j for j in range(n)if j not in B];out=[]
 for f in free:
  v={f:Q(1)}
  for j,row in sorted(B.items(),reverse=True):
   z=-sum((c*v.get(k,Q())for k,c in row.items()if k!=j),Q())
   if z:v[j]=z
  assert all(not sum((c*v.get(k,Q())for k,c in row.items()),Q())for row in A.values());out.append(v)
 return len(B),selected,out
def main(out):
 began=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.signal(signal.SIGALRM,lambda *_:(_ for _ in()).throw(TimeoutError('single60s budget')));signal.alarm(60);out.mkdir(parents=True,exist_ok=False)
 opath=ROOT/'reports/S61_NONSQUARE_FULL_BIANCHI_OPERATOR_20261104_v63/operator.json';parent=opath.parent/'report.json';pr=read(parent);assert pr['all_full_raw_map_and_literal_identity_checks_passed']and sha(opath)==pr['operator_sha256']
 for rel,h in pr['input_sha256'].items():assert sha(ROOT/rel)==h,rel
 op=read(opath);cpath=ROOT/'reports/S61_BINARY_HYPERPLANE_WEIGHTED_CHART_20261102.json';chart=read(cpath);w=op['local_coordinate_weights'];ew=op['equation_weights'];assert w==chart['local_coordinate_weights']and ew==chart['equation_weights']
 adpath=ROOT/'reports/S61_NONSQUARE_ORACLE_FRAME_ADAPTER_20261103_v63/adapter.json';ad=read(adpath)
 anchorpath=ROOT/'reports/S61_NONSQUARE_UNIVERSAL_QUAD22_20261103_v63/certificate.json';anchor=read(anchorpath)
 gp=ROOT/'reports/S61_BINARY_HYPERPLANE_GENERIC_RANKS_20261102.json';pointdata=read(gp);assert pointdata['all_passed']and pointdata['point']==[1,2,3];point=tuple(map(Q,pointdata['point']))
 gpcode=ROOT/'code/s61_binary_hyperplane_generic_ranks_20261102.py';assert sha(gpcode)==pointdata['input_sha256'][str(gpcode.relative_to(ROOT))]
 assert "a,e,d=s.symbols('a e d');point={a:1,e:2,d:3}"in gpcode.read_text()
 T=chart['T34'];U=chart['U34'];normal=chart['normal34'];assert len(normal)==8 and len(T)==len(U)==45 and len(T[0])==11 and len(U[0])==34
 assert [i for i,c in enumerate(w)if c==0]==[19,23]
 for i in[19,23]:assert all(T[i][j]=='0'for j in normal+[1])and all(c=='0'for c in U[i])
 assert [[T[i][j]for j in[2,5]]for i in[19,23]]==[['1','0'],['0','1']]
 assert sorted(normal+[1,2,5])==list(range(11))
 # v8 -> (x0,x1,x2,x3,x4,Z,W,H), with v6=Z+x1-3*x0 and v3=W.
 M=[[0]*8 for _ in range(8)]
 for i,j in[(0,0),(1,1),(2,2),(3,6),(4,3),(5,7),(6,5),(7,4)]:M[i][j]=1
 M[6][1]=1;M[6][0]=-3
 Uweights=[]
 for j in range(34):
  weights={w[i]for i in range(45)if U[i][j]!='0'};assert len(weights)==1;Uweights.append(weights.pop())
 neww=[1]*7+[2,1]+Uweights;assert len(neww)==43
 for i,j in enumerate(normal):assert{w[r]for r in range(45)if T[r][j]!='0'}=={chart['normal_weights34'][i]}
 assert{w[r]for r in range(45)if T[r][1]!='0'}=={1}
 pieces=[(i,r)for i,r in enumerate(op['local_weight_pieces'])if r['weight']==5];assert len(pieces)==20
 # Save the EXACT full layout before the single old-point rank test. No source equation is discarded.
 preflight={'operator_sha256':sha(opath),'frame_source_sha256':{str(p.relative_to(ROOT)):sha(p)for p in[cpath,adpath,anchorpath,gp,gpcode]},'weight':5,'source_piece_indices':[i for i,_ in pieces],
 'full_equation_column_count':135,'unknown_core_combination_count':20,'new_independent_variables':43,'eliminated_ell_variable':1,
 'dense_remainder_constraint_count':135*43,'dense_remainder_keys':'g=0..134, multiplier=-1 or newvariable0..42 except1',
 'raw_cofactor_record_count':sum(len(r['coefficients'])for _,r in pieces),'raw_core_term_count':sum(len(ts)for _,r in pieces for _,_,ts in r['coefficients']),
 'normal_frame_recipe':'raw45 = T34[:,normal34]*D*M*normal8 + T34[:,1]*base1 + U34*normalizer34; nuisance19,23=0, no other variable or normalizer set0',
 'D':ad['Tden'],'M_oldv8_from_newnormal8':M,'new_variable_weights':neww,'retained_base_T_column':1,'retained_U_columns':list(range(34)),
 'ell':'the frozen x7 ell with H coefficient0; x1=-sum(ell_j*x_j for j!=1)/(ell_x1)',
 'known_unit_ell_x1':'-16*C*J/243','open':'FULLU=a*d*F*B*C*J!=0',
 'extra_tG6_piece_weights':[r['weight']for r in op['local_weight_pieces']if r['spatial_degree']=='extra_tG6'],
 'F5_F6_and_tG6_indices_retained_even_if_their_weight5_cofactors_zero':True,
 'source_identity':'Every row is a frozen true homogeneous-weight component of a full135 identity, in an independent LINEAR extended chart. No nonlinear normal graph or source319 transport.',
 'attempt_scope':'QQ(a,e,d)-linear combinations constant in all43 localvariables, division by THIS parameter-dependent ell, trueweight5 only'}
 (out/'preflight.json').write_text(json.dumps(preflight,indent=2)+'\n')
 # Exact Fraction evaluation of the SAME accepted point; AST avoids floating division.
 values=dict(zip(['a','e','d'],point));cache={}
 def evnode(n):
  if isinstance(n,ast.Constant):assert isinstance(n.value,int);return Q(n.value)
  if isinstance(n,ast.Name):return values[n.id]
  if isinstance(n,ast.UnaryOp):return -evnode(n.operand)if isinstance(n.op,ast.USub)else evnode(n.operand)
  assert isinstance(n,ast.BinOp);a,b=evnode(n.left),evnode(n.right)
  if isinstance(n.op,ast.Add):return a+b
  if isinstance(n.op,ast.Sub):return a-b
  if isinstance(n.op,ast.Mult):return a*b
  if isinstance(n.op,ast.Div):assert b;return a/b
  assert isinstance(n.op,ast.Pow)and b.denominator==1;return a**b.numerator
 def ev(st):
  if st not in cache:cache[st]=evnode(ast.parse(st,mode='eval').body)
  return cache[st]
 def poly(ts):return sum((Q(c)*point[0]**ex[0]*point[1]**ex[1]*point[2]**ex[2]for ex,c in ts),Q())
 D=list(map(ev,ad['Tden']));frame=[[Q()for _ in range(43)]for _ in range(45)]
 for v in range(45):
  for j in range(8):frame[v][j]=sum((ev(T[v][normal[i]])*D[i]*M[i][j]for i in range(8)),Q())
  frame[v][8]=ev(T[v][1])
  for j in range(34):frame[v][9+j]=ev(U[v][j])
 assert all(not c for i in[19,23]for c in frame[i])
 assert all(not c or w[v]==neww[j]for v,row in enumerate(frame)for j,c in enumerate(row))
 blockproof=[]
 for weight in sorted(set(neww)):
  ri=[i for i in range(45)if w[i]==weight];ci=[j for j in range(43)if neww[j]==weight];assert len(ri)==len(ci)
  square=[[frame[i][j]for j in ci]for i in ri];det=determinant(square);assert det
  blockproof.append({'weight':weight,'raw_local_rows':ri,'new_coordinate_columns':ci,'square_matrix':[[str(c)for c in row]for row in square],'exact_determinant':str(det)})
 ell=[Q()for _ in range(8)]
 for mon,ts in anchor['ell']:assert sum(mon)==1;ell[mon.index(1)]+=poly(ts)
 aa,ee,dd=point;C=100*aa+18*ee**2-9*ee-2;J=5*aa*dd-6*ee-4
 assert ell[1]==-Q(16,243)*C*J and ell[1]and not ell[7]
 matrices={};before={};sourceg=collections.Counter()
 for col,(pieceid,row)in enumerate(pieces):
  fwd={}
  for g,var,ts in row['coefficients']:
   assert 0<=g<135 and -1<=var<45;assert ew[g]+(0 if var==-1 else w[var])==5
   if var in[19,23]:continue
   c=poly(ts);sourceg[ew[g]]+=1
   if var==-1:add(fwd,(g,-1),c)
   else:
    for j,z in enumerate(frame[var]):add(fwd,(g,j),c*z)
  for key,c in fwd.items():before.setdefault(key,{})[col]=c
  reduced=fwd.copy()
  for g in range(135):
   coef=reduced.pop((g,1),Q())
   if coef:
    for j in range(7):
     if j!=1:add(reduced,(g,j),-coef*ell[j]/ell[1])
  for key,c in reduced.items():matrices.setdefault(key,{})[col]=c
 assert all(var==-1 and ew[g]==5 or var>=0 and ew[g]+neww[var]==5 for g,var in matrices)
 rank,selected,ker=kernel(matrices,20);square=[];det=None
 if rank==20:
  square=[[matrices[k].get(j,Q())for j in range(20)]for k in selected];det=determinant(square);assert det
 cert={'old_accepted_point':[1,2,3],'parameter_order':['a','e','d'],'frame_block_proofs':blockproof,
 'oldpoint_ell_newnormal8':[str(c)for c in ell],'oldpoint_denominators_nonzero_verified':True,
 'all135_remainder_constraint_layout':preflight['dense_remainder_keys'],'full_dense_shape':[5805,20],
 'nonzero_oldpoint_constraint_rows':len(matrices),'exact_oldpoint_rank':rank,
 'selected_constraint_keys':[list(k)for k in selected],'selected20_minor':[[str(c)for c in row]for row in square],
 'selected20_minor_determinant':str(det)if det is not None else None,
 'all_nonzero_oldpoint_constraints':[[list(k),[[j,str(c)]for j,c in sorted(row.items())]]for k,row in sorted(matrices.items())],
 'oldpoint_kernel_vectors_diagnostic_only':[[str(v.get(j,Q()))for j in range(20)]for v in ker],
 'equation_weight_record_counts_before_linear_chart':dict(sourceg),
 'equationweight5_constant_constraint_rows':sum(ew[g]==5 and var==-1 for g,var in matrices),
 'equationweight6_constraint_rows_explicitly_zero':True,'extra_tG6_weight6_not_eligible_for_weight5':True,
 'generic_negative_proved':rank==20,'generic_survivor_proved':False,'point_rankdrop_implies_no_generic_kernel_claim':rank<20,
 'scope':'Only the20 frozen trueweight5 rows, combination coefficients inQQ(a,e,d) constant inlocalvariables, division bythe actualell in the full43-variable extended linear chart. No general saturation/no allpoint nonexistence.',
 'no_nonlinear_normalgraph_source319_or_kernel_count_claim':True}
 (out/'certificate.json').write_text(json.dumps(cert,indent=2)+'\n')
 inputs=[Path(__file__),ROOT/'research/notes'/(ID+'.md'),opath,parent,ROOT/'reports/S61_NONSQUARE_FULL_BIANCHI_OPERATOR_20261104_v63_REVIEW/report.json',cpath,adpath,ROOT/'reports/S61_NONSQUARE_ORACLE_FRAME_ADAPTER_20261103_v63_REVIEW/report.json',anchorpath,gp,gpcode]
 report={'id':ID,'all_exact_source_frame_and_oldpoint_checks_passed':True,'source_rows':20,'full_equation_columns_retained':135,'full_dense_constraint_shape':[5805,20],
 'nonzero_oldpoint_rows':len(matrices),'exact_QQ_rank_at_frozen_oldpoint':rank,'no_generic_core_combination_divisible_by_this_ell_proved':rank==20,
 'no_generic_survivor_claim':True,'no_new_point_prime_GB_source_regeneration_normalgraph_or319_transport':True,
 'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in inputs},'evidence_sha256':{p.name:sha(p)for p in out.iterdir()if p.is_file()},
 'artifact_bytes':sum(p.stat().st_size for p in out.iterdir()if p.is_file()),'seconds':time.monotonic()-began,'frozen':True}
 assert report['artifact_bytes']<32*1024*1024
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=ROOT/'reports'/ID);a=p.parse_args()
 try:main(a.output.resolve())
 except(TimeoutError,MemoryError)as error:
  a.output.mkdir(parents=True,exist_ok=True);(a.output/'BLOCKED.json').write_text(json.dumps({'completed':False,'constraint':type(error).__name__,'message':str(error),'no_auto_retry':True},indent=2)+'\n');raise
