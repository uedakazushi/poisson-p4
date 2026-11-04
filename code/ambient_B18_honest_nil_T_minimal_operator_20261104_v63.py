#!/usr/bin/env python3
"""All-order typed Schur5x849 and d3-only actual linear subspace.
No new source/gauge jet, coefficient series or X2/m1 evaluation.
"""
import hashlib,itertools,json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
from ambient_B18_honest_nil_T_first_values_20261104_v63 import pivots,invq,mv
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_HONEST_NIL_T_MINIMAL_OPERATOR_20261104_v63'
PATHS={
 'honest':'reports/AMBIENT_B18_HONEST_NIL_T_OPERATOR_20261104_v63/actual-honest-Tnil128-972-formal-operator.json',
 'honest_peer':'reports/AMBIENT_B18_HONEST_NIL_T_OPERATOR_20261104_v63_REVIEW/report.json',
 'values_report':'reports/AMBIENT_B18_HONEST_NIL_T_FIRST_VALUES_20261104_v63/report.json',
 'values_peer':'reports/AMBIENT_B18_HONEST_NIL_T_FIRST_VALUES_20261104_v63_REVIEW/report.json',
 'm0':'reports/AMBIENT_B18_HONEST_NIL_T_FIRST_VALUES_20261104_v63/phase1-actual-m0-and-section.json',
 'X1':'reports/AMBIENT_B18_HONEST_NIL_T_FIRST_VALUES_20261104_v63/phase2-actual-radical-X-firsts.json',
 'W0':'reports/AMBIENT_B18_HONEST_NIL_T_FIRST_VALUES_20261104_v63/phase3-actual-W84-rank-and-mu.json',
 'Koszul':'reports/AMBIENT_B18_REDUCED_CBAR_KOSZUL_INPUT_20261104_v63/actual-reduced-Cbar-three-sector-Koszul.json',
 'oldBfirst':'reports/AMBIENT_B18_COMPARISON_X_FIRST640_LINEAR_20261104_v63/actual-first640-X-mod-m2.json'}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def enc(x):
 if isinstance(x,Q):return str(x)
 if isinstance(x,dict):return {str(k):enc(v)for k,v in x.items()}
 if isinstance(x,(tuple,list)):return [enc(v)for v in x]
 return x
def count(x):
 if isinstance(x,Q):return 1
 if isinstance(x,dict):return sum(count(v)for v in x.values())
 if isinstance(x,(tuple,list)):return sum(count(v)for v in x)
 return 0
def sparse(v):return [[i,x]for i,x in enumerate(v)if x]
class Matrices:
 def __init__(self):self.nodes=[];self.shapes=[]
 def put(self,x,r,c):i=len(self.nodes);self.nodes.append(x);self.shapes.append([r,c]);return i
 def external(self,ns,key,r,c):return self.put(['actual_external_matrix',ns,key],r,c)
 def literal(self,A):return self.put(['QQ_matrix_rows',A],len(A),len(A[0])if A else 0)
 def eye(self,n):return self.put(['identity'],n,n)
 def zero(self,r,c):return self.put(['zero'],r,c)
 def select(self,A,I,J):
  r,c=self.shapes[A];assert all(0<=i<r for i in I)and all(0<=j<c for j in J)
  return self.put(['row_column_selection',A,I,J],len(I),len(J))
 def mul(self,*A):
  assert all(self.shapes[a][1]==self.shapes[b][0]for a,b in zip(A,A[1:]))
  return self.put(['matrix_product',list(A)],self.shapes[A[0]][0],self.shapes[A[-1]][1])
 def add(self,A,B):assert self.shapes[A]==self.shapes[B];return self.put(['matrix_sum',A,B],*self.shapes[A])
 def neg(self,A):return self.put(['matrix_neg',A],*self.shapes[A])
 def inv(self,A,origin,det):
  r,c=self.shapes[A];assert r==c and det
  return self.put(['origin_unit_matrix_inverse',A,{'literal_origin_matrix':origin,'determinant':det}],r,r)
 def hcat(self,*A):
  assert len({self.shapes[a][0]for a in A})==1
  return self.put(['horizontal_blocks',list(A)],self.shapes[A[0]][0],sum(self.shapes[a][1]for a in A))
 def vcat(self,*A):
  assert len({self.shapes[a][1]for a in A})==1
  return self.put(['vertical_blocks',list(A)],sum(self.shapes[a][0]for a in A),self.shapes[A[0]][1])
 def schur(self,A,I,J,origin,det):
  r,c=self.shapes[A];Ir=[i for i in range(r)if i not in I];Jr=[j for j in range(c)if j not in J]
  U=self.select(A,I,J);B=self.select(A,I,Jr);C=self.select(A,Ir,J);D=self.select(A,Ir,Jr)
  Ui=self.inv(U,origin,det);S=self.add(D,self.neg(self.mul(C,Ui,B)))
  # P maps the old target into the retained target; P*A removes the unit columns.
  PI=self.select(self.eye(r),I,list(range(r)));PR=self.select(self.eye(r),Ir,list(range(r)))
  P=self.add(PR,self.neg(self.mul(C,Ui,PI)))
  return S,P,{'input':A,'unit_rows':I,'unit_columns':J,'retained_rows':Ir,'retained_columns':Jr,
   'unit_block':U,'unit_inverse':Ui,'B':B,'C':C,'D':D,'Schur':S,'target_projection':P,
   'unit_origin_matrix':origin,'unit_origin_determinant':det,
   'exact_identity':'P*A has zero unit columns and Schur on retained columns; invertible source/target block operations give diag(U,Schur)'}
def main():
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30))
 out=ROOT/'reports'/STEM;out.mkdir(parents=True,exist_ok=False)
 D={k:read(ROOT/p)for k,p in PATHS.items()}
 for k in ['honest_peer','values_peer']:assert D[k]['all_passed']
 assert D['values_report']['mu_R0_Tnil']==5
 for p,h in D['values_report']['input_sha256'].items():assert sha(ROOT/p)==h
 H,M,X,W,K=(D[k]for k in ['honest','m0','X1','W0','Koszul'])
 m0=[[Q(q)for q in r]for r in M['m0_24_by44']];E=M['selected_source_columns24'];C=[i for i in range(44)if i not in E]
 Xi=X['X0_pivot_columns'];assert len(Xi)==18 and len(C)==20
 xd=[]
 for col in X['all780_source_columns44_sparse_dual9']:
  d={}
  for i,sp in col:d[i]={k:Q(v)for k,v in sp}
  xd.append(d)
 XC=[[xd[j].get(i,{}).get(0,Q(0))for j in Xi]for i in C]
 # Select an explicit compatible top18 minor after the24 section cancellation.
 rowchosen,_=pivots(XC);assert len(rowchosen)==18
 J0=[XC[i]for i in rowchosen];Jinv,jdet=invq(J0)
 coeffs=[[Q(q)for q in r]for r in X['all780_X0_factorization_coefficients']]
 assert all(mv(J0,c)==[xd[j].get(C[i],{}).get(0,Q(0))for i in rowchosen]for j,c in enumerate(coeffs))
 top2=[i for i in range(20)if i not in rowchosen];assert len(top2)==2
 Wr=W['W_pivot_rows'];Wp=W['W_pivot_residual_columns'];assert len(Wr)==len(Wp)==81
 Wcols=[]
 for col in W['all_W84_columns_sparse']:
  v=[Q(0)]*84
  for i,q in col:v[i]=Q(q)
  Wcols.append(v)
 Wi=[[Q(q)for q in r]for r in W['W_fixed_minor_inverse']];Wj=[[Wcols[j][i]for j in Wp]for i in Wr]
 assert all(sum((x*y for x,y in zip(row,col)),Q(0))==Q(i==j)for i,row in enumerate(Wj)for j,col in enumerate(zip(*Wi)))
 WC=[i for i in range(84)if i not in Wr];assert len(WC)==3
 # Exact constant quotient projection onto the three remaining wedge targets.
 projection=[]
 for t in WC:
  a=[sum((Wcols[Wp[k]][t]*Wi[k][j]for k in range(81)),Q(0))for j in range(81)]
  v=[Q(0)]*84;v[t]=Q(1)
  for i,c in zip(Wr,a):v[i]-=c
  projection.append(v)
 assert all(not q for c in Wcols for q in mv(projection,c))
 m=Matrices();mop=m.external('honest','new_m24_by44',24,44)
 xop=m.external('honest','actual_X44_by780',44,780);hop=m.external('honest','literal_h84_by780',84,780)
 d2=m.external('honest','d2_new',24,84);d3=m.external('honest','d3_new',84,168)
 EL=m.literal([[Q(i==j)for j in E]for i in range(44)])
 ME=m.mul(mop,EL);ME0=[[m0[i][j]for j in E]for i in range(24)]
 MEi=m.inv(ME,ME0,Q(M['selected_minor_determinant']));n=m.mul(EL,MEi)
 T=m.vcat(m.hcat(xop,n,m.zero(44,168)),m.hcat(hop,m.zero(84,24),d3))
 assert m.shapes[T]==[128,972]
 # Stage1: n is supported exactly on E, so its lower retained block is zero.
 nE0=[[Q(q)for q in row]for row in M['selected_minor_inverse']]
 T1,P1,S1=m.schur(T,E,list(range(780,804)),nE0,1/Q(M['selected_minor_determinant']))
 assert m.shapes[T1]==[104,948]
 T2,P2,S2=m.schur(T1,rowchosen,Xi,J0,jdet)
 assert m.shapes[T2]==[86,930]
 assert S2['retained_columns'][:762]==W['residual_original_X_columns']
 T3,P3,S3=m.schur(T2,[2+i for i in Wr],Wp,Wj,Q(W['W_minor_determinant']))
 assert m.shapes[T3]==[5,849]
 Ptotal=m.mul(P3,P2,P1)
 original=list(range(972));routes=[]
 for stage in [S1,S2,S3]:
  routes.append({'cancelled_original_relation_columns':[original[j]for j in stage['unit_columns']]})
  original=[original[j]for j in stage['retained_columns']]
 assert len(original)==849 and sum(j<780 for j in original)==681 and [j for j in original if j>=804]==list(range(804,972))
 assert sorted(original+sum((s['cancelled_original_relation_columns']for s in routes),[]))==list(range(972))
 # d3 has origin zero; after all Schur steps its linear target is (0_2, projection*d3_1).
 known=[]
 for col in K['all8_differentials_sparse_columns'][2]:
  v=[Q(0)]*40
  for row,mon,sgn in col:
   if sum(mon)==1:
    var=mon.index(1)
    for i in range(3):v[8*(i+2)+var]+=projection[i][row]*Q(sgn)
  known.append(v)
 kp,kr=pivots(known);kinv,kdet=invq([[known[j][i]for j in kp]for i in kr])if kp else([],Q(1))
 for col in known:
  c=mv(kinv,[col[i]for i in kr])if kp else[]
  assert all(sum((known[kp[j]][i]*c[j]for j in range(len(kp))),Q(0))==col[i]for i in range(40))
 # Pin whether source A2*Bold0 remains after constant X0 cancellation; do not assume.
 oldB0=[]
 for col in D['oldBfirst']['source_columns640_rows56_sparse_dual9']:
  v=[Q(0)]*48
  for i,sp in col:
   if i<48:
    for k,q in sp:
     if not k:v[i]=Q(q)
  oldB0.append(v)
 oldB0 += [[Q(0)]*48 for _ in range(140)]
 b0survive=[]
 for j in W['residual_original_X_columns']:
  v=[oldB0[j][i]-sum((coeffs[j][k]*oldB0[Xi[k]][i]for k in range(18)),Q(0))for i in range(48)]
  if any(v):b0survive.append([j,sparse(v)])
 support=sorted({i for j,sp in b0survive for i,q in sp});assert support==[32,36,40,44]
 bspan=[]
 for j,sp in b0survive:
  q=dict(sp);bspan.append([q.get(i,Q(0))for i in support])
 bp,br=pivots(bspan);assert len(bp)==4
 bminor=[[bspan[j][i]for j in bp]for i in br];bminv,bdet=invq(bminor)
 profile={
  'no_X2_m1_values_computed':True,
  'constant_residual_source_order':'c_j=e_j−sum(X0_factor_jk e_pivotk); X0*c_j=0',
  'required_h1_input':'(mX*c_j)_2=m0*(X2*c_j)+m1*(X1*c_j); m2*X0*c_j cancels, m1 does NOT cancel automatically',
  'section_pivot_h0_input':'h_pivot,0 from (mX_pivot)_1=m0*X1_pivot+m1*X0_pivot; needed with the first-order top18 elimination coefficients',
  'stage3_unit_schur_firsts':'S3_1=D1−C1 U0^-1 B0−C0 U0^-1 B1+C0 U0^-1 U1 U0^-1 B0; B0 generally nonzero for X columns',
  'm1_data':{'all_three_sectors_zero7_coefficients':'full quadratic coefficients of actual Vgamma after ordered-prefix restriction',
    'sector0_1_last_r_quotient':'pure r² coefficient of each Vgamma component',
    'sector2_last_r_squared_quotient':'pure r³ coefficient; its z first derivatives vanish by canonical definition',
    'known_cache':'C44/B11 supplies r² and z0r,z1r,z4r,z5r only; it does not supply full z_i*z_j or all z_i*r'},
  'X2_data':{'first640':'ResAT0*Bold2 + ResAT1*Bold1 + ResAT2*Bold0, plus actual P8 seconds',
   'last140':'full base-degree2 of restricted native Koszul6/HB35 coefficient functions and both monics',
   'constant_residual_Bold0_nonzero_count':len(b0survive),'constant_residual_Bold0_nonzero_list':b0survive,
   'surviving_Bold0_span_rank':4,'needed_ResAT2_old_columns_only':support,
   'needed_physical_zero_graph_rows_only':[64,68,73,74],
   'surviving_Bold0_span_minor':bminor,'surviving_Bold0_span_inverse':bminv,'surviving_Bold0_span_determinant':bdet,
   'ResAT2_Bold0_drop_claimed':not b0survive},
  'd3_only_linear_source_does_not_require_any_of_the_above':True,
  'source_or_gauge_jet_job_authorized_or_run':False,
  'full_first_Loewy_layer_or_nil_generators_or_cube_decided':False}
 data={'source_artifacts':{k:{'path':p,'sha256':sha(ROOT/p)}for k,p in PATHS.items()},
  'actual_S_frame':K['actual_S_variable_order'],'operator_nodes':m.nodes,'operator_shapes':m.shapes,
  'external_matrix_namespace':'honest maps named matrices to frozen HONEST_NIL_T_OPERATOR matrices, including original formal h; same ell0 pullback',
  'matrices':{'m':mop,'X':xop,'h':hop,'direct_unit_section_n':n,'base_T128_by972':T,
   'after_section104_by948':T1,'after_X0_86_by930':T2,'minimal_T5_by849':T3,'full_target_projection5_by128':Ptotal},
  'three_exact_unit_Schur_steps':[S1,S2,S3],
  'all972_original_relation_routes':routes,'retained_original_relation_columns849':original,
  'minimal_target_labels':[{'kind':'original_free44_source','index':C[i]}for i in top2]+[{'kind':'Koszul_wedge84','index':i,'label':K['basis_sector_wedge'][2][i]}for i in WC],
  'all_order_completeness_proof':'Replace the old section by E(mE)^-1; mn=I,Vn=d1 gives the same complete128x972 comparison. Three exact unit Schur cancellations give the same cokernel; no relation is dropped without its all-order target projection.',
  'minimality_proof':'all remaining origins are zero by m0 rank24, X0 factorization and W84 complete factorization; hence five targets are minimal',
  'fixed_minor_data':{'compatible_X0_top18_rows_in_C20':rowchosen,'compatible_X0_top18_minor':J0,'inverse':Jinv,'determinant':jdet},
  'd3_only_linear':{'remaining_wedge_indices3':WC,'constant_projection3_by84':projection,
   'all168_known_minimal_linear_columns40_sparse':[sparse(c)for c in known],
   'rank':len(kp),'pivot_columns':kp,'pivot_rows':kr,'minor_determinant':kdet,'minor_inverse':kinv,
   'full_column_factorization_verified':True,'all_first_two_target_linear_coefficients_zero':True,
   'known_first_Loewy_dimension_upper':40-len(kp),'known_length_mod_m_squared_upper':45-len(kp),
   'is_full_linear_relation_space':False},
  'full_first_layer_missing_data_profile':profile,
  'no_scalar_minimal849_coefficients_or_new_source_jet_expanded':True}
 assert count(data)<150000
 ev=out/'actual-minimal-Tnil5-849-and-d3-linear-operator.json';txt=json.dumps(enc(data),separators=(',',':'))+'\n';assert len(txt.encode())<32<<20;ev.write_text(txt)
 note=ROOT/'research/notes'/f'{STEM}.md'
 report={'id':STEM,'all_passed':True,'status':'PROVED_all_order_complete_minimal5x849_and_known_d3_linear_subspace',
  'actual_minimal_presentation_shape':[5,849],'unit_pair_counts':[24,18,81],
  'all780_original_X_and_all168_d3_relations_routed':True,'all_origin_determinant_units_proved':True,
  'd3_only_linear_rank':len(kp),'first_Loewy_dimension_upper_from_d3_only':40-len(kp),
  'full_first_Loewy_layer_or_nil_generators_cube_computed':False,'X2_m1_or_new_source_gauge_jet_computed':False,
  'constant_residual_Bold0_nonzero_count':len(b0survive),'visible_typed_QQ_support':count(data),'bytes':ev.stat().st_size,
  'input_sha256':{p:sha(ROOT/p)for p in PATHS.values()},'code_sha256':sha(Path(__file__)),'note_sha256':sha(note),
  'evidence_sha256':{ev.name:sha(ev)},'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
if __name__=='__main__':main()
