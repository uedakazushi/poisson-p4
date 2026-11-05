#!/usr/bin/env python3
"""Actual natural Told->Tnil map, complete source transport and cached origin.
No source/gauge/jet, no higher scalar coefficient or implicit solve.
"""
import hashlib,itertools,json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
import ambient_B18_radical_comparison_operator_20261104_v63 as P
from ambient_B18_honest_nil_T_minimal_operator_20261104_v63 import Matrices
from ambient_B18_honest_nil_T_first_values_20261104_v63 import pivots,invq,enc,nqq,polynomial,mv,sparse
from ambient_B18_native128_first_20261104_v63 import Algebra
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_TOLD_TNIL_NATURAL_MAP_20261104_v63'
PATHS={
 'old_unit':'reports/AMBIENT_B18_FULLV56_COMPARISON_UNIT_20261104_v63/actual-fullV56-m0-and-typed-comparison.json',
 'old_unit_peer':'reports/AMBIENT_B18_FULLV56_COMPARISON_UNIT_20261104_v63_REVIEW/report.json',
 'old_honest':'reports/AMBIENT_B18_HONEST_T_EULER_OPERATOR_20261104_v63/actual-honest-T116-936-Euler-operator.json',
 'old_minimal':'reports/AMBIENT_B18_T_MINIMAL5_SCHUR_OPERATOR_20261104_v63/actual-minimal5-by825-Schur-operator.json',
 'new_honest':'reports/AMBIENT_B18_HONEST_NIL_T_OPERATOR_20261104_v63/actual-honest-Tnil128-972-formal-operator.json',
 'new_honest_peer':'reports/AMBIENT_B18_HONEST_NIL_T_OPERATOR_20261104_v63_REVIEW/report.json',
 'new_minimal':'reports/AMBIENT_B18_HONEST_NIL_T_MINIMAL_OPERATOR_20261104_v63/actual-minimal-Tnil5-849-and-d3-linear-operator.json',
 'new_minimal_peer':'reports/AMBIENT_B18_HONEST_NIL_T_MINIMAL_OPERATOR_20261104_v63_REVIEW/report.json',
 'radical':'reports/AMBIENT_B18_RADICAL_COMPARISON_OPERATOR_20261104_v63/actual-radical-normalizer-comparison-operator.json',
 'original_native':'reports/AMBIENT_B18_COMPLETE_NORMALIZER_SYZYGY_OPERATOR_20261104_v63/actual-complete-Res4-syzygy-operator.json',
 'm0':'reports/AMBIENT_B18_HONEST_NIL_T_FIRST_VALUES_20261104_v63/phase1-actual-m0-and-section.json',
 'X1':'reports/AMBIENT_B18_HONEST_NIL_T_FIRST_VALUES_20261104_v63/phase2-actual-radical-X-firsts.json',
 'first':'reports/AMBIENT_B18_HONEST_NIL_T_FIRST_LINEAR_20261104_v63/phase3-all849-minimal-first-relations-and-rank.json',
 'first_peer':'reports/AMBIENT_B18_HONEST_NIL_T_FIRST_LINEAR_20261104_v63_REVIEW/report.json',
 'monics1':'reports/AMBIENT_B18_NATIVE128_FIRST_20261104_v63/phase1-source-firsts.json',
 'nilfinite':'reports/AMBIENT_B18_NIL_DVR_FINITENESS_20261104_v63/report.json',
 'nilfinite_peer':'reports/AMBIENT_B18_NIL_DVR_FINITENESS_20261104_v63_REVIEW/report.json',
 'old_Cbar':'reports/AMBIENT_B18_CBAR_LITERAL_RESOLUTION_20261104_v63/actual-Cbar-minimal-resolution-and-three-targets.json',
 'new_Cbar':'reports/AMBIENT_B18_REDUCED_COMPARISON_CBAR_VALUES_20261104_v63/actual-eight-nil-images-and-reduced-Cbar-cokernel.json'}
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def poly(p):return {tuple(m):Q(c)for m,c in p if Q(c)}
def pmul(a,b):return P.mul(a,b)
def padd(a,b):return P.add(a,b)
def pscale(a,c):return P.scale(a,c)
def pencode(p):return [[list(m),c]for m,c in sorted(p.items())]
def main():
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30));out=ROOT/'reports'/STEM;out.mkdir(parents=True,exist_ok=False)
 D={k:read(ROOT/p)for k,p in PATHS.items()}
 for k in ['old_unit_peer','new_honest_peer','new_minimal_peer','first_peer','nilfinite_peer']:assert D[k]['all_passed']
 O,H,M,R=(D[k]for k in ['old_unit','new_honest','new_minimal','radical'])
 tails=D['old_minimal']['actual_free_five_source_rows'];assert tails==[22,23,29,30,55]
 assert D['old_honest']['actual_source32_complement_rows']==[j for j in range(56)if j not in O['fixed_selected24_columns']]
 # Literal polynomial map of every36 Koszul(R9) column into complete35.
 R9=[poly(p)for p in R['R9_generators_polynomials9']];DD=[[poly(p)for p in row]for row in R['G8_polynomial_factor_D8_by9']]
 S=[{i:poly(p)for i,p in col}for col in R['native_tensor_Koszul6_HB2_resolution'][1]['columns']]
 pairs9=list(itertools.combinations(range(9),2));pairs8=list(itertools.combinations(range(8),2))
 def coladd(a,b,c=Q(1)):
  z=dict(a)
  for i,p in b.items():
   q=padd(z.get(i,{}),pscale(p,c))
   if q:z[i]=q
   elif i in z:del z[i]
  return z
 def key(c):return tuple((i,tuple(sorted(p.items())))for i,p in sorted(c.items()))
 SK={key(c):j for j,c in enumerate(S)}
 c=P.variable(6);r=P.variable(7);l=P.variable(8)
 hb0={6:l,7:pscale(r,-1)};hb1={6:l,8:pscale(c,-1)};b0=SK[key(hb0)];b1=SK[key(hb1)]
 routing=[]
 for i,j in pairs9:
  target={i:pscale(R9[j],-1),j:R9[i]}
  if key(target)in SK:q={SK[key(target)]:P.const(1)}
  elif key({k:pscale(p,-1)for k,p in target.items()})in SK:q={SK[key({k:pscale(p,-1)for k,p in target.items()})]:P.const(-1)}
  elif(i,j)==(6,7):q={b0:pscale(c,-1)}
  elif(i,j)==(6,8):q={b1:pscale(r,-1)}
  elif(i,j)==(7,8):q={b0:l,b1:pscale(l,-1)}
  else:raise AssertionError((i,j))
  check={}
  for k,a in q.items():
   for row,b in S[k].items():check=coladd(check,{row:pmul(a,b)})
  assert key(check)==key(target);routing.append(q)
 G=[]
 for row in DD:
  g={}
  for a,b in zip(row,R9):g=padd(g,pmul(a,b))
  G.append(g)
 oldsource=D['original_native']['native_S12_by32'];oldS=[{}for _ in range(32)]
 for row,col,gen,sign in oldsource['koszul28_in_first8_rows']:oldS[col][row]=pscale(G[gen],Q(sign))
 J=[]
 for i,j in pairs8:
  q={}
  for p,(a,b)in enumerate(pairs9):
   wedge=padd(pmul(DD[i][a],DD[j][b]),pscale(pmul(DD[i][b],DD[j][a]),-1))
   for k,v in routing[p].items():q=coladd(q,{k:pmul(wedge,v)},-1) # oldS is opposite standardKoszul
  J.append(q)
 J += [{}for _ in range(4)];assert len(J)==32
 checks=0
 for j,q in enumerate(J):
  lhs={}
  for k,a in q.items():
   for row,b in S[k].items():lhs=coladd(lhs,{row:pmul(a,b)})
  rhs={}
  for i,g in oldS[j].items():
   for row,a in enumerate(DD[i]):rhs=coladd(rhs,{row:pmul(a,g)})
  assert key(lhs)==key(rhs);checks+=1
 polynz=sum(len(p)for col in routing+J for p in col.values());assert polynz<150000
 (out/'preflight.json').write_text(json.dumps({'literal_native35_by32_polynomial_terms':polynz,'typed_node_upper':200,'origin_QQ_slots_upper':15000,'seconds':60,'bytes':32<<20,'new_source_or_jet':False},indent=2)+'\n')
 # All-order map defined on five ORIGINAL source-bound cycle functions.
 mat=Matrices();mo=mat.external('old_unit','typed_m24_by56',24,56);no=mat.external('old_unit','typed_n56_by24',56,24)
 E=mat.literal([[Q(i==j)for j in tails]for i in range(56)])
 oldcy=mat.mul(mat.add(mat.eye(56),mat.neg(mat.mul(no,mo))),E)
 Z=mat.external('new_honest','old_to_new_Z44_by56',44,56);cy=mat.mul(Z,oldcy)
 mn=mat.external('new_honest','new_m24_by44',24,44);mcy=mat.mul(mn,cy)
 vectors=[];hcols=[]
 def vec(op,n):i=len(vectors);vectors.append({'operation':op,'dimension':n});return i
 for j in range(5):
  a=vec(['ORIGINAL_matrix_column_functions',mcy,j],24);parts=[]
  for sec in range(3):
   b=vec(['contiguous_slice',a,8*sec,8],8)
   parts.append(vec(['ORIGINAL_complete_Koszul_Hadamard_Euler_lift',b,[1,1,2][sec],'new_honest.shared_homotopy_template_v1'],28))
  hcols.append(vec(['ordered_concatenation',parts],84))
 hh=mat.put(['matrix_from_original_columns',hcols],84,5)
 honestmap=mat.vcat(cy,hh);pi=mat.external('new_minimal','full_target_projection5_by128',5,128);natural=mat.mul(pi,honestmap)
 # Cached origin values require NO old cycle first derivative: Zrow20_0=0.
 m0=[[Q(c)for c in row]for row in O['all24_by56_m0']];Inv=[[Q(c)for c in row]for row in O['fixed24_minor_inverse']];sel=O['fixed_selected24_columns']
 n0=[[Q(0)]*24 for _ in range(56)]
 for i,row in zip(sel,Inv):n0[i]=row
 assert all(sum((m0[i][k]*n0[k][j]for k in range(56)),Q(0))==Q(i==j)for i in range(24)for j in range(24))
 old0=[]
 for j in tails:old0.append([Q(i==j)-sum((n0[i][k]*m0[k][j]for k in range(24)),Q(0))for i in range(56)])
 assert all(not v for c0 in old0 for v in mv(m0,c0))
 a36=Algebra(8,2,[[Q(c)for c in v]for v in D['monics1']['nu_first9']]);A=D['X1']['new_A12_by9_C36']
 Z0=[[Q(0)]*56 for _ in range(44)];Z20first=[[Q(0)]*56 for _ in range(3)]
 for i in range(12):
  for mask in range(4):
   b=[Q(0)]*36;b[mask]=Q(1)
   for gen in range(9):
    v=a36.mul(tuple(map(Q,A[i][gen])),tuple(b))
    for outmask in range(4):Z0[4*gen+outmask][4*i+mask]=v[outmask]
    if gen==5:
     assert not v[0]
     for k,z in enumerate([0,4,5]):Z20first[k][4*i+mask]=v[4*(z+1)]
 for i in range(8):Z0[36+i][48+i]=Q(1)
 new0=[mv(Z0,c0)for c0 in old0];new20first=[mv(Z20first,c0)for c0 in old0]
 assert all(not q for c0 in new0 for q in mv([[Q(c)for c in row]for row in D['m0']['m0_24_by44']],c0))
 Xv=[]
 for col in D['X1']['all780_source_columns44_sparse_dual9']:
  v={}
  for i,sp in col:v[i]={k:Q(c)for k,c in sp}
  Xv.append(v)
 xp=D['X1']['X0_pivot_columns'];free=[i for i in range(44)if i not in D['m0']['selected_source_columns24']]
 I=[free[i]for i in M['fixed_minor_data']['compatible_X0_top18_rows_in_C20']];Ji=[[Q(c)for c in row]for row in M['fixed_minor_data']['inverse']]
 hi=[[Q(c)for c in row]for row in D['first']['all54_actual_pivot_h0_readback']]
 mapped=[]
 for c0,der in zip(new0,new20first):
  q=mv(Ji,[c0[i]for i in I]);out0=[]
  for k in [42,43]:out0.append(c0[k]-sum((Xv[j].get(k,{}).get(0,Q(0))*q[p]for p,j in enumerate(xp)),Q(0)))
  for k in range(3):out0.append(Q(5,2)*der[k]-sum((hi[k][p]*q[p]for p in range(18)),Q(0)))
  mapped.append(out0)
 cp,cr=pivots(mapped);rank=len(cp)
 rankmat=[[mapped[j][i]for j in cp]for i in cr];rankinv,rankdet=invq(rankmat)
 fact=[]
 for col in mapped:
  q=mv(rankinv,[col[i]for i in cr]);assert all(sum((mapped[cp[k]][i]*q[k]for k in range(rank)),Q(0))==col[i]for i in range(5));fact.append(q)
 assert mapped[4][1]==1,'literal old gA generator must retain its unit component'
 data={'source_artifacts':{k:{'path':p,'sha256':sha(ROOT/p)}for k,p in PATHS.items()},
  'source_C_ring_and_frame':'same actual Cfree4 over R0, exact ell0 CA-F9; no abstract type substitution',
  'old_minimal_generator_free56_rows':tails,'old_minimal_honest116_lifts':'w84=0; free32 literal tail vectors; both old Schur codomain right inverses preserve these lifts exactly',
  'native_Koszul36_to_complete35_columns':[[[i,pencode(p)]for i,p in sorted(col.items())]for col in routing],
  'actual_native_old32_to_radical35_columns':[[[i,pencode(p)]for i,p in sorted(col.items())]for col in J],
  'all32_native_full_polynomial_chain_identities_checked':checks,
  'native_transport_orientation':'oldS columns G_j e_i−G_i e_j, hence negative wedge²D followed by standardKoszul36→HB/Koszul35',
  'actual_full_chain_map_old768_to_new780':{'first640':'literal identity','last128':'same actual C/R0 Res4 of saved native35x32 J','dimensions':[780,768],
   'identity':'Z*Xold=Xnew*diag(I640,Res4(J)); A^t N12^t S=(N12 A)^t S=(injection D)^t S; last four zero native columns map0'},
  'matrix_operator_nodes':mat.nodes,'matrix_operator_shapes':mat.shapes,'original_vector_lift_nodes':vectors,
  'matrices':{'old_five_actual_cycles56':oldcy,'actual_Z44_by56':Z,'new_five_actual_cycles44':cy,
   'five_ORIGINAL_honest128_lifts':honestmap,'actual_natural_map_minimal5_by5':natural},
  'linearity_and_descent':'Each of five original actual cycle functions is contracted before normal form, then extend the five columns S-linearly. Full source chain map sends every old boundary to a new boundary. Under accepted complete H1 identifications this is the same D0-kernel inclusion. Chosen honest lifts are fixed by both section-change maps, since Phi_n(c,Hmc)=c.',
  'actual_origin_data':{'old_m0_24_by56':m0,'old_n0_56_by24':n0,'all_five_old_cycle_origins56':old0,
   'actual_Z0_44_by56':Z0,'actual_Z_row20_zi_firsts3_by56':Z20first,'all_five_new_cycle_origins44':new0,
   'all_five_new_cycle_row20_zi_firsts3':new20first,'actual_origin_map_columns5':mapped,'rank':rank,
   'origin_rank_pivot_columns':cp,'origin_rank_pivot_rows':cr,
   'origin_rank_minor_matrix':rankmat,'origin_rank_minor_inverse':rankinv,'origin_rank_minor_determinant':rankdet,
   'all5_origin_column_factorizations':fact,
   'old_cycle_derivatives_not_needed':'actual Zrow20 origin0 kills them; canonical mnew row7 has zero z-derivatives',
   'literal_unit_entry':'new target2,old generator5 =1 (old free55 maps same C0 gA coordinate to newfree43)'},
  'exact_sequence_proof':{
   'sequence':'0→Told→Tnil→NE/ellNE→C/ellC→Cprime/ellCprime→0',
   'hypotheses':'ell regular on actual E4 and Ered; NE actual D-free9; both complete comparison H1 identify with the kernels of the same D0/ell→normalizer quotient maps',
   'derived_argument':'Short exact complex NE[0]→(D0→E4)→(D0→Ered) and Tor1_D(NE,k)=0 give injection. Common Tor1(D0) is in homological degree2, so no H2 is discarded.',
   'lengths':[9,8,4],'actual_cokernel_length':5,'actual_cokernel_mu':5-rank},
  'zero_origin_shift_shortcut_refuted':rank>0,'Told_equals_mTnil_proved':False,
  'cokernel_Q_consequence_if_rank2':{'actual_rank2':rank==2,'possible_Hilbert_functions':[[3,2],[3,1,1]],
   'm_cubed_Q_zero':True,'m_cubed_Tnil_contained_in_image_Told':True,
   'proof':'Q has length5 and mu3; mQ has length2. Its first layer has dimension1 or2, yielding precisely (3,1,1) or (3,2). Both have m³Q=0; no choice between them is inferred.'}if rank==2 else None,
  'no_actual_HF2_total_length_nil_generators_or_cube_inferred':True,
  'visible_polynomial_terms':polynz,'no_new_source_gauge_jet_or_higher_coefficient':True}
 assert nqq(data)<150000
 ev=out/'actual-natural-Told-to-Tnil-source-map-and-origin.json';blob=json.dumps(enc(data),separators=(',',':'))+'\n';assert len(blob.encode())<32<<20;ev.write_text(blob)
 report={'id':STEM,'all_passed':True,'status':'PROVED_actual_injective_natural_map_and_nonzero_origin_countergate',
  'actual_natural_map_origin_rank':rank,'actual_cokernel_length':5,'actual_cokernel_mu':5-rank,
  'old_gA_to_new_gA_origin_entry_one':True,'zero_origin_HF_shift_shortcut_false':rank>0,
  'complete_native35x32_chain_transport_polynomial_identities':32,'no_new_source_gauge_jet_or_higher_scalar_coefficient':True,
  'actual_HF2_Tnil_or_nil_cube_decided':False,'typed_QQ_support':nqq(data),'bytes':ev.stat().st_size,
  'input_sha256':{p:sha(ROOT/p)for p in PATHS.values()},'code_sha256':sha(Path(__file__)),
  'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ev.name:sha(ev)},'seconds':time.monotonic()-t}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
if __name__=='__main__':main()
