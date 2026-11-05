#!/usr/bin/env python3
"""Coherent original head23 edge as two unpaired cache polynomials."""
import hashlib,json,signal,time
from fractions import Fraction as Q
from pathlib import Path
import ambient_B18_C_genuine_projection7_source_v2_20261105_v63 as S
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_C_CONNECTING_HEAD23_SOURCE_20261105_v63'
PATHS={
 'adjoint':'reports/AMBIENT_B18_C_CONNECTING_SIGMA2_ADJOINT_SOURCE_20261105_v63/actual-three-row-sigma2-source-adjoint-and-coherent-fourth-frontier.json',
 'adjoint_peer':'reports/AMBIENT_B18_C_CONNECTING_SIGMA2_ADJOINT_SOURCE_20261105_v63_REVIEW/report.json',
 'q7':'reports/AMBIENT_B18_C_GENUINE_PROJECTION7_SOURCE_20261105_v63/sparse-actual-Dlinear-projector7-source-and-five-head-order-frontier.json',
 'q7_peer':'reports/AMBIENT_B18_C_GENUINE_PROJECTION7_SOURCE_20261105_v63_REVIEW/report.json',
 'unit':'reports/AMBIENT_B18_FULLV56_COMPARISON_UNIT_20261104_v63/actual-fullV56-m0-and-typed-comparison.json',
 'unit_peer':'reports/AMBIENT_B18_FULLV56_COMPARISON_UNIT_20261104_v63_REVIEW/report.json',
 'plane':'reports/AMBIENT_B18_C_CONNECTING_PLANE_ORDER_SOURCE_20261105_v63/all56-plane-orders-surviving-sigma2-adjoint-and-U9-namespace-correction.json',
 'C2':'reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63/actual-physical12-and-C2-right-exact-quotient.json',
 'C2_peer':'reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63_REVIEW/report.json',
 'raw33':'reports/AMBIENT_B18_ELL2_RAW33_QUOTIENT_20261105_v63/actual-ell2-product33-and-raw45-actions.json',
 'Q20':'reports/AMBIENT_B18_CPRIME_QUADRATIC20_VALUES_20261105_v63/phase1-common105-and-114-source-equations.json',
 'Q20_peer':'reports/AMBIENT_B18_CPRIME_QUADRATIC20_VALUES_20261105_v63_REVIEW/report.json',
 'lambda':'reports/AMBIENT_B18_C3_RIGHT2_QUARTIC_CACHED_VALUES_20261105_v63/phase1-six-root-existing-cache-readbacks.json',
 'lambda_peer':'reports/AMBIENT_B18_C3_RIGHT2_QUARTIC_CACHED_VALUES_20261105_v63_REVIEW/report.json',
 'mu':'reports/AMBIENT_B18_C_FULL_HB_PIVOT03_VALUES_20261105_v63/phase1-six-cached-pivot-pure-ell4-values.json',
 'mu_peer':'reports/AMBIENT_B18_C_FULL_HB_PIVOT03_VALUES_20261105_v63_REVIEW/report.json',
 'Cbar':'reports/AMBIENT_B18_CBAR_LITERAL_RESOLUTION_20261104_v63/actual-Cbar-minimal-resolution-and-three-targets.json',
}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def add(a,b):
 v=a.copy()
 for k,q in b.items():
  v[k]=v.get(k,Q(0))+q
  if not v[k]:del v[k]
 return v
def sc(a,q):return{k:v*q for k,v in a.items()if v*q}
def mul(a,b):
 v={}
 for x,c in a.items():
  for y,d in b.items():
   key=tuple(sorted(x+y));assert len(key)<=2
   v[key]=v.get(key,Q(0))+c*d
 return{k:q for k,q in v.items()if q}
def const(q):return{():Q(q)}if q else{}
def enc(p):return[[[list(label)for label in mon],str(q)]for mon,q in sorted(p.items())]
def main():
 started=time.monotonic();signal.alarm(60);pins={p:sha(ROOT/p)for p in PATHS.values()};receipts=[]
 fields={
  'unit':['all24_by56_m0','fixed_selected24_columns','fixed24_minor_inverse'],
  'adjoint':['actual_sigma0_56_by5_sparse','fixed_L_Cn0_row24_sparse'],
  'q7':['kappa_k_A_B_UNPAIRED','projector7_by660_sparse_columns','beta2_plane_after_actual_s_equals_r_plus_ell_over5'],
  'plane':['actual_plane2_B_three_monomials_rows3_columns56','actual_plane2_H_three_monomials_rows3_columns56'],
  'C2':['physical12_coefficients'],'raw33':['actual33_basis'],
  'Q20':['base27_monomials','nu_A_B27'],
  'lambda':['all3_lambda_readbacks'],'mu':['all3_pivot_mu03_readbacks'],
  'Cbar':['CB_shift_coefficients7'],
 }
 d={}
 for k,p in PATHS.items():
  d[k]=S.selective(ROOT/p,fields[k],receipts)if k in fields else json.loads((ROOT/p).read_text())
  if k.endswith('_peer'):assert d[k]['all_passed']
 assert d['adjoint']['actual_sigma0_56_by5_sparse'][1]==[[23,'1']]
 u=d['unit'];selected=u['fixed_selected24_columns'];iv=[[Q(x)for x in row]for row in u['fixed24_minor_inverse']]
 n0=[[Q(0)]*24 for _ in range(56)]
 for i,j in enumerate(selected):n0[j]=iv[i]
 B2=[[Q(q)for q in row]for row in d['plane']['actual_plane2_B_three_monomials_rows3_columns56']]
 H2=[[Q(q)for q in row]for row in d['plane']['actual_plane2_H_three_monomials_rows3_columns56']]
 beta=[[Q(q)for q in row]for row in d['q7']['beta2_plane_after_actual_s_equals_r_plus_ell_over5']]
 shifts=list(map(Q,d['Cbar']['CB_shift_coefficients7']));a1,b1=Q(8,25),Q(96,275)
 L=dict((i,Q(q))for i,q in d['adjoint']['fixed_L_Cn0_row24_sparse']);assert set(L)=={8,12,13}
 raw=d['raw33']['actual33_basis'];phys={r['physical_row']:r for r in d['C2']['physical12_coefficients']}
 mon27=d['Q20']['base27_monomials'];nus=d['Q20']['nu_A_B27'];inputs={};locators={}
 def var(label,value=None,locator=None):
  key=tuple(label)
  if value is not None:
   q=Q(value)
   if key in inputs:assert inputs[key]==q
   inputs[key]=q;locators[key]=locator
  return{(key,):Q(1)}
 def delta(row,ix):return var(('cache',row,ix),dict(phys[row]['delta']).get(ix,'0'),['C2','physical12_coefficients',row,'delta',ix])
 def b(row,ix):return var(('bB',row,ix),dict(phys[row]['bB']).get(ix,'0'),['C2','physical12_coefficients',row,'bB',ix])
 def nu(f,mon):
  ix=mon27.index(mon);return var(('nu',f,ix),nus[f][ix],['Q20','nu_A_B27',f,ix])
 lam={r['physical_row']:r['delta_r_ell2']for r in d['lambda']['all3_lambda_readbacks']}
 mu={r['physical_row']:r['delta_ell3']for r in d['mu']['all3_pivot_mu03_readbacks']}
 def dl(row):return var(('lambda12',row),lam[row],['lambda','all3_lambda_readbacks',row,'delta_r_ell2'])
 def dm(row):return var(('mu03',row),mu[row],['mu','all3_pivot_mu03_readbacks',row,'delta_ell3'])
 def affine(expr):
  p={}
  for key,q in expr:
   assert key[0]=='cache';p=add(p,sc(delta(key[1],key[2]),Q(q)))
  return p
 kap=[[affine(expr)for expr in row]for row in d['q7']['kappa_k_A_B_UNPAIRED']]
 qi={}
 for col in d['q7']['projector7_by660_sparse_columns']:
  exp=col['monomial']
  if col['raw_target']==2 and sum(exp)==1 and any(exp[:7]):
   i=exp.index(1);assert col['target7_affine'][0][0]==6;qi[i]=affine(col['target7_affine'][0][1])
 assert set(qi)=={0,4,5}
 rr,re,ee=24,25,26;assert raw[rr]==[0,2,0]and raw[re]==[0,1,1]
 # Treat newly needed lambda(bB63) as a formal symbol. No root is evaluated.
 new=var(('NEW_bB63_lambda12',));A=nu(0,[7,7]);ArE=nu(0,[7,8]);AE2=nu(0,[8,8])
 B=nu(1,[7,7]);BrE=nu(1,[7,8]);BE2=nu(1,[8,8])
 Az={i:nu(0,[i])for i in[0,4,5]}
 bB67=[-H2[k][21]for k in range(3)]
 w=[{}for _ in range(24)]
 w[7]=sc(B,-bB67[0]);w[15]=sc(add(A,sc(sum_poly(sc(Az[i],shifts[i])for i in[0,4,5]),Q(-1))),-shifts[5])
 for i in range(7):w[16+i]=delta(67,raw.index([i+1,1,0]))
 w[23]=delta(67,29)
 v=[sum_poly(sc(w[k],-n0[j][k])for k in range(24))for j in range(56)]
 assert all(not v[j]for j in[0,4,8,12,16,20,24])
 def dot_const(row,vec):return sum_poly(sc(p,q)for p,q in zip(vec,row))
 h11=dot_const(H2[1],v);h02=dot_const(H2[2],v)
 # NATIVE quartic tauAB67 term -nuB2*bB67,2; its other high is0.
 nativeA=add(sc(B,-bB67[2]+a1*bB67[1]),add(sc(BrE,-bB67[1]+a1*bB67[0]),sc(BE2,-bB67[0])))
 Ah=add(nativeA,sum_poly(mul(kap[0][k],p)for k,p in enumerate([delta(67,30),dl(67),dm(67)])))
 Bh=sum_poly(mul(kap[1][k],p)for k,p in enumerate([delta(67,30),dl(67),dm(67)]))
 # Complete sigma1_r A3/B3 terms. All first7 F_j unit columns have v=0.
 RAv={};RBv={}
 for row,col in[(18,1),(60,17),(67,21)]:RAv=add(RAv,mul(add(dl(row),sc(delta(row,30),-a1)),v[col]))
 RAv=add(RAv,sc(v[28],Q(1,6)-a1*Q(5,3)))
 for row,col in[(18,2),(60,18)]:RBv=add(RBv,mul(add(dl(row),sc(delta(row,30),-b1)),v[col]))
 RBv=add(RBv,mul(add(sc(new,Q(-1)),sc(b(63,30),b1)),v[28]))
 Ah=add(Ah,add(RAv,add(mul(kap[0][0],h11),mul(kap[0][1],h02))))
 Bh=add(Bh,add(RBv,add(mul(kap[1][0],h11),mul(kap[1][1],h02))))
 Rz=[];prefix_z_proofs=[]
 for i in[0,4,5]:
  Ai=b1*beta[i][0]-beta[i][1];Bi=-beta[i][0]
  # ALL sigma1_r mixed lower carry, and -L*m1_r*sigma1_r.
  mixed_ell=sum_poly(mul(p,v[j])for p,j in[(sc(b(63,raw.index([i+1,0,1])),Q(-1)),28),(delta(18,raw.index([i+1,0,1])),2),(delta(60,raw.index([i+1,0,1])),18)])
  mixed_r=sum_poly(mul(p,v[j])for p,j in[(sc(b(63,raw.index([i+1,1,0])),Q(-1)),28),(delta(18,raw.index([i+1,1,0])),2),(delta(60,raw.index([i+1,1,0])),18)])
  Bh=add(Bh,add(sc(mixed_ell,Ai),sc(mixed_r,-L[8+i])))
  # Primitive delta67 mixed term and primitive nuA*z5 contribution.
  Bh=add(Bh,mul(qi[i],delta(67,raw.index([i+1,0,1]))))
  if i==5:Bh=add(Bh,add(sc(ArE,Ai),sc(AE2,Bi)))
  # sigma1_zi=-n0*m1_zi*e23. Its leftadjoint has ONLYBz rows,
  # so all H-sector unknown delta2_zz is killed BEFORE a cache read.
  row=[Ai*B2[1][j]+Bi*B2[2][j]for j in range(56)]
  rz=[sum((row[j]*n0[j][k]for j in range(56)),Q(0))for k in range(24)]
  assert all(not q for k,q in enumerate(rz)if k not in[8,12,13])
  Rz.append([i,[[k,str(q)]for k,q in enumerate(rz)if q]])
  Bh=add(Bh,sc(Az[i],-rz[13]))
  prefix_z_proofs.append([i,'B H0/H4=0; B H5 first zi=nuA_zi; all H-delta2zz terms annihilated by rz n0'])
 # Retained -L*m2_rr*e23, with ORIGINAL-z ordered-prefix correction.
 m2={8:sc(Az[0],-shifts[5]),12:sc(Az[4],-shifts[5]),13:add(A,sc(Az[5],-shifts[5]))}
 Bh=add(Bh,sc(sum_poly(sc(m2[k],q)for k,q in L.items()),Q(-1)))
 canonical=[sc(Ah,Q(5,3)),add(Bh,sc(Ah,Q(-121,48)))]
 # Two old low rows remain UNPAIRED affine source operators for all5heads.
 c0=d['adjoint']['actual_sigma0_56_by5_sparse'];low=[]
 for h,col in enumerate(c0):
  vec=dict((j,Q(q))for j,q in col);pair=[]
  for order in[1,2]:
   poly={}
   for row,j in[(18,1),(60,17),(67,21)]:
    if vec.get(j):poly=add(poly,sc(var(('bB2',row,order),-H2[order][j],['plane','H2',order,j,'negative']),-vec[j]))
   pair.append(poly)
  low.append(pair)
 assert all(not p for h in[1,2,3]for p in low[h])
 polys=[Ah,Bh]+canonical+[p for pair in low for p in pair]
 all_labels={label for p in polys for mon in p for label in mon}
 assert all_labels-set(inputs)=={('NEW_bB63_lambda12',)}
 coefficients=[q for p in polys for q in p.values()]
 assert max(len(mon)for p in polys for mon in p)<=2
 owned=S.count(d)+S.count(v)+S.count(w)+S.count(polys)+S.count(inputs)+S.count(receipts)+20000;assert owned<100000
 data={
  'scope':'FULL coherent original oldhead1=free56row23=physical67*tauAB, plus two low rows onfiveheads. Pure source/cached-input assembly, no output polynomial or new root coefficient is paired.',
  'actual_cycle_origin':'sigma0=e23 exactly; m0e23=0, section fixed. All firstsigma=−n0*m1sigma0, including n1 term by exact fixedsection recurrence.',
  'post_polynomial_semantics':'Monomials are multisets of immutable cache/input labels, QQ coefficients. All output expressions degree<=2; NEW_bB63_lambda12 is a single unpaired source3 coefficient, not a new physical jet solve.',
  'all_input_values_already_cached':[[list(k),str(q),locators[k]]for k,q in sorted(inputs.items())],
  'new_single_cached_composition_required':{'label':['NEW_bB63_lambda12'],'original_physical_bB63_root':38648,'coefficient':'r*ell² afters=r+ell/5','same_given_joint105_degree3_suffices':True,'no_value_emission_or_newIFT_performed':True},
  'sigma1_r_source24':[[i,enc(p)]for i,p in enumerate(w)if p],
  'actual_sigma1_r56_as_cache_polynomials':[[i,enc(p)]for i,p in enumerate(v)if p],
  'first7_Funit_columns_sigma1r_identically0':[0,4,8,12,16,20,24],
  'sigma1_r_source_order_proof':'A curve r4 term is−nuB2_rr*bB67,2rr. B curve usesz=−c*r² andbA67=z5, giving−c5*(nuA2_rr−Σc_i nuA1_zi). H curve uses its OWNunshiftedz7/r²sector, givingdelta67_r³, NOTtheBshift. Allother m1r entries are H z_i Hadamard delta67_zi*r. ThusNOprimitivefourthsourceenterssigma1r forTHIShead.',
  'sigma1_z_complete_leftadjoints_Rz_n0':Rz,'sigma1_z_delta2zz_killed_before_read':prefix_z_proofs,
  'sigma2_full_retained_correction':'Use−L*m2sigma0−L*m1sigma1, with all3Hadamardrows. m2sigma0 isnuA*z5, evaluatedat ORIGINALz0, soitscontributionusesonlynuA2_rr andnuA1_z0/z4/z5.',
  'native_H_second_and_third_carry':'H column23=delta67. Thethree kappa_A/B plane-H corrections, all q(ziH)*delta67_ziell and all mixedB q(zi*rB/ziellB) corrections are explicit in thepolynomials. No isolated primitive projection is substituted for thisfullcoherent column.',
  'head1_low_two_coordinates_literal0':True,
  'head1_genuine_high_A_B_UNPAIRED':[enc(Ah),enc(Bh)],
  'head1_canonical_socle_coordinates_rg_k2_UNPAIRED':[enc(p)for p in canonical],
  'five_head_low_rows_rH_ellH_UNPAIRED':[[enc(p)for p in pair]for pair in low],
  'kernel_type_criterion':'M=J2²+k³ surjectsG=J2+k². Sourcek³basis(h0,h1,h2−h3/4). If the simple-target2×2 minor onh0,h1 is nonzero, ker(M→G)=J2+k; thisminorislow_ellH(h0)*k2(h1) sinceh1low0. No rank claimorvalueisgiven. WithoutthatconditionkerMcaninsteadbek³.',
  'all_other_heads_higher_values_retained_formal':True,
  'positive_worker_gate_still_needed':'Compilefaithfulsingle-root38648cacheclosure+oldloweroverlaps, all inputheightmajorants and postpolynomialbounds, wholecode/custody/sourcepeer/rootONEGO BEFORE any coefficient composition/pairing.',
  'source_polynomial_counts':[len(p)for p in polys],'highest_polynomial_degree':2,
  'whole_source_owned_numeric_upper':owned,'selective_source_receipts':receipts,
 }
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);ep=out/'full-coherent-head23-two-socle-cache-polynomials-and-low-two-row-source.json';ep.write_text(json.dumps(data,indent=2)+'\n')
 rep={'id':STEM,'all_passed':True,'status':'SOURCE_FULL_HEAD23_COHERENT_TWO_HIGHS_ONE_CACHED_THIRD_COMPOSITION_UNPAIRED',
  'primitive_source4_or_higher_needed_for_head23':False,'new_physical_source_solve_needed':False,
  'one_preexisting_joint105_third_cache_composition_needed':True,'head23_edge_values_computed':False,
  'input_sha256':pins,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/(STEM+'.md')),
  'evidence_sha256':{ep.name:sha(ep)},'whole_numeric_upper':owned,'seconds':time.monotonic()-started}
 (out/'report.json').write_text(json.dumps(rep,indent=2)+'\n');print(json.dumps(rep));signal.alarm(0)
def sum_poly(xs):
 p={}
 for q in xs:p=add(p,q)
 return p
if __name__=='__main__':main()
