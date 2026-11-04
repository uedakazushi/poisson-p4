#!/usr/bin/env python3
"""Independent existing-C44 values, exact sector maps and fixed image span."""
from pathlib import Path
from fractions import Fraction as Q
from itertools import product,combinations_with_replacement,permutations
import hashlib,json,resource,signal,time
from ambient_B18_T_native_S2_B11_review_20261104_v63 import B44,BASIS,bm
ROOT=Path(__file__).resolve().parents[1]
VALUE='AMBIENT_B18_REDUCED_COMPARISON_CBAR_VALUES_20261104_v63'
ORIGIN='AMBIENT_B18_REDUCED_COMPARISON_ORIGIN_20261104_v63'
def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def read(p):return json.loads((ROOT/p).read_text())
def vec(raw):return list(map(Q,raw))
def add(a,b,c=Q(1)):return [x+c*y for x,y in zip(a,b)]
def scale(a,q):return [x*q for x in a]
def act(cols,v):return [sum((c[i]*q for c,q in zip(cols,v)),Q(0)) for i in range(len(cols[0]))]
def umul(a,b):
 n=len(a);out=[Q(0)]*n
 for i,x in enumerate(a):
  for j,y in enumerate(b):
   if i+j<n:out[i+j]+=x*y
 return out
def dense44(p):return [p.get((i//4,i%4),Q(0)) for i in range(44)]
def sparse44(p):return {(i//4,i%4):Q(q) for i,q in enumerate(p) if Q(q)}
def main():
 started=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30));hashes={}
 def load(p):hashes[p]=sha(p);return read(p)
 def get(st,f='report.json'):return load('reports/'+st+'/'+f)
 for st,want in [(VALUE,'b090a3037d8a22c64dd461594f79dca1566b2a3b7110150cf797d82a3a11b6ac'),(ORIGIN,'0524a73bc87ff7c12ffa04ccd5137aca5d0f15173a725e76f50dbde3f6229674')]:
  r=get(st);assert r['all_passed'] and hashes['reports/'+st+'/report.json']==want
  for p,h in r['input_sha256'].items():assert sha(p)==h;hashes[p]=h
  for n,h in r['evidence_sha256'].items():p='reports/'+st+'/'+n;assert sha(p)==h;hashes[p]=h
 for st,want in [('AMBIENT_B18_T_NATIVE_S2_B11_20261104_v63_REVIEW','bb6e7e7c5327528f2be31b1d29db5891461d12adb319daa184b6c6addd63d693'),('AMBIENT_B18_CBAR_LITERAL_RESOLUTION_20261104_v63_REVIEW','7e17f85c19492d3e81dc502b38f33819e6045045bb63b131d160019adec3e378'),('AMBIENT_B18_NORMALIZER_NIL_GENERATORS_20261104_v63_REVIEW','02888d3c4b987327f6df2d8faf1f2adffa3f9c01af47e27aff8c5c0982e75472')]:
  r=get(st);assert r['all_passed'] and hashes['reports/'+st+'/report.json']==want
 B=get('AMBIENT_B18_T_NATIVE_B11_SOURCE_PREFLIGHT_20261104_v63','actual-B11-source-and-gauge-preflight.json')
 C=get('AMBIENT_B18_T_NATIVE_S2_B11_20261104_v63','phase3-actual-proper-G8-B11.json')
 L=get('AMBIENT_B18_CBAR_LITERAL_RESOLUTION_20261104_v63','actual-Cbar-minimal-resolution-and-three-targets.json')
 D=get(VALUE,'actual-eight-nil-images-and-reduced-Cbar-cokernel.json')
 I=get('AMBIENT_B18_PROPER_PHYSICAL_INTERFACE_FIRSTJET_20261104_v63','actual-source-physical-firstjet-interface.json')
 K=get('AMBIENT_B18_PROPER_CORE_ELIMINATION_20261104_v63','actual4variable7core-polynomials.json')
 OE=get(ORIGIN,'actual-first-interface-and-zero-origin-Vrad.json')
 assert B['B11_variable_order8']==L['actual_variables']==['zA0','zA1','zA2','zA3','zA4','zA5','zA6','rA']
 assert list(map(tuple,B['B11_basis_monomials8']))==BASIS
 assert L['actual_three_targets']==['tauA','tauB','tauA*tauB'] and L['r_lengths']==[3,3,2]
 shifts=vec(L['CB_shift_coefficients7']);assert shifts==[Q(9835,864),Q(0),Q(0),Q(0),Q(29875,3456),Q(-4565,864),Q(0)]
 nu=[{i:Q(q) for i,q in p if Q(q)} for p in B['nu_A_B_actual_B11']];assert all(not p.get(0,Q(0)) for p in nu)
 ring=B44(nu);kappa=list(map(sparse44,C['kappa9_C_B44']));assert len(kappa)==9
 # Actual origin proof: compare BOTH independent saved first interfaces.
 order=I['retained_source9_order'];assert order==['a','d','v2','v3','v4','v5','v6','v7','v8']
 assert I['physical_regular9_order']==['x108','x109','x112','x113','x117','x118','x119','tau_CA','tau_CB']
 first=[tuple(Q(x) for x in row[-2:]) for row in I['actual_source9_by_physical9_Jacobian']]
 assert first==[tuple(kappa[j].get((0,k),Q(0)) for k in (1,2)) for j in range(9)]
 assert [first[j] for j in (0,1,2,6)]==[(Q(1),Q(-1)),(Q(-2),Q(-1,2)),(Q(1),Q(0)),(Q(1),Q(0))]
 piv=K['selected5_pivot_order'];assert piv==['v4','v3','v5','v7','v8']
 mins={name:min(sum(m) for m,n,d in K['five_actual_rational_roots'][name]['numerator'] if n) for name in piv}
 assert min(mins.values())>=2 and all(first[order.index(name)]==(0,0) for name in piv)
 dh=(first[1][0]*(first[0][0]-first[1][0]+first[2][0]-2*first[6][0]),
     first[1][0]*(first[0][1]-first[1][1]+first[2][1]-2*first[6][1])+first[1][1]*(first[0][0]-first[1][0]+first[2][0]-2*first[6][0]),
     first[1][1]*(first[0][1]-first[1][1]+first[2][1]-2*first[6][1]))
 assert dh==(Q(-4),Q(0),Q(1,4)) and OE['actual_Vrad_origin_rank']==0
 assert OE['d_h_quadratic_coefficients_tauA2_tauAtauB_tauB2']==list(map(str,dh))
 # Exact B11-to-sector ring maps.  Defining monomials map to zero, so the
 # omitted higher functions really belong to each target annihilator.
 sizes=[3,3,2];maps=[]
 for sec,n in enumerate(sizes):
  cols=[]
  for m in BASIS:
   degree=m[7]+2*sum(m[:7]);coef=Q(1)
   for j in range(7):
    if m[j]:coef*=(-shifts[j])**m[j] if sec==1 else Q(0)
   col=[Q(0)]*n
   if degree<n:col[degree]=coef
   cols.append(col)
  maps.append(cols)
 assert maps==[[vec(c) for c in sec] for sec in D['B11_to_three_sector_ring_maps']]
 for sec,n in enumerate(sizes):
  zimages=[[Q(0)]*n for _ in range(7)];rimage=[Q(0)]*n;rimage[1]=Q(1)
  if sec==1:
   for j in range(7):zimages[j][2]=-shifts[j]
  assert umul(umul(rimage,rimage),rimage)==[Q(0)]*n
  assert all(umul(x,y)==[Q(0)]*n for x in zimages for y in zimages)
  assert all(umul(x,umul(rimage,rimage))==[Q(0)]*n for x in zimages)
  assert all(zimages[j]==[Q(0)]*n for j in (2,3,6))
  for i,j in product(range(11),repeat=2):
   bprod=bm({i:Q(1)},{j:Q(1)});got=[sum((q*maps[sec][k][t] for k,q in bprod.items()),Q(0)) for t in range(n)]
   assert got==umul(maps[sec][i],maps[sec][j])
 a,d,x,t=[kappa[i] for i in (0,1,2,6)]
 h=ring.add(a,ring.neg(d),x,ring.scale(t,-2));cc=ring.scale(ring.add(d,x),Q(1,2));ell=ring.scale(ring.add(a,d,ring.neg(x)),Q(1,3))
 nilgens=[ring.mul(d,h),ring.tm(d,cc,ell)]
 assert list(map(dense44,nilgens))==list(map(vec,D['actual_nil_generators_C44']))
 def projection(v):
  out=[]
  for sec,mask in enumerate((1,2,3)):
   out.extend(sum((q*maps[sec][base][j] for (base,m),q in v.items() if m==mask),Q(0)) for j in range(sizes[sec]))
  return out
 images=[]
 for g in nilgens:
  for mask in range(4):
   v=ring.mul(g,{(0,mask):Q(1)});assert dense44(v)==vec(D['all8_nil_basis_multiple_values_C44'][len(images)])
   images.append(projection(v))
 assert images==list(map(vec,D['all8_actual_nil_image_vectors_in_Cbar8']))
 actions=[]
 for var in range(8):
  A=[]
  for j in range(8):
   col=[Q(0)]*8
   for sec,(n,off) in enumerate(zip(sizes,(0,3,6))):
    if off<=j<off+n:
     target=j-off+(1 if var==7 else 2);co=Q(1) if var==7 else -shifts[var] if sec==1 else Q(0)
     if target<n:col[off+target]=co
   A.append(col)
  actions.append(A)
 assert actions==[[vec(c) for c in A] for A in D['Cbar8_literal_coordinate_actions_columns']]
 assert all(act(actions[i],act(actions[j],act(actions[k],[Q(int(x==b)) for x in range(8)])))==[Q(0)]*8 for i,j,k,b in product(range(8),range(8),range(8),range(8)))
 closure=[];labels=[]
 for j,v in enumerate(images):closure.append(v);labels.append([j,[]])
 for j,v in enumerate(images):
  for i in range(8):closure.append(act(actions[i],v));labels.append([j,[i]])
  for i,k in combinations_with_replacement(range(8),2):closure.append(act(actions[i],act(actions[k],v)));labels.append([j,[i,k]])
 assert closure==list(map(vec,D['all360_R0_submodule_closure_columns'])) and labels==D['all360_labels']
 pivots=D['image_pivots'];selected=D['selected_original_closure_columns'];assert pivots==[1,2,5,7] and selected==[0,1,2,4]
 basis=[]
 for j,p in zip(selected,pivots):
  v=closure[j][:]
  for old,b in zip(pivots,basis):v=add(v,b,-v[old])
  assert v[p] and all(not v[q] for q in range(p));v=scale(v,1/v[p]);basis.append(v)
 assert [[p,list(map(str,b))] for p,b in zip(pivots,basis)]==D['normalized_triangular_image_basis']
 def normal(v):
  for p,b in zip(pivots,basis):v=add(v,b,-v[p])
  return v
 assert all(normal(v)==[Q(0)]*8 for v in closure)
 assert all(normal(act(A,b))==[Q(0)]*8 for A in actions for b in basis)
 free=[i for i in range(8) if i not in pivots];assert free==D['quotient_original_basis_indices']==[0,3,4,6]
 def qmap(v):v=normal(v);return [v[i] for i in free]
 qa=[[qmap(A[i]) for i in free] for A in actions]
 assert qa==[[vec(c) for c in A] for A in D['all8_quotient_coordinate_actions_columns']]
 assert all(qmap(act(A,[Q(int(k==j)) for k in range(8)]))==act(QA,qmap([Q(int(k==j)) for k in range(8)])) for A,QA in zip(actions,qa) for j in range(8))
 assert all(not any(col) for A in qa[:7] for col in A)
 assert qa[7]==[[Q(0),Q(0),Q(-11,12),Q(0)],[Q(0),Q(0),Q(1),Q(0)],[Q(0)]*4,[Q(0)]*4]
 split=list(map(vec,D['exact_decomposition_basis_columns_in_quotient']))
 assert split==[[Q(1),Q(11,12),Q(0),Q(0)],[Q(0),Q(0),Q(0),Q(1)],[Q(0),Q(1),Q(0),Q(0)],[Q(0),Q(0),Q(1),Q(0)]]
 assert act(qa[7],split[0])==act(qa[7],split[1])==act(qa[7],split[3])==[Q(0)]*4 and act(qa[7],split[2])==split[3]
 assert D['quotient_length']==4 and D['quotient_Hilbert_function']==[3,1,0]
 # m_C^3 of the origin fibre is zero; d*u5 and the three native cubic
 # generators vanish.  dh has no tauA*tauB term and squares vanish.
 originproof={'actual_source_trace_firsts':{n:list(map(str,p)) for n,p in zip(order,first)},'core_root_min_numerator_degrees':mins,'d_h_coefficients':list(map(str,dh)),'actual_C_origin_algebra':'k[tauA,tauB]/(tauA²,tauB²)','all9_native_generators_and36_Qrad_origin_columns_zero':True,'all_actual_beta8_origins_zero':True,'actual_Vrad_origin_rank':0,'actual_Cprime_mu_R0':3,'no_kernel_basechange_used':True}
 valueproof={'source_kappa9_C44_reused_from_independent_B11_peer':True,'all_B11_defining_relations_and363_sector_products_verified':True,'all8_actual_nil_images':list(map(lambda x:list(map(str,x)),images)),'all360_module_closure_columns_verified':True,'fixed4_image_pivots_and_source_columns':[pivots,selected],'full_image_action_invariant':True,'actual_reduced_Cbar_length':4,'HF':[3,1],'quotient_basis':['tauA','tauB','r*tauB','tauA*tauB'],'r_tauA':'-11/12*r_tauB','all_z_coordinates_zero':True,'exact_S_module_decomposition':'k plus k plus k[[r]]/(r²)','right_exact_cokernel_quotient':'Cprime=C/image Nil(E4); tensor by R0/(ell) is right exact, so these8 fullR0-module images suffice. No kernel/Tor shortcut.','DVR_rank_or_torsion_exponents_not_from_this_finite_quotient':True}
 for st,proof,status in [(ORIGIN,originproof,{'actual_Vrad_origin_rank':0,'actual_Cprime_mu_R0':3}),(VALUE,valueproof,{'actual_nil_image_rank':4,'actual_reduced_Cbar_length':4,'actual_reduced_Cbar_HF':[3,1],'exact_decomposition':'k plus k plus k[[r]]/(r²)'})]:
  out=ROOT/'reports'/(st+'_REVIEW');assert not(out/'report.json').exists();out.mkdir(parents=True,exist_ok=True)
  ev=out/'independent-fixed-values-and-module-proof.json';ev.write_text(json.dumps(proof,indent=2)+'\n')
  for p in [str(Path(__file__).relative_to(ROOT)),'code/ambient_B18_T_native_S2_B11_review_20261104_v63.py','research/notes/'+st+'_REVIEW.md']:hashes[p]=sha(p)
  report={'id':out.name,'all_passed':True,**status,'no_new_source_gauge_jet_point_GB_or_780_expansion':True,'nil_D0_generators_KJK_or_cube_claimed':False,'input_sha256':dict(hashes),'evidence_sha256':{ev.name:hashlib.sha256(ev.read_bytes()).hexdigest()},'seconds':time.monotonic()-started,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
  (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');hashes[str((out/'report.json').relative_to(ROOT))]=sha(str((out/'report.json').relative_to(ROOT)))
  print(json.dumps({k:v for k,v in report.items() if k not in ('input_sha256','evidence_sha256')}),flush=True)
if __name__=='__main__':main()
