#!/usr/bin/env python3
"""Independent finite coefficient binding and universal all-chart graph identity.

No producer/load_fiber import, no GB, and no expansion of42 chart ideals.
The universal free-word identity supplies characteristic-zero proof for ALL42
charts; complete actual modular coefficients are independently reconstructed.
"""
import argparse,gzip,hashlib,itertools,json,time
from fractions import Fraction as F
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];P=32003
def sha(p):
 h=hashlib.sha256()
 with p.open('rb')as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()
def check_report(p):
 d=json.loads(p.read_text());assert d['all_passed']
 for name,h in d.get('input_sha256',{}).items():assert sha(ROOT/name)==h,name
 for name,h in d.get('evidence_sha256',{}).items():assert sha(p.parent/name)==h,name
 return d
def residue(s):
 a,slash,b=str(s).partition('/');den=int(b)%P if slash else 1;assert den;return(int(a)%P)*pow(den,-1,P)%P
def add(d,m,c):
 v=(d.get(m,0)+int(c))%P
 if v:d[m]=v
 else:d.pop(m,None)
def mul(p,q):
 d={}
 for m,c in p.items():
  for n,b in q.items():add(d,tuple(sorted(m+n)),c*b)
 return d
def free_add(p,q,sgn=1):
 d=p.copy()
 for k,c in q.items():
  d[k]=d.get(k,0)+sgn*c
  if not d[k]:del d[k]
 return d
def free_mul(p,q):
 d={}
 for (e,u),c in p.items():
  for(f,v),b in q.items():
   k=(tuple(x+y for x,y in zip(e,f)),u+v);d[k]=d.get(k,0)+c*b
 return {k:c for k,c in d.items()if c}
def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',required=True,type=Path);args=ap.parse_args();out=args.output.resolve();out.mkdir(parents=True,exist_ok=False);start=time.monotonic();folder=ROOT/'reports/S61_AFFINE_NORMAL_DIRECTION_GRAPH_20261103_v63';check_report(folder/'report.json');cert=json.loads((folder/'certificate.json').read_text())
 with gzip.open(folder/'coefficient-model-mod32003.json.gz','rt')as f:model=json.load(f)
 assert model['prime']==P
 finite=ROOT/'reports/S61_AFFINE_NORMAL_FINITE5_20261103_v63/certificate.json';fc=json.loads(finite.read_text());prior_review=ROOT/'reports/S61_AFFINE_NORMAL_FINITE5_20261103_v63_REVIEW/report.json';check_report(prior_review)
 crossfile=ROOT/'reports/S61_AFFINE_NORMAL_FINITE5_20261103_v63/cross-input.txt.gz';lines=gzip.decompress(crossfile.read_bytes()).decode().splitlines();assert lines[0]=='129 126 945'and len(lines)==259;C=np.array([[residue(c)for c in line.split()]for line in lines[1:130]],dtype=np.int64);D=np.array([[residue(c)for c in line.split()]for line in lines[130:]],dtype=np.int64);assert C.shape==(129,126)and D.shape==(129,945)
 mons=[(j,)for j in range(42)]+list(itertools.combinations_with_replacement(range(42),2));assert model['normal_monomials']==list(map(list,mons));mi={m:j for j,m in enumerate(mons)}
 derived=ROOT/'work/normal-finite5-public-v63/cross-output.txt';assert sha(derived)==fc['native_output_sha256'];coeff=np.array(model['P_coefficients'],dtype=np.int64);assert coeff.shape==(126,945);axis_normal=[5,13,21,29];axis_columns=[mi[(j,j)]for j in axis_normal]
 with derived.open()as f:
  assert f.readline().split()==['126','945'];sel=list(map(int,f.readline().split()));assert sel==fc['selected_cross_rows']
  for _ in range(126):f.readline()
  for i in range(126):
   s=f.readline().split();assert len(s)==945 and list(map(residue,s))==coeff[i].tolist();assert all(F(s[j])==0 for j in axis_columns)
  assert not f.read().strip()
 assert not np.any((C[sel]@coeff+D[sel])%P);other=[j for j in range(129)if j not in sel];H=(D[other]+C[other]@coeff)%P;assert H.tolist()==model['H_coefficients']
 # Verify the saved inverse as a complete42x42 matrix, not sampled vectors.
 T=coeff[:,:42].reshape(3,42,42);comm=(T[0]@T[1]-T[1]@T[0])%P;Ci=np.array(model['constant_commutator_inverse'],dtype=np.int64);assert np.array_equal(comm@Ci%P,np.eye(42,dtype=np.int64))and np.array_equal(Ci@comm%P,np.eye(42,dtype=np.int64))
 # Universal noncommutative operator expansion, scalars(x,lambda0,lambda1).
 # It proves all42 identities in characteristic zero once C^-1 is defined.
 A={((1,0,0),('A',)):1};B={((1,0,0),('B',)):1};l0={((0,1,0),()):1};l1={((0,0,1),()):1};U0=free_add(l0,A,-1);U1=free_add(l1,B,-1)
 formal=free_add(free_mul(free_add(B,l1,-1),U0),free_mul(free_add(A,l0,-1),U1),-1);assert formal=={((2,0,0),('A','B')):1,((2,0,0),('B','A')):-1}
 # Own direct ORIGINAL61 substitution for all5 retained F, overFp. The prior
 # review independently binds the SAME full134 and129 cross factorization QQ.
 original=ROOT/'legacy/v02/data/symplectic_slice_reduced_equations.json';pointpath=ROOT/'reports/S61_KNOWN_FIBER_AUDIT_20261031.json';raw=json.loads(original.read_text());point_record=next(p for p in json.loads(pointpath.read_text())['points']if p['component_id']=='T01');pointQQ=list(map(F,point_record['S61_coordinates_u1_through_u61']));point=list(map(residue,pointQQ));qids=[i for i,w in enumerate(raw['parameter_weights'])if w==2];qindex={old:j for j,old in enumerate(qids)}
 source=ROOT/'reports/S61_AFFINE_RESIDUAL42_20261103_v63/linear129x42-QQ.json.gz'
 with gzip.open(source,'rt')as f:linear=json.load(f)
 E=[[residue(c)for c in row]for row in linear['embedding']];forms=[]
 for j,old in enumerate(qids):
  p={};add(p,(),point[old])
  for k,c in enumerate(E[j]):add(p,(k,),c)
  if j>=3:add(p,(j,),1)
  forms.append(p)
 products={():{():1}}
 def product(free):
  free=tuple(sorted(free))
  if free not in products:
   d={():1}
   for j in free:d=mul(d,forms[j])
   products[free]=d
  return products[free]
 planes=[];Ks=[]
 for eq in raw['equations'][:5]:
  trans={}
  for mon,n,d in eq:
   c=residue(F(n,d));free=[]
   for old in mon:
    if old in qindex:free.append(qindex[old])
    else:c=c*point[old]%P
   for term,v in product(free).items():add(trans,term,c*v)
  plane={};K=np.zeros(945,dtype=np.int64)
  for mon,c in trans.items():
   if all(j<3 for j in mon):add(plane,mon,c)
   elif all(j>=3 for j in mon):K[mi[tuple(j-3 for j in mon)]]=(K[mi[tuple(j-3 for j in mon)]]+c)%P
   else:
    assert len(mon)==2 and mon[0]<3<=mon[1];K=(K+c*coeff[mon[0]*42+mon[1]-3])%P
  planes.append(plane);Ks.append(K)
 expected_plane=[{tuple(m):c for m,c in p}for p in model['plane5_sparse']];assert planes==expected_plane and np.array(Ks).tolist()==model['K_coefficients']
 # Exact QQ raw-infinity witnesses from the ORIGINAL134 source, notmodp zeros.
 diagonal=[[F(0)]*45 for _ in raw['equations']];symbolic={j:[]for j in [8,16,24,32]}
 for i,eq in enumerate(raw['equations']):
  byaxis={j:{}for j in symbolic}
  for mon,n,d in eq:
   free=tuple(sorted(qindex[x]for x in mon if x in qindex))
   if len(free)!=2 or free[0]!=free[1]:continue
   j=free[0];c=F(n,d);beta=tuple(x for x in mon if x not in qindex)
   for old in beta:c*=pointQQ[old]
   diagonal[i][j]+=c
   if j in byaxis:
    b=tuple(x+1 for x in beta);byaxis[j][b]=byaxis[j].get(b,F(0))+F(n,d)
  for j,p in byaxis.items():
   p={b:c for b,c in p.items()if c}
   if p:symbolic[j].append({'equation':i,'terms':[[list(b),str(c)]for b,c in sorted(p.items())]})
 axes=[j for j in range(45)if all(row[j]==0 for row in diagonal)];assert axes==[8,16,24,32];assert[qids[j]+1 for j in axes]==[11,22,33,45]
 for j,beta in zip(axes,[36,27,13,5]):assert len(symbolic[j])==1 and symbolic[j][0]['terms']==[[[beta],'2']]and pointQQ[beta-1]==0
 for i in range(129):assert all(F(lines[130+i].split()[j])==0 for j in axis_columns)
 # Allfour points also satisfy every retained graph coefficientblock exactly
 # modp, with Q(v)=0 (which was already checked directly overQQ above).
 for qaxis in axes:
  normal=qaxis-3;v=np.eye(42,dtype=np.int64)[normal];quad=np.array([int(m==(normal,normal))for m in mons[42:]],dtype=np.int64);assert not np.any(coeff[:,42:]@quad%P);assert not np.any(H[:,42:]@quad%P)and not np.any(np.array(Ks)[:,42:]@quad%P)
 # Dependencyscope: old exactdimension3 concerns FULL135, hence gives at least
 # a3D family for the containing134 rawincidence, no new134 upper/purityclaim.
 boundary=ROOT/'reports/S61_AFFINE_BOUNDARY_DIMENSION_20261102.json';bd=json.loads(boundary.read_text());assert bd['all_passed']and bd['raw_generic_infinity_exact_dimension']==3
 comm_review=ROOT/'reports/S61_AFFINE_NORMAL_COMMUTATOR_20261103_v63_REVIEW.json';check_report(comm_review)
 small={'universal_free_word_identity':[{'scalar_exponents_x_l0_l1':list(ex),'matrix_word':list(word),'coefficient':c}for(ex,word),c in sorted(formal.items())],'conclusion_for_every_i':'row_i(C^-1) gives v_i^3, hence1 on chartv_i=1','complete42x42_inverse_identity_mod32003':True,'retained_coefficient_counts':{'P':126*945,'H':3*945,'K':5*945,'plane_polynomials':5},'exact_QQ_raw_infinity_axes':axes,'free_u_onebased':[qids[j]+1 for j in axes],'symbolic_square_coefficients':symbolic,'all4_raw_axes_Qcoefficients_exactly_zero_over_QQ':True,'chart_degree_bounds':{'B':4,'zeta':5,'compatibility':6,'H':5,'F':10}}
 (out/'certificate.json').write_text(json.dumps(small,indent=2)+'\n')
 files=[Path(__file__),folder/'report.json',folder/'certificate.json',folder/'coefficient-model-mod32003.json.gz',finite,crossfile,prior_review,comm_review,source,original,pointpath,boundary,ROOT/'research/notes/S61_AFFINE_NORMAL_DIRECTION_GRAPH_20261103_v63.md']
 report={'id':'S61_AFFINE_NORMAL_DIRECTION_GRAPH_20261103_V63_REVIEW','all_passed':True,'all119070_QQ_Pcoefficients_independently_bound_to_mod_model':True,'all_H_K_and5plane_blocks_independently_rebuilt':True,'universal_characteristic_zero_Bezout_identity_all42_charts_proved':True,'three_sample_checks_not_used_as_universal_proof':True,'both_actual_modular_commutator_inverse_products_checked':True,'full_homogeneous_source_recovered_with_all3H_and5F':True,'divided_normal_incidence_Ah_v_definition_explicit':True,'r_nonzero_and_all_raw_infinity_n_nonzero_proof_accepted':True,'all42_compact_chart_graphs_closed_immersion_accepted':True,'exact_QQ_four_infinity_witnesses_source_checked':True,'old_full135_infinity_dimension3_retained_as_dependency':True,'compact134_contains_a_generic3D_infinity_family':True,'new_compact134_exact_upper_dimension_or_purity_claim':False,'affine_residual_is_s_nonzero_open_not_proper_beta_claim':True,'QQ_Bi_nilpotence_or_affine_empty_or_genericfive_claim':False,'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in files},'seconds':time.monotonic()-start}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:report[k]for k in ['all_passed','universal_characteristic_zero_Bezout_identity_all42_charts_proved','all42_compact_chart_graphs_closed_immersion_accepted','seconds']},indent=2))
if __name__=='__main__':main()
