#!/usr/bin/env python3
"""Independent corrected-lift genuine R0 action proof, no new action dots."""
import hashlib,json,pathlib,resource,signal,time
from fractions import Fraction as Q
ROOT=pathlib.Path(__file__).resolve().parents[1]
PRI='AMBIENT_B18_C_GENUINE_R0_TRANSPORT_SOURCE_20261105_v63';ID=PRI+'_REVIEW'
PRIMARY='555070524491511c804a2a2c8edc38f1dceea0d92333ac08887f566f634b6e07'
def sha(p):return hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()
def save(p,d):pathlib.Path(p).write_text(json.dumps(d,indent=2)+'\n')
def decode(a):return {tuple(k):Q(q) for k,q in a}
def add(a,b):
 d=dict(a)
 for k,v in b.items():
  d[k]=d.get(k,Q(0))+v
  if not d[k]:del d[k]
 return d
def scale(a,c):return {k:v*c for k,v in a.items() if v*c}
def mm(A,B):return [[sum((A[i][k]*B[k][j] for k in range(len(B))),Q(0)) for j in range(len(B[0]))] for i in range(len(A))]
def main():
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30))
 out=ROOT/'reports'/ID;out.mkdir(exist_ok=False);pins={}
 def pin(p,h=None):
  p=pathlib.Path(p);v=sha(p)
  if h:assert v==h,(str(p),v,h)
  pins[str(p.relative_to(ROOT))]=v
  return json.loads(p.read_text()) if p.suffix=='.json' else None
 pr=pin(ROOT/'reports'/PRI/'report.json',PRIMARY);docs={}
 for p,h in pr['input_sha256'].items():
  d=pin(ROOT/p,h);docs[pathlib.Path(p).name]=d
  if p.endswith('_REVIEW/report.json'):assert d['all_passed']
 D=pin(ROOT/'reports'/PRI/'genuine-four-R0-action-transport-and-cached-HB-operator-identities.json',
  pr['evidence_sha256']['genuine-four-R0-action-transport-and-cached-HB-operator-identities.json'])
 pin(ROOT/'code/ambient_B18_C_genuine_R0_transport_source_20261105_v63.py',pr['code_sha256'])
 pin(ROOT/'research/notes'/f'{PRI}.md',pr['note_sha256'])
 formal=docs['full-HB6-by2-affine-source-operator-and-old27-leaves.json'];P=[formal['P3_A'],formal['P3_B']]
 # Reconstruct the normalized cubic curve coefficients independently from
 # the full quartic reduction. No coefficient-leaf pairing is performed.
 N4=[[Q(q) for q in row] for row in formal['full_NF4']];ops=[];corrections=[]
 for side in range(2):
  a1=-N4[3*side][5*side+1];a2=a1*a1-N4[3*side][5*side];a3=-N4[3*side+2][5*side+1]
  n0=[];n1=[];n2=[]
  for k in range(3):
   p=[decode(x) for x in P[side][k]]
   n0.append(add(scale(p[0],-a1),p[1]))
   n1.append(add(scale(p[0],-a2),p[2]))
   n2.append(add(scale(p[0],-a3),p[3]))
  for j in [1,2]:
   # q(r*k) primitive contribution, then the actual Y2 free correction.
   alpha=scale(n0[j],-Q(1))
   response=add(add(scale(n0[j-1],-Q(1)),scale(alpha,a1)),n1[j])
   expected=scale(decode(formal['full_HB_6_by2_affine_operators'][3*side][j-1]),-Q(1))
   assert response==expected and all(k[0]!='mu03' for k in response)
   record=next(x for x in D['four_r_action_operator_identities'] if x['side']==side and x['low_generator']==j-1)
   assert response==decode(record['operator'])
   y=next(x for x in D['free_degree2_correction_operators_UNPAIRED'] if x['side']==side and x['low_generator']==j-1)
   assert scale(n1[j],-Q(1))==decode(y['r_ell_free_coefficient_c'])
   assert scale(n2[j],-Q(1))==decode(y['ell2_free_coefficient_d'])
   ops.append([side,j-1,'exactly minus accepted HB operator'])
   corrections.append([side,j-1,'Y2 includes both r ell and ell² terms; mu retained in ell² term'])
 frame=docs['actual-first21-quadratic-frame-and-twelve-cubic-zero-sensitivity-proof.json']
 raw=docs['actual-ell2-product33-and-raw45-actions.json']['actual33_basis']
 C={x['physical_row']:x for x in docs['actual-physical12-and-C2-right-exact-quotient.json']['physical12_coefficients']}
 ell={x['physical_row']:x for x in docs['all164-cached-pure-ell-readbacks-and-physical12-second-values.json']['all12_actual_physical_coefficients_pure_ell3']}
 zeros=[]
 for i,record in enumerate(frame['all12_verified_source_records'][:7]):
  p=record['physical_row'];assert record['delta_first9_zero']
  assert list(map(Q,record['bA_first_z7']))==[Q(i==j) for j in range(7)]
  for kind,ellkind in [('delta','delta_jet3'),('bA','bA_jet3')]:
   values=dict(C[p][kind])
   for mon in [[0,2,0],[0,1,1]]:assert not Q(values.get(raw.index(mon),'0'));zeros.append([p,kind,mon,'0'])
   assert not Q(ell[p][ellkind][2]);zeros.append([p,kind,[0,0,2],'0'])
 assert len(zeros)==42
 HB=docs['phase2-full-six-coordinate-HB-matrix-and-fixed-rank-certificate.json']['actual_full_HB6_by2_matrix']
 B=[[-Q(q) for q in HB[i]] for i in [0,3]]
 assert B==[[Q(3,5),-Q(6)],[Q(121,80),-Q(121,8)]]
 R=[[Q(0)]*4 for _ in range(4)];R[2][:2]=B[0];R[3][:2]=B[1]
 assert [[str(q) for q in row] for row in R]==D['actual_r_matrix4']
 assert mm(R,R)==[[0]*4 for _ in range(4)] and all(row[1]==-10*row[0] for row in B)
 T=[[Q(q) for q in row] for row in D['basis_change4']]
 I=[[Q(q) for q in row] for row in D['basis_change_inverse4']]
 ident=[[Q(i==j) for j in range(4)] for i in range(4)]
 assert mm(T,I)==ident==mm(I,T)
 J=[[Q(0)]*4 for _ in range(4)];J[1][0]=1;assert mm(mm(I,R),T)==J
 proof={'actual_genuine_target_not_naive_Dquotient_action':True,
  'corrected_lift_proof':[
   'The low kernel-coordinate representatives have ell image inF16, starting in degree3 because H quadratics vanish in gr2 and the T7 cubic coordinates cancel by alpha/beta.',
   'The degree3 free image has c*r ell²+d*ell³ on eachside. Subtract the actual degree2 free Y2=c*r ell+d*ell², retaining BOTH coefficients.',
   'The remaining ell image lies in m4C4, whose literalbasis is r ell³A/B. It has a free degree3 preimageY3 underell. SubtractingY3 gives an actual finitekerE23 representative. Its coefficients stay UNKNOWN.',
   'mY3 lies in m4C4 subset imE³. Thus Y3 affects no genuineR0action. Highkernelclasses have order3 corrected representatives and are similarly killed bym.',
   'When computing q(r*the complete corrected representative), r*dell² terms lie inF and q kills them; they need not be boundaries individually. The whole corrected expression is a cycle. This avoids falsely treating q asR0-linear.'
  ],
  'four_affine_action_identities_checked_coefficientwise':ops,
  'all_Y2_components_retained':corrections,
  'all42_literal_plane2_zero_bindings':zeros,
  'z_annihilation_proof':[
   'Allz targetclasses lie inm2 by the first21unitframe. Verifieddelta_i2 andbAi2 have no pureplane term; RHi=(0,delta_i,bAi),bAi1=zi thenimpliesziH in m3.',
   'For low corrected heads, zi*(rH/ellH) lies inm4. Their r²A/B terms also land inm4 sinceziA=0 andziB in m2.',
   'For the SPECIFIC free degree2 Y2 terms ell²A/B,r ellA/B, ziY2 lies inm4 byziA=0,ziB in m2. This is not a claim about an arbitrary raw degree2 vector.',
   'UnknownY3 terms giveziY3 in m4. All these contributions die in imE³; therefore eachzi acts0 on genuineC[ell]. Ell acts0 by definition.'
  ],
  'actual_r_low_to_high_matrix':[[str(q) for q in row] for row in B],
  'r_rank1_and_squarezero':True,'all_z_and_ell_zero':True,
  'Jordan_basis_and_both_inverse_products_verified':True,
  'actual_R0_type':'R0/(z0,...,z6,ell,r²) plus k²','actual_mu_R0':3,'socle_dimension':3,'m2_annihilator':True,
  'no_new_r_action_leaf_pairing':True,
  'unknown_Y3_full_C4_action_tables_and_connecting_values_uncomputed':True,
  'no_original_H0_nil_cube_or_index_transfer':True}
 ep=out/'independent-corrected-lift-R0-action-and-all-affine-HB-identities.json';save(ep,proof)
 report={'id':ID,'all_passed':True,'frozen':True,'source_only':True,'input_sha256':pins,'code_sha256':sha(__file__),
  'evidence_sha256':{ep.name:sha(ep)},'actual_genuine_R0_J2_plus_k2_transport_accepted':True,
  'no_new_action_dots_or_connecting_values':True,'elapsed_seconds':time.monotonic()-start}
 save(out/'report.json',report);signal.alarm(0)
 print(json.dumps({'all_passed':True,'report':str(out/'report.json'),'sha256':sha(out/'report.json'),'seconds':report['elapsed_seconds']}))
if __name__=='__main__':main()
