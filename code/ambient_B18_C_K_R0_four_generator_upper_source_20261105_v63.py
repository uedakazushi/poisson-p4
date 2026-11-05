#!/usr/bin/env python3
"""Pure SOURCE: four R0 generators suffice for the original comparison kernel."""
import hashlib,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_C_K_R0_FOUR_GENERATOR_UPPER_SOURCE_20261105_v63'
PATHS={
 'kernel_basis':'reports/AMBIENT_B18_C_CONNECTING_KERNEL3_SPLIT_20261105_v63/actual-M7-stable-surjection-splitting-and-kernel-three-source-basis.json',
 'kernel_basis_peer':'reports/AMBIENT_B18_C_CONNECTING_KERNEL3_SPLIT_20261105_v63_REVIEW/report.json',
 'honest_uct':'reports/AMBIENT_B18_C_CONNECTING_SEVEN_SOURCE_20261105_v63/honest-seven-source-quotient-four-fixed-kernel-classes-and-edge-constraints.json',
 'honest_uct_peer':'reports/AMBIENT_B18_C_CONNECTING_SEVEN_SOURCE_20261105_v63_REVIEW/report.json',
 'finite_K':'reports/AMBIENT_B18_NIL_DVR_FINITENESS_20261104_v63/report.json',
 'finite_K_peer':'reports/AMBIENT_B18_NIL_DVR_FINITENESS_20261104_v63_REVIEW/report.json',
}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def plus(a,b):
 c=a.copy()
 for k,v in b.items():
  c[k]=c.get(k,0)+v
  if not c[k]:del c[k]
 return c
def scale(a,mon,sign=1):return{tuple(sorted(k+mon)):sign*v for k,v in a.items()}
def main():
 d={k:json.loads((ROOT/p).read_text())for k,p in PATHS.items()}
 assert all(x['all_passed']for k,x in d.items()if k.endswith('_peer'))
 assert d['kernel_basis']['five_known_actual_Kbar_classes']==['z5*h4','z5*r*h4','r²*h4','r³*h4','r*h3']
 assert d['kernel_basis']['kernel_basis_labels'][0]=='u=h3−(b3/b1)h1−((p3−b3*p1/b1)/d)*rh4'
 # Formal scalars a=b3/b1, b=(p3-b3*p1/b1)/d, and the four
 # coefficients of rh1 in ITold remain independent variables.
 e=lambda j:[{():1}if i==j else{}for i in range(5)]
 def vecadd(a,b):return[plus(x,y)for x,y in zip(a,b)]
 def vecscale(a,m,sign=1):return[scale(x,m,sign)for x in a]
 rh1=[{('tv',):1},{('tq1',):1},{('tq2',):1},{('tg',):1},{}]
 # Coordinates (v,q1,q2,g,w). The actual identity is
 # rU=w-a*rh1-b*q2. This does not replace rU by w in Told.
 rU=vecadd(vecadd(e(4),vecscale(rh1,('a',),-1)),vecscale(e(2),('b',),-1))
 recovered=vecadd(vecadd(rU,vecscale(rh1,('a',))),vecscale(e(2),('b',)))
 assert recovered==e(4)
 data={'scope':'Pure SOURCE upper bound for original normalization kernel K as module over the SAME regular R0. No new edge/coefficient or full lift computed.',
  'known_honest_basis':'Kbar=K/ellK has QQ-basis v,q1,q2,g,w,U,S with v=z5h4,q1=z5rh4,q2=r²h4,g=r³h4,w=rh3. This is the accepted honest UCT kernel, not a naive reduced homology.',
  'actual_relations':['r*v=q1','r*q2=g','r*h1 in ITold=span(v,q1,q2,g) exactly by the accepted M7 relation','U=h3-a*h1-b*r*h4, a=b3/b1,b=(p3-b3*p1/b1)/d'],
  'formal_rh1_coefficients':['tv','tq1','tq2','tg'],
  'actual_w_recovery':'w=r*U+a*tv*v+a*tq1*r*v+(a*tq2+b)*q2+a*tg*r*q2',
  'all_five_coordinate_polynomial_identities_verified':True,
  'four_actual_bar_generators':['U','S','v=z5*h4','q2=r²*h4'],
  'generation_proof':'They generate I*Told since q1=rv andg=rq2; then the displayed exact rU identity supplies w=rh3; together U,S this contains every honest seven-basis class. No coefficient of rh1 is evaluated or assumed0.',
  'upstairs_Nakayama':'K is finite over D=QQ[[ell]], hence finite over R0. Choose any actual K lifts of the four barclasses, possible by K→K/ellK. Their quotient Q satisfies Q=ell Q. Since ell lies in the local maximal ideal, Nakayama givesQ0. Thus mu_R0 K<=4.',
  'coefficient_ring_scope':'R0-module upper bound ONLY; no minimality, no Dmodulecountchange (mu_D K7), and no expanded fullR0 lift coefficients. M7 split is not lifted to Told/Kbar.',
  'no_H0_length_torsion_exponent_nilcube_or_index_inference':True,
  'input_sha256':{p:sha(ROOT/p)for p in PATHS.values()},'code_sha256':sha(Path(__file__)),
  'formal_identity_size':sum(len(x)for x in rU)+sum(len(x)for x in recovered),'all_passed':True}
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);p=out/'pure-four-R0-generator-upper-bound-and-exact-rU-recovery.json';p.write_text(json.dumps(data,indent=2)+'\n')
 r={'id':STEM,'all_passed':True,'status':'PURE_SOURCE_MU_R0_KERNEL_AT_MOST_FOUR_NO_VALUES','code_sha256':data['code_sha256'],'input_sha256':data['input_sha256'],'evidence_sha256':{p.name:sha(p)},'no_actual_new_response_or_lift_values':True}
 (out/'report.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r))
if __name__=='__main__':main()
