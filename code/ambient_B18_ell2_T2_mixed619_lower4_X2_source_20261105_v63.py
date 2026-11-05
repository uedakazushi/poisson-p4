#!/usr/bin/env python3
"""Known Gamma source rows only; no X2 or mixed619 coefficient evaluation."""
import json,math,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
from ambient_B18_ell2_T2_mu_review_20261105_v63 import fields
from ambient_B18_ell2_T2_nonactive171_source_review_20261105_v63 import sha,sp,count,add
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_ELL2_T2_MIXED619_LOWER4_X2_SOURCE_20261105_v63'
def main():
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024));inputs={}
 def read(p,keys=None):
  inputs[p]={'path':p,'sha256':sha(ROOT/p)}
  return fields(ROOT/p,keys) if keys else json.loads((ROOT/p).read_text())
 pp='reports/AMBIENT_B18_ELL2_T2_MIXED619_LOWER2_ADJOINT_20261105_v63';r=read(pp+'/report.json');assert r['all_passed'] and r['both_adjoint_identities_exist']
 ep=pp+'/actual-two-mixed-output-adjoint-certificates-or-obstruction.json';v=read(ep);assert sha(ROOT/ep)==r['evidence_sha256'][Path(ep).name]
 kp='reports/AMBIENT_B18_ELL2_T2_MIXED619_LOWER4_KNOWN_INPUT_GATE_20261105_v63/pruned-known619-lower4-input-and-prior-resource-gate.json';d=read(kp)
 m=read('reports/AMBIENT_B18_ELL2_T2_MU_INPUT_GATE_20261105_v63/actual-T2-mu-input-source-and-cache-gate.json',['Mnew0_columns'])
 M=list(map(sp,m['Mnew0_columns']));P=list(map(sp,d['Pi3_source56_columns']));assert len(M)==len(P)==56
 rho=[{},{}]+list(map(sp,v['solutions_on136_literal_outputs']));chi=[sp(x) for x in d['lower_duals2_rows3']];mono=[(7,7),(7,8),(0,7),(0,8),(4,7),(4,8),(5,7),(5,8)]
 maximum=0
 def mark(a):
  nonlocal maximum
  maximum=max(maximum,a[0].bit_length(),a[1].bit_length());return a
 def ex(x):x=Q(x);return mark((x.denominator,abs(x.numerator)))
 def times(a,b):return mark((a[0]*b[0],a[1]*b[1])) if a[1] and b[1] else (1,0)
 def plus(a,b):
  z=math.lcm(a[0],b[0]);return mark((z,a[1]*(z//a[0])+b[1]*(z//b[0])))
 def accum(a,k,b):
  if b[1]:a[k]=plus(a.get(k,(1,0)),b)
 bounds=[]
 for target in range(4):
  a={}
  for source,col in enumerate(P):
   for i,c in chi[target//2].items():accum(a,8*source+target%2,times(ex(c),ex(col.get(i,Q(0)))))
  for out,c in rho[target].items():
   row,j=divmod(out,8)
   for source,col in enumerate(M):accum(a,8*source+j,times(ex(c),ex(col.get(row,Q(0)))))
  bounds.append(a)
 prior=maximum;assert prior<=384
 owned=count([d,M,P,rho,bounds])+12000;assert owned<150000
 G=[]
 for target,bound in enumerate(bounds):
  a={}
  for source,col in enumerate(P):
   x=sum((c*col.get(i,Q(0)) for i,c in chi[target//2].items()),Q(0))
   if x:a[8*source+target%2]=x
  for out,c in rho[target].items():
   row,j=divmod(out,8)
   for source,col in enumerate(M):
    if row in col:a=add(a,{8*source+j:c*col[row]},Q(-1))
  for k,x in a.items():
   D,U=bound[k];assert D%x.denominator==0 and abs(x.numerator)*(D//x.denominator)<=U
  G.append(a)
 source='reports/AMBIENT_B18_ELL2_T2_KNOWN_FIRST_PROFILE_20261105_v63/actual-known-bottom-first-rank-and-corrected-source-kernel.json'
 c=read(source,['all_zero_bottom_and_Y0_source_combinations'])['all_zero_bottom_and_Y0_source_combinations'];assert len(c)==619 and all(x['full56_Xgamma1_equals_Xp0a1'] for x in c)
 source_records=sum(len(x['original768_combination']) for x in c);a1records=sum(len(x['actual_a1_27_times9']) for x in c)
 assert source_records==992 and a1records==1887
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);e=out/'actual-four-Gamma-X2-source-rows-and-lazy619-recipe.json'
 ev={'ordinary_second_pair_order8':mono,'Gamma4_on_source56_times8':[[[i,str(x)] for i,x in sorted(a.items())] for a in G],
 'Gamma_slot_interpretation':'8*source56row+ordinarysecondpairindex;nofactorials','all_four_lower_channels_retained':True,
 'zero_alpha_A_does_not_mean_zero_Gamma_A':True,'needed_source_rows56':sorted({i//8 for a in G for i in a}),
 'needed_second_pairs':sorted({tuple(mono[i%8]) for a in G for i in a}),
 'source619_count':619,'original768_records':source_records,'actual_a1_records':a1records,'lazy_recipe':'For eachc,a1:Gamma_t(X_c2-Xp1*a1); useordinarySym2convolutionandtheSAMEactual27pivotcolumns. No coefficient queried or replaced0.',
 'exact_identity':'G2_1Ybar1=M0Z2, Z2=X_c2-Xp1*a1. alpha=rho*G2_1 onall918slots. Correctedlower_t=chi_tPi3Z2-alpha_tYbar1=Gamma_tZ2 moduloallknownlinearproductsandQ2. M1andG2_2Y0dropbecauseZ0=Z1=Yc0=0.',
 'known_operator_prior_bits':prior,'whole_known_source_owned':owned,'actual_X2_callback_or_619_values_evaluated':False,
 'next_required_gates':['independentpostpeerofadjointidentity','sourcefoldoldB/P/nativenamespaceandfixedXpfirstcarry','newweightedhonestsecondLINsensitivities','nativeS0cNseparatezero','literalraw33quadraticfunctionaldescentforprimitiveandcanonicalgaugeinputs','wholeheight/ownershipbeforeactual619values']}
 e.write_text(json.dumps(ev,separators=(',',':'))+'\n')
 report={'id':STEM,'all_passed':True,'source_only':True,'Gamma_support4':[len(a) for a in G],'known_operator_prior_bits':prior,'independent_adjoint_review_pending':True,'input_sha256':inputs,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{e.name:sha(e)},'seconds':time.monotonic()-t};(out/'report.json').write_text(json.dumps(report,separators=(',',':'))+'\n');print(json.dumps({'all_passed':True,'support':[len(a) for a in G],'prior':prior,'neededpairs':ev['needed_second_pairs'],'neededrows':ev['needed_source_rows56'],'report_sha256':sha(out/'report.json')}))
if __name__=='__main__':main()
