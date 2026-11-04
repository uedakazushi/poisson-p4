#!/usr/bin/env python3
"""Bounded EXACT rational recovery/cancellation only; no elimination or new points."""
import argparse,gzip,hashlib,json,resource,signal,time
from pathlib import Path
import sympy as S
ROOT=Path(__file__).resolve().parents[1]
STEM='Q_NEW_FAMILY_TWO_SHAPE_RECOVERY_20261104_v63'
def read(p):
 with(gzip.open(p,'rt')if p.suffix=='.gz'else p.open())as f:return json.load(f)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main(out):
 started=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));signal.alarm(60);out.mkdir(parents=True,exist_ok=False)
 r,s,q,B,R,Z0,Z1,Z2,ZN=S.symbols('r s q B R Z0 Z1 Z2 ZN');h=s*s-25*r*q*q;F=100*r*(r-1)*q*q-4*r*s+s*s+4*r
 base=ROOT/'reports/Q_NEW_FAMILY_LOW_INVARIANT_BIRATIONAL_COMPRESSION_20261104_v63';review=ROOT/'reports/Q_NEW_FAMILY_LOW_INVARIANT_BIRATIONAL_COMPRESSION_20261104_v63_REVIEW/report.json'
 assert read(review)['all_passed'];n={};maps={};paths=[];frontier=[]
 def savefront(stage):
  frontier.append(dict(stage=stage,seconds=time.monotonic()-started));(out/'frontier.json').write_text(json.dumps(frontier,indent=2)+'\n')
 for name in['F4','p0','p1','p2','p3','p4','N12']:
  path=base/(name+'.json.gz');paths.append(path);d=read(path)
  n[name]=S.Poly.from_dict({tuple(ex):S.Rational(z)for z,*ex in d['terms']},(r,s,q,B)).as_expr()/S.Integer(d['denominator_integer'])
  maps[name]=n[name]/h**d['denominator_h_power']
 assert S.expand(maps['F4']-S.Rational(2,405)*F)==0
 p=[maps['p'+str(i)]for i in range(5)];f=maps['F4'];L0=S.Poly(n['p0'],B).coeff_monomial(B);D=2359296000*r**4*(r-1)*q*q*(-s+5*q+1)
 assert S.expand(L0-S.Rational(2,405)**2*D)==0
 assert S.Poly(n['p0'],B).coeff_monomial(1)==0
 assert S.expand(n['p3'].subs(B,0)+4*n['p4'].subs(B,0))==0
 Rmap=S.cancel((p[3]+4*p[4])/p[0]);Rexpected=(80*q*r*r-6*r*r*s+6*r*r-11*r*s+10*r+s)/(8*r*(r-1)*(5*q-s+1));assert S.cancel(Rmap-Rexpected)==0
 T=(r-1)*R-2*r;L=8*R*r*(r-1)-6*r*r-11*r+1;Nq=(s-1)*L-(r-1);qback=Nq/(40*r*T)
 assert S.cancel(Rmap.subs(q,qback)-R)==0 and S.cancel(qback.subs(R,Rmap)-q)==0
 J=(10*r-1)*s-10*r
 assert S.cancel((5*q-s+1).subs(q,qback)-(r-1)*J/(8*r*T))==0
 Bback=Z0*h*F**2/D;assert S.cancel((p[0]/f**2).subs(B,Bback)-Z0)==0
 b0=S.cancel(Bback.subs(Z0,p[0]/f**2));assert S.cancel(b0-B)==0
 savefront('B_and_Mobius_q_two_way_recovery_verified')
 # Constant-in-B linear relations are metadata: no implicit ideal or GB is built.
 consts=[n['p'+str(i)].subs(B,0)for i in range(1,5)]+[S.expand(f*f*h)]
 mon=sorted(set().union(*(S.Poly(a,r,s,q).monoms()for a in consts)));M=S.Matrix([[S.Poly(a,r,s,q).coeff_monomial(m)for a in consts]for m in mon]);null=M.nullspace();assert M.rank()==3
 const_kernel=[list(map(str,v))for v in null]
 def rational_record(expr,variables):
  value=S.cancel(expr);num,den=S.fraction(value);pn,pd=S.Poly(num,*variables),S.Poly(den,*variables)
  assert S.gcd(pn,pd).total_degree()==0
  return value,dict(variables=list(map(str,variables)),numerator=str(num),denominator=str(den),
   numerator_terms=len(pn.terms()),denominator_terms=len(pd.terms()),numerator_multidegree=list(pn.degree_list()),denominator_multidegree=list(pd.degree_list()))
 # Cancel3-variable coefficient functions BEFORE combiningZ0, keeping the job small.
 coefficients={};values={};tables={}
 for i in[1,2]:
  Ai=S.Poly(n['p'+str(i)],B).coeff_monomial(1);Li=S.Poly(n['p'+str(i)],B).coeff_monomial(B)
  a,ar=rational_record((Ai/(h*f*f)).subs(q,qback),(r,s,R));b,br=rational_record((Li/L0).subs(q,qback),(r,s,R))
  coefficients['Z'+str(i)]=dict(constant=ar,Z0_coefficient=br)
  value,record=rational_record(a+Z0*b,(r,s,R,Z0));values['Z'+str(i)]=value;tables['Z'+str(i)]=record
  # Exact pre-q-substitution identity binds the whole short formula to the source.
  assert S.cancel((p[i]/f**2).subs(B,Bback)-(Ai/(h*f*f)+Z0*Li/L0))==0
  savefront('normalized_p'+str(i)+'_two_shape_function_completed')
 # p4 determined by the universal identity; p3 byR andp4. No extra shape equation.
 assert S.cancel((p[4]/f**2).subs(B,Bback)-(S.Rational(1,1600)-Z0/4))==0
 assert S.cancel((p[3]/f**2).subs(B,Bback)-(Z0*(Rmap+1)-S.Rational(1,400)))==0
 M0=S.Poly(n['N12'],B).coeff_monomial(B);M1=S.Poly(n['N12'],B).coeff_monomial(B**2);assert S.Poly(n['N12'],B).coeff_monomial(1)==0
 a,ar=rational_record((M0/(h*f*L0)).subs(q,qback),(r,s,R));b,br=rational_record((f*M1/L0**2).subs(q,qback),(r,s,R))
 coefficients['ZN']=dict(Z0_coefficient=ar,Z0_squared_coefficient=br)
 value,record=rational_record(Z0*a+Z0**2*b,(r,s,R,Z0));values['ZN']=value;tables['ZN']=record
 assert S.cancel((maps['N12']/f**3).subs(B,Bback)-(Z0*M0/(h*f*L0)+Z0**2*f*M1/L0**2))==0
 savefront('normalized_N12_two_shape_function_completed')
 constraints={};summary=[]
 for name,target in[('Z1',Z1),('Z2',Z2),('ZN',ZN)]:
  num,den=S.fraction(values[name]);poly=S.Poly(S.expand(num-target*den),r,s,R,Z0,target)
  constraints[name]=dict(variables=list(map(str,poly.gens)),polynomial=str(poly.as_expr()),terms=len(poly.terms()),multidegree=list(poly.degree_list()),
   shape_bidegree=list(poly.degree_list()[:2]),denominator_localizer=str(den))
  summary.append(dict(name=name,terms=len(poly.terms()),shape_bidegree=list(poly.degree_list()[:2]),parameter_degrees=dict(zip(map(str,poly.gens[2:]),poly.degree_list()[2:]))))
 # All additional recovery localizers are displayed; no assertion of extension on theirzeros.
 hR=S.cancel(h.subs(q,qback));FR=S.cancel(F.subs(q,qback));DR=S.cancel(D.subs(q,qback));DRexpected=184320*r*(r-1)**2*Nq**2*J/T**3;assert S.cancel(DR-DRexpected)==0
 recipe=dict(normalized_target=dict(Z0='p0/F4^2',R='(p3+4p4)/p0',Z1='p1/F4^2',Z2='p2/F4^2',ZN='N12/F4^3'),
  R_expression=str(Rexpected),T=str(T),L=str(L),Nq=str(Nq),q_recovery=str(qback),B_recovery_before_q=str(Bback),
  recovered_minus_s_plus5q_plus1=str((r-1)*J/(8*r*T)),D_after_q=str(DRexpected),h_after_q=str(hR),F_after_q=str(FR),
  known_tau_source_field_reparametrized='QQ(xi,r,s,q,B)=QQ(xi,r,s,R,Z0) on displayed dense unitopen',
  base_localizers=['xi','r','s','Z0','T','D_after_q','h_after_q','F_after_q'],
  function_denominator_localizers={name:rec['denominator']for name,rec in tables.items()},
  scope='rational dense source chart only; retain ALL displayed nonzero denominators, no global extension/properness/index1 assertion')
 result=dict(recipe=recipe,coefficient_functions=coefficients,normalized_two_shape_functions=tables,cleared_constraints=constraints,
  constant_B0_relation_columns=['p1','p2','p3','p4','F4^2'],constant_B0_rank=3,constant_B0_kernel=const_kernel)
 (out/'certificate.json').write_text(json.dumps(result,indent=2)+'\n');savefront('all_three_constraints_saved_no_elimination')
 report=dict(id=STEM,all_passed=True,exact_B_cancellation_and_two_way_Mobius_q_recovery=True,remaining_unknown_shape_variables=['r','s'],base_target_parameters=['Z0','R'],
  normalized_constraint_summary=summary,no_Groebner_resultant_elimination_or_source26_or_new_points=True,
  known_tau_source_reparametrization_only_not_actual_target_field_index1=True,
  generic_finiteness_for_two_spectral_constraints_follows_existing_spectral_base_dominance=True,
  no_actual_degree_genus_or_whole_chart_properness_claim=True,
  limits=dict(cpu_and_wall_seconds=60,address_space_bytes=2_000_000_000),
  input_sha256={str(p.relative_to(ROOT)):sha(p)for p in[Path(__file__),review,*paths]},
  evidence_sha256={p.name:sha(p)for p in out.iterdir()if p.is_file()},seconds=time.monotonic()-started)
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');signal.alarm(0)
 print(json.dumps({k:report[k]for k in['all_passed','normalized_constraint_summary','seconds']},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
