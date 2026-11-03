#!/usr/bin/env python3
"""Exact QQ finiteness of pair rank failure; no Groebner basis or field lift.

The source replay reconstructs all 804 actual coefficients and the pair chart.
Projection finiteness is certified by a literal degree33 resultant identity and
an H-content Bezout identity. t is never inverted unless already a known unit.
"""
from pathlib import Path
from fractions import Fraction
import json,hashlib,subprocess,sys,time,itertools
import sympy as s
ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'work/nonsquare-v62/pair-natural-finite-v76'
OUT=Path(sys.argv[1]).resolve()if len(sys.argv)>1 else ROOT/'work/nonsquare-v62/pair-finite-proof-v76'
OUT.mkdir(parents=True,exist_ok=True)

def main():
 start=time.monotonic();a,e,P,t,v,H=s.symbols('a e P t v H')
 with(OUT/'actual-source-replay.log').open('w')as log:
  subprocess.run([sys.executable,str(ROOT/'code/s61_nonsquare_pair_gate_frontier_20261103_v63.py'),str(OUT/'actual-source')],check=True,stdout=log)
 source=OUT/'actual-source/pair-after.txt';it=iter(source.read_text().split());n,m,k=[int(next(it))for _ in range(3)];assert(n,m,k)==(2,6,0);labels=[next(it)for _ in range(n+k)];old=[]
 for i in range(n):
  row=[]
  for j in range(m):
   nt=int(next(it));q=0
   for z in range(nt):
    c=s.Rational(next(it));ex=[int(next(it))for _ in range(3)];q+=c*a**ex[0]*e**ex[1]*P**ex[2]
   row.append(q)
  old.append(row)
 natural=[[s.Poly(s.cancel(q.subs(a,t*P).subs(e,(v-20*t)/3)/P**4).subs(P,H-8*t),t,v,H,domain=s.QQ)for q in row]for row in old]
 h0=(160*t-15*H-2*v+31)/36;rows=[natural[0],[s.Poly(q3.as_expr()-h0*q0.as_expr(),t,v,H,domain=s.QQ)for q0,q3 in zip(*natural)]]
 f=(v+3)*(3*v+4);E=3*H*(v+4)+5*v-12
 assert s.expand(rows[0][5].as_expr()-27*f/40000)==0
 assert s.expand(rows[1][5].as_expr()+3*E/32000)==0
 # Three necessary proportionality equations are linear in t after removing t-units.
 equations=[];normalizations=[];ab=[]
 for j in [2,3,4]:
  q=rows[0][5]*rows[1][j]-rows[1][5]*rows[0][j];power=0
  while not q.is_zero:
   quo,rem=q.div(s.Poly(t,t,v,H,domain=s.QQ))
   if rem:break
   q=quo;power+=1
  content,q=q.primitive();pt=s.Poly(q.as_expr(),t);assert pt.degree()<=1
  alpha=s.Poly(pt.coeff_monomial(t),H,v,domain=s.QQ);beta=s.Poly(pt.coeff_monomial(1),H,v,domain=s.QQ)
  ab.append((alpha,beta));equations.append(q);normalizations.append(dict(pair_coefficient_index=j,t_unit_power=power,QQcontent=str(content)))
 elim=[];elim_normalizations=[]
 for i,j in itertools.combinations(range(3),2):
  q=ab[i][0]*ab[j][1]-ab[j][0]*ab[i][1];content,q=q.primitive();fpower=0;fp=s.Poly(f,H,v,domain=s.QQ)
  while not q.is_zero:
   quo,rem=q.div(fp)
   if rem:break
   q=quo;fpower+=1
  elim.append(q);elim_normalizations.append(dict(indices=[i,j],QQcontent=str(content),f_chart_power=fpower))
 # Every simultaneous zero of all t-coefficients on the necessary equations
 # would satisfy beta3/f=0, -beta4/f=0, alpha4=0; a literal QQ identity excludes it.
 alpha_cert=json.loads((DATA/'alpha-nonzero-certificate.json').read_text())
 parsed={key:s.Poly(s.sympify(alpha_cert[key]),H,v,domain=s.QQ)for key in ['C','beta3_div_f','alpha4','Bezout_coefficient_C','Bezout_coefficient_beta3_div_f','Bezout_coefficient_alpha4']}
 assert parsed['C']==s.Poly(s.cancel(-ab[2][1].as_expr()/f),H,v,domain=s.QQ)
 assert parsed['beta3_div_f']==s.Poly(s.cancel(ab[1][1].as_expr()/f),H,v,domain=s.QQ)
 assert parsed['alpha4']==ab[2][0]
 assert sum((parsed[c]*parsed[q]for c,q in [('Bezout_coefficient_C','C'),('Bezout_coefficient_beta3_div_f','beta3_div_f'),('Bezout_coefficient_alpha4','alpha4')]),s.Poly(0,H,v,domain=s.QQ))==s.Poly(1,H,v,domain=s.QQ)
 # Nonzero R(v) and content(A in H)=1 show that the projection is finite.
 projection=json.loads((DATA/'finite-projection-certificate.json').read_text());load=lambda terms:s.Poly(sum(s.Rational(c)*H**ex[0]*v**ex[1]for ex,c in terms),H,v,domain=s.QQ)
 AA,BB,RR,CA,CB=[load(projection[key])for key in ['A','B','R','Bezout_A','Bezout_B']]
 assert AA==elim[1] and BB==elim[2] and not RR.is_zero and RR.degree(H)==0 and RR.degree(v)==33
 assert CA*AA+CB*BB==RR
 loadv=lambda terms:s.Poly(sum(s.Rational(c)*v**ex[0]for ex,c in terms),v,domain=s.QQ)
 actual_coeff=[s.Poly(c,v,domain=s.QQ)for c in s.Poly(AA.as_expr(),H).all_coeffs()]
 saved_coeff=list(map(loadv,projection['A_H_coefficients']));weights=list(map(loadv,projection['A_content_Bezout']));assert actual_coeff==saved_coeff
 assert sum((c*q for c,q in zip(weights,actual_coeff)),s.Poly(0,v,domain=s.QQ))==s.Poly(1,v,domain=s.QQ)
 # On f=0, a nonzero first row forces E=0, leaving two explicit H values.
 branches=[]
 for vv,hh in [(s.Integer(-3),s.Integer(9)),(-s.Rational(4,3),s.Rational(7,3))]:
  assert s.Poly(rows[0][5].as_expr().subs(v,vv),H,t).is_zero
  assert s.solve(E.subs(v,vv),H)==[hh]
  rr=[[s.Poly(q.as_expr().subs({v:vv,H:hh}),t,domain=s.QQ)for q in row]for row in rows];minors=[];normalized=[];factor_records=[]
  for i,j in itertools.combinations(range(5),2):
   q=rr[0][i]*rr[1][j]-rr[0][j]*rr[1][i];minors.append(q);powers=[]
   for unit in [t,hh-8*t]:
    power=0;up=s.Poly(unit,t,domain=s.QQ)
    while not q.is_zero:
     quo,rem=q.div(up)
     if rem:break
     q=quo;power+=1
    powers.append(power)
   content,q=q.primitive();normalized.append(q);factor_records.append(dict(powers_t_P=powers,QQcontent=str(content)))
  g=s.Poly(0,t,domain=s.QQ);cofs=[]
  for q in normalized:
   if g.is_zero and q.is_zero:cofs.append(s.Poly(0,t,domain=s.QQ));continue
   u,z,h=s.gcdex(g,q);cofs=[s.Poly(c*u,t,domain=s.QQ)for c in cofs]+[s.Poly(z,t,domain=s.QQ)];g=s.Poly(h,t,domain=s.QQ)
  assert g==s.Poly(1,t,domain=s.QQ)
  assert sum((c*q for c,q in zip(cofs,normalized)),s.Poly(0,t,domain=s.QQ))==g
  branches.append(dict(v=str(vv),H=str(hh),raw_minors=[str(q.as_expr())for q in minors],normalized_minors=[str(q.as_expr())for q in normalized],normalization_records=factor_records,Bezout_coefficients=[str(q.as_expr())for q in cofs]))
 enc=lambda q,vs:[[list(ex),str(c)]for ex,c in s.Poly(q,*vs,domain=s.QQ).terms()]
 certificate=dict(actual_pair_source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),natural_coordinate_substitution={'a':'t*(H-8*t)','e':'(v-20*t)/3','P':'H-8*t','J':'5H-2v-4'},row_unit='P^4',Q3_row_subtraction=str(h0),natural_rows=[[enc(q.as_expr(),(t,v,H))for q in row]for row in rows],linear_t_equations=[enc(q.as_expr(),(t,v,H))for q in equations],linear_equation_normalizations=normalizations,elimination_polynomials=[enc(q.as_expr(),(H,v))for q in elim],elimination_normalizations=elim_normalizations,alpha_nonzero_Bezout=alpha_cert,projection_certificate=projection,special_branches=branches)
 (OUT/'certificate.json').write_text(json.dumps(certificate,indent=2)+'\n')
 result=dict(all_exact_QQ_actual_source_checks_passed=True,pair_rank_failure_dimension_at_most_zero=True,open='a*P*J!=0',special_f_zero_branches_excluded=True,general_projection_finite_degree33_resultant=True,H_content1_exact_Bezout=True,no_all_t_coefficients_zero_exact_Bezout=True,at_most_one_t_per_projected_pair=True,raw24_defect_finite_using_prior_gcd_and_binary_packets=True,raw24_consequence_open='a*P*F*B*C*J!=0',whole_open_pair_unit_or_raw24_unit_not_proved=True,original17_kernel_bridge_not_reaudited_here=True,no_original_product41_or_full57_or_ambient_claim=True,seconds=time.monotonic()-start)
 (OUT/'report.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)
if __name__=='__main__':main()
