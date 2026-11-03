#!/usr/bin/env python3
"""Exact all-aPJ rank2 of the actual pair after omitting its B^2 column."""
from pathlib import Path
import sympy as s
import json,time,sys,subprocess,itertools
ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'work/nonsquare-v62/quad78-small-gates-v83/omit-B2/certificate.json'
OUT=Path(sys.argv[1]).resolve()if len(sys.argv)>1 else ROOT/'work/nonsquare-v62/pair-omit-B2-replay-v83'
OUT.mkdir(parents=True,exist_ok=True)

def main():
 start=time.monotonic();t,v,H=s.symbols('t v H')
 with(OUT/'whole-pair-source.log').open('w')as log:subprocess.run([sys.executable,str(ROOT/'code/s61_nonsquare_pair_empty_QQ_20261103_v63.py'),str(OUT/'whole-pair-source')],check=True,stdout=log)
 source=json.loads((OUT/'whole-pair-source/finite-source/certificate.json').read_text());whole=json.loads((ROOT/'work/nonsquare-v62/pair-finite-subresultant-v78/whole-pair-cleanup-certificate.json').read_text());extra=json.loads(DATA.read_text())
 load=lambda terms,vars:s.Poly(sum(s.Rational(c)*s.prod(z**k for z,k in zip(vars,ex))for ex,c in terms),*vars,domain=s.QQ)
 rows=[[load(q,(t,v,H)).as_expr()for q in row]for row in source['natural_rows']];Q=load(whole['degree23_branch']['Q'],(v,));assert Q==load(extra['degree23_Q'],(v,));Hb=load(whole['degree23_branch']['Hbar'],(v,));Tb=load(whole['degree23_branch']['Tbar'],(v,))
 def evaluate(expr):
  q=s.Poly(expr,t,H,v,domain=s.QQ);hp=[s.Poly(1,v,domain=s.QQ)];tp=[s.Poly(1,v,domain=s.QQ)];vp=[s.Poly(1,v,domain=s.QQ)]
  for _ in range(q.degree(H)):hp.append((hp[-1]*Hb).rem(Q))
  for _ in range(q.degree(t)):tp.append((tp[-1]*Tb).rem(Q))
  for _ in range(q.degree(v)):vp.append((vp[-1]*s.Poly(v,v,domain=s.QQ)).rem(Q))
  return sum(((tp[ex[0]]*hp[ex[1]]*vp[ex[2]]).rem(Q)*c for ex,c in q.terms()),s.Poly(0,v,domain=s.QQ)).rem(Q)
 actual=evaluate(rows[0][5]*rows[1][1]-rows[1][5]*rows[0][1]);assert actual==load(extra['remaining_actual_pair_minor_col1'],(v,));inv=load(extra['remaining_minor_inverse_mod_Q'],(v,));assert(actual*inv).rem(Q)==s.Poly(1,v,domain=s.QQ)
 for branch in extra['special_branches']:
  vv=s.Rational(branch['v']);hh=s.Rational(branch['H']);unit=s.Rational(branch['q0_cross_unit_coefficient']);assert unit
  assert s.expand((rows[0][3]-s.Rational(4,15)*rows[0][4]).subs(v,vv)-unit*t)==0
  rr=[[s.Poly(q.subs({v:vv,H:hh}),t,domain=s.QQ)for q in row]for row in rows];norm=[]
  for record,saved_raw,saved_primitive in zip(branch['normalization'],branch['raw_minors'],branch['primitive_minors']):
   i,j=record['columns'];raw=rr[0][i]*rr[1][j]-rr[0][j]*rr[1][i];assert raw==load(saved_raw,(t,));primitive=load(saved_primitive,(t,));p,q=record['unitpowers_t_P'];assert raw==primitive*s.Poly(s.Rational(record['QQcontent'])*t**p*(hh-8*t)**q,t,domain=s.QQ);norm.append(primitive)
  cofs=[load(q,(t,))for q in branch['Bezout_coefficients']];assert sum((c*q for c,q in zip(cofs,norm)),s.Poly(0,t,domain=s.QQ))==s.Poly(1,t,domain=s.QQ)
 result=dict(all_literal_QQ_checks_passed=True,actual_pair_retained_columns=[1,2,3,4,5],omitted_column='B^2',retained_matrix_rank2_on_every_a_P_J_point=True,open='a*P*J!=0',special_f_zero_retained_first_row_nonzero=True,degree23_actual_col1_minor_inverse_verified=True,conditional_five_quadric_linear_syzygies_exactly2_if_ternary_pair_coprime=True,ternary_pair_commonfactor_gate_not_proved=True,no_quad78_or_fullquartic_claim=True,seconds=time.monotonic()-start)
 (OUT/'report.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)
if __name__=='__main__':main()
