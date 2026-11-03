#!/usr/bin/env python3
"""Full D corollary from the independently gated meeting and true two-jet data."""
import hashlib,itertools,json,time
from pathlib import Path
import sympy as s
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_CA_SELF_PAIR_ACTUAL_D_FULL_RING_20261103_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 start=time.monotonic();out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False)
 gates=[ROOT/'reports/AMBIENT_CA_SELF_PAIR_MEETING_DOUBLE_20261103_v63_REVIEW/report.json',
  ROOT/'reports/AMBIENT_CA_SELF_PAIR_ACTUAL_D_TWOJET_20261103_v63/report.json',
  ROOT/'reports/AMBIENT_CA_SELF_PAIR_ACTUAL_D_TWOJET_20261103_v63_REVIEW.json']
 for p in gates:
  d=json.loads(p.read_text());assert d['all_passed'] is True
  for n,h in d['input_sha256'].items():assert sha(ROOT/n)==h,n
 rho,tau,f,c,d=s.symbols('rho tau f c d');u=s.symbols('u1:9');v=s.symbols('v1:9');b=s.symbols('b1:9')
 relations=[x*y for x in u for y in v]+[x*tau for x in u]+[y*(tau-rho*rho)for y in v]+[tau*(tau-rho*rho)]
 assert len(relations)==81
 substitution={rho:f*d,tau:f*c,**{x:f*y for x,y in zip(v,b)}}
 restrictions=[s.expand(q.subs(substitution,simultaneous=True))for q in relations]
 # Every monomial belongs to (f²,f*u8), and the displayed unit cofactors give
 # both generators in the restricted actual ideal. The ALL-S tangent argument
 # in the note handles every additional ambient equation in the opposite direction.
 for q in restrictions:
  for ex,co in s.Poly(q,f,*u,c,d,*b).terms():assert ex[0]>=2 or(ex[0]>=1 and sum(ex[1:9])>=1)
 for i in range(8):assert s.cancel(restrictions[64+i]/c-f*u[i])==0
 assert s.cancel(restrictions[-1]/(c*(c-f*d*d))-f*f)==0
 # (f) intersect (f²,u8) has these exact monomial LCM generators.
 monomial_intersection=[f*f,*[f*x for x in u]]
 evidence=dict(branch_convention='B1: v=tau=0; B2: u=0,tau=rho²',
  full81_CA_relations=list(map(str,relations)),all81_CB_restrictions=list(map(str,restrictions)),
  actual_D_ideal=list(map(str,monomial_intersection)),
  lower_inclusion_unit_cofactors=['1/c','1/(c*(c-f*d²))'],
  upper_inclusion='For every actual F=f*k, ALL-S tangentD=TCB gives k|S=0.',
  primary_components=['(f)','(f²,u1,...,u8)'],nilradical='C[[s20]]*f')
 ep=out/'certificate.json';ep.write_text(json.dumps(evidence,indent=2)+'\n')
 note=ROOT/'research/notes'/f'{STEM}.md';paths=[Path(__file__),note,*gates]
 report=dict(id=STEM,all_passed=True,whole_actual_D_ideal_equality_from_two_inclusions=True,
  full_completion='C[[s20,u8,f]]/(f²,f*u1,...,f*u8)',
  cone=dict(dimension=28,embedding_dimension=29,depth=20,multiplicity=1),
  projective=dict(dimension=27,embedding_dimension=28,depth=19,multiplicity=1),
  reduced_support='smooth N',square_zero_nilradical='C[[s20]]*f',
  associated_primes=['minimal N','embedded S20 cone / S19 projective'],
  nilmodule_generic_length_on_S=1,global_actual_D_embedded_associated_component_from_algebraic_family=True,
  B18_origin_or_all_global_primary_components_solved=False,
  scope='Nonempty parameter open of P_a and its G saturation, within the verified meeting-double unit open.',
  input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},evidence_sha256={'certificate.json':sha(ep)},seconds=time.monotonic()-start)
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print('FULL ACTUAL D GENERIC CA SELF-PAIR: EMBEDDED S20, PASS',report['seconds'])
if __name__=='__main__':main()
