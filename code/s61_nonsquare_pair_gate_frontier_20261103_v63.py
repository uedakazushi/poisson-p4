#!/usr/bin/env python3
"""Verify actual pair-gate coordinates and first safe QQ scalar pivot; no GB."""
from pathlib import Path
import sympy as s
import hashlib,json,subprocess,shlex,sys,time
ROOT=Path(__file__).resolve().parents[1]
OUT=Path(sys.argv[1]).resolve()if len(sys.argv)>1 else ROOT/'work/nonsquare-v62/pair-gate-replay-v64'
OUT.mkdir(parents=True,exist_ok=True)
def main():
 start=time.monotonic();a,e,P=s.symbols('a e P');x=s.symbols('x0:6');u,v,z=s.symbols('u v z');J=5*P-6*e-4
 source=ROOT/'reports/S61_BINARY_NONSQUARE_MOD_ELL_COORDINATES_20261103_v62.json';assert hashlib.sha256(source.read_bytes()).hexdigest()=='80f7c46b1bd4ffef4a78f492ad4f7a2bb40b2ffa27cddcd23d5b5d660db4a7e0';data=json.loads(source.read_text());qs=[]
 for i,row in enumerate(data['cleared_qfour_sparse']):
  q=0
  for mon,terms in row:q+=sum(s.Rational(co)*a**(ex[0]-ex[2]+[0,2,2,1][i])*e**ex[1]*P**ex[2]for ex,co in terms)*s.prod(y**k for y,k in zip(x,mon[:6]))
  qs.append(q)
 Rc=15*P**2-6*P*e+6*P+80*a;Tc=-5*P+6*e-2;Sc=21*P*e+16*P+140*a;Vc=9*P*e+4*P+60*a;Kc=85*P**2-12*P*e+58*P+600*a;Ka=10*P**2+3*P*e+13*P+100*a
 B=30*a*Rc*x[0]-9*P*Sc*x[1]+90*P*a*Tc*x[2]-18*P*Vc*x[3]-2*a*Kc*x[4];A=(B+10*a*P*J*x[4])/3
 t=(30*P*a*Tc*v-6*P*Vc*u-4*a*Ka*z)/(60*P**2);images={x[0]:0,x[1]:3*t,x[2]:v,x[3]:-s.Rational(7,2)*t+u,x[4]:z,x[5]:0};assert s.cancel(A.subs(images,simultaneous=True))==0;assert s.cancel(B.subs(images,simultaneous=True)+10*a*P*J*z)==0
 Q2=s.cancel((qs[2]-e*qs[1])/(a*a*J));assert s.cancel(qs[1].subs(images,simultaneous=True))==0 and s.cancel(Q2.subs(images,simultaneous=True))==0
 rows=[s.Poly(s.cancel(P**4*qs[i].subs(images,simultaneous=True)),u,v,z)for i in [0,3]];mons=[(i,j,2-i-j)for i in range(3)for j in range(3-i)]
 with(OUT/'actual-pair-native.txt').open('w')as f:
  f.write('2 6 0\n')
  for q in rows:
   for mon in mons:
    pp=s.Poly(q.coeff_monomial(s.prod(y**k for y,k in zip((u,v,z),mon))),a,e,P,domain=s.QQ);terms=[]if pp.is_zero else pp.terms();f.write(str(len(terms))+''.join(' '+str(c)+' '+' '.join(map(str,ex))for ex,c in terms)+'\n')
 sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest();assert sha(OUT/'actual-pair-native.txt')==sha(ROOT/'work/nonsquare-v62/four-quadric-pair-gate-v52/pair-native.txt')
 flags=shlex.split(subprocess.check_output(['pkg-config','--cflags','--libs','flint'],text=True));bins={}
 for stem in ['polynomial_unit_span_pivots_20261103_v62','polynomial_exact_scalar_unit_check_20261103_v62']:
  exe=OUT/(stem+'-bin');subprocess.run(['g++','-O2','-std=c++17',str(ROOT/'code'/(stem+'.cpp')),'-o',str(exe),*flags],check=True);bins[stem]=exe
 units=ROOT/'work/nonsquare-v62/raw-fourquad-structure-v49/aPJ-units.txt';pivot=bins['polynomial_unit_span_pivots_20261103_v62'];W=ROOT/'work/nonsquare-v62/pair-gate-small-units-v64'
 with(OUT/'prefix.log').open('w')as f:subprocess.run([str(pivot),str(OUT/'actual-pair-native.txt'),str(OUT/'pair-after.txt'),str(units),'0'],stdout=f,stderr=subprocess.STDOUT,check=True)
 assert sha(OUT/'pair-after.txt')==sha(ROOT/'work/nonsquare-v62/four-quadric-pair-gate-v52/pair-after.txt')
 lines=(OUT/'pair-after.txt').read_text().splitlines();(OUT/'row0.txt').write_text('1 6 0\n'+'\n'.join(lines[3:9])+'\n')
 with(OUT/'unit.log').open('w')as f:subprocess.run([str(bins['polynomial_exact_scalar_unit_check_20261103_v62']),str(OUT/'row0.txt'),str(W/'D2-row0.witness'),str(units)],stdout=f,check=True)
 assert 'EXPONENTS 1 4 0' in(OUT/'unit.log').read_text()
 with(OUT/'pivot.log').open('w')as f:subprocess.run([str(pivot),str(OUT/'pair-after.txt'),str(OUT/'last-row.txt'),str(units),'1',str(W/'D2-row0.witness'),'0'],stdout=f,stderr=subprocess.STDOUT,check=True)
 assert sha(OUT/'last-row.txt')==sha(W/'last-row.txt')
 r=dict(all_exact_actual_source_and_QQ_pivot_checks_passed=True,pair_space='A=x0=W=0, B free',pair_coordinate_denominators_only_P=True,first_row_known_unit='a*P^4',cofactor_terms=10,auxiliary_pair_safe_pivots=1,remaining_shape=[1,6],pair_rank2_over_QQ_not_proved=True,GF32003_rank2_diagnostic_only=True,no_global_raw24_or_product41_or_ambient_claim=True,seconds=time.monotonic()-start);(OUT/'replay-report.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r),flush=True)
if __name__=='__main__':main()
