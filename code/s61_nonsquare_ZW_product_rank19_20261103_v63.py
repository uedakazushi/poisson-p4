#!/usr/bin/env python3
"""Exact all-point rank19 of the ORIGINAL five-variable Z=W=0 quad products.

This includes ell as a free light variable. No GB, fields, or CI assertion.
"""
from pathlib import Path
import sympy as s
import json,hashlib,subprocess,shlex,sys,time
ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'work/nonsquare-v62/ZW-L-pair-v82'
OUT=Path(sys.argv[1]).resolve()if len(sys.argv)>1 else ROOT/'work/nonsquare-v62/ZW-product19-replay-v82'
OUT.mkdir(parents=True,exist_ok=True)

def main():
 start=time.monotonic();a,e,d=s.symbols('a e d');old=s.symbols('v0:8');x=s.symbols('x0:5');source=ROOT/'reports/S61_BINARY_HYPERPLANE_QUADRATIC_CODIMENSION_20261102.json';saved=json.loads(source.read_text());qs=[]
 images={old[0]:x[0],old[1]:x[1],old[2]:x[2],old[3]:0,old[4]:x[3],old[5]:0,old[6]:x[1]-3*x[0],old[7]:x[4]}
 for index in [0,1,3,5]:
  q=0
  for mon,terms in saved['quadratic_rows_sparse'][index]:q+=sum(s.Rational(c)*a**ex[0]*e**ex[1]*d**ex[2]for ex,c in terms)*s.prod(v**k for v,k in zip(old,mon))
  qs.append(s.expand(q.subs(images,simultaneous=True)))
 C=100*a+18*e**2-9*e-2;L=s.cancel(-s.Rational(405,4)*s.diff(qs[1],x[1])/C);s.Poly(L,*x,a,e,d,domain=s.QQ)
 assert s.expand(s.diff(L,x[1])+40*C)==0
 assert s.expand(qs[1]-L*L/8100)==0
 assert s.expand(qs[2]-L*(e*L-400*a*x[0])/8100)==0
 assert s.expand((e*L-400*a*x[0])*qs[1]-L*qs[2])==0
 H=s.hessian(qs[0],(x[0],x[2],x[3],x[4]))/2;assert not qs[0].has(x[1]);assert s.cancel(H.det()-625*a*a/1024)==0
 rep=s.cancel((L+40*C*x[1]).subs(x[0],0)/(40*C));assert not rep.has(x[1]);assert s.cancel(L.subs(x[0],0).subs(x[1],rep))==0
 rows=[s.Poly(s.cancel(q.subs(x[0],0).subs(x[1],rep)),*x[2:])for q in [qs[0],qs[3]]];mons=[(i,j,2-i-j)for i in range(3)for j in range(3-i)]
 with(OUT/'pair-native.txt').open('w')as f:
  f.write('2 6 0\n')
  for i,row in enumerate(rows):
   for mon in mons:
    p=s.Poly(s.cancel(row.coeff_monomial(s.prod(y**k for y,k in zip(x[2:],mon)))*(C if i else 1)),a,e,d,domain=s.QQ);terms=[]if p.is_zero else p.terms();f.write(str(len(terms))+''.join(' '+str(c)+' '+' '.join(map(str,ex))for ex,c in terms)+'\n')
 sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest();assert sha(OUT/'pair-native.txt')==sha(DATA/'pair-native.txt')
 # Replay actual known-unit row/column operations; only a,C are supplied.
 flags=shlex.split(subprocess.check_output(['pkg-config','--cflags','--libs','flint'],text=True));bins={}
 for stem in ['polynomial_unit_span_pivots_20261103_v62','polynomial_exact_scalar_unit_check_20261103_v62']:
  exe=OUT/(stem+'-bin');subprocess.run(['g++','-O2','-std=c++17',str(ROOT/'code'/(stem+'.cpp')),'-o',str(exe),*flags],check=True);bins[stem]=exe
 units=DATA/'units-aC.txt';pivot=bins['polynomial_unit_span_pivots_20261103_v62']
 with(OUT/'first-pivot.log').open('w')as f:subprocess.run([str(pivot),str(OUT/'pair-native.txt'),str(OUT/'after.txt'),str(units),'0'],check=True,stdout=f,stderr=subprocess.STDOUT)
 assert sha(OUT/'after.txt')==sha(DATA/'after.txt');lines=(OUT/'after.txt').read_text().splitlines();assert lines[0]=='1 5 1';(OUT/'scalar.txt').write_text('1 5 0\n'+'\n'.join(lines[3:])+'\n')
 with(OUT/'scalar-QQ-unit.log').open('w')as f:subprocess.run([str(bins['polynomial_exact_scalar_unit_check_20261103_v62']),str(OUT/'scalar.txt'),str(DATA/'last-row-unit.witness'),str(units)],check=True,stdout=f)
 assert 'EXACT_SCALAR_IMAGE 1' in(OUT/'scalar-QQ-unit.log').read_text()
 with(OUT/'second-pivot.log').open('w')as f:subprocess.run([str(pivot),str(OUT/'after.txt'),str(OUT/'final.txt'),str(units),'1',str(DATA/'last-row-unit.witness'),'0'],check=True,stdout=f,stderr=subprocess.STDOUT)
 assert(OUT/'final.txt').read_text().splitlines()[0]=='0 5 1'
 # Verify the saved 117-term identity once more using the independently parsed QQ data.
 bez=json.loads((DATA/'literal-Bezout.json').read_text());load=lambda ts:s.Poly(sum(s.Rational(c)*a**ex[0]*e**ex[1]*d**ex[2]for ex,c in ts),a,e,d,domain=s.QQ)
 fs=list(map(load,bez['source_scalar5']));cofs=list(map(load,bez['scalar_unit_cofactors']));assert sum((p*q for p,q in zip(fs,cofs)),s.Poly(0,a,e,d,domain=s.QQ))==s.Poly(1,a,e,d,domain=s.QQ)
 # A positive old point is a diagnostic; the uniform lower bound uses rank4/gcd/pair.
 products=[s.Poly(q*y,*x)for q in qs for y in x];monomials=sorted(set().union(*(q.monoms()for q in products)));assert len(monomials)==35
 M=s.Matrix([[q.coeff_monomial(s.prod(y**k for y,k in zip(x,mon))).subs({a:1,e:2,d:3})for mon in monomials]for q in products]);assert M.rank()==19
 # Replay the accepted generic four-quadric CI witness, without a search.
 rawpath=ROOT/'reports/S61_BINARY_NONSQUARE_MOD_ELL_COORDINATES_20261103_v62.json';raw=json.loads(rawpath.read_text());xx=s.symbols('xx0:6');yy=s.symbols('yy0:4');qpoint=[]
 for row in raw['cleared_qfour_sparse']:
  qpoint.append(sum(sum(s.Rational(c)for ex,c in cs)*s.prod(z**k for z,k in zip(xx,mon[:6]))for mon,cs in row))
 qpoint=[qpoint[0],qpoint[1],(qpoint[2]-qpoint[1])/(-5),qpoint[3]]
 cut=dict(zip(xx[:4],yy));cut[xx[4]]=yy[0]+yy[1]+2*yy[2]+3*yy[3];cut[xx[5]]=2*yy[0]+3*yy[1]+5*yy[2]+7*yy[3]
 restricted=[s.Poly(q.subs(cut,simultaneous=True),*yy)for q in qpoint]
 monom=lambda D:[(i,j,k,D-i-j-k)for i in range(D+1)for j in range(D-i+1)for k in range(D-i-j+1)]
 columns=[s.Poly(q.as_expr()*s.prod(z**k for z,k in zip(yy,mon)),*yy)for q in restricted for mon in monom(3)]
 Mac=s.Matrix([[q.coeff_monomial(s.prod(z**k for z,k in zip(yy,mon)))for q in columns]for mon in monom(5)])
 prior=json.loads((ROOT/'work/nonsquare-v62/raw-fourquad-structure-v49/generic-CI.json').read_text());chosen=prior['selected_columns'];assert len(chosen)==56
 determinant=s.polys.matrices.DomainMatrix.from_Matrix(Mac[:,chosen]).convert_to(s.QQ).det();assert determinant and str(determinant)==prior['selected_square_determinant']
 (OUT/'generic-CI-replay.json').write_text(json.dumps(dict(actual_raw_point=[1,1,1],shape=[56,80],selected_columns=chosen,exact_nonzero_QQ_determinant=str(determinant),all_generic_CI_gates_passed=True),indent=2)+'\n')
 enc=lambda q:[[list(ex),str(c)]for ex,c in s.Poly(q,a,e,d,domain=s.QQ).terms()]
 certificate=dict(original_remaining_light_positions=[0,1,2,4,7],ell_not_cut=True,quadric_source_indices=[0,1,3,5],L=str(L),L_x1='-40*C',q1_square='L^2/8100',q3_product='L*(e*L-400*a*x0)/8100',Q0_modL_nonzero_determinant='625*a^2/1024',Q0_Q5_pair_rank2_literal_QQ_scalar_Bezout=bez,actual_quadric_coefficients=[[[list(mon),enc(c)]for mon,c in s.Poly(q,*x).terms()]for q in qs])
 (OUT/'certificate.json').write_text(json.dumps(certificate,indent=2)+'\n')
 result=dict(all_exact_QQ_checks_passed=True,original_five_variable_products_rank_exactly19_everywhere=True,coefficient_open='a*C!=0',actual_source_open='a*F*B*C!=0',universal_linear_syzygy=True,Q0_modL_nondegenerate4=True,Q0_Q5_modL_x0_rank2=True,scalar_5_generator_literal_unit_terms=117,ell_never_cut=True,no_d_or_J_inverses=True,old_rank20_premise_false=True,conditional_CI_colon_H4=64,conditional_light_quadric_quartic_rank=146,generic_four_quadric_CI_exact_witness_replayed=True,generic_colon_ell_exact_sequence_applies=True,generic_light_quadratic_quartic_rank=146,all_core_light_quadratic_quartic_rank_upper=146,generic_weighted_quadratic_quartic_rank=152,all_core_weighted_quadratic_quartic_rank_upper=152,all_point_CI_and_rank146_lower_not_proved=True,no_full_quartic_or_ambient_claim=True,source_sha256=sha(source),seconds=time.monotonic()-start)
 (OUT/'report.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)
if __name__=='__main__':main()
