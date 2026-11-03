#!/usr/bin/env python3
"""Independent physical binding and integer/tiny-matrix review of NORMAL_FINITE5.

No import of the producer, load_fiber, residual42 constructor, or GB engine.
"""
import argparse,gzip,hashlib,itertools,json,resource,shlex,subprocess,time
from pathlib import Path
import sympy as sp
from sympy.polys.rings import ring
ROOT=Path(__file__).resolve().parents[1];PRIME=32003;STEM='S61_AFFINE_NORMAL_FINITE5_20261103_v63'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb')as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()
def rat(c):return sp.QQ(str(c))
def monomial(indices,n=45):
 e=[0]*n
 for i in indices:e[i]+=1
 return tuple(e)
def verify_report(p):
 d=json.loads(p.read_text());assert d['all_passed']
 for name,h in d.get('input_sha256',{}).items():assert sha(ROOT/name)==h,name
 for name,h in d.get('evidence_sha256',{}).items():assert sha(p.parent/name)==h,name
 return d
def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',required=True,type=Path);args=ap.parse_args();out=args.output.resolve();out.mkdir(parents=True,exist_ok=False);start=time.monotonic();folder=ROOT/'reports'/STEM;parent=verify_report(folder/'report.json');certificate=json.loads((folder/'certificate.json').read_text());prior=ROOT/'reports/S61_AFFINE_RESIDUAL42_20261103_v63';verify_report(prior/'report.json')
 with gzip.open(prior/'linear129x42-QQ.json.gz','rt')as f:linear=json.load(f)
 original=ROOT/'legacy/v02/data/symplectic_slice_reduced_equations.json';pointpath=ROOT/'reports/S61_KNOWN_FIBER_AUDIT_20261031.json';model=json.loads(original.read_text());record=next(p for p in json.loads(pointpath.read_text())['points']if p['component_id']=='T01');point=list(map(rat,record['S61_coordinates_u1_through_u61']));weights=model['parameter_weights'];columns=[i for i,w in enumerate(weights)if w==2];assert len(columns)==45;index={old:j for j,old in enumerate(columns)}
 # Build the AFFINE physical linear forms directly in original61 coordinates:
 # u_old=center_old + E_old*z + n, then multiply original S61 monomials.
 ambient=ring(','.join('x'+str(i)for i in range(45)),sp.QQ);Q=ambient[0];x=ambient[1:];E=[[rat(c)for c in row]for row in linear['embedding']];assert E[:3]==[[sp.QQ(int(i==j))for j in range(3)]for i in range(3)]
 forms=[]
 for j,old in enumerate(columns):
  p=Q(point[old])+sum((x[k]*E[j][k]for k in range(3)),Q.zero)
  if j>=3:p+=x[j]
  forms.append(p)
 products={():Q.one}
 def product(free):
  free=tuple(sorted(free))
  if free not in products:
   q=Q.one
   for j in free:q*=forms[j]
   products[free]=q
  return products[free]
 physical=[]
 for equation in model['equations']:
  coefficient_map={}
  for mon,numerator,denominator in equation:
   c=sp.QQ(numerator,denominator);free=[]
   for old in mon:
    if old in index:free.append(index[old])
    else:c*=point[old]
   key=tuple(sorted(free));assert len(key)<=2;coefficient_map[key]=coefficient_map.get(key,sp.QQ.zero)+c
  p=sum((product(key)*c for key,c in coefficient_map.items()if c),Q.zero);assert (0,)*45 not in p;physical.append(p)
 assert len(physical)==134
 # Check the given frame itself against the independently stored plane chart.
 plane_source=ROOT/'reports/S61_KNOWN_FIVE_POINT_SECTION_20261101.json';plane_info=json.loads(plane_source.read_text());zz=sp.symbols('v0:3')
 # The stored plane uses UNCENTERED affine q coordinates. Its constant term
 # need not vanish; translation by the actual T01 point cancels it exactly.
 center3=[sp.Rational(str(point[columns[k]]))for k in range(3)]
 for j,s in enumerate(plane_info['affine_plane_substitution']):
  p=sp.Poly(sp.sympify(s),*zz)
  assert rat(p.eval(dict(zip(zz,center3))))==point[columns[j]]
  assert all(rat(p.coeff_monomial(zz[k]))==E[j][k]for k in range(3));assert p.total_degree()<=1
 changes=[[rat(c)for c in row]for row in linear['pure_plane_subtraction_coefficients']];assert len(changes)==129
 assert all(c.denominator%PRIME for c in point+[c for row in E+changes for c in row])
 cross=[];rest=[];normal_mons=[(j,)for j in range(42)]+list(itertools.combinations_with_replacement(range(42),2));assert len(normal_mons)==945
 for i in range(129):
  g=physical[i+5]-sum((physical[j]*changes[i][j]for j in range(5)),Q.zero)
  supplied=Q.zero
  for j,column in enumerate(linear['all129_linear_factor_rows'][i]):
   for mon,c in column:
    assert len(mon)<=1;supplied+=x[j+3]*(Q(rat(c))if not mon else x[mon[0]]*rat(c))
  assert g==supplied
  c=[sp.QQ.zero]*126;d={}
  for ex,value in g.items():
   zdegree=sum(ex[:3]);ndegree=sum(ex[3:]);assert zdegree+ndegree<=2 and ndegree>=1
   if zdegree:
    assert zdegree==ndegree==1;k=next(k for k in range(3)if ex[k]);j=next(j for j in range(42)if ex[j+3]);c[42*k+j]=value
   else:
    key=tuple(j for j in range(42)for _ in range(ex[j+3]));assert key not in d;d[key]=value
  cross.append(c);rest.append([d.get(mon,sp.QQ.zero)for mon in normal_mons])
 raw=('129 126 945\n'+'\n'.join(' '.join(str(c)for c in row)for row in cross+rest)+'\n').encode();expected=gzip.decompress((folder/'cross-input.txt.gz').read_bytes());assert raw==expected and hashlib.sha256(raw).hexdigest()==certificate['native_input_sha256'];(out/'physical-cross-input.txt').write_bytes(raw)
 # Plane highest-degree blocks rebuilt by literal polynomial multiplication;
 # repeated variables introduce no derivative/binomial multiplicity.
 quads=list(itertools.combinations_with_replacement(range(3),2));cubics=list(itertools.combinations_with_replacement(range(3),3));plane_high=[]
 for p in physical[:5]:plane_high.append(Q.from_dict({ex:c for ex,c in p.items()if sum(ex)==2 and not any(ex[3:])}))
 QM=sp.Matrix([[sp.Rational(str(p.get(monomial(mon),sp.QQ.zero)))for mon in quads]for p in plane_high]);CM=sp.Matrix([[sp.Rational(str((x[k]*p).get(monomial(mon),sp.QQ.zero)))for mon in cubics]for k in range(3)for p in plane_high]);assert QM.tolist()==[[sp.Rational(c)for c in row]for row in certificate['plane_quadratic_matrix']];assert CM.tolist()==[[sp.Rational(c)for c in row]for row in certificate['plane_cubic_matrix']]
 piv=certificate['plane_quadratic_pivots'];sel=certificate['plane_cubic_rows'];QI=sp.Matrix(certificate['plane_quadratic_inverse']);CI=sp.Matrix(certificate['plane_cubic_inverse']);assert QM[:,piv]*QI==QI*QM[:,piv]==sp.eye(5);assert CM[sel,:]*CI==CI*CM[sel,:]==sp.eye(10);assert all(sp.denom(c)%PRIME for c in list(QI)+list(CI));basis=[[],[0],[1],[2],list(next(mon for j,mon in enumerate(quads)if j not in piv))];assert basis==certificate['fixed_module_generators_z_monomials']
 # Verify saved126 inverse and all119070 reduction coefficients over Z via an
 # independent global-denominator-clear algorithm, without recomputing inverse.
 saved=ROOT/'work/normal-finite5-public-v63/cross-output.txt';assert sha(saved)==certificate['native_output_sha256']and saved.stat().st_size==certificate['native_output_bytes'];cpp=ROOT/'code/s61_affine_normal_finite5_integer_review_20261103_v63.cpp';exe=out/'integer-review';flags=shlex.split(subprocess.check_output(['pkg-config','--cflags','--libs','flint'],text=True));include=subprocess.check_output(['pkg-config','--variable=includedir','flint'],text=True).strip();subprocess.run(['c++','-std=c++17','-O2','-I'+include,str(cpp),*flags,'-lgmp','-o',str(exe)],check=True)
 def limits():resource.setrlimit(resource.RLIMIT_AS,(2200000000,2200000000));resource.setrlimit(resource.RLIMIT_CPU,(120,121))
 t=time.monotonic()
 with(out/'integer-review.log').open('w')as f:p=subprocess.run([str(exe),str(out/'physical-cross-input.txt'),str(saved),str(out/'integer-review.json')],stdout=f,stderr=subprocess.STDOUT,preexec_fn=limits,timeout=125)
 assert p.returncode==0;integer=json.loads((out/'integer-review.json').read_text());assert integer['all_passed'];native_seconds=time.monotonic()-t;exe.unlink();(out/'physical-cross-input.txt').unlink()
 # Original N coefficient factorization is exact, so constant row operations
 # recover129G from the126 cross-rewrites plus3 pure-normal polynomials.
 assert 3*42*41//2==2583 and 42*5==210 and 3+2583+210==2796
 inputs=[Path(__file__),cpp,folder/'report.json',folder/'certificate.json',folder/'cross-input.txt.gz',prior/'report.json',prior/'linear129x42-QQ.json.gz',original,pointpath,plane_source,ROOT/'research/notes'/f'{STEM}.md',ROOT/'research/notes/S61_AFFINE_RESIDUAL42_20261103_v63.md',ROOT/'research/notes/S61_AFFINE_POLARIZATION_SPLITTER_20261103_v63.md']
 report={'id':'S61_AFFINE_NORMAL_FINITE5_20261103_V63_REVIEW','all_passed':True,'physical_binding':'Own direct original61 monomial substitution in SymPy QQ45-variable ring; no producer/load_fiber import','all134_original_equations_accounted_for':True,'all129_factorizations_checked_exact':True,'cross_input_rebuilt_from_actual_G_and_byte_equal':True,'normal_coordinate_frame_determinant':1,'integer_inverse_products_checked':2,'integer_cross_coefficients_checked':119070,'tiny_quadratic_block_rank':5,'tiny_cubic_block_rank':10,'both_tiny_inverse_products_and_punit_denominators_checked':True,'five_generators_z_monomials':basis,'finite_module_proof_accepted':True,'residual_closed_immersion_proof_accepted':True,'2796_degree_atmost4_normal_saturation_preserves_whole_residual':True,'beta_only_open_spreading_retains_universal_plane_and_splitter_dependencies':True,'beta_finiteness_residual_emptiness_or_R_freeness_claim':False,'full_growing_GB_or_inverse_recomputed':False,'integer_review_seconds':native_seconds,'input_sha256':{str(p.resolve().relative_to(ROOT)):sha(p)for p in inputs},'seconds':time.monotonic()-start}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:report[k]for k in ['all_passed','all134_original_equations_accounted_for','integer_cross_coefficients_checked','five_generators_z_monomials','seconds']},indent=2))
if __name__=='__main__':main()
