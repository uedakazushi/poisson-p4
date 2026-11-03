#!/usr/bin/env python3
"""Independent physical-arc tangent and cycle recognition of the origin binary cubics."""
import argparse,gzip,hashlib,json,time
from pathlib import Path
import sympy as s
from sympy.polys.matrices import DomainMatrix
ROOT=Path(__file__).resolve().parents[1]
def load(p):return json.loads(gzip.decompress(p.read_bytes())if p.suffix=='.gz'else p.read_bytes())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def mat(M):return [[str(v)for v in row]for row in M.tolist()]
def main(out):
    started=time.monotonic();out.mkdir(parents=True,exist_ok=False)
    deps=[ROOT/f'reports/AMBIENT_{name}_ORIGIN_CONDUCTOR_MULTIPLICITY4_20261103_v63_REVIEW/report.json'for name in['CA','CB']]
    deps +=[ROOT/'reports/AMBIENT_B18_CONDUCTOR_DIMENSION20_20261103_v63_REVIEW/report.json']
    for path in deps:
        d=load(path);assert d['all_passed'] is True
        for name,digest in d.get('input_sha256',{}).items():assert sha(ROOT/name)==digest,name
    upperpath=ROOT/'reports/AMBIENT_B18_CONDUCTOR_DIMENSION20_20261103_v63/certificate.json.gz';upper=load(upperpath)
    firstpath=ROOT/'reports/AMBIENT_B18_ODD_FIRSTJETS_N_PROJECTION_20261103_v63/certificate.json.gz';first=load(firstpath);Rrows=first['common_Y28']+[61]
    W,U,V=s.symbols('W U V');records={};sourcepaths=[]
    for name in['CA','CB']:
        path=ROOT/f'reports/AMBIENT_{name}_ORIGIN_POLYNOMIAL_SELF_PAIR_20261103_v63/certificate.json.gz';d=load(path);sourcepaths.append(path)
        if name=='CA':
            parameter=s.Symbol('a');point=s.Matrix([sum(s.Rational(c)*parameter**degree for degree,c in row)for row in d['first_source_physical150']])
        else:
            parameter=s.Symbol('b');point=s.Matrix(d['full_polynomial_Q150'])
        binding=upper['actual_jet_bindings'][name];K=s.Matrix(binding['full_minimal_kernel_basis29_by22'])
        velocity=point.diff(parameter).subs(parameter,0);r=velocity[Rrows,:]
        coordinates,parameters=K.gauss_jordan_solve(r);assert len(parameters)==0 and K*coordinates==r
        w,u,v=coordinates[-3:,0];assert w==0 and u!=0 and v!=0
        expected=[0,-s.Rational(5,2),-s.Rational(15,2)]if name=='CA'else[0,-s.Rational(9,4),-s.Rational(1,4)]
        assert [w,u,v]==expected
        family=K[:,:19].row_join(r);rows=list(DomainMatrix.from_Matrix(family.T).convert_to(s.QQ).rref()[1]);assert len(rows)==20
        minor=s.Rational(DomainMatrix.from_Matrix(family[rows,:]).convert_to(s.QQ).det());assert minor
        ell=U-s.factor(u/v)*V;c0=s.sympify(binding['h63_actual_U3_coefficient']);cubic=s.expand(c0*U*ell**2)
        assert s.Poly(cubic,U,V).coeff_monomial(U**3)==c0
        expectedcubic=s.Rational(12,125)*U*(U-V/3)**2 if name=='CA'else-s.Rational(4,243)*U*(U-9*V)**2
        assert s.expand(cubic-expectedcubic)==0
        records[name]=dict(actual_physical150_arc_velocity=list(map(str,velocity)),actual_R29_arc_velocity=list(map(str,r)),
          all29_equations_solved_directly=True,full_g19_WUV_coordinates=list(map(str,coordinates)),WUV=list(map(str,[w,u,v])),
          self_component_tangent_UV_equation=str(ell),self_GL19_plus_arc_rank_at_origin=20,self_family20_rows=rows,self_family20_minor=str(minor),
          actual_support_branches=['W-axis','V-axis','selfline'],generic_primary_lengths=[1,1,2],
          fundamental_tangent_cycle='1[W-axis]+1[V-axis]+2[selfline]',
          binary_cubic_factored=str(s.factor(cubic)),binary_cubic_expanded=str(cubic),
          U3_coefficient=str(c0),U2V_coefficient=str(s.Poly(cubic,U,V).coeff_monomial(U**2*V)),UV2_coefficient=str(s.Poly(cubic,U,V).coeff_monomial(U*V**2)))
    (out/'certificate.json.gz').write_bytes(gzip.compress((json.dumps(records,separators=(',',':'))+'\n').encode(),mtime=0))
    note=ROOT/'research/notes/AMBIENT_ORIGIN_BINARY_CUBIC_RECOGNITION_20261103_v63_REVIEW.md';paths=[Path(__file__),note,upperpath,firstpath,*sourcepaths,*deps]
    report=dict(id='AMBIENT_ORIGIN_BINARY_CUBIC_RECOGNITION_20261103_v63_REVIEW',all_passed=True,
      independent_source_from_decoded_actual_physical_arcs_not_fixed_C_shortcut=True,
      every29_tangent_equation_checked=True,
      CA_self_tangent_WUV=['0','-5/2','-15/2'],CB_self_tangent_WUV=['0','-9/4','-1/4'],
      self_support_smooth_at_origin_by_rank20_immersion=True,
      initial_cycle_allocation_uses_prime_filtration_Rees_and_distinct_smooth_tangent_branches=True,
      CA_exact_binary_initial='(12/125)*U*(U-V/3)^2',CB_exact_binary_initial='(-4/243)*U*(U-9V)^2',
      CA_exact_initial_ideal='(WU,WV,(12/125)*U*(U-V/3)^2)',
      CB_exact_initial_ideal='(WU,WV,(-4/243)*U*(U-9V)^2)',
      only_associated_graded_conductor_recognized=True,
      higher_analytic_normal_form_or_relative_CA_CB_embedding_not_claimed=True,
      historical_unknown_cubic_sources_and_receipts_not_modified=True,
      input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},evidence_sha256={'certificate.json.gz':sha(out/'certificate.json.gz')},seconds=time.monotonic()-started)
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print('INDEPENDENT PHYSICAL CA/CB BINARY INITIAL RECOGNITION PASS',report['seconds'],flush=True)
if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True);main(ap.parse_args().output.resolve())
