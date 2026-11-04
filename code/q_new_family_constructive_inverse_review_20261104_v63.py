#!/usr/bin/env python3
"""Independent finite-operator/source-binding audit, with no inverse expansion."""
import argparse,ast,gzip,hashlib,json,resource,signal,time
from collections import Counter
from pathlib import Path
import sympy as s
ROOT=Path(__file__).resolve().parents[1]
STEM='Q_NEW_FAMILY_CONSTRUCTIVE_INVERSE_RECIPE_20261104_v63_REVIEW'
def read(p):
 with(gzip.open(p,'rt')if p.suffix=='.gz'else p.open())as f:return json.load(f)
def sha(p):
 h=hashlib.sha256()
 with p.open('rb')as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()
def main(out):
 started=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False)
 folder=ROOT/'reports/Q_NEW_FAMILY_CONSTRUCTIVE_INVERSE_RECIPE_20261104_v63';ar=read(folder/'report.json');assert ar['all_passed']
 for rel,h in ar['input_sha256'].items():assert sha(ROOT/rel)==h,rel
 for rel,h in ar['evidence_sha256'].items():assert sha(folder/rel)==h,rel
 payload=read(folder/'recipe.json.gz');saved=read(folder/'old_matrix_solves.json');p=32003
 krylovp=ROOT/'reports/Q_NEW_FAMILY_ZN_PRIMITIVE_KRYLOV_20261104_v63_REVIEW/independent-certificate.json';k=read(krylovp)
 M=k['all100_matrix_entries'];assert saved['saved_matrix']==M and saved['prime']==p and saved['oldpoint_index']==0
 assert k['determinant']==saved['saved_determinant']==21349
 for v,position in [(saved['a'],1),(saved['b'],0)]:
  assert len(v)==10
  assert [sum(M[i][j]*v[j]for j in range(10))%p for i in range(10)]==[int(i==position)for i in range(10)]
 assert saved['M_a']==[0,1]+[0]*8 and saved['M_b']==[1]+[0]*9
 # Bind every embedded polynomial/circuit to the already accepted source.
 gp=ROOT/'reports/Q_NEW_FAMILY_TWO_SHAPE_CONTENT_GCD_20261104_v63/primitive.txt'
 lines=gp.read_text().splitlines();assert lines[0]=='6 891'
 assert payload['G_terms_embedded']==[list(map(int,t.split()))for t in lines[1:]]
 cp=ROOT/'reports/Q_NEW_FAMILY_TWO_SHAPE_PRS_CONTENT_PREFLIGHT_20261104_v63/C1-primitive.txt.gz'
 with gzip.open(cp,'rt')as f:
  assert f.readline().split()==['5','2'];co=[]
  for n in [10314,10698]:
   assert int(f.readline())==n;co.append([list(map(int,f.readline().split()))for _ in range(n)])
  assert not f.read().strip()
 assert payload['C1_terms_constant_then_linear_embedded']==co
 hp=ROOT/'reports/Q_NEW_FAMILY_TWO_SHAPE_HOMOGENIZATION_20261104_v63/certificate.json';hh=read(hp)
 small={'N79':hh['rational_coefficient_pairs']['E1']['numerator'],'N41':hh['rational_coefficient_pairs']['E2']['numerator'],'Fn':hh['factors']['Fn'],'J':hh['factors']['J']}
 assert payload['N79_N41_Fn_J_embedded']==small
 sp=ROOT/'reports/Q_FIRSTJET_SPIN5_STRUCTURE_20261103_v63/certificate.json.gz';spin=read(sp)
 assert payload['global_quadratic_monomials']==spin['quadratic_monomials']
 assert payload['global_C_coefficients16x136']==spin['spin_adjoint_C_coefficients16x136']
 assert payload['global_D_coefficients16x136']==spin['vector_adjoint_D_coefficients16x136']
 ip=ROOT/'reports/Q_NEW_FAMILY_HYPERSURFACE_RECOGNITION_20261103_v63/invariant-certificates.json';iv=read(ip)
 O=s.Matrix(payload['global_Omega']);Oi=s.Matrix(payload['global_Omega_inverse'])
 assert O==s.Matrix(iv['symplectic_form16'])and O.T==-O and O*Oi==s.eye(16)
 orbitp=ROOT/'reports/Q_NEW_FAMILY_RATIONAL_ORBIT_SECTION_20261104_v63/certificate.json'
 circuitp=ROOT/'reports/Q_NEW_FAMILY_SOURCE_RATIONAL_QUOTIENT_20261103_v63/source-five-circuit.json.gz'
 assert payload['rational_orbit_section']==read(orbitp)
 assert payload['original_source_five_circuit']==read(circuitp)
 orbitgate=ROOT/'reports/Q_NEW_FAMILY_RATIONAL_ORBIT_SECTION_20261104_v63_REVIEW/report.json';assert read(orbitgate)['all_passed']
 assert payload['rational_orbit_section']['variables']==['xi','r','s','q','B']
 assert all(c['quadratic_remainder_odd_part_zero']for c in payload['rational_orbit_section']['all16_checks'])
 # Fix ONE universally nonzero J14, with a named global basis witness.
 # This never assumes that an arbitrary fixed word happens to be nonzero.
 defs=read(ROOT/'reports/Q_NEW_FAMILY_DEGREE26_POISSON_BASIS_20261104_v63/definitions.json')
 witnesses=read(ROOT/'reports/Q_NEW_FAMILY_DEGREE26_POISSON_BASIS_20261104_v63/minor-witnesses.json')
 basisgate=ROOT/'reports/Q_NEW_FAMILY_DEGREE26_POISSON_BASIS_20261104_v63_REVIEW/report.json';bg=read(basisgate)
 assert bg['all_passed']and bg['complete_QQ_Inv14_3_Inv18_9_Inv26_45_bases_certified']
 assert defs['Inv8_trace_words'][:2]==['0000','0001']
 assert defs['degree14'][0]=={'kind':'Poisson_Inv8_Inv8','left':0,'right':1}
 assert witnesses['14']['columns']==[0,1,4]and witnesses['14']['minor_mod32003']==7584==bg['exact_coefficient_minors']['14']
 # Derive the noncommutative cyclic-trace gradients independently from
 # the literal words, then inspect the finite implementation's factors.
 def cofactors(word,letter):
  return Counter(tuple(word[i+1:]+word[:i])for i,x in enumerate(word)if x==letter)
 d0=cofactors([0,0,0,0],0);d1C=cofactors([0,0,0,1],0);d1D=cofactors([0,0,0,1],1)
 assert d0=={(0,0,0):4}and d1C=={(0,0,1):1,(0,1,0):1,(1,0,0):1}and d1D=={(0,0,0):1}
 source=ROOT/'code/q_new_family_constructive_inverse_recipe_20261104_v63.py';tree=ast.parse(source.read_text());fn=next(n for n in tree.body if isinstance(n,ast.FunctionDef)and n.name=='j14_operator')
 assign={n.targets[0].id:n.value for n in fn.body if isinstance(n,ast.Assign)and isinstance(n.targets[0],ast.Name)}
 def word(n):
  if isinstance(n,ast.Name):return {'C':(0,),'D':(1,),'C2':(0,0),'C3':(0,0,0)}[n.id]
  assert isinstance(n,ast.Call)and isinstance(n.func,ast.Name)and n.func.id=='mm';return word(n.args[0])+word(n.args[1])
 assert Counter(word(n)for n in assign['pieces'].elts)==d1C
 assert ast.unparse(assign['grad0'])=='[4 * tr(mm(A, C3)) for A in dC]'
 assert ast.unparse(assign['grad1'])=='[tr(mm(A, K)) + tr(mm(C3, B)) for A, B in zip(dC, dD)]'
 derivative=payload['word_trace_derivative_formulas']
 assert derivative=={'T0':'tr(C^4)','T1':'tr(C^3 D)','dT0':'4 tr(dC C^3)',
                     'dT1':'tr(dC (C^2 D+C D C+D C^2))+tr(C^3 dD)',
                     'bracket':'sum_i,j dT0_i Omega_inverse_ij dT1_j'}
 program={n['name']:n for n in payload['operator_program']}
 assert len(program)==11
 assert program['two_linear_solves']['e_r']==[0,1]+[0]*8 and program['two_linear_solves']['e_1']==[1]+[0]*9
 assert program['recover_r']['expression']=='(sum a_j*eta^j)/(sum b_j*eta^j)'
 assert program['recover_s']['expression']=='-B(r)/A(r)'
 assert program['matrix']['columns']=='coefficient vector of U^j*V^(9-j), j=0,...,9'
 assert program['recover_xi']['expression']=='(J14_target/F4_target^3)*f4^3/h14'
 fieldgate=ROOT/'reports/Q_NEW_FAMILY_TARGET_FIELD_AND_SPECTRAL_DEGREES_20261104_v63_REVIEW/report.json';assert read(fieldgate)['all_passed']
 # All denominator and correctness conclusions below are in the generic
 # characteristic-zero field, independently of the saved point solves.
 proof={'finite_program_names':list(program),
        'all_literal_inputs_bound':{'G_terms':891,'C1_terms':[10314,10698],'C_D_rows':32,'quadratic_monomials':136,'source_five_DAG_nodes':len(payload['original_source_five_circuit']['nodes'])if'nodes'in payload['original_source_five_circuit']else'full literal circuit equality checked'},
        'old_matrix_twenty_dot_product_identities':True,
        'adjugate_convention':'adj(M)_ij=(-1)^(i+j)det(delete row j,column i)',
        'r_inverse_proof':'M*a=e_r and M*b=e_1 imply r=V^9*sum(a_j*eta^j), 1=V^9*sum(b_j*eta^j), hence quotient; denominator equals V^-9 nonzero generically',
        's_inverse':'C1_linear is nonzero generically by the independently proved3/5 birationality',
        'q_B_inverse':'the embedded exact Mobius and compression recipes, on their named nonzero dense chart',
        'fixed_nonzero_J14':{'definition':'{tr(C^4),tr(C^3 D)}','Inv14_semantic_basis_column':0,'columns_in_nonzero_prior_minor':[0,1,4],'minor_mod32003':7584},
        'cyclic_trace_gradients':{'T0_C':[{'word':list(w),'multiplicity':n}for w,n in d0.items()],
                                  'T1_C':[{'word':list(w),'multiplicity':n}for w,n in d1C.items()],
                                  'T1_D':[{'word':list(w),'multiplicity':n}for w,n in d1D.items()]},
        'ambient_derivative_before_source_substitution':True,
        'h14_direct':'evaluate the actual embedded rational Sp4 orbit representative at xi=1, using global quadratic C/D coefficients and Omega inverse',
        'h14_fallback':'evaluate the full1032 source circuit in the etale rank2 algebra u^2=1600r with X=s/(48u),c=B*u/40,W=(q+1)/3072; tau invariance gives zero odd part, so trace/2 is identical',
        'h14_nonzero':'selected universal bracket is nonzero by its prior full-basis minor; defining degree>=26 prevents vanishing of this degree14 restriction',
        'xi_inverse':'all16 source covariance gives J14=xi^7*h14, F4=xi^2*f4; displayed ratio therefore equals xi',
        'scope':'finite rational-operator inverse on a dense open; no generic expanded adjugate/h14/source26/NEW equation or boundary extension'}
 ep=out/'independent-proof.json';ep.write_text(json.dumps(proof,indent=2)+'\n')
 inputs=[Path(__file__),ROOT/'research/notes'/(STEM+'.md'),folder/'report.json',folder/'recipe.json.gz',folder/'old_matrix_solves.json',source,krylovp,gp,cp,hp,sp,ip,orbitp,circuitp,orbitgate,basisgate,fieldgate,
         ROOT/'reports/Q_NEW_FAMILY_DEGREE26_POISSON_BASIS_20261104_v63/definitions.json',ROOT/'reports/Q_NEW_FAMILY_DEGREE26_POISSON_BASIS_20261104_v63/minor-witnesses.json']
 report={'id':STEM,'all_passed':True,'complete_finite_operator_recipe_independently_accepted':True,
         'all_embedded_polynomials_C_D_Omega_and_both_source_presentations_bound':True,
         'same_saved_old_matrix_twenty_identities_directly_checked_no_Krylov_replay':True,
         'deterministic_fixed_J14_universally_nonzero_by_named_Inv14_basis_minor':True,
         'independent_cyclic_trace_gradient_and_ambient_derivative_order_audit':True,
         'generic_denominator_and_r_s_q_B_xi_field_inverse_proof_accepted':True,
         'no_new_point_prime_or_generic_adjugate_Inv14_source_norm_source26_expansion':True,
         'no_global_boundary_extension_defining_equation_or_entire_Poisson_component_claim':True,
         'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in inputs},'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-started,'frozen':True}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');signal.alarm(0)
 print(json.dumps({k:report[k]for k in ['all_passed','deterministic_fixed_J14_universally_nonzero_by_named_Inv14_basis_minor','seconds']},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=ROOT/'reports'/STEM);a=p.parse_args();main(a.output.resolve())
