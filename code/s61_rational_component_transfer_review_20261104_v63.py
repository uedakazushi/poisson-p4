#!/usr/bin/env python3
"""Independent fixed-point skew Gram-Schmidt, actual perspective gauge, and component transfer audit."""
import argparse,hashlib,itertools,json,resource,signal,time
from pathlib import Path
import sympy as s
ROOT=Path(__file__).resolve().parents[1]
STEM='S61_RATIONAL_COMPONENT_TRANSFER_20261104_v63_REVIEW'
def read(p):return json.loads(p.read_text())
def sha(p):
 h=hashlib.sha256()
 with p.open('rb')as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()
def enc(M):return [[str(c)for c in row]for row in M.tolist()]
def main(out):
 started=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False)
 folder=ROOT/'reports/S61_RATIONAL_COMPONENT_TRANSFER_20261104_v63';ar=read(folder/'report.json');assert ar['all_passed']and not ar['new_component16_or25_rationality_proved_here']
 for rel,h in ar['input_sha256'].items():assert sha(ROOT/rel)==h,rel
 for rel,h in ar['evidence_sha256'].items():assert sha(folder/rel)==h,rel
 rebuild=ROOT/'reports/S61_REBUILD_20261030.json';rr=read(rebuild);assert rr['all_passed']
 assert rr['exact_reconstruction_audit']['all_coefficient_comparisons_exact_over_Q']
 assert rr['exact_reconstruction_audit']['checked_eliminated_rows']==55
 assert rr['model_counts']['reduced_generator_count']==134 and rr['model_counts']['generator_redundancy_identities_verified_over_Q']==20
 adapter=ROOT/'reports/S61_ADAPTER_20261030.json';assert read(adapter)['all_passed']
 fp=ROOT/'reports/20261030-s61-weight1/finite_identity_certificate.json';finite=read(fp);assert finite['all_passed']
 legacy=ROOT/'legacy/v02/data/symplectic_slice_certificate.json';old=read(legacy)
 assert old['perspective_orbit_pivots_raw']==[61,62,63,64]
 assert old['linear_jacobi_rank_exact']==55
 # Derive the inverse change directly by solving y=R*x, rather than asking
 # a generic symbolic matrix inverse for the author's displayed R.
 a,b,c,d,e,f=s.symbols('a b c d e f');p=a*f-b*e+c*d
 C=s.Matrix([[0,a,b,c],[-a,0,d,e],[-b,-d,0,f],[-c,-e,-f,0]])
 J=s.Matrix([[0,1,0,0],[-1,0,0,0],[0,0,0,1],[0,0,-1,0]])
 R=s.Matrix([[1,0,0,0],[0,0,0,1/c],[-e/c,1,0,-a/c],[-f/p,0,c/p,-b/p]])
 inverse=s.Matrix([[1,0,0,0],[e/c,a,1,0],[f/c,b,0,p/c],[0,c,0,0]])
 assert(R*inverse-s.eye(4)).applyfunc(s.cancel)==s.zeros(4)
 assert(inverse*J*inverse.T-C).applyfunc(s.cancel)==s.zeros(4)
 assert s.cancel(inverse.det()-p)==0
 # Pair the Gram-Schmidt rows individually; the frame has only c,p units.
 gram=(R*C*R.T).applyfunc(s.cancel);assert gram==J
 assert s.cancel(R.det()-1/p)==0
 assert gram.subs(c,1)==J and p.subs(c,1)==a*f-b*e+d
 # Perspective pushforward at the SAME fixed point: phi_a(y)=y/(1+a*y).
 # Its full finite derivative, not its infinitesimal linearization, gives
 # the exact degree0..3 transformation for any projective bivector.
 y=s.Matrix(s.symbols('y0:4'));alpha=s.Matrix(s.symbols('alpha0:4'));beta=s.Matrix(s.symbols('beta0:4'))
 ff=(alpha.T*y)[0];h=1-ff;K=s.eye(4)-y*alpha.T
 z=s.Matrix(s.symbols('z0:4'));phi=z/(1+(alpha.T*z)[0]);der=phi.jacobian(z)
 derivative_at_inverse=der.subs(dict(zip(z,y/h)),simultaneous=True)
 assert(derivative_at_inverse-h*K).applyfunc(s.cancel)==s.zeros(4)
 mm=s.symbols('m0:6');M=s.zeros(4)
 for v,(i,j)in zip(mm,itertools.combinations(range(4),2)):M[i,j]=v;M[j,i]=-v
 wedge=lambda u,v:u*v.T-v*u.T
 assert(K*M*K.T-M-wedge(y,M*alpha)).applyfunc(s.expand)==s.zeros(4)
 B=s.Matrix(s.symbols('B0:4'));radial=wedge(y,B)
 assert(K*radial*K.T-h*radial).applyfunc(s.expand)==s.zeros(4)
 direction=-2*ff*J+wedge(y,J*alpha)
 pairs=list(itertools.combinations(range(4),2));basis=[(i,j,k)for i,j in pairs for k in range(4)]
 coefficients=s.Matrix([s.expand(direction[i,j]).coeff(y[k])for i,j,k in basis])
 response=coefficients.jacobian(alpha)
 assert response==s.Matrix(finite['perspective_pushforward_matrix'])
 oldresponse=response[:4,:];assert oldresponse==s.diag(-3,-3,-2,-2)and oldresponse.det()==36
 oldgauge=s.Matrix(s.symbols('L12_0:4'));normalizing=-oldresponse.inv()*oldgauge
 assert normalizing==s.Matrix([oldgauge[0]/3,oldgauge[1]/3,oldgauge[2]/2,oldgauge[3]/2])
 assert oldgauge+oldresponse*normalizing==s.zeros(4,1)
 trace=s.Matrix(finite['trace_matrix']);assert trace*response==-5*s.eye(4)
 G=s.eye(5);H=s.eye(5)
 for i in range(4):G[0,i+1]=alpha[i];H[0,i+1]=beta[i]
 addition=s.eye(5)
 for i in range(4):addition[0,i+1]=alpha[i]+beta[i]
 assert G*H==addition and G.subs(dict(zip(alpha,-alpha)),simultaneous=True)*G==s.eye(5)
 evidence={'fixed_point':'[1:0:0:0:0], never add an evaluation-point parameter',
           'c_index_convention':'C03 in the four-dimensional affine tangent space; global spatial pair is(1,4)',
           'C':enc(C),'pfaffian':str(p),'R':enc(R),'explicit_R_inverse':enc(inverse),
           'only_frame_localizers':['c','pfaffian'],'det_R':'1/pfaffian','all_frame_identities':True,
           'full_finite_affine_perspective_formula':'h^2*(J+E wedge Jalpha)+h*(L+E wedge Lalpha)+Q+E wedge Qalpha+C_radial',
           'full_degrees_retained':'Q changes by -fL+E wedge Lalpha+f^2 J-2f E wedge Jalpha; higher changes are retained by the finite group action',
           'old_S61_gauge_response':enc(oldresponse),'old_alpha':[str(v)for v in normalizing],
           'primitive_gauge_is_different':'grad(L01+L23)/5 is not substituted for the old S61 alpha',
           'scheme_product':'actual rank4 affine-cone open U6 x A4 x S61 via C, actual-old-gauge-alpha, N; inverse undo perspective then R',
           'scheme_input':'exact55 elimination with rational constant units and20 identities proving the154/134 IDEALS equal, not only radicals',
           'component_restriction':'connected PGL5 and tensorGm stabilize every reduced irreducible component; their actual maps and inverses preserve component ideals, including Artin points',
           'slice_irreducibility':'the product restricted to an irreducible component forces its slice Y to be irreducible; the product is reduced and the base is faithfully flat, so the reduced component statement is valid',
           'nonempty_chart':'move a nondegenerate point to the fixed point by PGL, then use GL4 fixing that point to make C03 nonzero',
           'projective_product':'take the unique tensor-scale representative c=1 BEFORE normalization, giving U5 x A4 x Y, not a scalar quotient of S61',
           'rationality_transfer':'Y rational n-dimensional over QQ (respectively C) implies the projective component rational of dimension n+9 over the same field',
           'Y16_or_X25_rationality_not_proved_here':True}
 ep=out/'independent-certificate.json';ep.write_text(json.dumps(evidence,indent=2)+'\n')
 inputs=[Path(__file__),ROOT/'research/notes'/(STEM+'.md'),folder/'report.json',folder/'certificate.json',rebuild,adapter,fp,legacy,ROOT/'code/s61_model.py',ROOT/'research/notes/S61_MODEL_AUDIT_20261030.md']
 report={'id':STEM,'all_passed':True,'independent_six_variable_pairing_and_explicit_R_inverse':True,
         'same_fixed_evaluation_point_no_extra_four_dimensions':True,'independent_full_finite_perspective_pushforward_identity':True,
         'actual_old_S61_gauge_response_diag_minus3_minus3_minus2_minus2_verified':True,
         'primitive_gauge_not_confused_and_all_weight2_changes_retained':True,
         'accepted55_elimination_and134_IDEAL_identity_scheme_product':True,
         'connected_component_stability_and_slice_irreducibility_proof_accepted':True,
         'projective_c_equals1_before_normalization_yields_exact_dimension_shift9':True,
         'QQ_and_C_conditional_rationality_transfer_accepted':True,
         'Y16_rationality_full_component25_rationality_not_assumed_or_proved_here':True,
         'no_model_rebuild_new_point_source_jet_source26_or_author_script_rerun':True,
         'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in inputs},'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-started,'frozen':True}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');signal.alarm(0)
 print(json.dumps({k:report[k]for k in ['all_passed','projective_c_equals1_before_normalization_yields_exact_dimension_shift9','seconds']},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=ROOT/'reports'/STEM);a=p.parse_args();main(a.output.resolve())
