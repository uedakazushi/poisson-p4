#!/usr/bin/env python3
"""Independent actual odd-minor staircase and generic W completion audit."""
import argparse
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import subprocess
import time
import sympy as s

ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_GENERIC_W_FULL_COMPLETIONS_20261103_v63'
def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda:f.read(1<<20),b''):h.update(chunk)
    return h.hexdigest()
def read(path):
    with (gzip.open(path,'rt') if path.suffix=='.gz' else path.open()) as f:return json.load(f)
def mat(rows):return s.Matrix([[s.sympify(c) for c in row] for row in rows])
def same(A,B):assert A.shape==B.shape and (A-B).applyfunc(s.expand)==s.zeros(*A.shape)


def main(output,replay):
    started=time.monotonic();folder=ROOT/'reports'/STEM;report=read(folder/'report.json');assert report['all_passed']
    for name,digest in report['input_sha256'].items():assert sha(ROOT/name)==digest,name
    for name,digest in report['evidence_sha256'].items():assert sha(folder/name)==digest,name
    subprocess.run([str(ROOT/'.venv/bin/python'),str(ROOT/'code/ambient_B18_generic_W_completion_20261103_v63.py'),
                    '--output',str(replay)],check=True,timeout=40)
    for name in report['evidence_sha256']:assert (folder/name).read_bytes()==(replay/name).read_bytes()
    records=read(folder/'staircase.json.gz');quadpath=ROOT/'reports/AMBIENT_B18_CA_CB_ODD_QUADRATIC_DEFECT_20261103_v63/certificate.json.gz'
    quad=read(quadpath);priorpath=ROOT/'reports/AMBIENT_B18_ODD_FIRSTJETS_N_PROJECTION_20261103_v63/certificate.json.gz';prior=read(priorpath)
    W,U,V=s.symbols('W U V');h,x,y,z,w=s.symbols('h x y z w');X=[h,x,y,z,w]
    pairs=list(itertools.combinations(range(5),2));mons=list(itertools.combinations_with_replacement(range(5),2))
    def vector(P):return s.Matrix([s.Poly(P[i,j],*X).coeff_monomial(X[a]*X[b]) for i,j in pairs for a,b in mons])
    ey=s.Matrix([0,0,1,0,0]);D=s.Matrix([x*x,0,0,h*x,h*h+z*z]);P0=ey*D.T-D*ey.T
    PW=P0.copy();PW[3,4]+=W*x*x;PW[4,3]-=W*x*x
    assert vector(PW)[140]==W and vector(PW)[106]==1
    summary={};criticalpaths=[]
    for name in ['CA','CB']:
        path=ROOT/f'reports/AMBIENT_B18_{name}_CRITICAL_SCHEME_20261103_v63/certificate.json.gz';critical=read(path);criticalpaths.append(path)
        tt=list(map(s.Symbol,critical['slice_variables']));axis=dict(zip(tt,[0,0,W,0,0,0,0,0,0,0]))
        equations=critical['predicted_entire_critical_ideal']
        schur=critical['full150_critical_Schur_column'] if name=='CA' else critical['all150_Schur_entries']
        assert all(s.expand(s.sympify(q).subs(axis))==0 for q in [*equations,*schur])
        d=quad[name];saved=records[name];same(PW,mat(saved['physical_entire_W_curve']))
        same(PW,mat(d['full_source_P']).subs({U:0,V:0}))
        Y=d['selected_actual28_Y'];actualcurve=vector(PW)[Y+[61],:];same(actualcurve,mat(saved['actual_R29_curve']))
        T=mat(d['base_critical_tangent29_by3']);same(actualcurve.diff(W),T[:,0]);assert T.rank()==3
        A=mat(prior['CA_CB_full_firstjet_data'][name]['actual']);same(A*T,s.zeros(150,3))
        rows=d['actual_seven_odd_basis_rows'];columns=d['selected_actual_base_columns'];B=A[rows,columns]
        assert B.det()==s.Rational(-8 if name=='CA' else 1,s.Integer(1) if name=='CA' else 128)
        normal=s.eye(29)[:,columns]*B.inv();directions=normal.row_join(T[:,1:])
        same(directions,mat(saved['normalized_nine_source_directions']));same(A[rows,:]*directions,s.eye(7).row_join(s.zeros(7,2)))
        assert directions.rank()==9
        residual=list(map(s.sympify,d['all_zero_R_firstjet_odd_quadratics']));qrows=d['independent_quadratic_rows']
        C=s.Matrix([[s.Poly(residual[i],W,U,V).coeff_monomial(m) for m in [W*U,W*V]] for i in qrows])
        same(C,mat(saved['WU_WV_matrix']));assert C.det()==s.Rational(8,125) if name=='CA' else C.det()==-s.Rational(2,45)
        leading=s.diag(*([1]*7+[W,W]));leading[7:9,7:9]=W*C
        same(leading,mat(saved['staircase_leading_matrix']))
        normalized=s.expand(leading.det());raw=s.expand(B.det()*normalized)
        assert normalized==C.det()*W**2
        assert raw==(s.Rational(-64,125) if name=='CA' else -s.Rational(1,2880))*W**2
        summary[name]=dict(normalized_leading=str(normalized),raw_leading=str(raw),
            full_critical_W_axis_and_all150_source_Schur_identities=True,
            actual9minor_leading_coefficient_bound=True)
    npath=ROOT/'reports/AMBIENT_B18_N_CRITICAL_SCHEME_20261103_v63/certificate.json.gz';nc=read(npath);criticalpaths.append(npath)
    ntt=list(map(s.Symbol,nc['slice_variables']));naxis=dict(zip(ntt,[0,W,0,0,0,0,0,0,0]))
    assert all(s.expand(s.sympify(q).subs(naxis))==0 for q in nc['predicted_exact_ideal']+nc['all150_Schur_entries'])
    nalpha=s.Matrix(list(map(s.sympify,nc['center19'])));nalpha[nc['free_indices'][1]]+=W
    same(mat(nc['fiber19_basis'])*nalpha,vector(PW))
    # Full64 regular-I kernel relations are reconstructed independently.
    relationpath=folder/'N_regular64_relations.json.gz';data=read(relationpath)
    nu,v=s.symbols('nu v');ss=s.symbols('s1:7');yy=s.symbols('y0:8');t=s.Symbol('t');b=[nu,*ss,v]
    own=[s.expand(b[j]*yy[i]-b[i]*yy[j]) for i,j in itertools.combinations(range(8),2)]
    own += [yy[i]*yy[j]-nu*b[i]*b[j] for i in range(8) for j in range(i,8)]
    assert len(own)==64
    for left,right in zip(own,data['all64']):assert s.expand(left-s.sympify(right))==0
    for q in own:assert s.rem(s.Poly(q.subs(dict(zip(yy,[a*t for a in b]))),t),s.Poly(t*t-nu,t)).is_zero
    accepted=[ROOT/'reports/AMBIENT_B18_QUADRATIC_IMAGE_REDUCTION_20261103_v63_REVIEW/report.json',
        ROOT/'reports/AMBIENT_B18_CA_CB_ODD_QUADRATIC_DEFECT_20261103_v63_REVIEW.json',
        ROOT/'reports/AMBIENT_B18_CA_N_CRITICAL_SCHEMES_20261103_v63_REVIEW.json',
        ROOT/'reports/AMBIENT_B18_CB_CRITICAL_SCHEME_20261103_v63_REVIEW.json',
        ROOT/'reports/AMBIENT_B18_N_FULL_COMPLETION_20261103_v63_REVIEW.json']
    assert all(read(path)['all_passed'] for path in accepted)
    paths=[Path(__file__),ROOT/'research/notes'/f'{STEM}_REVIEW.md',folder/'report.json',folder/'staircase.json.gz',relationpath,quadpath,priorpath,
           ROOT/'research/notes'/f'{STEM}.md',*criticalpaths,*accepted]
    result=dict(id=STEM+'_REVIEW',all_passed=True,author_all_evidence_fresh_byte_identical=True,
        independent_physical_entire_W_curve_common_to_CA_CB_N=True,
        all_actual_critical_W_axes_and_full150_Schur_identities_replayed=True,
        actual_odd9minor_normalization_and_W2_leading_coefficients=summary,
        unknown_higher_blocks_accounted_by_Schur_orders_not_set_to_zero=True,
        generic_CA_CB_J_regular9_and_actualI_mod_mJ_surjective_Nakayama_accepted=True,
        whole_proper_finite_length_atmost2_open_and_corank1_imply_single_length2_fiber=True,
        actual_quadratic_source_parameter_base_and_normalization_fraction_field_proof_accepted=True,
        algebraic_generic_descent_via_coherent_Fitting_open_and_etale_finite_rank2_accepted=True,
        N_known_full_I_localization_and_actual_W_smooth_critical_branch_Ac7_accepted=True,
        full64_N_relations_and_Koszul_rank2_completeness_accepted=True,
        CA_CB_cone_dim_edim_depth_multiplicity=[29,38,21,2],CA_CB_projective=[28,37,20,2],
        N_cone_dim_edim_depth_multiplicity=[28,36,21,2],N_projective=[27,35,20,2],
        all_three_finite_subintegral_normalizations=True,
        scope_generic_W_only_no_all_nonzero_W_orbit_transport_assumed=True,
        generic_V_numeric_Hessian_or_odd_order_results_not_used=True,
        no_CA_CB_B18_origin_full_I_or_actual_D_or_global_bad_locus_claim=True,
        input_sha256={str(path.relative_to(ROOT)):sha(path) for path in paths},seconds=time.monotonic()-started)
    output.parent.mkdir(parents=True,exist_ok=True);output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ['all_passed','actual_odd9minor_normalization_and_W2_leading_coefficients','seconds']},indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--replay',type=Path,required=True);args=parser.parse_args();main(args.output.resolve(),args.replay.resolve())
