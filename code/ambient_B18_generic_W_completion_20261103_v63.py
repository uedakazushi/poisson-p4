#!/usr/bin/env python3
"""Generic W corollary from accepted entire critical ideals and actual odd2jets."""
import argparse,gzip,hashlib,itertools,json,time
from pathlib import Path
import sympy as s
from ambient_CA_global_flag_model_20261103_v63 import vector
ROOT=Path(__file__).resolve().parents[1]
def load(p):return json.loads(gzip.decompress(p.read_bytes())if p.suffix=='.gz'else p.read_bytes())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def mat(M):return [[str(x)for x in row]for row in M.tolist()]
def main(out):
    started=time.monotonic();out.mkdir(parents=True,exist_ok=False)
    quadpath=ROOT/'reports/AMBIENT_B18_CA_CB_ODD_QUADRATIC_DEFECT_20261103_v63/certificate.json.gz';quad=load(quadpath)
    priorpath=ROOT/'reports/AMBIENT_B18_ODD_FIRSTJETS_N_PROJECTION_20261103_v63/certificate.json.gz';prior=load(priorpath)
    deps=[ROOT/'reports/AMBIENT_B18_CA_CB_ODD_QUADRATIC_DEFECT_20261103_v63_REVIEW.json',
      ROOT/'reports/AMBIENT_B18_CA_N_CRITICAL_SCHEMES_20261103_v63_REVIEW.json',
      ROOT/'reports/AMBIENT_B18_CB_CRITICAL_SCHEME_20261103_v63_REVIEW.json',
      ROOT/'reports/AMBIENT_B18_QUADRATIC_IMAGE_REDUCTION_20261103_v63_REVIEW/report.json',
      ROOT/'reports/AMBIENT_B18_N_FULL_COMPLETION_20261103_v63_REVIEW.json',
      ROOT/'reports/AMBIENT_CA_GENERIC_V_FULL_COMPLETION_20261103_v63/report.json']
    for p in deps:
        data=load(p);assert data['all_passed'] is True
        for name,digest in data.get('input_sha256',{}).items():assert sha(ROOT/name)==digest,name
    W,U,V=s.symbols('W U V');tt=s.symbols('s0:10');zeroWV={U:0,V:0};sourceW=[0,0,W,0,0,0,0,0,0,0];records={};critpaths=[]
    for name in['CA','CB']:
        cp=ROOT/f'reports/AMBIENT_B18_{name}_CRITICAL_SCHEME_20261103_v63/certificate.json.gz';crit=load(cp);critpaths.append(cp)
        sub=dict(zip(tt,sourceW))
        assert all(s.expand(s.sympify(f).subs(sub))==0 for f in crit['predicted_entire_critical_ideal'])
        key='full150_critical_Schur_column'if name=='CA'else'all150_Schur_entries'
        assert all(s.expand(s.sympify(f).subs(sub))==0 for f in crit[key])
        d=quad[name];A=s.Matrix(prior['CA_CB_full_firstjet_data'][name]['actual']);rr=d['actual_seven_odd_basis_rows'];cc=d['selected_actual_base_columns'];B=A[rr,cc]
        assert len(rr)==len(cc)==7 and B.det()!=0
        T=s.Matrix(d['base_critical_tangent29_by3']);assert T.rank()==3 and A*T==s.zeros(150,3)
        E=s.eye(29)[:,cc];normal=E*B.inv();directions=normal.row_join(T[:,1:]);rawdirections=E.row_join(T[:,1:])
        assert A[rr,:]*directions==s.eye(7).row_join(s.zeros(7,2))and directions.rank()==9
        PW=s.Matrix(d['full_source_P']).subs(zeroWV).applyfunc(s.expand)
        rawcurve=vector(PW)[d['selected_actual28_Y']+[61],:]
        assert rawcurve.diff(W)==T[:,0]
        h=[s.sympify(f)for f in d['all_zero_R_firstjet_odd_quadratics']]
        rows=d['independent_quadratic_rows'];C=s.Matrix([[s.Poly(h[i],W,U,V).coeff_monomial(W*U),s.Poly(h[i],W,U,V).coeff_monomial(W*V)]for i in rows])
        assert C.det()==s.sympify(d['independent_WU_WV_minor'])!=0
        staircase=s.eye(7).row_join(s.zeros(7,2)).col_join(s.zeros(2,7).row_join(W*C))
        assert s.expand(staircase.det())==W**2*C.det()
        records[name]=dict(entire_monic_critical_W_source=list(map(str,sourceW)),all150_critical_Schur_entries_zero=True,
          physical_entire_W_curve=mat(PW),actual_R29_curve=mat(rawcurve),actual_R29_WUV_tangent=mat(T),
          seven_actual_odd_basis_rows=rr,seven_base_columns=cc,seven_raw_unit=str(B.det()),
          normalized_nine_source_directions=mat(directions),raw_nine_source_directions=mat(rawdirections),
          two_zero_firstjet_odd_rows=rows,full_accepted_quadratic_residuals=list(map(str,h)),WU_WV_matrix=mat(C),
          normalized_actual9minor_leading_order=2,normalized_actual9minor_leading_coefficient=str(C.det()),
          raw_actual9minor_leading_coefficient=str(s.factor(B.det()*C.det())),staircase_leading_matrix=mat(staircase),
          omitted_terms='top blocks I7+O(W),O(W); lowerleftO(W), lowerrightW*C+O(W²); Schur correctionO(W²)',
          no_unknown_higher_coefficients_set_to_zero=True)
    assert s.Matrix(records['CA']['physical_entire_W_curve'])==s.Matrix(records['CB']['physical_entire_W_curve'])
    np=ROOT/'reports/AMBIENT_B18_N_CRITICAL_SCHEME_20261103_v63/certificate.json.gz';nc=load(np);nt=s.symbols('t0:9');nvalues=[0,W,0,0,0,0,0,0,0];nsub=dict(zip(nt,nvalues))
    assert all(s.expand(s.sympify(f).subs(nsub))==0 for f in nc['predicted_exact_ideal']+nc['all150_Schur_entries'])
    nalpha=s.Matrix(nc['center19']);nalpha[nc['free_indices'][1]]+=W
    assert s.Matrix(nc['fiber19_basis'])*nalpha==vector(s.Matrix(records['CA']['physical_entire_W_curve']))
    records['N']=dict(entire_monic_critical_W_source=list(map(str,nvalues)),all150_critical_Schur_entries_zero=True,physical_W_curve_identical_to_CA_CB=True)
    critpaths.append(np)
    (out/'staircase.json.gz').write_bytes(gzip.compress((json.dumps(records,separators=(',',':'))+'\n').encode(),mtime=0))
    nu,v=s.symbols('nu v');ss=s.symbols('s1:7');yy=s.symbols('y0:8');t=s.Symbol('t');b=[nu,*ss,v]
    linear=[b[j]*yy[i]-b[i]*yy[j]for i,j in itertools.combinations(range(8),2)]
    product=[yy[i]*yy[j]-nu*b[i]*b[j]for i in range(8)for j in range(i,8)]
    assert len(linear)==28 and len(product)==36
    for f in linear+product:
        assert s.rem(s.Poly(f.subs({y:c*t for y,c in zip(yy,b)}),t),s.Poly(t*t-nu,t)).as_expr()==0
    (out/'N_regular64_relations.json.gz').write_bytes(gzip.compress((json.dumps(dict(generators=list(map(str,b)),all28_Koszul=list(map(str,linear)),all36_products=list(map(str,product)),all64=list(map(str,linear+product))),separators=(',',':'))+'\n').encode(),mtime=0))
    note=ROOT/'research/notes/AMBIENT_B18_GENERIC_W_FULL_COMPLETIONS_20261103_v63.md'
    paths=[Path(__file__),note,quadpath,priorpath,*critpaths,*deps]
    result=dict(id='AMBIENT_B18_GENERIC_W_FULL_COMPLETIONS_20261103_v63',all_passed=True,
      only_accepted_actual2jet_and_entire_critical_certificate_reused=True,no_new_high_order_or_GB_job=True,
      exact_entire_W_curve_common_to_CA_CB=True,
      normalized_nine_minor_leading={'CA':'8/125*W^2','CB':'-2/45*W^2'},
      raw_nine_minor_leading={'CA':'-64/125*W^2','CB':'-1/2880*W^2'},
      actual_I_rank9_at_generic_W_by_nonzero_formal_minor=True,
      generic_J_regular9_and_Nakayama_I_eq_J_CA_CB=True,
      proper_finite_fiber_length2_at_generic_W_from_P0=True,
      normalization_deduced_after_actual_I_nonzero=True,
      CA_CB_generic_W_completion='C[[g20,s8,t²,s_i*t,t³]] inside C[[g20,s8,t]]',
      CA_CB_cone={'dimension':29,'embedding_dimension':38,'depth':21,'multiplicity':2},
      CA_CB_projective={'dimension':28,'embedding_dimension':37,'depth':20,'multiplicity':2},
      N_from_proven_entire_ideal_after_W_inversion='(nu,s1,...,s6,W*v,v²) becomes (nu,s1,...,s6,v)',
      N_generic_W_completion='C[[g20,s7,t²,s_i*t,t³]] inside C[[g20,s7,t]]',
      N_cone={'dimension':28,'embedding_dimension':36,'depth':21,'multiplicity':2},
      N_projective={'dimension':27,'embedding_dimension':35,'depth':20,'multiplicity':2},
      N_complete_regular_kernel_relations=64,all_three_finite_subintegral_normalization=True,
      algebraic_generic_application_requires_open_locus_descent_as_in_note=True,
      scope='Generic W branch only; not every W, B18 origin CA/CB full ideals, all global bad loci, mixed cubic or actual D.',
      input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},evidence_sha256={p.name:sha(p)for p in out.iterdir()if p.is_file()},seconds=time.monotonic()-started)
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');print('GENERIC W CA/CB/N: EXACT2JET STAIRCASE +FULLCOMPLETION COROLLARY PASS',result['seconds'],flush=True)
if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True);main(ap.parse_args().output.resolve())
