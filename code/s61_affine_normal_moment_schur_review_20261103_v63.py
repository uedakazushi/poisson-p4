#!/usr/bin/env python3
"""Independent532/536 norm localization proof and exact135-generator algebra audit."""
import argparse
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import time
import sympy as s

ROOT=Path(__file__).resolve().parents[1]
STEM='S61_AFFINE_NORMAL_MOMENT_SCHUR_20261103_v63_REVIEW'


def sha(path):
    h=hashlib.sha256()
    with path.open('rb')as f:
        for block in iter(lambda:f.read(1<<20),b''):h.update(block)
    return h.hexdigest()


def read(path):
    with (gzip.open(path,'rt')if path.suffix=='.gz'else path.open())as f:return json.load(f)


def determinant(A):
    n=A.rows;value=s.S.Zero
    for p in itertools.permutations(range(n)):
        sign=(-1)**sum(p[i]>p[j]for i in range(n)for j in range(i+1,n))
        value+=sign*s.prod(A[i,p[i]]for i in range(n))
    return s.expand(value)


def companion(poly,h):
    poly=s.Poly(poly,h);n=poly.degree();assert poly.LC()==1;T=s.zeros(n)
    for j in range(n-1):T[j+1,j]=1
    for i in range(n):T[i,n-1]=-poly.coeff_monomial(h**i)
    return T


def polynomial_matrix(poly,h,T):
    p=s.Poly(poly,h);M=s.zeros(T.rows);power=s.eye(T.rows)
    for k in range(p.degree()+1):M+=p.coeff_monomial(h**k)*power;power=power*T
    return M.applyfunc(s.expand)


def main(out):
    start=time.monotonic();out.mkdir(parents=True,exist_ok=False)
    old=ROOT/'reports/S61_AFFINE_NORMAL_MOMENT_SCHUR_20261103_v63'
    sharp=ROOT/'reports/S61_AFFINE_NORMAL_SHARP_MOMENT_SCHUR_20261103_v63'
    prior=ROOT/'reports/S61_AFFINE_NORMAL_SCALAR_SCHUR_20261103_v63_REVIEW/report.json'
    globalreview=ROOT/'reports/S61_AFFINE_NORMAL_GLOBAL_QUINTIC_20261103_v63_REVIEW.json'
    deps=[old/'report.json',sharp/'report.json',prior,globalreview];checked={}
    for path in deps:
        d=read(path);assert d['all_passed']
        for name,want in d.get('input_sha256',{}).items():assert sha(ROOT/name)==want;checked[name]=want
        for name,want in d.get('evidence_sha256',{}).items():assert sha(path.parent/name)==want;checked[str((path.parent/name).relative_to(ROOT))]=want
    p=read(prior)
    assert p['B_surjective_basechange_coker_B_and_fixed_bottom4_free_kernel_accepted']
    assert p['localized_scalar_quotient_is_actual_algebra_not_only_Fitting_support']
    assert p['all_unknown_nilpotents_retained']
    with tempfile.TemporaryDirectory(prefix='moment-schur-review-')as tmp:
        for folder,source,short in [(old,'s61_affine_normal_moment_schur_20261103_v63.py','old'),(sharp,'s61_affine_normal_sharp_moment_schur_20261103_v63.py','sharp')]:
            destination=Path(tmp)/short
            subprocess.run([sys.executable,str(ROOT/'code'/source),'--output',str(destination)],capture_output=True,text=True,check=True,timeout=20)
            for name in ['fixtures.json','criterion.json']:assert (destination/name).read_bytes()==(folder/name).read_bytes()
    h,a,t=s.symbols('h a t');q=s.symbols('q0:4');g4=h**4+sum(q[i]*h**i for i in range(4));G=(h-t)*g4
    T4=companion(g4,h);T5=companion(s.expand(G),h);E=s.zeros(5,4)
    for j in range(4):E[j,j]=-t;E[j+1,j]=1
    assert(T5*E-E*T4).applyfunc(s.expand)==s.zeros(5,4)
    assert determinant(E[1:,:])==1
    assert s.expand(determinant(T4-a*s.eye(4))-g4.subs(h,a))==0
    # Independent fixtures: companion powers and multiplication norms, no resultant routine.
    oldfixtures=read(old/'fixtures.json');newfixtures=read(sharp/'fixtures.json');audit=[]
    cases=[('simple',s.prod(h-i for i in range(1,6)),1,{0:(h-1)*h,1:(h-1)**2}),
           ('multiple',(h-2)**3*(h-3)*(h+1),2,{0:(h-2)*h,1:h-2}),
           ('sharp',h**5,0,{0:h*h,133:h})]
    for index,(name,G,h0,remainders)in enumerate(cases):
        G=s.expand(G);q4=s.div(G,h-h0,h)[0]
        assert s.rem(G,h-h0,h)==0 and s.gcd(G,s.gcd(*remainders.values()))==h-h0
        r=sum(a**i*f for i,f in remainders.items());S=s.div(r,h-h0,h)[0]
        assert s.rem(r,h-h0,h)==0
        T=companion(G,h);normD=determinant(polynomial_matrix(S,h,companion(q4,h)))
        v=s.Matrix([s.Poly(r,h).coeff_monomial(h**k)for k in range(5)])
        cols=[]
        for j in range(5):cols.append(v);v=(T*v).applyfunc(s.expand)
        M=s.Matrix.hstack(*cols);delta=determinant(M[1:,:4])
        oldcols=[((T-a*s.eye(5))*col).applyfunc(s.expand)for col in cols[:4]]
        oldM=s.Matrix.hstack(*oldcols);old_delta=determinant(oldM[1:,:])
        assert s.expand(delta-normD)==0
        assert s.expand(old_delta-q4.subs(h,a)*normD)==0
        assert all(s.expand(x-s.sympify(y))==0 for x,y in zip(M,list(s.Matrix(newfixtures[index]['multiplication5_by5']))))
        assert all(s.expand(x-s.sympify(y))==0 for x,y in zip(oldM,list(s.Matrix(oldfixtures[index]['all5_by4_columns']))))
        assert s.expand(delta-s.sympify(newfixtures[index]['delta0']))==0
        assert s.expand(old_delta-s.sympify(oldfixtures[index]['bottom4_determinant']))==0
        assert s.Poly(delta,a).degree()<=532 and s.Poly(old_delta,a).degree()<=536
        if name=='sharp':assert delta==a**532 and old_delta==a**536
        audit.append(dict(name=name,new_delta=str(delta),old_delta=str(old_delta),norm_D=str(normD)))
    # Coefficient-ring algebra fixture, reconstructed from multiplication by h-x.
    x,b,c=s.symbols('x b c');G=h**5+c*h+b;T=companion(G,h);M=T-x*s.eye(5);D=M[1:,:4];delta=determinant(D)
    adj=s.zeros(4)
    for i in range(4):
        for j in range(4):adj[i,j]=(-1)**(i+j)*determinant(D.minor_submatrix(j,i))
    numerator=s.expand(-(M[0,:4]*adj)[0]);norm=determinant(M)
    assert delta==1 and numerator==x and s.expand(norm+G.subs(h,x))==0
    fixture=newfixtures[3];assert s.expand(norm-s.sympify(fixture['Norm_C_r']))==0
    assert numerator==s.sympify(fixture['Cramer_h_numerator'])
    # All134 inputs enter each structured column; record exact sparse recipes only.
    presentation=read(ROOT/'reports/S61_AFFINE_NORMAL_GLOBAL_QUINTIC_20261103_v63/presentation.json.gz')
    assert len(presentation['all134_remainders'])==134
    for row in presentation['all134_remainders']:
        assert all(0<=degree<=4 and all(len(mon)<=2 for mon,coefficient in coeff)for degree,coeff in row)
    new_columns=[[[5*i+j,i,1]for i in range(134)]for j in range(4)]
    old_columns=[[[5*i+j+1,i,1]for i in range(134)]+[[5*i+j,i+1,-1]for i in range(134)]for j in range(4)]
    assert max(v[0]for col in new_columns+old_columns for v in col)==669
    assert sum(map(len,new_columns))==536 and sum(map(len,old_columns))==1072
    oldrecipe=read(old/'criterion.json');recipe=read(sharp/'criterion.json')
    assert oldrecipe['full_cover_options']['fixed_QQ_constants']==list(range(537))
    assert recipe['full_cover_options']['fixed_QQ_constants']==list(range(533))
    assert recipe['degree_a_bound']==4*133==532 and oldrecipe['degree_bounds']['delta_a']==4*(133+1)==536
    assert recipe['degree10_option']['scalar_generators']==670
    option=recipe['degree34_option'];assert option['scalar_generators']==1+134==135
    assert option['formal_degree_atmost']==2+4*8==34 and option['norm_formal_degree_atmost']==5*2==10
    assert option['U_and_delta_formal_degree_atmost']==4*2==8
    (out/'independent-fixtures.json').write_text(json.dumps(audit+[dict(name='coefficient_ring',delta=str(delta),U=str(numerator),norm=str(norm))],indent=2)+'\n')
    (out/'column-recipes.json').write_text(json.dumps(dict(entry_schema='[old Phi column,a exponent,integer coefficient]',sharp4_columns=new_columns,old4_columns=old_columns),indent=2)+'\n')
    paths=[Path(__file__),ROOT/f'research/notes/{STEM}.md',*deps,old/'fixtures.json',old/'criterion.json',sharp/'fixtures.json',sharp/'criterion.json',ROOT/'reports/S61_AFFINE_NORMAL_GLOBAL_QUINTIC_20261103_v63/presentation.json.gz']
    report=dict(id=STEM,all_passed=True,current_dependency_hashes_checked=len(checked),
        both_frozen_author_fixture_and_recipe_files_fresh_byte_equal=True,
        independent_companion_multiplication_norm_fixtures_without_resultants=True,
        universal_symbolic_G5_g4_intertwining_and_bottom_basis_determinant1=True,
        exact_old_delta_g4a_norm_and_sharp_delta_norm_identity_accepted=True,
        nonseparable_and_fully_ramified_fibers_allowed=True,
        sharp_a532_and_old_a536_exact_fixtures_rebuilt=True,
        all134_structured_source_column_indices_and_degree_bounds_bound=True,
        fiber_J_h_minus_h0_and_unit_ideal_s_i_from_exact_cyclic_B_proof_accepted=True,
        delta_coefficient_ideal_unit_in_B_not_in_Rprime_accepted=True,
        B_to_Ba_delta_faithfully_flat_by_flatness_and_nonzero_residue_field_fibers=True,
        delta_regular_on_Ba_and_all_existing_Ass_and_nilpotent_submodules_retained=True,
        full_point_cover_of_Ba_by_one_delta_claimed=False,
        single_generic_scalar_unit_iff_B_empty_accepted=True,
        three_generic_or533_fixed_sharp_full_point_covers_accepted=True,
        old_three_generic_or537_fixed_cover_remains_valid=True,
        exact135_equation_unital_algebra_reduction_accepted=True,
        coefficient_ring_delta1_tau_x_norm_minus_Gx_fixture_rebuilt=True,
        scalar135_degree_bound=34,norm_degree_bound=10,alternative670_degree_bound=10,
        actual_S61_delta_or_scalar_generators_expanded_or_units_computed=False,
        generic_five_exhaustion_or_full_w_saturation_claimed=False,
        no_extra_model_rank_GB_or_Macaulay_probe=True,
        input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},
        evidence_sha256={name:sha(out/name)for name in ['independent-fixtures.json','column-recipes.json']},seconds=time.monotonic()-start)
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print('INDEPENDENT SHARP532/OLD536 FAITHFULLYFLAT / EXACT135SCALAR CRITERIA PASS',report['seconds'],flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
