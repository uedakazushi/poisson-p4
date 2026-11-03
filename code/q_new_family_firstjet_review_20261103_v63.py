#!/usr/bin/env python3
"""Independent full seven-variable raw-jet rank25 / primitive beta-image15 audit."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import time
import itertools

import sympy as s
from sympy.polys.matrices import DomainMatrix
from sympy.polys.rings import ring

from q_positive_weight6_small_model_review_20261103_v63 import read,qm,rank

ROOT=Path(__file__).resolve().parents[1]


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=False);started=time.monotonic()
    packet=ROOT/'reports/Q_NEW_FAMILY_FIRSTJET_IMAGE_20261103_v63';old=read(packet/'report.json');assert old['all_passed']
    for p,digest in old['input_sha256'].items():assert sha(ROOT/p)==digest
    for p,digest in old['evidence_sha256'].items():assert sha(packet/p)==digest
    data=read(packet/'certificate.json.gz')
    familypath=ROOT/'reports/Q_POSITIVE_WEIGHT6_NEW_FAMILY_20261103_v63/certificate.json.gz'
    family=read(familypath);n=s.symbols('n0:3');p=s.symbols('p0:4');variables=[*n,*p]
    x=s.Matrix(s.symbols('x y z w'));t=s.Matrix(s.symbols('t0:4'));pairs=list(itertools.combinations(range(4),2))
    P=s.Matrix([[s.sympify(v)for v in row]for row in family['full_affine_P']])
    C=P.subs(dict(zip(x,p)),simultaneous=True)
    L=sum((P.diff(x[i]).subs(dict(zip(x,p)),simultaneous=True)*t[i]for i in range(4)),s.zeros(4))
    def extract(polynomial):
        zero=dict.fromkeys(t,0)
        return s.Matrix([polynomial[i,j].subs(zero)for i,j in pairs]+[
            s.diff(polynomial[i,j],tt).subs(zero)for i,j in pairs for tt in t])
    jet=C+L;columns=[extract(jet.diff(v))for v in variables]
    labels=list(map(str,variables))
    # Derive both linear and perspective actions from vector fields directly.
    # The perspective field is quadratic; extracting only its first jet is
    # independent of copying the producer's closed formula.
    for i,j in itertools.product(range(4),repeat=2):
        vector=s.zeros(4,1);vector[i]=t[j];derivative=vector.jacobian(t)
        delta=derivative*jet+jet*derivative.T-sum((vector[k]*jet.diff(t[k])for k in range(4)),s.zeros(4))
        columns.append(extract(delta));labels.append(f'GL4_{i}_{j}')
    for i in range(4):
        vector=-t[i]*t;derivative=vector.jacobian(t)
        delta=derivative*jet+jet*derivative.T-sum((vector[k]*jet.diff(t[k])for k in range(4)),s.zeros(4))
        columns.append(extract(delta));labels.append(f'perspective_{i}')
    columns.append(extract(jet));labels.append('P_scalar')
    M=s.Matrix.hstack(*columns).applyfunc(s.expand)
    assert labels==data['source_column_labels']and M.shape==(30,28)
    expected=s.Matrix([[s.sympify(v)for v in row]for row in data['full30_by28_polynomial_matrix']]);assert M==expected
    R,*rv=ring('n0,n1,n2,p0,p1,p2,p3',s.QQ)
    A=[[R.from_expr(v)for v in row]for row in M.tolist()];rowids=list(range(30));colids=list(range(28))
    for step in data['twenty_rational_unit_pivots']:
        i,j=step['current_row'],step['current_column'];c=s.Rational(step['pivot'])
        assert step['original_row']==rowids[i]and step['original_column']==colids[j]
        assert A[i][j]==R.from_expr(c)and c!=0
        pivot=A[i];A=[[v-row[j]/c*pivot[k]for k,v in enumerate(row)if k!=j]for ii,row in enumerate(A)if ii!=i]
        del rowids[i];del colids[j]
    assert len(data['twenty_rational_unit_pivots'])==20
    assert rowids==data['remaining_original_rows']and colids==data['remaining_original_columns']
    assert len(A)==10 and len(A[0])==8
    assert all(v==R.from_expr(s.sympify(data['residual10_by8_polynomial_matrix'][i][j]))for i,row in enumerate(A)for j,v in enumerate(row))
    # Different exact backend: full rational function field seven variables.
    K=s.QQ.frac_field(*variables);fieldmatrix=[[K.from_sympy(v.as_expr())for v in row]for row in A]
    reduced,pivots=DomainMatrix(fieldmatrix,(10,8),K).rref();assert len(pivots)==5
    assert all(v==K.zero for row in reduced.to_list()[5:]for v in row)
    (args.output/'own-field-rref.json').write_text(json.dumps(dict(field=str(K),pivot_columns=list(pivots),
        rref=[[str(v)for v in row]for row in reduced.to_list()],
        all_last_five_rows_exactly_zero=True),indent=2)+'\n')
    # Also replay the audited native function-field engine on our own input.
    native=ROOT/'code/rational_function_matrix_rank_20261103_v62.cpp';binary=ROOT/'work/q-new-family-firstjet-review-native'
    flags=subprocess.check_output(['pkg-config','--cflags','--libs','flint'],text=True).split()
    subprocess.run(['g++','-O2','-std=c++17',str(native),*flags,'-o',str(binary)],check=True,timeout=30)
    inputpath=args.output/'own-native-input.txt'
    with inputpath.open('w')as stream:
        stream.write('7 10 8\n')
        for row in A:
            for value in row:
                den=s.ilcm(1,1,*(s.denom(s.Rational(c))for c in value.values()))
                stream.write(str(len(value))+'\n')
                for exponent,c in sorted(value.items()):stream.write(str(int(s.Rational(c)*den))+' '+' '.join(map(str,exponent))+'\n')
                stream.write('1\n'+str(den)+' 0 0 0 0 0 0 0\n')
    run=subprocess.run([str(binary),str(inputpath)],capture_output=True,text=True,timeout=30,check=True)
    assert 'EXACT_FUNCTION_FIELD_RANK 5 ROWS 10 COLUMNS 8'in run.stdout
    (args.output/'own-native.log').write_text(run.stdout)
    sample={v:s.Rational(data['sample'][str(v)])for v in variables};SM=M.subs(sample)
    minor=SM[data['lower25_minor_rows'],data['lower25_minor_columns']]
    det=DomainMatrix.from_Matrix(minor).convert_to(s.QQ).det()
    assert det==s.Rational(data['lower25_minor_determinant'])!=0
    assert rank(SM[:6,7:23])==6 and SM[:6,23:27]==s.zeros(6,4)and rank(SM[:,23:27])==4
    cv=C.subs(sample);pf=cv[0,1]*cv[2,3]-cv[0,2]*cv[1,3]+cv[0,3]*cv[1,2]
    assert pf==-s.Rational(16,81)
    (args.output/'dimension-certificate.json').write_text(json.dumps(dict(
        own_raw30_by28_polynomial_matrix_rebuilt=True,all7source_variables_retained=True,
        own_perspective_directions_from_quadratic_vectorfields=True,
        all20_rational_unit_Schur_steps_checked=True,own_seven_variable_field_rank=5,
        raw_jet_image_dimension=25,value_nuisance_dimension=6,perspective_nuisance_dimension=4,
        primitive_beta_image_dimension=15,irreducible_hypersurface_image=True,
        full_new_component_S61_dimension=16,generic_fiber_over_its_beta_image_dimension=1,
        general_W16_dominance=False,generic_five_point_counterexample=False,
        unfolded_hypersurface_polynomial_or_defining_prime_not_computed=True),indent=2)+'\n')
    paths=[Path(__file__),native,packet/'report.json',packet/'certificate.json.gz',familypath,
        ROOT/'reports/Q_POSITIVE_WEIGHT6_NEW_FAMILY_20261103_v63_REVIEW/report.json',
        ROOT/'research/notes/Q_NEW_FAMILY_FIRSTJET_IMAGE_20261103_v63_REVIEW.md']
    result=dict(id='Q_NEW_FAMILY_FIRSTJET_IMAGE_20261103_v63_REVIEW',all_passed=True,
        own_raw30_matrix_and_full_group_directions=True,all20_unit_operations_exact=True,
        two_exact_seven_variable_function_field_rank_engines_agree=True,
        exact_primitive_image_dimension15_not_dominant=True,
        new_component_generic_fiber_dimension1_over_its_image=True,
        no_generic_upper_from_samples=True,hypersurface_recognition_not_yet_computed=True,
        input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},
        evidence_sha256={p.name:sha(p)for p in args.output.iterdir()if p.is_file()and p.name!='report.json'},seconds=time.monotonic()-started)
    (args.output/'report.json').write_text(json.dumps(result,indent=2)+'\n')
    print('INDEPENDENT RAWJET25 / PRIMITIVE BETA15 REVIEW PASS',result['seconds'],flush=True)


if __name__=='__main__':main()
