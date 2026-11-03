#!/usr/bin/env python3
"""Independent exact degree22 restriction review on old points only."""
import argparse
from fractions import Fraction as Q
import hashlib
import itertools
import json
from pathlib import Path
import resource
import signal
import subprocess
import time

from q_firstjet_spin4_mixed_discriminant_review_20261103_v63 import mat,mm,tp,det,inv

ROOT=Path(__file__).resolve().parents[1]
ID='Q_NEW_FAMILY_DEGREE22_INVARIANT_EXCLUSION_20261103_v63_REVIEW'


def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def necklace(n):
    return sorted({min(word[i:]+word[:i]for i in range(n))for word in itertools.product([0,1],repeat=n)})


def main(output):
    started=time.monotonic();signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
    output.mkdir(parents=True,exist_ok=False)
    author=ROOT/'reports/Q_NEW_FAMILY_DEGREE22_INVARIANT_EXCLUSION_20261103_v63/report.json'
    probe=ROOT/'reports/Q_NEW_FAMILY_DEGREE22_INVARIANT_PROBE_20261103_v63/report.json'
    primary_exact=author.parent/'exact-evaluation.json'
    source36=ROOT/'reports/Q_NEW_FAMILY_DEGREE20_INVARIANT_EXCLUSION_20261103_v63/samples.json'
    source36review=ROOT/'reports/Q_NEW_FAMILY_DEGREE20_INVARIANT_EXCLUSION_20261103_v63_REVIEW/report.json'
    old18=ROOT/'reports/Q_NEW_FAMILY_DEGREE18_INVARIANT_EXCLUSION_20261103_v63/samples.json'
    old18review=ROOT/'reports/Q_NEW_FAMILY_DEGREE18_INVARIANT_EXCLUSION_20261103_v63_REVIEW/report.json'
    old14=ROOT/'reports/Q_NEW_FAMILY_DEGREE12_14_INVARIANT_EXCLUSION_20261103_v63/certificate.json'
    invreview=ROOT/'reports/Q_NEW_FAMILY_INVARIANT_DEGREE16_20261103_v63_REVIEW.json'
    charreview=ROOT/'reports/Q_FIRSTJET_TWO_ADJOINT_RELATIONS_20261103_v63_REVIEW/report.json'
    invpath=ROOT/'reports/Q_NEW_FAMILY_HYPERSURFACE_RECOGNITION_20261103_v63/invariant-certificates.json'
    ar,pr,exact,oldpoints=map(read,[author,probe,primary_exact,source36])
    assert ar['all_passed']and pr['all_passed']and exact['all_passed']
    checked={}
    for p in [author,probe,source36review,old18review,invreview,charreview]:
        data=read(p);assert data['all_passed'],p
        for key in ['input_sha256','evidence_sha256']:
            for name,want in data.get(key,{}).items():
                found=next((x for x in [ROOT/name,p.parent/name]if x.is_file()and sha(x)==want),None)
                assert found is not None,(p,name)
                checked[str(found.relative_to(ROOT))]=want
    assert read(charreview)['character_only_Inv22']==21
    assert read(invreview)['theta_and_H_same_conjugation_and_constant_Poisson_bracket_invariance_proved']
    assert read(source36review)['all_passed']and read(old18review)['all_passed']
    assert len(oldpoints['actual_QQ_source_points'])==36 and pr['existing_point_count']==36 and pr['new_beta_points_generated']==0
    assert pr['same_accepted_QQ_samples']==oldpoints['actual_QQ_source_points']
    rows,columns=pr['selected_rows'],pr['selected_columns'];assert len(rows)==len(columns)==len(set(rows))==len(set(columns))==21
    assert all(0<=i<36 for i in rows)and all(0<=j<24 for j in columns)
    selected=[oldpoints['actual_QQ_source_points'][i]for i in rows]
    saved_samples=read(author.parent/'samples.json')
    assert saved_samples['existing_accepted_QQ_source_points']==selected
    assert saved_samples['probe_rows']==rows and saved_samples['selected_columns']==columns
    assert all(Q(p['value_pivot'])and Q(p['value_pfaffian'])for p in selected)
    w8=necklace(4)[:5];assert len(w8)==5
    d18=read(old18);w18=[tuple(w)for w in d18['selected_words']]
    assert w18==[necklace(9)[i]for i in [6,10,12,18,19,21,23,26,32]]
    d14=read(old14);pairs=[d14['degree14_bracket_pairs'][i]for i in d14['degree14_basis_columns']]
    assert pairs==[[1,2],[1,3],[2,3]]
    definitions=[{'kind':'F4_times_Inv18_basis','basis_index':i,'trace9_word':list(w)}for i,w in enumerate(w18)]
    definitions += [{'kind':'Inv8_times_Inv14','basis8':i,'basis14':j,'trace4_word':list(w8[i]),'bracket_pair':pairs[j]}for i in range(5)for j in range(3)]
    assert definitions==pr['definitions']
    info=read(invpath);omega=mat(info['symplectic_form16']);oi=inv(omega)
    reps=list(map(mat,info['all_ten_representation_matrices']))
    J=mat([[0,1,0,0],[-1,0,0,0],[0,0,0,1],[0,0,-1,0]])
    generators=[]
    for i in range(4):
        for j in range(i,4):
            symmetric=[[Q(0)for _ in range(4)]for _ in range(4)];symmetric[i][j]=symmetric[j][i]=1
            generators.append(mm(J,symmetric))
    gram=[[sum(mm(A,B)[i][i]for i in range(4))for B in generators]for A in generators];gi=inv(gram)
    assert mm(omega,oi)==[[Q(i==j)for j in range(16)]for i in range(16)]
    assert gi==tp(gi)
    assert all(mm(omega,A)==tp(mm(omega,A))for A in reps)
    # Rebuild the complete exact-native input semantically before giving it to
    # an independent forward-dual evaluator.
    inp=output/'own-native-input.txt'
    with inp.open('w')as stream:
        def write(A):
            stream.write(f'{len(A)} {len(A[0])}\n')
            for row in A:stream.write(' '.join(map(str,row))+'\n')
        stream.write('21 24\n')
        for A in [omega,oi,gi,*reps]:write(A)
        for w in [*w8,*w18]:stream.write(''.join(map(str,w))+'\n')
        for pair in pairs:stream.write(' '.join(map(str,pair))+'\n')
        stream.write(' '.join(map(str,columns))+'\n')
        for point in selected:write([[Q(x)]for x in point['beta16']])
    assert inp.read_text().split()==(author.parent/'native-input.txt').read_text().split()
    cpp=ROOT/'code/q_new_family_degree22_review_20261103_v63.cpp';exe=output/'own-native'
    flags=subprocess.check_output(['pkg-config','--cflags','--libs','flint'],text=True).split()
    compile_run=subprocess.run(['c++','-O2','-std=c++17',str(cpp),*flags,'-o',str(exe)],check=True,capture_output=True,text=True)
    (output/'compile.log').write_text(compile_run.stdout+compile_run.stderr)
    resultpath=output/'own-exact-evaluation.json'
    cp=subprocess.run([str(exe),str(inp),str(resultpath)],check=True,capture_output=True,text=True,timeout=50)
    (output/'own-native.log').write_text(cp.stdout+cp.stderr)
    assert 'OWN_FORWARD_DUAL_QQ_ALL378_FACTORS_ALL504_PRODUCTS_21_MINOR_PASS' in cp.stdout
    own=read(resultpath);assert own['all_passed']
    assert own['evaluation']==exact['evaluation']
    assert own['factors']==exact['factors_F4_Inv8five_Inv14three_Inv18nine']
    assert len(own['evaluation'])==21 and all(len(row)==24 for row in own['evaluation'])
    assert len(own['factors'])==21 and all(len(row)==18 for row in own['factors'])
    assert len(own['gradients_of_Inv8_indices1_2_3'])==21
    for factor,value in zip(own['factors'],own['evaluation']):
        numbers=list(map(Q,factor))
        expected=[numbers[0]*x for x in numbers[9:]]+[x*y for x in numbers[1:6]for y in numbers[6:9]]
        assert expected==list(map(Q,value))
    # A second exact determinant algorithm, using Python Fraction elimination.
    determinant=det([[Q(row[j])for j in columns]for row in own['evaluation']])
    assert determinant==Q(own['exact_selected21_minor'])==Q(exact['exact_selected21_minor'])and determinant
    assert abs(determinant.numerator).bit_length()==4922 and determinant.denominator.bit_length()==3186
    assert determinant.denominator%32003 and determinant.numerator*pow(determinant.denominator,-1,32003)%32003==3746
    evidence={
        'old_source_point_count':36,'used_old_point_indices':rows,'candidate_definitions':definitions,
        'no_point_records_added':True,'all_native_input_entries_rebuilt':True,
        'independent_forward_dual_trace_gradients_all16_directions':True,
        'all378_factor_values_exact_QQ_equal':True,'all504_product_values_exact_QQ_equal':True,
        'own21_minor_equals_author_and_separate_Fraction_determinant':str(determinant),
        'numerator_bits':4922,'denominator_bits':3186,'auxiliary_mod32003':3746,
        'Inv14_are_Poisson_brackets_not_trace7':True,
        'complete_Inv22_dimension21_from_prior_independent_character_review':True,
    }
    (output/'independent-certificate.json').write_text(json.dumps(evidence,indent=2)+'\n')
    exe.unlink()
    note=ROOT/'research/notes'/(ID+'.md')
    files=[Path(__file__),cpp,note,author,probe,primary_exact,source36,source36review,old18,old18review,old14,invreview,charreview,invpath,
           ROOT/'code/q_firstjet_spin4_mixed_discriminant_review_20261103_v63.py',
           ROOT/'code/q_new_family_degree22_invariant_exact_20261103_v63.cpp',
           ROOT/'code/q_new_family_degree22_invariant_exclusion_20261103_v63.py',
           ROOT/'research/notes/Q_NEW_FAMILY_DEGREE22_INVARIANT_EXCLUSION_20261103_v63.md']
    report={
        'id':ID,'all_passed':True,
        'old36_QQ_records_only_selected21_used_no_new_points':True,
        'all24_universal_degree22_product_definitions_verified':True,
        'Inv14_basis_is_constant_Poisson_brackets_not_trace7':True,
        'own_forward_dual_all16_dtheta_dH_trace_derivatives':True,
        'all378_factor_and504_product_values_QQ_equal_to_primary':True,
        'QQ_selected21_minor_nonzero_independently_recomputed_twice':True,
        'full_Inv22_dimension21_and_restriction_injective':True,
        'NEW25_firstjet_hypersurface_equation_degree_lower_bound':24,
        'degree24_probe_existence_equation_or_global_exhaustion_claim':False,
        'independent_native_verification_runs':1,
        'author_native_or_probe_rerun':False,
        'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in files},
        'checked_dependency_sha256':checked,
        'evidence_sha256':{p.name:sha(p)for p in output.iterdir()if p.is_file()},
        'peak_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        'child_peak_RSS_KiB':resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss,
        'seconds':time.monotonic()-started,
    }
    (output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    signal.alarm(0);print(json.dumps({'all_passed':True,'degree_lower_bound':24,'seconds':report['seconds']}))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=ROOT/'reports'/ID)
    main(parser.parse_args().output.resolve())
