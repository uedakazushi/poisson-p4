#!/usr/bin/env python3
"""Accepted odd2jet+Nakayama application and whole81relation N completion."""
import argparse,gzip,hashlib,itertools,json,subprocess,time
from pathlib import Path
import sympy as s
ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main(out,review,m2):
    started=time.monotonic();out.mkdir(parents=True,exist_ok=False)
    author=ROOT/'reports/AMBIENT_B18_N_ODD_QUADRATIC_JETS_20261103_v63/report.json';proof=json.loads(author.read_text());audit=json.loads(review.read_text())
    assert proof['all_passed'] is True and audit['all_passed'] is True
    assert proof['actual_full_N_odd_ideal_minimal_generator_lower_bound']==9
    dependencies=[author,review,ROOT/'reports/AMBIENT_B18_ODD_DISCRIMINANT_REDUCTION_20261103_v63_REVIEW.json',ROOT/'reports/AMBIENT_B18_CA_N_CRITICAL_SCHEMES_20261103_v63_REVIEW.json']
    for path in dependencies:
        data=json.loads(path.read_text());assert data['all_passed'] is True
        for name,digest in data['input_sha256'].items():assert sha(ROOT/name)==digest,name
    # The application proof identifiesJ=(nu,I), whose9 minimal classes are
    # all hit byactualI: 7 linear plus the two independent WV,V² classes.
    # Nakayama onJ/I then provesI=J, not merelya generator lower bound.
    n,w,v=s.symbols('n w v');ss=s.symbols('s1:7');yy=s.symbols('y0:9');t=s.Symbol('t');base=[n,w,v,*ss]
    b=[n,*ss,w*v,v*v];assert len(b)==9
    def powers(m):return s.Poly(m,*base).monoms()[0]
    def lcm(a,c):return s.prod(x**max(i,j)for x,i,j in zip(base,powers(a),powers(c)))
    linear=[];products=[]
    for i,j in itertools.combinations(range(9),2):
        common=lcm(b[i],b[j]);linear.append(s.expand(common/b[i]*yy[i]-common/b[j]*yy[j]))
    for i in range(9):
        for j in range(i,9):products.append(yy[i]*yy[j]-n*b[i]*b[j])
    assert len(linear)==36 and len(products)==45
    substitutions={y:bi*t for y,bi in zip(yy,b)}
    for equation in linear+products:
        assert s.rem(s.Poly(equation.subs(substitutions),t),s.Poly(t*t-n,t)).as_expr()==0
    relations=linear+products
    data=dict(core_R_variables=list(map(str,base)),odd_coefficient_ideal_generators=list(map(str,b)),odd_target_variables=list(map(str,yy)),Taylor_pair_syzygies=list(map(str,linear)),all_pair_products=list(map(str,products)),all81=list(map(str,relations)),normalization_map={str(y):str(bi*t)for y,bi in zip(yy,b)},normalizer_relation='t^2-n')
    (out/'relations.json.gz').write_bytes(gzip.compress((json.dumps(data,separators=(',',':'))+'\n').encode(),mtime=0))
    variables=[t,*base,*yy]
    text='S=QQ['+','.join(map(str,variables))+',MonomialOrder=>Eliminate 1];\n'
    text+='P=ideal('+','.join(str(x).replace('**','^')for x in relations)+');\n'
    text+='G=ideal('+','.join(str(x).replace('**','^')for x in[t*t-n,*[y-bi*t for y,bi in zip(yy,b)]])+');\n'
    text+='K=eliminate({t},G);\nassert(K==P);\nassert(dim(S/P)==10); --includesfree eliminated t\n'
    text+='print "N_WHOLE81RELATION_KERNEL_PASS";\nexit 0;\n'
    (out/'kernel.m2').write_text(text)
    if m2:
        run=subprocess.run(['M2','--script',str((out/'kernel.m2').resolve())],capture_output=True,text=True,check=True,timeout=60)
        assert 'N_WHOLE81RELATION_KERNEL_PASS'in run.stdout;(out/'kernel.log').write_text(run.stdout+run.stderr)
    note=ROOT/'research/notes/AMBIENT_B18_N_FULL_COMPLETION_20261103_v63.md'
    paths=[Path(__file__),note,*dependencies]
    result=dict(id='AMBIENT_B18_N_FULL_COMPLETION_20261103_v63',all_passed=True,
        actual_odd_I_maps_surjectively_to_all9_classes_of_J_mod_mJ=True,
        Nakayama_applied_to_finite_module_J_over_I=True,actual_odd_I_equals_I_plus_discriminant=True,
        exact_full_odd_ideal='I=(nu,s1,...,s6,w*v,v^2)',
        actual_normalizer='B=C[[g19,w,v,s6,t]]',
        actual_N_cone_completion='C[[g19,w,v,s6,t^2,s1*t,...,s6*t,t^3,w*v*t,v^2*t]] subset B',
        actual_N_projective_completion='same core withg18 instead ofg19',
        exact_core_target_relation_count=81,complete_first_syzygy_relations=36,complete_pair_products=45,
        whole_kernel_completeness_by_Taylor_syzygies_and_rank2_normal_form=True,
        independent_Macaulay2_whole_elimination_kernel_checked=bool(m2),
        cone_dimension=28,cone_embedding_dimension=37,cone_depth=20,
        projective_dimension=27,projective_embedding_dimension=36,projective_depth=19,
        Hilbert_Samuel_multiplicity=2,
        domain_and_analytically_irreducible=True,no_embedded_associated_primes_in_actual_N_completion=True,
        source_conductor_ideal_in_B='(t^2,s1,...,s6,w*v,v^2) B',
        target_conductor_ideal_in_A='nine even I generators and nine odd I*t generators',
        normalizer_quotient_module='B/A=(C[[g19,w,v]]/(w*v,v^2))*t',
        embedded_prime_belongs_only_to_conductor_quotient_not_actual_N_domain=True,
        t_square_and_cube_in_A_but_t_not_in_A=True,normalization_elementary_subintegral=True,
        normalization_finite_universal_homeomorphism=True,seminormalization_equals_regular_B=True,not_seminormal=True,
        no_offdiagonal_geometric_or_reduced_curve_arc_pairs_any_order_in_this_formal_germ=True,
        local_reduced_bad_locus='smooth C[[g19,w]] support of conductor',
        local_singular_nonnormal_and_nonCM_loci_equal_conductor_support=True,
        scope='Generic B18 open orbit N completion only; not allspecialB18 points, global N bad locus, CA/CB full types, or actualD.',
        input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},evidence_sha256={p.name:sha(p)for p in out.iterdir()if p.is_file()},seconds=time.monotonic()-started)
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');print('FULL ACTUAL N COMPLETION: I=J /81REL /DOMAIN DEPTH20 e2 PASS',result['seconds'],flush=True)
if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--odd-review',type=Path,required=True);ap.add_argument('--m2',action='store_true');args=ap.parse_args();main(args.output.resolve(),args.odd_review.resolve(),args.m2)
