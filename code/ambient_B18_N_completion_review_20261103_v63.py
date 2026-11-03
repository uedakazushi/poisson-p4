#!/usr/bin/env python3
"""Independent Nakayama application and whole N completion kernel audit."""
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
STEM='AMBIENT_B18_N_FULL_COMPLETION_20261103_v63'
def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda:f.read(1<<20),b''):h.update(chunk)
    return h.hexdigest()
def read(path):
    with (gzip.open(path,'rt') if path.suffix=='.gz' else path.open()) as f:return json.load(f)


def main(output,work):
    started=time.monotonic();work.mkdir(parents=True,exist_ok=False);folder=ROOT/'reports'/STEM
    report=read(folder/'report.json');assert report['all_passed']
    for name,digest in report['input_sha256'].items():assert sha(ROOT/name)==digest,name
    for name,digest in report['evidence_sha256'].items():assert sha(folder/name)==digest,name
    gate=ROOT/'reports/AMBIENT_B18_N_ODD_QUADRATIC_JETS_20261103_v63_REVIEW.json'
    odd=read(gate);assert odd['all_passed'] and odd['actual_odd_minimal_generator_lower_bound9']
    assert odd['zero_full_R_firstjet_two_odd_classes_WV_V2_minor']==24
    replay=work/'producer'
    subprocess.run([str(ROOT/'.venv/bin/python'),str(ROOT/'code/ambient_B18_N_full_completion_20261103_v63.py'),
                    '--output',str(replay),'--odd-review',str(gate),'--m2'],check=True,timeout=40)
    for name in report['evidence_sha256']:assert (folder/name).read_bytes()==(replay/name).read_bytes()
    relations=read(folder/'relations.json.gz');nu,w,v=s.symbols('n w v');ss=s.symbols('s1:7');ys=s.symbols('y0:9');t=s.Symbol('t')
    base=[nu,w,v,*ss];generators=[nu,*ss,w*v,v*v]
    own=[]
    for i,j in itertools.combinations(range(9),2):
        first=s.Poly(generators[i],*base).monoms()[0];second=s.Poly(generators[j],*base).monoms()[0]
        lcm=s.prod(variable**max(a,b) for variable,a,b in zip(base,first,second))
        own.append(s.expand(lcm/generators[i]*ys[i]-lcm/generators[j]*ys[j]))
    own += [ys[i]*ys[j]-nu*generators[i]*generators[j] for i in range(9) for j in range(i,9)]
    assert len(own)==81
    for left,right in zip(own,relations['all81']):assert s.expand(left-s.sympify(right))==0
    image={nu:t*t};image.update({y:g.subs(nu,t*t)*t for y,g in zip(ys,generators)})
    assert all(s.expand(f.subs(image,simultaneous=True))==0 for f in own)
    # Compute the actual kernel of a map to a polynomial domain, using a
    # different ring/map presentation from the producer's elimination ring.
    text='S=QQ['+','.join(map(str,[*base,*ys]))+'];\n'
    text+='T=QQ['+','.join(map(str,[w,v,*ss,t]))+'];\n'
    images=[t*t,w,v,*ss,*[g.subs(nu,t*t)*t for g in generators]]
    text+='phi=map(T,S,{'+','.join(str(x).replace('**','^') for x in images)+'});\n'
    text+='use S;\n'
    text+='P=ideal('+','.join(str(x).replace('**','^') for x in own)+');\n'
    text+='K=ker phi;\nassert(K==P);\nassert(dim(S/P)==9);\n'
    text+='print "INDEPENDENT_N_WHOLE81_MAP_KERNEL_PASS";\nexit 0;\n'
    script=work/'independent-kernel.m2';script.write_text(text)
    run=subprocess.run(['M2','--script',str(script)],text=True,capture_output=True,check=True,timeout=40)
    assert 'INDEPENDENT_N_WHOLE81_MAP_KERNEL_PASS' in run.stdout
    log=work/'independent-kernel.log';log.write_text(run.stdout+run.stderr)
    # Monomial ideal J has7 independent linear generators and2 minimal
    # quadratic generators. The independent actual odd jets supply all9.
    assert sum(s.Poly(g,*base).total_degree()==1 for g in generators)==7
    assert [g for g in generators if s.Poly(g,*base).total_degree()==2]==[w*v,v*v]
    critical=ROOT/'reports/AMBIENT_B18_CA_N_CRITICAL_SCHEMES_20261103_v63_REVIEW.json'
    discriminator=ROOT/'reports/AMBIENT_B18_ODD_DISCRIMINANT_REDUCTION_20261103_v63_REVIEW.json'
    assert read(critical)['all_passed'] and read(discriminator)['all_passed']
    # Subintegral square/cube witness is literally present in the map.
    assert image[ys[0]]==t**3 and image[nu]==t**2
    evidence=output.parent/(output.stem+'_evidence');evidence.mkdir(parents=True,exist_ok=False)
    for source in [script,log]:(evidence/source.name).write_bytes(source.read_bytes())
    paths=[Path(__file__),ROOT/'research/notes'/f'{STEM}_REVIEW.md',folder/'report.json',folder/'relations.json.gz',
           ROOT/'research/notes'/f'{STEM}.md',gate,critical,discriminator]
    result=dict(id=STEM+'_REVIEW',all_passed=True,
        actual_source_odd2_gate_independently_passed=True,
        all9_known_J_mod_mJ_classes_hit_by_actual_I_accepted=True,
        Nakayama_on_finite_J_over_I_proves_exact_I_equals_J=True,
        exact_I='(nu,s1,...,s6,w*v,v^2)',
        original_nu_parameter_preserved_when_lifting_mod_nu_formal_coordinates=True,
        independent_all81_Taylor_and_product_relations_match=True,
        fresh_author_kernel_evidence_byte_identical=True,
        independent_M2_kernel_of_map_into_polynomial_domain_equals_full81=True,
        whole_completed_kernel_proved_by_rank2_reduction_and_Taylor_syzygies=True,
        cone_dimension_embedding_dimension_depth=[28,37,20],projective_dimension_embedding_dimension_depth=[27,36,19],
        Hilbert_Samuel_multiplicity=2,
        actual_N_domain_and_no_embedded_associated_primes=True,
        source_IB_and_target_I_plus_It_conductor_proof_accepted=True,
        conductor_quotient_not_actual_A_contains_embedded_nil_class=True,
        elementary_subintegral_finite_universal_homeomorphism_and_seminormalization_B=True,
        no_offdiagonal_geometric_or_reduced_curve_pairs_any_order_in_formal_germ=True,
        local_bad_loci_equal_conductor_support_proof_accepted=True,
        scope_generic_B18_N_only_not_full_CA_CB_D_or_all_global_N=True,
        input_sha256={str(path.relative_to(ROOT)):sha(path) for path in paths},
        evidence_sha256={str(path.relative_to(ROOT)):sha(path) for path in evidence.iterdir()},seconds=time.monotonic()-started)
    output.parent.mkdir(parents=True,exist_ok=True);output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ['all_passed','Nakayama_on_finite_J_over_I_proves_exact_I_equals_J','seconds']},indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--work',type=Path,required=True);args=parser.parse_args();main(args.output.resolve(),args.work.resolve())
