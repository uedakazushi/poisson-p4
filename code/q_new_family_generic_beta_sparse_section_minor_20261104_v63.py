#!/usr/bin/env python3
"""ONE capped symbolic transversality minor and its geometric field-cover gate."""
import argparse,gzip,json,resource,shlex,subprocess,time
from fractions import Fraction as Q
from pathlib import Path
from q_new_family_v_zero_degree26_support_preflight_20261104_v63 import read,sha
ROOT=Path(__file__).resolve().parents[1]
ID='Q_NEW_FAMILY_GENERIC_BETA_SPARSE_SECTION_MINOR_20261104_v63'
def encodecontent(poly):return {tuple(t[1:]):Q(poly['content'])*int(t[0])for t in poly['terms']}
def main(out,work):
    start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False);work.mkdir(parents=True,exist_ok=False)
    pref=ROOT/'reports/Q_NEW_FAMILY_GENERIC_BETA_SPARSE_SECTION_PREFLIGHT_20261104_v63';gate=read(pref/'report.json');assert gate['all_passed']
    for n,h in gate['input_sha256'].items():assert sha(ROOT/n)==h,n
    data=read(pref/'certificate.json');inp=work/'input.txt'
    with inp.open('w')as f:
        f.write('11 10\n'+' '.join(map(str,data['remaining_rows']))+'\n')
        for row in data['Schur11x10']:
            for p in row:
                f.write(str(len(p))+'\n')
                for t in p:f.write(' '.join(map(str,t))+'\n')
    with gzip.open(out/'native-input.txt.gz','wb')as f:f.write(inp.read_bytes())
    native=ROOT/'code/q_new_family_generic_beta_sparse_section_minor_20261104_v63.cpp';helper=ROOT/'code/q_new_family_v_zero_degree26_restriction_20261104_v63.cpp';binary=work/'symbolic-minor';raw=work/'native-certificate.json'
    flags=shlex.split(subprocess.check_output(['pkg-config','--cflags','--libs','flint'],text=True));subprocess.run(['c++','-O2','-std=c++17',str(native),*flags,'-lgmp','-o',str(binary)],check=True,timeout=max(1,60-(time.monotonic()-start)))
    result=None;error=None;found=None
    try:
        with(out/'native.log').open('w')as f:subprocess.run([str(binary),str(inp),str(raw)],stdout=f,stderr=subprocess.STDOUT,check=True,timeout=max(1,60-(time.monotonic()-start)))
        assert raw.stat().st_size<32*1024*1024;result=read(raw);assert result['all_passed']and result['all_Bareiss_divisions_exact_and_remultiplied']
        with gzip.open(out/'native-certificate.json.gz','wb')as f:f.write(raw.read_bytes())
        if result['nonzero_deleted_row']>=0:
            found=result['nonzero_deleted_row'];attempt=result['attempts'][-1];det=encodecontent(attempt['determinant']);assert det
            # Row and column permutation signs for the actual sorted15-row minor.
            sortedrows=[i for i in range(16)if i!=found];newrows=data['source_pivot_rows']+[i for i in data['remaining_rows']if i!=found];permutation=[sortedrows.index(i)for i in newrows];sign=(-1)**sum(permutation[i]>permutation[j]for i in range(15)for j in range(i+1,15))
            full15={tuple([e[0]+1,*e[1:]]):c*sign*45000 for e,c in det.items()};full16={e:c*((-1)**(found+1))for e,c in full15.items()}
            rows=lambda p:[[str(c),*e]for e,c in sorted(p.items())]
            section=dict(variables=['lambda','H','L','s0','q0','epsilon'],deleted_row=found,
                sparse_beta6=[z+(f'+epsilon'if i==found else'')for i,z in enumerate(data['new_beta16'])],
                coordinate_convention='original beta16 basis, zero-based',localizers=['lambda','H','L'],
                Schur_determinant_terms=rows(det),row_permutation_sign=sign,source_block_determinant='45000*H',
                full15_minor=dict(lambda_power=14,HLsq_polynomial_terms=rows(full15)),full16_Jacobian_determinant=dict(lambda_power=14,HLsq_polynomial_terms=rows(full16)),
                determinant_epsilon_specialization=0,nonzero_coefficient_witness=rows(full16)[0],
                actual_pulledback_Lie10_and_source6=True,dominant_generically_finite_to_W16=True,
                inference='nonzero16Jac at epsilon0 proves dominance inchar0; equal16dimensions gives finite functionfield extension; extension is faithfully flat, permitting generic-beta exhaustion after basechange',
                no_birationality_or_generic_exhaustion_solution_claim=True)
            assert found==4
            gamma={(0,0,0,0):Q(1),(0,0,0,1):Q(5),(0,0,1,0):Q(8),(2,0,0,0):Q(15),(2,0,0,1):Q(115)}
            assert det=={(e[0]+2,e[1]+2,e[2],e[3]):c*Q(2**30,5605041796875)for e,c in gamma.items()}
            assert full16=={(e[0]+3,e[1]+2,e[2],e[3]):c*Q(2**33,996451875)for e,c in gamma.items()}
            section['Gamma']='1+5*q0+8*s0+15*H^2+115*H^2*q0'
            section['Schur_determinant_factored']='2^30/5605041796875*H^2*L^2*Gamma'
            section['actual16_Jacobian_factored']='2^33/996451875*lambda^14*H^3*L^2*Gamma'
            section['all_factorization_coefficients_remultiplied']=True
            (out/'section6.json').write_text(json.dumps(section,indent=2)+'\n')
        assert sum(p.stat().st_size for p in out.iterdir()if p.is_file())<32*1024*1024
    except Exception as e:error=repr(e)
    paths=[Path(__file__),native,helper,ROOT/f'research/notes/{ID}.md',pref/'report.json',pref/'certificate.json']
    report=dict(id=ID,all_passed=error is None,error=error,nonzero_deleted_row=found,
        actual_symbolic_transversality_minor_nonzero=found is not None,dominant_generically_finite_six_parameter_section=found is not None,
        one_symbolic_job_max11_minors=True,no_point_prime_or_modular_rank=True,
        no_birationality_or_J0_substitution_or_generic_exhaustion_claim=True,
        input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},evidence_sha256={p.name:sha(p)for p in out.iterdir()if p.is_file()},seconds=time.monotonic()-start)
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(dict(all_passed=report['all_passed'],deleted_row=found,seconds=report['seconds'],error=error)))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--work',type=Path,required=True);a=p.parse_args();main(a.output,a.work)
