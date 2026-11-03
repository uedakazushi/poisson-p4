#!/usr/bin/env python3
"""Finite structural bounds and costs; no degree4 expansion or GB is launched."""
import argparse
import hashlib
import json
from math import comb
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--ranks',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=False)
    rankpath=args.ranks/'report.json';r=json.loads(rankpath.read_text());assert r['all_passed']
    n=23;m=57;c=17;linear_syzygies=n*m-r['linear_quadratic_multiple_rank'];assert linear_syzygies==101
    # A101-dimensional subspace of S_1^57 has101 distinct leading module
    # monomials.  If componenti containsr_i such axes, its linear multiples
    # have at leastg(r_i)=23*r_i-r_i*(r_i-1)/2 distinct degree2 leading terms.
    g=lambda t:n*t-t*(t-1)//2
    infinity=10**9;dp=[infinity]*(linear_syzygies+1);dp[0]=0
    for _ in range(m):
        nxt=[infinity]*(linear_syzygies+1)
        for s,value in enumerate(dp):
            if value==infinity:continue
            for t in range(min(n,linear_syzygies-s)+1):nxt[s+t]=min(nxt[s+t],value+g(t))
        dp=nxt
    shadow=dp[linear_syzygies];full,remainder=divmod(linear_syzygies,n)
    assert shadow==full*g(n)+g(remainder)==1275
    dimensions={degree:comb(n+degree-1,degree) for degree in range(6)}
    rows4=m*dimensions[2]+c*dimensions[1];rows5=m*dimensions[3]+c*dimensions[2]
    h4=dimensions[4]-rows4+shadow;assert h4==102
    rows={4:rows4,5:rows5};memory={}
    for degree in [4,5]:
        nnz=m*dimensions[degree-2]*dimensions[2]+c*dimensions[degree-3]*dimensions[3]
        memory[degree]=dict(rows=rows[degree],columns=dimensions[degree],
                            dense_limb8_bytes=rows[degree]*dimensions[degree]*8,
                            sparse_nnz_upper_bound=nnz,
                            sparse_col32_coeff32_plus_rowptr64_bytes_upper_bound=8*nnz+8*(rows[degree]+1))
    report=dict(id='Q_THROUGH3_SECTION_MACAULAY_ASSESSMENT_20261103_V63',all_passed=True,
                coefficient_field='GF(32003)',section_variables=n,quadrics=m,cubics=c,
                exact_linear_quadratic_syzygy_dimension=linear_syzygies,
                independent_leading_module_axes_count=linear_syzygies,
                minimum_linear_shadow_dimension_by_exhaustive_allocation_DP=shadow,
                minimizing_distribution_nonzero=[23,23,23,23,9],
                degree4_Macaulay_rows=rows4,degree4_polynomial_space_dimension=dimensions[4],
                degree4_ideal_rank_upper_bound=rows4-shadow,
                section_Hilbert_degree4_lower_bound=h4,
                Artinian_length_if_finite_at_least=sum(r['exact_homogeneous_section_Hilbert_values_mod_p_through3'])+h4,
                degree4_FULL_Macaulay_rank_Artinian_certificate_ruled_out=True,
                Artinian_GB_certificate_with_generators_of_degree_at_most4_NOT_ruled_out=True,
                dense_and_sparse_input_storage_estimates=memory,
                sparse_input_fitting_memory_does_not_certify_feasible_elimination_time_or_fill=True,
                actual_rank_degree4_or_higher_not_computed=True,
                no_new_GB_or_higher_DAG_truncation_performed=True,
                QQ_Hilbert_degree4_or_Krull_dimension_not_inferred=True,
                input_sha256={str(path.resolve().relative_to(ROOT)):sha(path) for path in [Path(__file__),rankpath]})
    (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print('PASS structuralGF h4>=',h4,'minimumshadow',shadow,'denseD5bytes',memory[5]['dense_limb8_bytes'])


if __name__=='__main__':main()
