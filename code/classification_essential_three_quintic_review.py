#!/usr/bin/env python3
"""Finite arithmetic receipt for the separate essential-three quintic proof.

This checks Hilbert/genus/Milnor counts and the complete finite binomial type
superset. BV, determinant, module, analytic-stratum and closure arguments are
mathematical proofs in the separate note, not claims made by a plan checker.
"""
from pathlib import Path
from math import comb,gcd
import hashlib,json,time,resource,signal
ROOT=Path(__file__).resolve().parents[1]
JOB='reports/CLASSIFICATION_ESSENTIAL_THREE_QUINTIC_REVIEW_JOB_20261107.json'
OUT='reports/CLASSIFICATION_ESSENTIAL_THREE_QUINTIC_REVIEW_20261107.json'

def main():
    started=time.monotonic();signal.alarm(30)
    resource.setrlimit(resource.RLIMIT_AS,(256*1024**2,)*2)
    path=ROOT/OUT
    if path.exists():raise FileExistsError(path)
    job=json.loads((ROOT/JOB).read_text())
    for b in job['exact_inputs_and_hashes']:assert hashlib.sha256((ROOT/b['path']).read_bytes()).hexdigest()==b['sha256'],b['path']
    dimension=5;coefficient_degree=2;jacobi_degree=2*coefficient_degree-1
    assert jacobi_degree+dimension-3==5
    assert 2*coefficient_degree+dimension-4==5
    numerator={0:1,4:-3,6:2}
    derivatives=[sum(c for d,c in numerator.items()),sum(d*c for d,c in numerator.items()),sum(d*(d-1)*c for d,c in numerator.items())]
    assert derivatives==[0,0,24];tau=derivatives[2]//2;assert tau==12
    hS=lambda n:comb(n+2,2) if n>=0 else 0
    Hilbert=[sum(c*hS(n-d) for d,c in numerator.items()) for n in range(13)]
    assert all(v==12 for v in Hilbert[4:])
    hI5=hS(5)-Hilbert[5];assert hI5==9
    pa=(5-1)*(5-2)//2;assert pa==6
    budget=[(g,B,2*pa-2*g-B) for g in range(pa+1) for B in range(2*pa+1) if 2*pa-2*g-B>=tau]
    assert budget==[(0,0,12)]
    local=[]
    for a in range(2,14):
        for b in range(a,14):
            mu=(a-1)*(b-1)
            if gcd(a,b)==1 and mu<=12:
                assert mu%2==0
                local.append(dict(a=a,b=b,mu=mu,delta=mu//2))
    assert [(z['a'],z['b']) for z in local]==[(2,3),(2,5),(2,7),(2,9),(2,11),(2,13),(3,4),(3,5),(3,7),(4,5)]
    configurations=[]
    def enumerate_types(first,remaining,indices):
        if not remaining:configurations.append(indices);return
        for i in range(first,len(local)):
            mu=local[i]['mu']
            if mu<=remaining:enumerate_types(i,remaining-mu,indices+[i])
    enumerate_types(0,12,[])
    assert configurations and all(sum(local[i]['mu'] for i in row)==12 and len(row)<=6 for row in configurations)
    tangent=hI5-1;fixedF=4+comb(4,2)-1;grassmann=3*(5-3)
    assert tangent==8 and fixedF==9 and grassmann+tangent+fixedF==23
    m=4;ordinary_delta=m*(m-1)//2;ordinary_mu=(m-1)**2
    assert ordinary_delta==pa==6 and ordinary_mu==9<tau
    out=dict(task_id='G',all_passed=True,
        theorem='All rank4 actual irreducible essential-three quintic divisors have projective source image dimension<=23; any rank-one squarefree m4 leading binaryquartic is impossible for rank4 with actual irreducible quintic divisor.',
        mathematical_proof_audit='Canonical divergence-zero strict lift/BV-Cartan, primitive-gradient base bracket vanishing, full algebraic determinant criterion, AR basis and all fiber dependence, full cotangent6mag source, saturated Jacobian, finite contact-type strata and moving-point tangent bounds independently checked in the separate note.',
        numerical_certificate=dict(canonical_weight_checks=[5,5],Jacobian_resolution_degrees=[[-6,-6],[-4,-4,-4]],Hilbert_numerator=numerator,Hilbert_numerator_derivatives_at1=derivatives,Milnor_algebra_Hilbert_function_0to12=Hilbert,total_Tjurina=tau,h0_I_Sigma5=hI5,projective_equianalytic_tangent_upper=tangent,
                                  all_genus_branch_budget_solutions=budget,all_binomial_local_type_superset=local,
                                  all_unordered_type_configurations_sum_mu12=configurations,
                                  global_realizability_of_list_not_claimed=True,
                                  fixed_projective_F_source_upper=fixedF,essential_variable_Grassmann_dimension=grassmann,whole_projective_source_upper=23,
                                  ordinary_quadruple_point=dict(delta=ordinary_delta,branches=4,mu=ordinary_mu,total_tau_upper=ordinary_mu,required_tau=tau,contradiction=True)),
        external_theorem_scope='Saito1971 isolated germ f in J=>analytically quasihomogeneous, read as the exact indexed statement of Maynadier translation Introduction; original German/full proof not retrieved. Determinant criterion is reproved algebraically in review note.',
        complements=['special-j rank-two leading pencils','generic repeated leading binaryroots','rank-two plane m2/m3','other nonnormal quintic supports','whole E membership'],
        no_CAS=True,inputs=job['exact_inputs_and_hashes'],reviewer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),seconds=time.monotonic()-started,maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    path.write_text(json.dumps(out,indent=2)+'\n')
    print('PASS',OUT,'localtypes',len(local),'typeconfigurations',len(configurations),'wholeG23',flush=True)
if __name__=='__main__':main()
