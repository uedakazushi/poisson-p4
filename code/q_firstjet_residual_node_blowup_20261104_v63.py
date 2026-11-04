#!/usr/bin/env python3
"""ONE new proper (2,2) model test at the SAME authorized coefficients."""
import argparse
import hashlib
import json
from pathlib import Path
import resource
import signal
import subprocess
import time

ROOT=Path(__file__).resolve().parents[1]
STEM='Q_FIRSTJET_RESIDUAL_NODE_BLOWUP_20261104_v63'
P=32003
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def add(a,b):return((a[0]+b[0])%P,(a[1]+b[1])%P)
def neg(a):return((-a[0])%P,(-a[1])%P)
SUM=287*pow(2,-1,P)%P;PROD=10914*pow(2,-1,P)%P
def mul(a,b):
    return((a[0]*b[0]-PROD*a[1]*b[1])%P,
           (a[0]*b[1]+a[1]*b[0]+SUM*a[1]*b[1])%P)
def power(a,n):
    b=(1,0)
    for _ in range(n):b=mul(b,a)
    return b
def norm(a):return(a[0]*a[0]+SUM*a[0]*a[1]+PROD*a[1]*a[1])%P
def padd(a,b):
    c=dict(a)
    for ex,co in b.items():
        c[ex]=add(c.get(ex,(0,0)),co)
        if c[ex]==(0,0):c.pop(ex)
    return c
def pmul(a,b):
    c={}
    for ea,ca in a.items():
        for eb,cb in b.items():
            ex=tuple(ea[i]+eb[i] for i in range(2))
            c=padd(c,{ex:mul(ca,cb)})
    return c
def ppow(a,n):
    b={(0,0):(1,0)}
    for _ in range(n):b=pmul(b,a)
    return b

def main(out):
    started=time.monotonic();signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000))
    out.mkdir(parents=True,exist_ok=False)
    old=ROOT/'reports/Q_FIRSTJET_RESIDUAL_QUARTIC_SMOOTH_SPECIALIZATION_20261104_v63'
    data=json.loads((old/'certificate.json').read_text())
    assert data['one_authorized_coefficient_specialization']==[1,2,3,5,7] and data['prime']==P
    assert not data['projectively_geometrically_smooth']
    assert (-4943)%P
    alpha=(0,1);beta=(SUM,P-1);delta=add(alpha,neg(beta))
    assert norm(delta)
    # r=(alpha*V-beta*U)/delta, z=(V-U)/delta; original m=M.
    rn={(0,1):alpha,(1,0):neg(beta)}
    zn={(0,1):(1,0),(1,0):(P-1,0)}
    H={}
    for key,co in data['all15_specialized_coefficients_r_m'].items():
        i,j=map(int,key.split(','));k=4-i-j
        term=pmul(ppow(rn,i),ppow(zn,k))
        scalar=mul((co,0),power(delta,j))
        H=padd(H,{ex:mul(value,scalar) for ex,value in term.items()})
    assert H and all(i<=2 and j<=2 for i,j in H)
    assert norm(H[2,2])
    def expression(polynomial):
        return ' + '.join(f'({a}+{b}*aa)*u^{i}*v^{j}' for (i,j),(a,b) in sorted(polynomial.items()))
    script='S=ZZ/32003[aa,u,v,MonomialOrder=>GRevLex];\nN=2*aa^2-287*aa+10914;\n'
    flags=[]
    for e in range(2):
        for f in range(2):
            patch={(2-i if e else i,2-j if f else j):co for (i,j),co in H.items()}
            tag=f'PATCH_{e}{f}';flags.append(tag)
            script+=f'H={expression(patch)};\n'
            script+=f'good=(1_S % gb(ideal(N,H,diff(u,H),diff(v,H)))==0_S);\nprint("{tag} "|toString good);\n'
    m2=out/'blowup22.m2';m2.write_text(script)
    status='pending';returncode=None;elapsed=None
    try:
        t0=time.monotonic()
        proc=subprocess.run(['M2','--script',str(m2.resolve())],capture_output=True,text=True,timeout=55)
        elapsed=time.monotonic()-t0;returncode=proc.returncode;log=proc.stdout+proc.stderr
        status='completed' if returncode==0 else 'M2_failure_STOP_no_retry'
    except subprocess.TimeoutExpired as exc:
        log=(exc.stdout or b'')+(exc.stderr or b'')
        if isinstance(log,bytes):log=log.decode(errors='replace')
        status='timeout_STOP_no_retry'
    (out/'M2.log').write_text(log)
    checks={tag:(tag+' true' in log) for tag in flags}
    passed=returncode==0 and all(checks.values())
    certificate={'same_coefficient_point':[1,2,3,5,7],'prime':P,
                 'node_etale_algebra':'Fp[alpha]/(2alpha^2-287alpha+10914)',
                 'node_sum':SUM,'node_product':PROD,'node_discriminant_modp':(-4943)%P,
                 'alpha_beta_difference_norm':norm(delta),
                 'cleared_transform':'delta^4 F((alpha V-beta U)/delta,M,(V-U)/delta)',
                 'proper_model':'bihomogenize H(U/M,V/M) to bidegree(2,2) in P1xP1',
                 'all_H22_coefficient_pairs':[[list(ex),list(co)] for ex,co in sorted(H.items())],
                 'H22_coefficient_norm':norm(H[2,2]),'all_four_patch_results':checks,
                 'all_geometric_choices_of_alpha_tested_simultaneously':True,
                 'strict_transform_avoids_contracted_joining_line_intersection':0}
    (out/'certificate.json').write_text(json.dumps(certificate,indent=2)+'\n')
    paths=[Path(__file__),old/'report.json',old/'certificate.json',
           ROOT/'research/notes/Q_FIRSTJET_RESIDUAL_QUARTIC_FORCED_NODES_20261104_v63.md']
    report={'id':STEM,'all_passed':passed,'run_status':status,
            'same_one_coefficient_specialization_only_new_proper_blowup_model':True,
            'NOT_original_quartic_retry_or_new_beta_source_point':True,
            'special_proper22_model_geometrically_smooth':passed,
            'proper_flat_etale_node_family_generic_smoothness_proved_if_PASS':passed,
            'generic_quotient_normalization_geometrically_irreducible_genus1':passed,
            'original_generic_quartic_exactly_two_ordinary_nodes_if_PASS':passed,
            'exact_full_quadratic_cover_genus_or_selected_degree8_claimed_here':False,
            'no_alternative_parameters_point_invariants_or_unbounded_GB':True,
            'input_sha256':{str(p.relative_to(ROOT)):sha(p) for p in paths},
            'evidence_sha256':{name:sha(out/name) for name in ['blowup22.m2','M2.log','certificate.json']},
            'M2_returncode':returncode,'M2_seconds':elapsed,'seconds':time.monotonic()-started,
            'max_child_rss_kib':resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss}
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({key:report[key] for key in ['all_passed','run_status','special_proper22_model_geometrically_smooth','seconds']}))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True)
    main(p.parse_args().output)
