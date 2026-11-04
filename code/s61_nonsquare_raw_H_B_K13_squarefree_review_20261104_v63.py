#!/usr/bin/env python3
"""Independent literal integer convolution, no gcd or factor search."""
import argparse,hashlib,json,math,resource,signal,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
STEM='S61_NONSQUARE_RAW_H_B_K13_SQUAREFREE_20261104_v63_REVIEW'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def dense(p):
    x=[0]*(max(int(i) for i,c in p)+1)
    for i,c in p:x[int(i)]=int(c)
    return x
def convolution(a,b):
    x=[0]*(len(a)+len(b)-1)
    for i,c in enumerate(a):
        for j,d in enumerate(b):x[i+j]+=c*d
    return x
def total(a,b):
    x=[0]*max(len(a),len(b))
    for i,c in enumerate(a):x[i]+=c
    for i,c in enumerate(b):x[i]+=c
    while len(x)>1 and not x[-1]:x.pop()
    return x
def main(out):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));out.mkdir(parents=True,exist_ok=False)
    b=ROOT/'reports/S61_NONSQUARE_RAW_H_B_K13_SQUAREFREE_20261104_v63'
    p=b/'report.json';r=read(p);assert sha(p)=='5944045153ac0e15f67c990565971702a6577fed71e9750499f57c8c23a59526' and r['all_passed']
    for q,h in r['input_sha256'].items():assert sha(ROOT/q)==h,q
    for q,h in r['evidence_sha256'].items():assert sha(b/q)==h,q
    c=read(b/'certificate.json');native=read(b/'native-certificate.json')
    source=ROOT/'reports/S61_NONSQUARE_RAW_H_B_K13_D_EUCLID_20261104_v63/checkpoint.json'
    actual=ROOT/'reports/S61_NONSQUARE_RAW_H_B_K13_ACTUAL_SOURCE_SCHEME_20261104_v63_REVIEW/report.json'
    assert read(actual)['all_passed']
    g=dense(c['g13']);assert g==dense(read(source)['g13']) and len(g)==14 and g[-1]
    dg=[(i+1)*g[i+1] for i in range(13)];assert dg==dense(c['derivative'])
    content=math.gcd(*map(abs,dg));assert content==int(c['derivative_content']) and content>0
    dp=[a//content for a in dg];assert dp==dense(c['primitive_derivative'])
    S,T=dense(c['derivative_cofactor']),dense(c['g13_cofactor']);R=int(c['nonzero_QQ_integer_constant']);assert R
    result=total(convolution(S,dg),convolution(T,g));assert result==[R]
    primitive=total(convolution(dense(native['S']),dp),convolution(dense(native['T']),g));assert primitive==[int(native['R'])] and primitive[0]
    assert S==dense(native['S']) and T==[content*a for a in dense(native['T'])] and R==content*int(native['R'])
    proof={'actual_source_g13_degree':13,'derivative_degree':12,'literal_full_integer_Bezout_remainder':[str(R)],
      'all_other_coefficients_exactly_zero':True,'native_primitive_and_original_derivative_identities_both_replayed':True,
      'own_engine':'dense Python integer convolution, no producer import or gcd/factor search',
      'cofactor_nonzero_term_counts':[sum(bool(a) for a in S),sum(bool(a) for a in T)],
      'cofactor_max_bits':max(abs(a).bit_length() for a in S+T),'constant_bits':abs(R).bit_length(),
      'derivative_content':str(content),'finite_etale_degree13_QQ_proved':True,
      'actual_Bbad_geometrically_reduced13points_by_independent_scheme_gate':True,
      'QQ_rational_points_or_QQ_irreducibility_not_claimed':True,
      'no_full_cubic_rankbad_or_Poisson_nil_extension':True}
    ev=out/'independent-integer-identity.json';ev.write_text(json.dumps(proof,indent=2)+'\n')
    paths=[Path(__file__),ROOT/'research/notes'/f'{STEM}.md',p,b/'certificate.json',b/'native-certificate.json',source,actual,
      ROOT/'code/s61_nonsquare_raw_H_B_K13_squarefree_20261104_v63.py',ROOT/'research/notes/S61_NONSQUARE_RAW_H_B_K13_SQUAREFREE_20261104_v63.md']
    rep={'id':STEM,'all_passed':True,'g13_squarefree_and_actual_geometric_reduced13_accepted':True,
      'literal_integer_cofactor_identity_independently_checked':True,
      'no_new_gcd_factor_source_point_prime_GB_or_producer_run':True,
      'input_sha256':{str(q.relative_to(ROOT)):sha(q) for q in paths},'evidence_sha256':{ev.name:sha(ev)},
      'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    (out/'report.json').write_text(json.dumps(rep,indent=2)+'\n');print(json.dumps({k:rep[k] for k in ['id','all_passed','seconds','max_rss_kib']}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
