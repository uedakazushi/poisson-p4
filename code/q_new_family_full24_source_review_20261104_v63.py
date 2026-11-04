#!/usr/bin/env python3
"""One independent full67 source coefficient minor; no source point evaluations."""
import argparse,gzip,hashlib,itertools,json,resource,shlex,signal,subprocess,time
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
STEM='Q_NEW_FAMILY_FULL_DEGREE24_SOURCE_RELATION_20261104_v63_REVIEW'
PRIME=32003
def read(p):
    with(gzip.open(p,'rt')if p.suffix=='.gz'else p.open())as f:return json.load(f)
def sha(p):
    h=hashlib.sha256()
    with p.open('rb')as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()
def detmod(A):
    A=[r.copy()for r in A];det=1
    for j in range(len(A)):
        i=next(i for i in range(j,len(A))if A[i][j])
        if i!=j:A[i],A[j]=A[j],A[i];det=-det%PRIME
        det=det*A[j][j]%PRIME;iv=pow(A[j][j],-1,PRIME)
        for i in range(j+1,len(A)):
            c=A[i][j]*iv%PRIME
            if c:A[i]=[(x-c*y)%PRIME for x,y in zip(A[i],A[j])]
    return det
def main(out,work):
    started=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
    out.mkdir(parents=True,exist_ok=False);work.mkdir(parents=True,exist_ok=False)
    folder=ROOT/'reports/Q_NEW_FAMILY_FULL_DEGREE24_SOURCE_RELATION_20261104_v63'
    ar,generators,cert=[read(folder/p)for p in ['report.json','generators.json','native-certificate.json.gz']]
    assert ar['all_passed']and ar['status']=='all_Inv24_QQ_restrictions_injective'
    assert cert['coefficient_rank_modp']==67 and cert['selected_columns']==list(range(67))
    gates=[ROOT/'reports/Q_NEW_FAMILY_FULL_DEGREE24_HESSIAN_BASIS_20261104_v63_REVIEW/report.json',
        ROOT/'reports/Q_NEW_FAMILY_HESSIAN_LOW_INVARIANT_EXPANSION_20261104_v63_REVIEW/report.json',
        ROOT/'reports/Q_NEW_FAMILY_LOW_INVARIANT_EXPANSION_20261104_v63_REVIEW/report.json',
        ROOT/'reports/Q_NEW_FAMILY_FULL_DEGREE24_LOGIC_20261104_v63_REVIEW/report.json',
        ROOT/'reports/Q_NEW_FAMILY_DEGREE22_INVARIANT_EXCLUSION_20261103_v63_REVIEW/report.json']
    assert all(read(p)['all_passed']for p in gates)
    names=['F4','p1','p2','p3','p4','N12','CD000101','CD000111','theta5H','theta7H']
    assert generators['generator_order']==names and generators['weights_beta']==[4,8,8,8,8,12,12,12,12,16]
    paths=[ROOT/'reports'/('Q_NEW_FAMILY_LOW_INVARIANT_EXPANSION_20261104_v63'if i<6 else 'Q_NEW_FAMILY_HESSIAN_LOW_INVARIANT_EXPANSION_20261104_v63')/(n+'.json.gz')for i,n in enumerate(names)]
    polys=list(map(read,paths))
    # Independently translate the complete universal70 definitions into the
    # actual10 source generators, rather than trusting the author's exponent file.
    all70=[]
    for r in range(4):
        for indices in itertools.combinations_with_replacement(range(1,5),r):
            e=[0]*10;e[0]=6-2*r
            for i in indices:e[i]+=1
            all70.append(e)
    for i in range(4):
        for j in range(5):
            e=[0]*10;e[5+i]=1;e[0]=3 if j==0 else 1
            if j:e[j]=1
            all70.append(e)
    for i,j in itertools.combinations_with_replacement(range(4),2):
        e=[0]*10;e[5+i]+=1;e[5+j]+=1;all70.append(e)
    for j in range(5):
        e=[0]*10;e[9]=1;e[0 if j==0 else j]=2 if j==0 else 1;all70.append(e)
    selected=[i for i in range(70)if i not in [62,68,69]]
    assert len(all70)==70 and all70==generators['all70_generator_exponents']
    assert selected==generators['selected_columns_from70']
    exponents=[all70[i]for i in selected];assert exponents==generators['selected67_generator_exponents']
    for j,e in enumerate(exponents):
        assert sum(a*b for a,b in zip(e,generators['weights_beta']))==24
        assert sum(a*p['denominator_D_power']for a,p in zip(e,polys))==24
        assert sum(a*p['denominator_B_power']for a,p in zip(e,polys))<=12
        d=1
        for a,p in zip(e,polys):d*=int(p['denominator_integer'])**a
        assert d==int(generators['column_constant_denominators'][j])and d%PRIME
    inp=out/'own-input.txt'
    with inp.open('w')as f:
        f.write('32003 4 10 67 67\n')
        for p in polys:
            terms=[[int(c)%PRIME,b,w,y,a]for c,a,y,w,b in p['terms']if int(c)%PRIME]
            f.write(f"{p['denominator_B_power']} {len(terms)}\n")
            for row in terms:f.write(' '.join(map(str,row))+'\n')
        for e in exponents:f.write(' '.join(map(str,e))+'\n')
        for e in cert['selected_rows']:f.write(' '.join(map(str,reversed(e)))+'\n')
    native=ROOT/'code/q_new_family_full24_source_review_20261104_v63.cpp';binary=work/'full67-coefficient-review'
    flags=shlex.split(subprocess.check_output(['pkg-config','--cflags','--libs','flint'],text=True))
    with(out/'compile.log').open('w')as f:subprocess.run(['c++','-O2','-std=c++17',str(native),*flags,'-lgmp','-o',str(binary)],check=True,stdout=f,stderr=subprocess.STDOUT,timeout=20)
    try:
        with(out/'native.log').open('w')as f:subprocess.run([str(binary),str(inp),str(out/'own-coefficient-matrix.json')],check=True,stdout=f,stderr=subprocess.STDOUT,timeout=50)
    except Exception as exc:
        (out/'blocked.json').write_text(json.dumps(dict(error=repr(exc),input_sha256=sha(inp),one_worker_only=True),indent=2)+'\n');raise
    own=read(out/'own-coefficient-matrix.json');assert own['all_passed']
    matrix=[[int(a)%PRIME for a in row]for row in cert['selected_integer_matrix']]
    assert own['matrix67']==matrix
    determinant=detmod(matrix);assert determinant==cert['selected_minor_modp']==9857!=0
    deps=[Path(__file__),native,ROOT/f'research/notes/{STEM}.md',*gates,*paths,
        *[folder/p for p in ['report.json','generators.json','native-certificate.json.gz']],
        ROOT/'code/q_new_family_full_degree24_source_relation_20261104_v63.py',
        ROOT/'code/q_new_family_full_degree24_source_relation_20261104_v63.cpp',
        ROOT/'research/notes/Q_NEW_FAMILY_FULL_DEGREE24_SOURCE_RELATION_20261104_v63.md']
    report=dict(id=STEM,all_passed=True,
        actual_source10_generator_maps_and_all70_semantic_definitions_bound=True,
        selected67_is_complete_universal_Inv24_QQ_basis=True,
        independent_reversed_LEX_nmod_coefficient_engine_no_author_integer_replay=True,
        all4489_selected_coefficients_and_minor9857_verified=True,
        full_Inv24_QQ_and_C_restriction_injective=True,
        NEW25_firstjet_defining_equation_degree_lower_bound=26,
        odd_degree25_absent_by_central_minus_identity=True,
        no_new_source_point_other_prime_higher_degree_probe_or_new_equation=True,
        input_sha256={str(p.relative_to(ROOT)):sha(p)for p in deps},
        evidence_sha256={p.name:sha(p)for p in out.iterdir()if p.is_file()},seconds=time.monotonic()-started)
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');signal.alarm(0)
    print(json.dumps({k:report[k]for k in ['all_passed','all4489_selected_coefficients_and_minor9857_verified','NEW25_firstjet_defining_equation_degree_lower_bound','seconds']},indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--work',type=Path,required=True);a=p.parse_args();main(a.output.resolve(),a.work.resolve())
