#!/usr/bin/env python3
"""One authorized117-product same-plane universal Inv28 basis job."""
import argparse,gzip,json,resource,shlex,subprocess,time
from pathlib import Path
from q_new_family_v_zero_degree26_support_preflight_20261104_v63 import read,sha
from q_new_family_full_degree24_universal_basis_20261104_v63 import echelon,det_mod
ROOT=Path(__file__).resolve().parents[1]
ID='Q_NEW_FAMILY_DEGREE28_PRODUCT_BASIS_20261104_v63'
def main(out,work):
    start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False);work.mkdir(parents=True,exist_ok=False)
    pre=ROOT/'reports/Q_NEW_FAMILY_DEGREE28_PRODUCT_PREFLIGHT_20261104_v63';gate=read(pre/'report.json');assert gate['all_passed']
    for n,h in gate['input_sha256'].items():assert sha(ROOT/n)==h,n
    formal=read(pre/'saved-formal-low-columns.json.gz');certificate=read(pre/'certificate.json');chars=read(pre/'character-records.json');assert chars[28]['invariant_dimension']==108
    polys=[(4,formal['F4'])]+[(8,z)for z in formal['p1_p2_p3_p4']]+[(12,z)for z in formal['primitive12']]+[(16,formal['primitive16'])]+[(14,z)for z in formal['Inv14']]+[(24,z)for z in formal['Inv24']];assert len(polys)==80
    inp=work/'formal-input.txt'
    with inp.open('w')as f:
        f.write('32003 80\n')
        for d,z in polys:
            assert len(z)==(d+1)*(d+2)//2;f.write(f'{d} {len(z)}\n'+' '.join(map(str,z))+'\n')
    with gzip.open(out/'formal-input.txt.gz','wb')as f:f.write(inp.read_bytes())
    native=ROOT/'code/q_new_family_degree28_product_basis_20261104_v63.cpp';helper=ROOT/'code/q_new_family_full_degree24_universal_basis_20261104_v63.cpp';binary=work/'universal28-native'
    flags=shlex.split(subprocess.check_output(['pkg-config','--cflags','--libs','flint'],text=True));subprocess.run(['c++','-O2','-std=c++17',str(native),*flags,'-lgmp','-o',str(binary)],check=True,timeout=max(1,60-(time.monotonic()-start)))
    error=None;data=None;minor=None;complete=False
    try:
        raw=work/'native-certificate.json'
        with(out/'native.log').open('w')as f:subprocess.run([str(binary),str(inp),str(raw)],stdout=f,stderr=subprocess.STDOUT,check=True,timeout=max(1,60-(time.monotonic()-start)))
        assert raw.stat().st_size<32*1024*1024;data=read(raw);assert data['all_passed']and len(data['columns117'])==117
        columns=data['columns117'];assert all(len(z)==435 for z in columns)
        selected,rows=echelon(columns);assert selected==data['selected_columns']and rows==data['selected_coefficient_rows'];assert len(selected)==data['rank']<=108
        minor=det_mod([[columns[j][i]for j in selected]for i in rows]);assert minor
        complete=data['rank']==108
        with gzip.open(out/'native-certificate.json.gz','wb')as f:f.write(raw.read_bytes())
        definitions=dict(candidates=certificate['candidates'],degree24_basis_definitions=certificate['degree24_basis_definitions'],selected_basis_indices=selected,basis=[certificate['candidates'][i]for i in selected],coefficient_order='i=0..28,j=0..28-i: S^i*T^j')
        (out/'definitions.json').write_text(json.dumps(definitions,indent=2)+'\n')
        (out/'minor.json').write_text(json.dumps(dict(prime=32003,dimension=108,rank=len(selected),minor=minor,columns=selected,coefficient_rows=rows,matrix=[[columns[j][i]for j in selected]for i in rows]),indent=2)+'\n')
    except Exception as e:error=repr(e)
    paths=[Path(__file__),native,helper,ROOT/f'research/notes/{ID}.md',pre/'report.json',pre/'certificate.json',pre/'saved-formal-low-columns.json.gz',pre/'character-records.json']
    report=dict(id=ID,all_passed=error is None,status=data['status']if data else'bounded_failure',error=error,
        complete_QQ_Inv28_basis=complete,character_dimension28=108,candidate_count=117,coefficient_rank=data['rank']if data else None,minor_mod32003=minor,
        SAME_fixed_formal_plane=True,no_new_point_or_source28_restriction_or_rank=True,
        no_candidate_enlargement_or_higher_degree=True,
        input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},evidence_sha256={p.name:sha(p)for p in out.iterdir()if p.is_file()},seconds=time.monotonic()-start)
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(dict(all_passed=report['all_passed'],rank=report['coefficient_rank'],minor=minor,complete=complete,error=error,seconds=report['seconds'])))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--work',type=Path,required=True);a=p.parse_args();main(a.output,a.work)
