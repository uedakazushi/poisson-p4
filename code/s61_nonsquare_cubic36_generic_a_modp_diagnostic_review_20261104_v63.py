#!/usr/bin/env python3
"""READ ONLY exact source/modp-input comparison, no Singular or GB invocation."""
import argparse
import gzip
import hashlib
import json
import math
from pathlib import Path
import re
import resource
import signal
import time

ROOT=Path(__file__).resolve().parents[1]
STEM='S61_NONSQUARE_CUBIC36_GENERIC_A_MODP_DIAGNOSTIC_20261104_v63_REVIEW'
P=32003
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def bind_inputs(r):
    for n,h in r.get('input_sha256',{}).items():assert sha(ROOT/n)==h,n

def main(out):
    start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60)
    assert not out.exists();out.mkdir(parents=True)
    base=ROOT/'reports/S61_NONSQUARE_CUBIC36_GENERIC_A_MODP_DIAGNOSTIC_20261104_v63'
    ap=base/'report.json';author=read(ap);bind_inputs(author)
    pp=ROOT/'reports/S61_NONSQUARE_CUBIC36_FUNCTION_FIELD_PREFLIGHT_20261104_v63'
    prpath=pp/'report.json';pre=read(prpath);bind_inputs(pre);assert pre['all_read_only_format_provenance_and_statistics_checks_passed']
    statpath=pp/'statistics.json';stats=read(statpath);assert sha(statpath)==pre['statistics_sha256']
    assert author['ring']=='F32003(a)[e,d,z]' and author['column_map']=='R^106 -> R^36'
    for n,rec in author['artifacts'].items():
        p=base/n;assert sha(p)==rec['sha256'] and p.stat().st_size==rec['bytes']
    obpath=base/'OBSERVED_RECEIPT.json';obs=read(obpath);assert obs['report_sha256']==sha(ap)
    zp=base/obs['reproducible_compressed_input']['path'];zrec=obs['reproducible_compressed_input']
    assert sha(zp)==zrec['sha256'] and zp.stat().st_size==zrec['bytes']
    data=gzip.decompress(zp.read_bytes());assert hashlib.sha256(data).hexdigest()==zrec['decompressed_sha256']==sha(base/'input.sing')
    source=ROOT/'work/nonsquare-v62/cubic27-v6/after21';assert sha(source)==author['source_sha256']==stats['source_sha256']
    phasepath=ROOT/'work/nonsquare-v62/cubic27-v6/phase21-check.json';phase=read(phasepath)
    inverse={};total=zero=kept=0;degrees=[0,0,0];td=ed=0;M={}
    with source.open() as f:
        assert next(f).split()==['36','70','21']
        used=list(map(int,next(f).split()));retained=list(map(int,next(f).split()))
        assert used==phase['used_original_columns']==stats['used_original_columns']
        assert retained==phase['remaining_original_rows']==stats['retained_original_rows']
        for i in range(36):
            for j in range(70):
                tokens=next(f).split();n=int(tokens[0]);assert len(tokens)==1+4*n
                row={};previous=None
                for k in range(n):
                    q=tokens[1+4*k];num,den=map(int,q.split('/')) if '/' in q else (int(q),1)
                    ex=tuple(map(int,tokens[2+4*k:5+4*k]));assert min(ex)>=0 and (previous is None or previous>ex);previous=ex
                    assert num and den>0 and math.gcd(abs(num),den)==1 and den%P
                    total+=1
                    for t in range(3):degrees[t]=max(degrees[t],ex[t])
                    td=max(td,sum(ex));ed=max(ed,sum(ex[1:]))
                    if den not in inverse:inverse[den]=pow(den,P-2,P)
                    c=num*inverse[den]%P
                    if c:row[ex]=c;kept+=1
                    else:zero+=1
                if row:M[i+1,j+1]=row
        assert not f.read().strip()
    assert (total,zero,kept,len(M))==(421667,14,421653,1631)
    assert degrees==pre['max_degree_a_e_d']==[12,18,7] and td==25 and ed==22
    # Parse the saved grouped input with a standalone restricted grammar.
    text=data.decode();prefix='ring R=(32003,a),(e,d,z),dp;\nmatrix M[36][106];\n'
    assert text.startswith(prefix);matrixpart,tail=text[len(prefix):].split('poly F=',1)
    parsed={};group=re.compile(r'\(([^()]*)\)((?:\*[ed](?:\^[0-9]+)?)*)')
    term=re.compile(r'(-?[0-9]+)(?:\*a(?:\^([0-9]+))?)?')
    for line in matrixpart.splitlines():
        mm=re.fullmatch(r'M\[([0-9]+),([0-9]+)\]=(.*);',line);assert mm
        key=tuple(map(int,mm.group(1,2)));assert key not in parsed
        expression=mm.group(3);position=0;poly={}
        while position<len(expression):
            if position:assert expression[position]=='+';position+=1
            gm=group.match(expression,position);assert gm;position=gm.end();cofs,mon=gm.groups();powers={'e':0,'d':0}
            for v,n in re.findall(r'\*([ed])(?:\^([0-9]+))?',mon):powers[v]=int(n) if n else 1
            pos=0
            while pos<len(cofs):
                if pos:assert cofs[pos]=='+';pos+=1
                tm=term.match(cofs,pos);assert tm;pos=tm.end();c=int(tm.group(1))%P
                raw=tm.group(0);aexp=int(tm.group(2)) if tm.group(2) else (1 if '*a' in raw else 0)
                ex=(aexp,powers['e'],powers['d']);assert ex not in poly and c;poly[ex]=c
        parsed[key]=poly
    assert parsed==M
    expected='500*a+90*e^2+3*e-24;\npoly B=50*a*d+18*d*e^2-9*d*e-2*d+60*e-40;\npoly C=100*a+18*e^2-9*e-2;\npoly rel=z*d*F*B*C-1;\nfor(int i=1;i<=36;i++){M[i,70+i]=rel;}\nmodule I=M;'
    assert tail.startswith(expected) and 'poly J=' not in text and not re.search(r'\bJ\b',tail.split('module G=std(I);')[0])
    sanitypath=base/'field-sanity.sing';sanity=sanitypath.read_text()
    assert 'matrix T[2][3]=a*e,1-a*e,0,0,0,a;' in sanity and 'matrix X[3][2]=1,0,1,0,0,1/a;' in sanity
    # Explicit product is [[ae+1-ae,0],[0,a/a]]=I2 over the field.
    assert author['sanity']['returncode']==0 and author['exact_field_std_sanity_passed']
    assert (base/'field-sanity.log').read_text()=='EXACT_FIELD_AND_STD_SANITY_PASS\n'
    assert (base/'main.log').read_text()=='EXACT_INPUT36_BY106_READY_NO_J\n'
    assert author['main']['returncode']==-9 and author['main']['stop_reason']=='shared_total_60s_deadline'
    assert not author['std_completed_marker_observed'] and not author['profile_completed_marker_observed'] and not author['profile_markers']
    assert 60<=author['seconds']<61 and author['hard_shared_seconds']==60 and author['child_address_space_limit_bytes']==2000000000
    assert sum(p.stat().st_size for p in base.rglob('*') if p.is_file())<32*1024*1024
    producer=ROOT/'code/s61_nonsquare_cubic36_generic_a_modp_diagnostic_20261104_v63.py';code=producer.read_text()
    for literal in ['deadline=started+60','preexec_fn=child_limit','start_new_session=True','os.killpg(proc.pid,signal.SIGKILL)',"'--cpus=0'",'2000000000']:
        assert literal in code
    assert author['no_auto_retry_or_extension'] and not author['QQ_promotion']
    codimpath=ROOT/'reports/S61_BINARY_HYPERPLANE_CUBIC_CODIMENSION_REVIEW_20261102.json'
    quarticpath=ROOT/'reports/S61_NONSQUARE_FIXED_QUARTIC_CODIMENSION_TRANSFER_20261103_v63_REVIEW.json'
    # These are different accepted geometric theorems; do not identify their cokernels.
    cd,qd=read(codimpath),read(quarticpath)
    assert cd['all_passed'] and qd['all_passed']
    paths=[Path(__file__),ROOT/f'research/notes/{STEM}.md',ap,prpath,statpath,obpath,zp,source,phasepath,producer,
           base/'input.sing',sanitypath,base/'field-sanity.log',base/'main.log',base/'OBSERVED_SCOPE.md',codimpath,quarticpath]
    report=dict(id=STEM,all_passed=True,review_kind='READ ONLY exact finite source/input comparison and timeout frontier audit',
        all2520_entries_and421667_QQ_terms_independently_parsed=True,all_source_denominators_p_units=True,
        terms_vanishing_mod32003=14,nonzero_reduced_terms=421653,nonzero_entries=1631,
        own_restricted_parser_compares_every_reduced_coefficient=True,source36_by70_column_orientation_and36_localizer_columns=True,
        ring='F32003(a)[e,d,z]',localizer='d*F*B*C',J_inverted=False,
        recorded_sanity_only_and_main_single_input_ready_then_deadline=True,recorded_wall_seconds=author['seconds'],
        no_STD_completion_quotient_unit_nonunit_dimension_vdim_annihilator_or_QQ_claim=True,
        prior_schur_prefix_p_integrality_or_physical_modp_equivalence_reaudited=False,
        cubic_and_quartic_codimension_theorems_distinguished=True,
        conditional_generic_QQ_a_finite_support_and_required_monic_V_no_escape_proof_recorded=True,
        Singular_solver_or_producer_invoked=False,new_point_field_GB_retry_or_expanded_source_regeneration=False,
        input_sha256={str(p.relative_to(ROOT)):sha(p) for p in paths},seconds=time.monotonic()-start)
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');signal.alarm(0)
    print(json.dumps({k:report[k] for k in ['all_passed','terms_vanishing_mod32003','Singular_solver_or_producer_invoked','seconds']},indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
