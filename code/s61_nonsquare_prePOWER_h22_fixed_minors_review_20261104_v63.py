#!/usr/bin/env python3
"""Independent fixed pre-POWER source/minor and named-strip witness review."""
import argparse,hashlib,json,resource,signal,subprocess,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_PREPOWER_H22_FIXED_MINORS_20261104_v63_REVIEW'
def read(p):
    with p.open()as f:return json.load(f)
def sha(p):
    with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def dump(p,d):p.write_text(json.dumps(d,indent=2)+'\n')
def main(out,resume=False):
    start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(58)
    if resume:assert out.exists()and(out/'source-preflight-frontier.json').exists()and not(out/'native.json').exists()and not(out/'report.json').exists()
    else:out.mkdir(parents=True,exist_ok=False)
    R=ROOT/'reports';pins={}
    def load(stem,file='certificate.json'):
        p=R/stem/file;pins[str(p.relative_to(ROOT))]=sha(p);return read(p)
    p3='S61_NONSQUARE_PREPOWER_H22_PURE3_FIXED_MINOR_20261104_v63';p4='S61_NONSQUARE_PREPOWER_H22_PURE4_NAMED_STRIP_20261104_v63'
    r3=load(p3,'report.json');r4=load(p4,'report.json');assert pins[f'reports/{p3}/report.json']=='56140959c7a6b5852f1da7dbaa2cbf9cb2306c6c59d84fb4c7fb13cc7e99b6a9'and pins[f'reports/{p4}/report.json']=='164ed080bee0c89f86c5fe7a1ad6d2533e45d46aa2e4d547da964e72223f2b3d'
    assert r3['all_passed']and r4['all_passed']
    for report in [r3,r4]:
        for p,h in report['input_sha256'].items():assert sha(ROOT/p)==h
    n3=load(p3,'native.json');n4=load(p4,'native.json');assert pins[f'reports/{p3}/native.json']==r3['evidence_sha256']['native.json']and pins[f'reports/{p4}/native.json']==r4['evidence_sha256']['native.json']
    q=load('S61_NONSQUARE_CUBIC_Q2_L_UNIT_COMPLETION_20261104_v63');cp=load('S61_NONSQUARE_CUBIC_ZW25_UNITS_20261104_v63','checkpoint.json');pd=load('S61_NONSQUARE_CUBIC_ZW25_PAIR_DUAL_PREFLIGHT_20261104_v63');h=load('S61_NONSQUARE_PURE_KERNEL12_SOURCE_PREFLIGHT_20261104_v63','prePOWER-h9-source.json');height=load('S61_NONSQUARE_PREPOWER_H22_PURE3_FRAGMENT_PREFLIGHT_20261104_v63')
    for stem in ['S61_NONSQUARE_CUBIC_Q2_L_UNIT_COMPLETION_20261104_v63_REVIEW','S61_NONSQUARE_CUBIC_ZW25_UNITS_20261104_v63_REVIEW','S61_NONSQUARE_CUBIC_LAST_Q9_POWER_BEZOUT_20261104_v63_REVIEW']:
        r=load(stem,'report.json');assert r['all_passed']
    nodes=q['scalar_DAG']['nodes'];purecols=q['new22_target_original32_columns'][:17];rows=[1,3,6];targets=[11,13,25];labels=h['original_target_labels'];assert labels==[11,13,14,15,17,18,19,22,25]
    assert [q['new15_source_original25_labels'][i]for i in rows]==[10,13,16]
    roots=[[q['remaining15x22_scalar_roots'][i][purecols.index(j)]for j in targets+[22]]for i in rows]
    assert [q['remaining15x22_scalar_roots'][0][q['new22_target_original32_columns'].index(j)]for j in labels]==h['source_prePOWER_DAG_roots']
    def atom(node):
        op=node[0]
        if op=='checkpoint':return cp['matrix25x32'][node[1]][node[2]]
        if op=='accepted_pair_second':return pd['Lframe_dual6x2'][node[1]][1]
        assert op=='new_unit_atom';return q['scalar_DAG']['new_unit_atoms_num_den'][node[1]]
    zcache={}
    def structuralzero(i):
        if i in zcache:return zcache[i]
        n=nodes[i];op=n[0]
        if op=='z':v=True
        elif op=='one':v=False
        elif op in ['checkpoint','accepted_pair_second','new_unit_atom']:v=not atom(n)[0]
        elif op=='m':v=structuralzero(n[1])or structuralzero(n[2])
        elif op=='s':v=n[1]==n[2]or(structuralzero(n[1])and structuralzero(n[2]))
        else:assert op=='a';v=structuralzero(n[1])and structuralzero(n[2])
        zcache[i]=v;return v
    for i in rows:
        for j in q['all5_remaining_heavy_quotient_columns_unchanged']:assert structuralzero(q['remaining15x22_scalar_roots'][i][q['new22_target_original32_columns'].index(j)])
    extra_tail_labels=[j for j in purecols if j not in labels];assert len(extra_tail_labels)==8
    extra_tail_roots=[q['remaining15x22_scalar_roots'][0][q['new22_target_original32_columns'].index(j)]for j in extra_tail_labels]
    for j in q['all5_remaining_heavy_quotient_columns_unchanged']:assert structuralzero(q['remaining15x22_scalar_roots'][0][q['new22_target_original32_columns'].index(j)])
    original=height['needed_scalar_DAG_roots'];assert len(original)==108 and set(original)=={z['old_DAG_root']for z in n3['node_values']}
    need=set();todo=original+h['source_prePOWER_DAG_roots']
    while todo:
        i=todo.pop()
        if i in need:continue
        need.add(i);n=nodes[i]
        if n[0]in ['a','s','m']:todo.extend(n[1:])
    ids=sorted(need);assert len(ids)==195;local={i:j for j,i in enumerate(ids)};text=[];literal=0;atom_count=0;maxatom=0
    def poly(p):
        nonlocal literal
        literal+=len(p);text.append(str(len(p)))
        for m,c in p:assert len(m)==3 and m[1]>=0 and m[2]>=0;text.append(' '.join(map(str,m))+' '+c)
    for i in ids:
        n=nodes[i];op=n[0]
        if op=='z':text.append('z')
        elif op=='one':text.append('o')
        elif op in ['a','s','m']:text.append(f'{op} {local[n[1]]} {local[n[2]]}')
        else:
            num,den=atom(n);assert len(den)==1 and den[0][0][1:]==[0,0]and den[0][0][0]>=0 and Q(den[0][1])
            atom_count+=1
            for pol in [num,den]:
                for m,c in pol:
                    v=Q(c);maxatom=max(maxatom,abs(v.numerator).bit_length(),v.denominator.bit_length())
            text.append(f'p {den[0][0][0]} {den[0][1]}');poly(num)
    text.append('108')
    for row in n3['node_values']:text.append(str(local[row['old_DAG_root']]));poly(row['Laurent_polynomial'])
    for root,p in zip(h['source_prePOWER_DAG_roots'],h['source_h9']):text.append(str(local[root]));poly(p)
    for row in roots:text.append(' '.join(str(local[i])for i in row))
    text.append(' '.join(map(str,[labels.index(22),*[labels.index(j)for j in targets]])))
    assert n3['fixed_source_actual135']==[18,39,68]and n3['fixed_columns_original32']==targets
    for row in n3['nine_actual_entries']:
        for p in row:poly(p)
    for p in n3['three_signed_cofactor_terms']:poly(p)
    poly(n3['exact_fixed3_determinant']);poly(n4['actual_h22']);poly(n4['exact_unscaled_pure4_minor'])
    assert n4['actual_h22']==h['source_h9'][labels.index(22)]and n4['original_fixed3']==n3['exact_fixed3_determinant']
    factors=[[[[1,0,0],'1']],[[[0,0,1],'1']],[[[1,0,0],'500'],[[0,2,0],'90'],[[0,1,0],'3'],[[0,0,0],'-24']],[[[1,0,1],'50'],[[0,2,1],'18'],[[0,1,1],'-9'],[[0,0,1],'-2'],[[0,1,0],'60'],[[0,0,0],'-40']],[[[1,0,0],'100'],[[0,2,0],'18'],[[0,1,0],'-9'],[[0,0,0],'-2']],[[[1,0,1],'5'],[[0,1,0],'-6'],[[0,0,0],'-4']]]
    for p in factors:poly(p)
    names=['a','d','F','B','C','J'];text.append(str(len(n4['steps'])))
    for s in n4['steps']:
        text.append(f"{names.index(s['factor'])} {s['power_before']} {int(s['exact'])}")
        for k in ['before','division_quotient','remainder']:poly(s[k])
    poly(n4['remaining_polynomial']);assert literal<60000;ip=out/'review-input.txt';ip.write_text(f'{literal} 195\n'+'\n'.join(text)+'\n');del text
    # Independently validate the saved rational-height envelope.
    for i,row in enumerate(n3['nine_actual_entries']):
        for j,p in enumerate(row):
            den,L1=map(int,height['nine_entry_integer_common_denominator_and_l1_majorants'][i][j]);vals=[Q(c)*den for m,c in p];assert all(v.denominator==1 for v in vals)and sum(abs(v.numerator)for v in vals)<=L1
    det=n3['exact_fixed3_determinant'];assert len(det)==1895 and max(sum(m)for m,c in det)==40
    den=int(height['determinant_common_QQ_denominator']);vals=[Q(c)*den for m,c in det];assert all(v.denominator==1 for v in vals)and max(abs(v.numerator)for v in vals)<=int(height['determinant_integer_coefficient_abs_upper'])
    assert len(n4['exact_unscaled_pure4_minor'])==818 and len(n4['remaining_polynomial'])==434 and max(sum(m)for m,c in n4['remaining_polynomial'])==20
    dump(out/'source-audit.json',{'all108_refs_source_bound':True,'extra_cached_Q9_ancestors':87,'all9_selected_Q9_coefficients_bound':True,'all8_additional_pure_Q9_tails_retained_lazily':True,'additional8_original32_labels':extra_tail_labels,'additional8_actual_prePOWER_roots':extra_tail_roots,'additional8_NOT_setzero_or_newly_evaluated':True,'all3_selected_rows_and_Q9_heavy5_structurally_zero':True,'source_rows_current15':rows,'actual135_labels':[18,39,68],'original25_source_labels':[10,13,16],'fixed_original32_columns':targets,'h22_column':22,'no_Hadj_or_B_denominator':True,'all_atom_denominators_QQ_times_a_power':True,'atom_count':atom_count,'max_atom_QQ_bits':maxatom,'all9_and_determinant_height_majorants_checked':True,'input_polynomial_terms':literal,'no_POWER_207coefficient_evaluation_or_author_import':True})
    cpp=Path(__file__).with_suffix('.cpp');exe=out/'review-native';flags=subprocess.check_output(['pkg-config','--cflags','--libs','flint'],text=True).split();subprocess.run(['g++','-O2','-std=c++17',str(cpp),'-o',str(exe),*flags],check=True,timeout=max(1,52-(time.monotonic()-start)))
    with(out/'native.log').open('w')as log:ret=subprocess.run([str(exe),str(ip),str(out/'native.json')],stdout=log,stderr=subprocess.STDOUT,timeout=max(1,56-(time.monotonic()-start)))
    audit=read(out/'native.json');assert ret.returncode==0 and audit['all_passed']
    proof={'Schur_Pluecker':'Over universal polynomial domain, localize at h=Q[q] and clear firstcolumn: det4=h*det(M/h)=det3/h^2. Multiplying yields a polynomial identity, so it holds before h inversion and at h=0.','exact_actual_result':'Actualpure4minor=a^4*C^3*R434, R nonzero nonconstant of totaldegree20.','chart_scope':'h22 is only a chartunit; no inverse of it is used to assert globalFULLU unit status. a,d,F,B,C,J are the only stripped named factors.','conclusion':'The fixed pure3 is generically independent and actualpure4 minor is source-bound. UnknownR remains; its zero set is not the fullrank defect locus.','not_claimed':['Fullkernel9 independence','WholeFULLUcubic57 orquartic239','A globalpure4unit','A newfactor/minor/point/prime/GB/search']}
    dump(out/'proof.json',proof);note=ROOT/f'research/notes/{ID}.md';assert note.exists();result={'id':ID,'all_passed':True,'frozen':True,'primary3_report_sha256':pins[f'reports/{p3}/report.json'],'primary4_report_sha256':pins[f'reports/{p4}/report.json'],'actualpure3_nonzero_terms':1895,'actualpure4_terms':818,'unknown_remaining_terms':434,'unknown_remaining_totaldegree':20,'named_unit_powers':[4,0,0,0,3,0],'h22_not_declared_globalunit':True,'native_arithmetic':audit,'input_sha256':pins,'review_source_sha256':{str(p.relative_to(ROOT)):sha(p)for p in [Path(__file__),cpp,note]},'evidence_sha256':{p.name:sha(p)for p in out.iterdir()if p.is_file()and p!=exe},'no_full57_or_kernel9_or_globalunit_claim':True,'seconds':time.monotonic()-start}
    assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024;signal.alarm(0);dump(out/'report.json',result);print(json.dumps({'id':ID,'all_passed':True,'native_seconds':audit['seconds'],'seconds':result['seconds'],'report_sha256':sha(out/'report.json'),'input_terms':literal},indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);p.add_argument('--resume-preflight',action='store_true');a=p.parse_args();main(a.output.resolve(),a.resume_preflight)
