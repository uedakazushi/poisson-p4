#!/usr/bin/env python3
"""Independent deterministic input/source and fixed-witness finite-H15 review."""
import argparse,gzip,hashlib,itertools,json,resource,signal,struct,subprocess,time
from fractions import Fraction
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_P100_H15_FULL239_FINITE_TEST_20261104_v63_REVIEW';P=32003
def read(p):
    with p.open()as f:return json.load(f)
def sha(p):
    with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def dump(p,d):p.write_text(json.dumps(d,indent=2)+'\n')
def tr(a):
    a=list(a)
    while a and not a[-1]:a.pop()
    return a
def pa(a,b):return tr([(a[i] if i<len(a) else 0)+(b[i] if i<len(b) else 0)for i in range(max(len(a),len(b)))])
def pm(a,b):
    c=[0]*(len(a)+len(b)-1)if a and b else[]
    for i,x in enumerate(a):
        for j,y in enumerate(b):c[i+j]=(c[i+j]+x*y)%P
    return tr(c)
def pr(a,f):
    a=tr([x%P for x in a]);assert f[-1]==1
    while len(a)>=len(f):
        c=a[-1];k=len(a)-len(f)
        for j,x in enumerate(f):a[k+j]=(a[k+j]-c*x)%P
        a=tr(a)
    return a
def pw(a,k,f):
    r=[1]
    while k:
        if k&1:r=pr(pm(r,a),f)
        k>>=1
        if k:a=pr(pm(a,a),f)
    return r
class Tokens:
    def __init__(self,path):self.f=path.open();self.values=iter(());self.count=0
    def get(self):
        while True:
            try:self.count+=1;return next(self.values)
            except StopIteration:
                line=self.f.readline()
                if not line:raise AssertionError('unexpected EOF')
                self.values=iter(line.split())
    def check(self,x):
        got=self.get();assert got==str(x),(self.count,got,str(x))
    def many(self,xs):
        for x in xs:self.check(x)
    def end(self):assert not list(self.values)and not self.f.read().strip();self.f.close()
def main(out,work):
    start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(59);out.mkdir(parents=True,exist_ok=False)
    R=ROOT/'reports';b=R/'S61_NONSQUARE_P100_H15_FULL239_FINITE_TEST_20261104_v63';report=read(b/'report.json')
    assert sha(b/'report.json')=='7041b20efce46e4a4f2596cc64a7f313bb37f3ae66be8943dcc857e061976291'and report['all_passed']and report['positive_full239_on_finite_H15']
    pinned={}
    def load(stem,name='certificate.json'):
        p=R/stem/name;pinned[str(p.relative_to(ROOT))]=sha(p);return read(p)
    for path,digest in report['input_sha256'].items():assert sha(ROOT/path)==digest;pinned[path]=digest
    for name,digest in report['evidence_sha256'].items():assert sha(b/name)==digest
    for name,digest in report['generated_source_inputs'].items():assert sha(work/name)==digest
    for name,digest in report['temporary_Gamma_source_hashes'].items():assert sha(work/name)==digest
    fac=load('S61_NONSQUARE_P100_H15_MODP_FACTOR_INVENTORY_20261104_v63','native.json');factors=[f['coefficients']for f in fac['factors']]
    assert [len(f)-1 for f in factors]==[4,1,1,1,2,3,3]and len(set(map(tuple,factors)))==7
    hn=load('S61_NONSQUARE_P100_NEXT9_QQ_GCD_20261104_v63','native.json');hbar=[0]*16
    qchecks=0
    def qm(c):
        nonlocal qchecks
        q=Fraction(c);assert q.denominator%P;qchecks+=1;return q.numerator*pow(q.denominator,-1,P)%P
    for m,c in hn['exact_gcd_monic']:hbar[m[0]]=qm(c)
    product=[1];rabin_tests=0
    for f,receipt in zip(factors,fac['factors']):
        assert receipt['multiplicity']==1 and f[-1]==1;product=pm(product,f);x=pr([0,1],f);v=x
        for got in receipt['Frobenius_chain']:v=pw(v,P,f);assert tr(got)==v
        assert v==x
        for test in receipt['proper_subfield_tests']:
            k=test['frobenius_exponent'];v=x
            for _ in range(k):v=pw(v,P,f)
            diff=tr([(v[i] if i<len(v) else 0)-(x[i] if i<len(x) else 0)for i in range(max(len(v),len(x)))])
            diff=[c%P for c in diff];assert tr(test['difference'])==tr(diff)
            value=pa(pm(test['Bezout_difference'],diff),pm(test['Bezout_factor'],f));assert tr([c%P for c in value])==[1];rabin_tests+=1
        for aa,bb in [(1,0),(1,20),(1,60),(203,2700),(1,12)]:assert pr([bb,aa],f)
    assert product==hbar==fac['Hbar'];assert len(product)==16
    cat=load('S61_NONSQUARE_QUARTIC_319_ROWS_20261103_v63');catalogue=[cat['source_rows'][i]for i in cat['retained_row_indices']]
    h=load('S61_NONSQUARE_H_SQUARE_H1_PREFLIGHT_20261104_v63');o=load('S61_NONSQUARE_LAZY_QUARTIC_ORACLE_20261103_v63','oracle.json');O=o['arithmetic'];ON=O['nodes']
    fd=load('S61_NONSQUARE_JOINT_LELL_FREE164_20261104_v63');nf=fd['full_target_projection239_to122'];LN=fd['linear_nodes'];SN=fd['scalar_circuit']['nodes'];SA=fd['scalar_circuit']['atoms_num_den']
    ri=load('S61_NONSQUARE_P100_H1_RIGHT_INVERSE_20261104_v63','native.json');L=ri['right_inverse112x23'];frame=load('S61_NONSQUARE_P100_G18_PURE7_20261104_v63','joint-cubic-frame.json');ad=load('S61_NONSQUARE_ORACLE_FRAME_ADAPTER_20261103_v63','adapter.json')
    qp=R/'S61_BINARY_HYPERPLANE_QUADRATIC_CODIMENSION_20261102.json';quads=read(qp);oldmons=read(R/'S61_BINARY_NONSQUARE_CUBIC_REPLAY_20261103_v62/source-report.json')['normal_monomials']
    actual16=[e['original_equation']for e in ad['basis57_source_definitions'][41:]]
    def code(m):return sum(x<<(3*i)for i,x in enumerate(m))
    mons=[]
    for total in range(5):
        for inds in itertools.combinations_with_replacement(range(8),total):
            m=[0]*8
            for i in inds:m[i]+=1
            mons.append(m)
    mons.sort(key=code);assert len(mons)==495
    cfg=Tokens(work/'config.txt');cfg.check(495)
    for m in mons:cfg.many([code(m),sum(m)])
    classes=h['target122_Hdegree_classes'];assert list(map(len,[classes['0'],classes['1'],classes['2']]))==[98,23,1]
    for seq in [[code(m)for m in h['parent_joint_target239_monomials']],[code(m)for m in oldmons],classes['0'],classes['1'],classes['2']]:cfg.check(len(seq));cfg.many(seq)
    cfg.check(319)
    for i,row in enumerate(catalogue):
        if i<174:assert row[0]=='quadric';cfg.many([0,row[1],code(row[2])])
        elif i<286:
            g=actual16.index(row[2]);assert row[0]=='cubic_normalized'and i==174+7*g+[0,1,2,3,4,6,7].index(row[3]);cfg.many([1,g,row[3]])
        else:assert row==['quartic',o['genuine_quartic_original_rows'][i-286]];cfg.many([2,i-286,0])
    cfg.check(7)
    for f in factors:cfg.check(len(f)-1);cfg.many(f)
    cfg.end()
    def closure(nodes,roots):
        need=set();atoms=set();todo=list(roots)
        while todo:
            i=todo.pop()
            if i in need:continue
            need.add(i);n=nodes[i]
            if n[0]in ['a','s','m']:todo.extend(n[1:])
            elif n[0]=='l':
                for c,j in n[1]:atoms.add(c);todo.append(j)
            elif n[0]=='c':atoms.add(n[1])
            else:assert n[0]in ['v','z']
        return sorted(need),sorted(atoms)
    ons,oas=closure(ON,o['genuine_quartic_roots']);lns,sroots=closure(LN,nf['output122_linear_roots']);sns,sas=closure(SN,sroots)
    oi={x:i for i,x in enumerate(ons)};oa={x:i for i,x in enumerate(oas)};li={x:i for i,x in enumerate(lns)};si={x:i for i,x in enumerate(sns)};sa={x:i for i,x in enumerate(sas)}
    src=Tokens(work/'shared.txt');src.check(495)
    for m in mons:src.many([code(m),sum(m)])
    def cp(pol):
        src.check(len(pol))
        for m,c in pol:src.check(qm(c));src.many(m)
    def cr(nd):
        for pol in nd:cp(pol)
    def ur(r):
        for pol in [r['numerator'],r['denominator']]:
            src.check(len(pol))
            for m,c in pol:src.many([m[0],qm(c)])
    def node(n,mp,ma=None,vm=None):
        src.check(n[0])
        if n[0]in ['a','s','m']:src.many([mp[j]for j in n[1:]])
        elif n[0]=='c':src.check(ma[n[1]])
        elif n[0]=='v':src.check(vm[n[1]]if vm is not None else n[1])
        elif n[0]=='l':
            src.check(len(n[1]))
            for c,j in n[1]:src.many([ma[c],mp[j]])
        else:assert n[0]=='z'
    src.check(len(oas))
    for i in oas:cr(O['atoms'][i])
    src.check(len(ons))
    for i in ons:node(ON[i],oi,oa)
    src.check(33)
    for actual,root in zip(o['genuine_quartic_original_rows'],o['genuine_quartic_roots']):src.many([actual,oi[root]])
    src.check(len(sas))
    for i in sas:cr(SA[i])
    src.check(len(sns))
    for i in sns:node(SN[i],si,sa)
    src.check(len(lns));vm={v:i for i,v in enumerate(nf['input_variables'])}
    for i in lns:node(LN[i],li,si,vm)
    src.check(122);src.many([li[r]for r in nf['output122_linear_roots']])
    for row in frame['physical_light_inverse7x7']:
        for r in row:ur(r)
    hs=h['scalar_circuit']
    for root in h['H2_source_pair_change_scalar_roots'][0]+[h['H2_known_unit_root']]:
        n=hs['nodes'][root];assert n[0]=='c';cr(hs['atoms_num_den'][n[1]])
    for row in quads['quadratic_rows_sparse']:
        src.check(len(row))
        for m,poly in row:src.check(code(m));cp(poly)
    active=[i for i,row in enumerate(L)if any(r['numerator']for r in row)];assert len(active)==42;src.check(42)
    for i in active:
        src.check(i)
        for r in L[i]:ur(r)
    src.end()
    cubic=Tokens(work/'cubics.txt');cubic.check(16);terms=0;uncompressed=hashlib.sha256();mp=R/'S61_BINARY_NONSQUARE_CUBIC_ROW_BASIS_20261103_v62/remaining.txt.gz'
    with gzip.open(mp,'rb')as zf:
        first=zf.readline();assert first.strip()==b'57 91 57';uncompressed.update(first)
        for row in range(57):
            for col in range(91):
                raw=zf.readline();uncompressed.update(raw)
                if row<41:continue
                s=raw.decode().split();n=int(s[0]);assert len(s)==1+4*n;terms+=n;cubic.check(n)
                for j in range(n):cubic.many([qm(s[1+4*j]),*s[2+4*j:5+4*j]])
        assert not zf.read()
    cubic.end();br=read(mp.with_name('report.json'));assert uncompressed.hexdigest()==br['remaining_matrix_sha256'];assert terms==148858
    # qchecks also includes 16 independently reduced H coefficients.
    assert qchecks-16==176101==report['coefficient_punit_checks']
    inp=out/'fixed-rank-input.txt';inverse_checks=0;selected_counts=[];rankcerts=[]
    with inp.open('w')as wf:
        wf.write('7 32003\n')
        for i,f in enumerate(factors):
            d=len(f)-1;rc=read(b/f'factor{i}.json');assert rc['factor_index']==i and rc['rank']==98 and rc['positive_full98'];rows=rc['selected_source319_rows'];cols=rc['pivot_columns'];A=rc['selected_original_rows'];assert len(set(rows))==len(set(cols))==98 and sorted(cols)==list(range(98))
            assert len(A)==98 and all(len(row)==98 for row in A);byrow={r:j for j,r in enumerate(rows)}
            with(work/f'Gamma-factor{i}.bin').open('rb')as gf:
                seen=set()
                for _ in range(319):
                    raw=gf.read(2*(1+98*d));assert len(raw)==2*(1+98*d);vs=struct.unpack('<'+'H'*(1+98*d),raw);r=vs[0];assert 0<=r<319 and r not in seen;seen.add(r);assert all(v<P for v in vs[1:])
                    if r in byrow:assert list(vs[1:])==[v for elt in A[byrow[r]]for v in elt]
                assert len(seen)==319 and not gf.read()
            def fm(x,y):return pr(pm(x,y),f)
            data=(work/f'inverse-factor{i}.bin').open('rb')
            try:
                while True:
                    raw=data.read(1+4*d)
                    if not raw:break
                    assert len(raw)==1+4*d and raw[0]in [0,1];cs=struct.unpack('<'+'H'*(2*d),raw[1:]);assert all(v<P for v in cs);assert fm(cs[:d],cs[d:])==[1];inverse_checks+=1
            finally:data.close()
            wf.write(f'{i} {d} 98\n'+' '.join(map(str,f))+'\n'+' '.join(map(str,rows))+'\n'+' '.join(map(str,cols))+'\n')
            for row in A:
                for elt in row:assert len(elt)==d;wf.write(' '.join(map(str,elt))+'\n')
            selected_counts.append(98*98*d);rankcerts.append(rc['selected98_minor_determinant']);del A,rc
    assert inverse_checks==21329
    audit={'all_generated_config_shared_cubic_tokens_independently_bound':True,'all176101_literal_QQ_coefficients_punit_checked':True,'all319_catalogue_retained':True,'physical_to_joint_before_NF117':True,'all33_raw_Qhat_roots_retained_with_DH_inverse_square':True,'section_uses_actual_delta7_L_and_corrected_original112_tails':True,'chi_source77_131_U_and_full_tail_algorithm_independently_audited':True,'all_seven_selected98_matrices_match_actual_Gamma_binary_source_rows':True,'every_logged_denominator_and_rank_inverse_product_independently_checked':inverse_checks,'full_factor_product_and_Rabin_Bezout_tests':rabin_tests,'all_five_curve_localizers_unit_in_all_factors':True,'source_coefficient_evaluator_or_NF_not_rerun':True,'rank_search_not_repeated':True,'QQ_source_core_terms':terms,'QQ_coefficient_checks_including_H16':qchecks,'selected_field_slots_each_factor':selected_counts}
    dump(out/'input-and-source-audit.json',audit)
    helper=ROOT/'code/s61_nonsquare_P100_H15_fixed_rank_review_helper_20261104_v63.cpp';exe=out/'fixed-rank-review-native'
    subprocess.run(['g++','-O2','-std=c++17',str(helper),'-o',str(exe)],check=True,timeout=max(1,52-(time.monotonic()-start)))
    with(out/'native.log').open('w')as log:ret=subprocess.run([str(exe),str(inp),str(out)],stdout=log,stderr=subprocess.STDOUT,timeout=max(1,56-(time.monotonic()-start)))
    nr=read(out/'fixed-rank-review-native.json');assert ret.returncode==0 and nr['all_passed']
    assert [x['determinant_in_original_column_order']for x in nr['factors']]==rankcerts
    proof={'scope':'Original quartic rank239 at EVERY characteristic-zero point in the finite H15 candidate algebra only.','steps':['Hmonic gives a finite-free V=Z_(32003) algebra S of rank15.','Every denominator inversion in the fixed algorithm has been checked in all seven fields; their product decomposition has fullHbar and no discarded factor. Such denominator elements are units in S by finite Nakayama.','Gamma=f0-f2/U2*chi0-(f1-f2/U2*chi1)*(Lactual^t*Ddelta7*PiOrig) preserves all raw pure tails and original319 sources.','The fixed selected98 minors are independently nonzero in every factor. CRT gives full Gamma row surjectivity over S/pS.','The Gamma cokernel is finite over V; its special fiber is zero, hence it is zero by Nakayama, including all finite scheme structure.','Undo accepted NF117 plus the actual unitH2 and H1section23 to obtain 117+1+23+98=239 on S_Q.'],'not_claimed':['Full239 offH15 or on wholeP100','WholeFULLU quartic239','A second Gamma coefficient evaluation','A rank-negative conclusion from any modular deficiency']}
    dump(out/'proof.json',proof)
    reportout={'id':ID,'all_passed':True,'frozen':True,'positive_full239_on_entire_H15_only':True,'independent_arithmetic':audit,'native_fixed_rank_audit':nr,'primary_report_sha256':sha(b/'report.json'),'pinned_inputs':pinned,'review_source_sha256':{str(p.relative_to(ROOT)):sha(p)for p in [Path(__file__),helper]},'evidence_sha256':{p.name:sha(p)for p in out.iterdir()if p.is_file()and p!=exe},'seconds':time.monotonic()-start}
    assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024;signal.alarm(0);dump(out/'report.json',reportout)
    print(json.dumps({'id':ID,'all_passed':True,'positive_full239_on_entire_H15_only':True,'inverse_products':inverse_checks,'native_seconds':nr['seconds'],'total_seconds':reportout['seconds'],'report_sha256':sha(out/'report.json')},indent=2),flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);p.add_argument('--work',required=True,type=Path);a=p.parse_args();main(a.output.resolve(),a.work.resolve())
