#!/usr/bin/env python3
"""Exact source-bound ell² raw product33 and raw cokernel45 base/actions.

Uses only frozen conductor cubic coefficients. No source/implicit solve,
physical-difference evaluation, or comparison homology is performed.
"""
import argparse
import hashlib
import itertools
import json
import resource
import signal
import time
from fractions import Fraction as F
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_ELL2_RAW33_QUOTIENT_20261105_v63'

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def add(a,b,c=F(1)):
    o=dict(a)
    for m,v in b.items():
        w=o.get(m,F(0))+c*v
        if w:o[m]=w
        else:o.pop(m,None)
    return o
def scale(a,c):return {m:v*c for m,v in a.items() if v*c}
def pmul(a,b,cap=5):
    o={}
    for (r,l),v in a.items():
        for (s,k),w in b.items():
            m=(r+s,l+k)
            if m[1]<2 and sum(m)<=cap:
                x=o.get(m,F(0))+v*w
                if x:o[m]=x
                else:o.pop(m,None)
    return o
def compose2(records,a,b):
    out={};one={(0,0):F(1)}
    for mon,c in records:
        v=one
        for _ in range(mon[0]):v=pmul(v,a)
        for _ in range(mon[1]):v=pmul(v,b)
        out=add(out,v,F(c))
    return out
def rec(a):return [[list(m) if isinstance(m,tuple) else m,str(c)] for m,c in sorted(a.items()) if c]
def srec(a):return [[i,str(v)] for i,v in sorted(a.items()) if v]
def sparse_rank_columns(columns):
    basis={};receipts=[]
    for j,orig in enumerate(columns):
        v=dict(orig);steps=[]
        while v:
            p=min(v)
            if p not in basis:
                d=v[p];basis[p]=scale(v,1/d)
                receipts.append({'source_column':j,'pivot':p,'subtract':[[k,str(c)] for k,c in steps],'divide_by':str(d),'normalized':srec(basis[p])})
                break
            c=v[p];steps.append((p,c));v=add(v,basis[p],-c)
    return len(basis),receipts
def mv(M,v):
    out={}
    for j,c in v.items():out=add(out,M[j],c)
    return out
def product_columns(A,B):return [mv(A,v) for v in B]
def eye(n):return [{i:F(1)} for i in range(n)]
def blockdiag(blocks):
    out=[];off=0
    for M in blocks:
        out.extend([{i+off:c for i,c in v.items()} for v in M]);off+=len(M)
    return out

def main(out,peer):
    start=time.monotonic();signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
    resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024))
    out=out.resolve();out.mkdir(parents=True,exist_ok=False);inputs={}
    def read(p):
        p=ROOT/p;inputs[str(p.relative_to(ROOT))]=sha(p);return json.loads(p.read_text())
    def write(name,obj):
        p=out/name;p.write_text(json.dumps(obj,separators=(',',':'))+'\n')
        assert sum(q.stat().st_size for q in out.iterdir() if q.is_file())<32*1024*1024
        return sha(p)
    progress={'phase':'source gate','numerical_rank_started':False}
    write('phase.json',progress)
    proof=read('reports/AMBIENT_B18_ELL2_RAW45_TRUNCATION_PREFLIGHT_20261105_v63/report.json')
    assert proof['all_passed']
    pp=ROOT/peer;pr=json.loads(pp.read_text());assert pr['all_passed'];inputs[peer]=sha(pp)
    assert pr.get('primary_sha256',proof.get('unused'))==sha(ROOT/'reports/AMBIENT_B18_ELL2_RAW45_TRUNCATION_PREFLIGHT_20261105_v63/report.json') or any(
      p.endswith('AMBIENT_B18_ELL2_RAW45_TRUNCATION_PREFLIGHT_20261105_v63/report.json') and h==sha(ROOT/p)
      for p,h in pr.get('input_sha256',{}).items()),'independent peer must pin exact repaired primary'
    src=read('reports/AMBIENT_B18_CONDUCTOR_PRODUCT_SPECIAL_FIBER20_20261104_v63/actual-special-fiber20-coefficients.json')
    source_path='reports/AMBIENT_B18_CONDUCTOR_PRODUCT_SPECIAL_FIBER20_20261104_v63/actual-special-fiber20-coefficients.json'
    assert inputs[source_path]==proof['input_sha256'][source_path]
    beta2=[];beta3=[];M2=[[{} for _ in range(7)] for _ in range(7)]
    for i,row in enumerate(src['beta2_before_ell']):
        assert all(not any(m[:7]) and sum(m)==2 for m,c in row)
        beta2.append([[m[7:],c] for m,c in row])
    for i,row in enumerate(src['beta3_before_ell']):
        plane=[]
        for m,c in row:
            assert sum(m)==3
            if any(m[:7]):
                assert sum(m[:7])==1
                j=m[:7].index(1);M2[i][j]=add(M2[i][j],{tuple(m[7:]):F(c)})
            else:plane.append([m[7:],c])
        beta3.append(plane)
    # Literal normal-form identity beta3=plane+M2*z; (I+M2)^-1 is unit.
    rebuilt=[]
    for i in range(7):
        terms={tuple([0]*7+m):F(c) for m,c in beta3[i]}
        for j in range(7):
            for rs,c in M2[i][j].items():
                m=[0]*7+list(rs);m[j]=1;terms=add(terms,{tuple(m):c})
        assert terms=={tuple(m):F(c) for m,c in src['beta3_before_ell'][i]}
        rebuilt.append(rec(terms))
    flatten={'M2_plane_coefficients7_by7':[[rec(p) for p in row] for row in M2],
      'beta3_source_reconstructed_full9':rebuilt,'unit_matrix':'I7+M2; M2 order2; determinant constant1',
      'through_degree3':'(I7-M2)*(z+beta2+beta3_plane+M2*z)=z+beta2+beta3_plane modulo m4',
      'no_unknown_fourth_coefficients_set_zero':True}
    r={(1,0):F(1)};ell={(0,1):F(1)};s=add(r,ell,F(1,5))
    beta=[add(compose2(a,r,s),compose2(b,r,s)) for a,b in zip(beta2,beta3)]
    gA=compose2(src['gA_plane3'],r,s)
    bot=[row[7:] for row in src['YB_in_YA_linear'][7:]]
    rb=add(scale(r,F(bot[0][0])),scale(s,F(bot[0][1])))
    sb=add(scale(r,F(bot[1][0])),scale(s,F(bot[1][1])))
    gB=compose2(src['gB_own_plane3'],rb,sb)
    assert gA.get((3,0))==F(125,24) and gB.get((3,0))==F(-15125,144)
    assert all(sum(m)==3 for m in gA) and all(sum(m)==3 for m in gB)
    # Ambient module coordinates (sector0=1, sectori+1=zi).
    ambient=[]
    for t in range(8):
        for l in range(2):
            for a in range(6-l-(t!=0)):
                ambient.append((t,a,l))
    ambient.sort(key=lambda m:(sum(m[1:])+(m[0]!=0),-m[0],m[2],m[1]))
    assert len(ambient)==74;index={m:i for i,m in enumerate(ambient)}
    def module_term(t,p):return {index[(t,a,l)]:c for (a,l),c in p.items() if (t,a,l) in index}
    relations=[]
    for i in range(7):
        for j in range(i+1,7):
            v=add(module_term(i+1,beta[j]),module_term(j+1,beta[i]),-1)
            relations.append((['pair',i,j],v))
    for i in range(7):relations.append((['gA',i],add(module_term(i+1,gA),module_term(0,pmul(gA,beta[i])))))
    for i in range(7):relations.append((['gB',i],module_term(i+1,gB)))
    assert len(relations)==35
    def plane_times(v,p):
        outv={}
        for k,c in v.items():
            t,a,l=ambient[k]
            for (b,h),d in p.items():
                q=(t,a+b,l+h)
                if q in index:outv=add(outv,{index[q]:c*d})
        return outv
    def raw_z(v,i):
        o={}
        for k,c in v.items():
            t,a,l=ambient[k]
            if not t:q=module_term(i+1,{(a,l):c})
            else:q=module_term(i+1,scale(pmul(beta[t-1],{(a,l):F(1)}),-c))
            o=add(o,q)
        return o
    rows=[];labels=[]
    for label,v in relations:
        if not v:continue
        order=min(sum(ambient[k][1:])+(ambient[k][0]!=0) for k in v)
        for l in range(2):
            for a in range(6-order-l):
                w=plane_times(v,{(a,l):F(1)})
                if w:rows.append(w);labels.append({'relation':label,'plane_multiplier':[a,l]})
    terms=sum(map(len,rows));assert len(rows)<=147 and terms<=1470
    assert terms+74*74+147*74+45*45*9<150000
    write('source-and-size-preflight.json',{'input_sha256':inputs,'basis74':ambient,
      'beta7_after_s_equals_r_plus_ell_over5':[rec(b) for b in beta],'gA':rec(gA),'gB':rec(gB),
      'relations35':[[label,srec(v)] for label,v in relations],
      'rows':len(rows),'columns':74,'input_nonzeros':terms,'module_shape_cap':[147,74],
      'coefficient_job_kind':'existing frozen polynomial finite quotient; no physical source evaluation',
      'limits':{'seconds':60,'bytes':2_000_000_000,'artifact_bytes':32*1024*1024,'total_coefficients':150000}})
    write('actual-unit-flattening-through3.json',flatten)
    progress={'phase':'quotient rank','numerical_rank_started':True};write('phase.json',progress)
    basis={};receipts=[]
    for j,orig in enumerate(rows):
        v=dict(orig);steps=[]
        while v:
            p=min(v)
            if p not in basis:
                c=v[p];basis[p]=scale(v,1/c)
                receipts.append({'pivot':p,'source_row':j,'subtract':[[k,str(a)] for k,a in steps],'divide_by':str(c),'normalized':srec(basis[p])})
                break
            c=v[p];steps.append((p,c));v=add(v,basis[p],-c)
    free=[i for i in range(74) if i not in basis];assert len(free)==33
    fi={p:i for i,p in enumerate(free)};reduced={}
    for i in reversed(range(74)):
        if i in fi:reduced[i]={fi[i]:F(1)}
        else:
            q={}
            for j,c in basis[i].items():
                if j!=i:q=add(q,reduced[j],-c)
            reduced[i]=q
    def nf(v):
        o={}
        for i,c in v.items():o=add(o,reduced[i],c)
        return o
    for v in rows:assert not nf(v)
    for label,v in relations:
        for i in range(7):assert not nf(raw_z(v,i))
    actions=[]
    for i in range(7):actions.append([nf(raw_z({k:F(1)},i)) for k in free])
    actions.append([nf(plane_times({k:F(1)},r)) for k in free])
    E=[nf(plane_times({k:F(1)},ell)) for k in free]
    actions.append(E)
    assert all(product_columns(A,B)==product_columns(B,A) for A in actions for B in actions)
    assert all(not c for c in product_columns(E,E))
    # Filtration certificate: normal forms never lower total degree.
    degree=lambda k:sum(ambient[k][1:])+(ambient[k][0]!=0)
    for k,v in reduced.items():assert all(degree(free[j])>=degree(k) for j in v)
    for A in actions:
        for j,v in enumerate(A):assert all(degree(free[k])>=degree(free[j])+1 for k in v)
    # Every truncated conductor product annihilates all33 basis vectors.
    def poly_action(p):
        ans=[{} for _ in range(33)]
        for (a,l),c in p.items():
            M=eye(33)
            for _ in range(a):M=product_columns(actions[7],M)
            for _ in range(l):M=product_columns(E,M)
            ans=[add(x,y,c) for x,y in zip(ans,M)]
        return ans
    A_g=poly_action(gA);B_g=poly_action(gB)
    IA=actions[:7]+[A_g]
    IB=[[add(x,y) for x,y in zip(actions[i],poly_action(beta[i]))] for i in range(7)]+[B_g]
    assert all(all(not col for col in product_columns(A,B)) for A in IA for B in IB)
    # Canonical6-dimensional quotients, same common r on both branches.
    b6=[(a,l) for l in range(2) for a in range(3)];i6={m:i for i,m in enumerate(b6)}
    branches=[];maps=[]
    for side,g in [('CA',gA),('CB',gB)]:
        c=-g.get((2,1),F(0))/g[(3,0)]
        def red6(p):
            o={}
            for (a,l),v in p.items():
                if l>=2:continue
                if a>=4 or (a>=3 and l):continue
                if a==3:a,l,v=2,1,v*c
                if v:o=add(o,{i6[(a,l)]:v})
            return o
        def mul6(a,b):
            return red6(pmul({b6[i]:v for i,v in a.items()},{b6[j]:w for j,w in b.items()}))
        er=red6(r);ee=red6(ell)
        ze=[{} for _ in range(7)] if side=='CA' else [red6(scale(p,-1)) for p in beta]
        Bact=[[mul6(z,{j:F(1)}) for j in range(6)] for z in ze+[er,ee]]
        Q=[]
        for k in free:
            t,a,l=ambient[k];v=red6({(a,l):F(1)})
            if t:v=mul6(ze[t-1],v)
            Q.append(v)
        assert all(product_columns(Q,A)==product_columns(B,Q) for A,B in zip(actions,Bact))
        qrank,qreceipts=sparse_rank_columns(Q);assert qrank==6
        branches.append({'side':side,'basis_plane_monomials':b6,'r3_equals_ell_r2_coefficient':str(c),
          'coordinate_actions9_sparse_columns':[[srec(v) for v in A] for A in Bact],
          'quotient_map6_by33_columns':[srec(v) for v in Q],'surjectivity_rank':qrank,'surjectivity_receipts':qreceipts})
        maps.append((Bact,Q))
    raw45=[blockdiag([maps[0][0][i],maps[1][0][i],actions[i]]) for i in range(9)]
    s_actions=[add(x,y,F(1,5)) for x,y in zip(actions[7],E)]
    evidence={'coordinate_order9':[f'zA_{i}' for i in range(7)]+['rA','ell'],
      'actual_CA_F9_order9':[f'zA_{i}' for i in range(7)]+['rA','sA'],
      'sA_equals_rA_plus_ell_over5':True,'basis74':ambient,'relation_rows_sparse':[srec(v) for v in rows],
      'relation_row_sources':labels,'rank41_receipts':receipts,'rank':len(basis),
      'actual_product_dimension':33,'free_column_indices':free,'actual33_basis':[ambient[k] for k in free],
      'all74_normal_forms':[srec(reduced[k]) for k in range(74)],
      'coordinate_actions9_sparse_columns':[[srec(v) for v in A] for A in actions],
      'sA_action_sparse_columns':[srec(v) for v in s_actions],
      'all9_actions_commute':True,'ell_squared_zero':True,'m6_zero_by_strict_degree_filtration':True,
      'all64_truncated_pair_products_annihilate_all33_basis_vectors':True,
      'actual_source_equality':'by repaired union-D6/truncate3 theorem and explicit unit flattening',
      'canonical_branch_quotients6':branches,
      'raw45_basis_sectors':{'tauA':list(range(6)),'tauB':list(range(6,12)),'tauA_tauB':list(range(12,45))},
      'raw45_coordinate_actions9_sparse_columns':[[srec(v) for v in A] for A in raw45],
      'physical12_difference_actions_not_evaluated':True,'comparison_T2_not_computed':True}
    write('actual-ell2-product33-and-raw45-actions.json',evidence)
    progress={'phase':'complete','rank41':len(basis),'dimension33':33,'raw45':45};write('phase.json',progress)
    evidence_hashes={p.name:sha(p) for p in out.iterdir() if p.is_file() and p.name!='phase.json'}
    slots=terms+sum(map(len,basis.values()))+sum(map(len,reduced.values()))+sum(len(c) for A in actions for c in A)+sum(len(c) for A in raw45 for c in A)
    assert slots<150000
    report={'id':STEM,'all_passed':True,'actual_ell2_product_dimension':33,'actual_ell2_raw_cokernel_dimension':45,
      'rank':41,'shape':[len(rows),74],'input_nonzeros':terms,'stored_QQ_slots':slots,
      'all9_coordinate_actions_exact':True,'all64_actual_source_pair_products_zero':True,
      'physical_source_or_gauge_solve_performed':False,'C2_or_T2_or_H0_computed':False,
      'input_sha256':inputs,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),
      'evidence_sha256':evidence_hashes,'seconds':time.monotonic()-start,
      'peak_RSS_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    write('report.json',report)
    print(json.dumps({k:report[k] for k in ['all_passed','actual_ell2_product_dimension','actual_ell2_raw_cokernel_dimension','shape','input_nonzeros','stored_QQ_slots','seconds']}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);p.add_argument('--lemma-peer',required=True);a=p.parse_args();main(a.output,a.lemma_peer)
