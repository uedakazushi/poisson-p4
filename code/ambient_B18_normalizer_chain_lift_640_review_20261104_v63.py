#!/usr/bin/env python3
"""Independent sparse-QQ/source-ID audit of the actual UP640=Q48 B operator.

No producer, source solve, scalar division, jet, GB or rank computation is run.
The polynomial alphabet consists of the SAME original C0 scalar function IDs.
"""
import argparse, hashlib, json, resource, signal, time
from fractions import Fraction as F
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
PRIMARY='AMBIENT_B18_NORMALIZER_CHAIN_LIFT_640_NEGATIVE_BOUNDARY_20261104_v63'
STEM='AMBIENT_B18_NORMALIZER_CHAIN_LIFT_640_20261104_v63_REVIEW'
Z={};O={():F(1)}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def add(*polys):
    r={}
    for p in polys:
        for m,c in p.items():
            q=r.get(m,F(0))+c
            if q:r[m]=q
            else:r.pop(m,None)
    return r
def sc(p,c):return {m:c*x for m,x in p.items() if c*x}
def mul(*polys):
    r=O
    for p in polys:
        q={}
        for m,c in r.items():
            for n,b in p.items():
                k=tuple(sorted(m+n));a=q.get(k,F(0))+c*b
                if a:q[k]=a
                else:q.pop(k,None)
        r=q;assert len(r)<150000
    return r
def atom(i):return {(i,):F(1)}
def va(*vs):return [add(*(v[i] for v in vs)) for i in range(4)]
def vs(v,c):return [mul(x,c) for x in v]
def scalar(p):return [p,Z,Z,Z]

def main(out):
    start=time.monotonic();signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
    out.mkdir(parents=True,exist_ok=False)
    folder=ROOT/'reports'/PRIMARY;rp=folder/'report.json';A=json.loads(rp.read_text())
    assert sha(rp)=='dbc64bb35c1d3708d90e45b5767c89ceb909e0b016165d753d9555a86e708df3'
    assert A['all_passed']
    ep=ROOT/A['evidence']['path'];assert sha(ep)==A['evidence']['sha256']=='d57cdec83150cde5f2f342a83fa817c71b60dfa36d158c9e8636b04d96fef379'
    E=json.loads(ep.read_text());D={}
    for label,r in A['input_hashes'].items():
        p=ROOT/r['path'];assert sha(p)==r['sha256'];D[label]=json.loads(p.read_text())
        if label.endswith('_peer'):assert D[label]['all_passed']
    assert E['scalar_source_artifacts']==A['input_hashes']
    S,G,C,N=[D[k] for k in ['source','relative','C0','normalizer']]
    full=C['DAG_scalar_nodes'];base=S['DAG_scalar_nodes']
    assert full[:len(base)]==base and full[:len(G['DAG_scalar_nodes'])]==G['DAG_scalar_nodes']
    rows=[18,19,28,29,60,67,72,63,64,68,73,74]
    assert E['physical_target_rows12']==rows
    ca,cb=[S['families'][k] for k in ['CA','CB']]
    aliases={};seed={};opaque=set();identity_records=[]
    lmemo={}
    def linear(i):
        if i in lmemo:return lmemo[i]
        n=base[i]
        if n[0]=='c':p={():F(n[1],n[2])} if n[1] else {}
        elif n[0]=='n':p=sc(linear(n[1]),-1)
        elif n[0]=='a':p=add(linear(n[1]),linear(n[2]))
        elif n[0]=='m':
            a,b=linear(n[1]),linear(n[2])
            p=mul(a,b) if (not a or not b or set(a)<={()} or set(b)<={()}) else atom(i)
        else:p=atom(i)
        lmemo[i]=p;return p
    for family in [ca,cb]:
        for r,i in family['actual_zero_odd_I0_generators']:
            j=family['complete_normalized_physical_w_coefficients'][0][r][1]
            assert linear(i)==linear(j);aliases[j]=i
            identity_records.append({'kind':'exact_QQ_linear_source_alias','odd_generator':i,'physical_odd_coefficient':j})
    def ref(i):
        if i not in seed:seed[i]=atom(i);opaque.add(i)
        return seed[i]
    def source(i):
        i=aliases.get(i,i)
        if i in seed:return seed[i]
        n=full[i];o=n[0]
        if o=='c':p={():F(n[1],n[2])} if n[1] else {}
        elif o=='a':p=add(source(n[1]),source(n[2]))
        elif o=='n':p=sc(source(n[1]),-1)
        elif o=='m':p=mul(source(n[1]),source(n[2]))
        elif o in ['r','g','i']:p=atom(i);opaque.add(i)
        elif o=='v':assert i in range(2,11);p=atom(i);opaque.add(i)
        else:raise AssertionError((i,n))
        seed[i]=p;return p
    cg=C['actual_C0_basis8_scalar_nodes'];assert len(cg)==8
    for i in cg:ref(i)
    ca_divs=C['all35_regular_coordinate_divisions']
    EA=[r['node'] for r in ca_divs[0]['all7_regular_coordinate_divisions']]
    for i in EA:ref(i)
    fA=ca['actual_zero_odd_I0_generators'][7][1]
    seed[fA]=add(source(cg[7]),*(mul(ref(e),source(z)) for e,z in zip(EA,cg[:7])))
    identity_records.append({'kind':'accepted_actual_CA_f63_tail','function':fA,'gA':cg[7],'coordinate_quotients':EA})
    for t,div in zip(ca_divs[1:],C['all4_known_regular_gA_divisions']):
        i=t['original_pure_CA_source_function'];coeff=[r['node'] for r in t['all7_regular_coordinate_divisions']]+[div['node']]
        seed[i]=add(*(mul(ref(j),source(z)) for j,z in zip(coeff,cg)))
        identity_records.append({'kind':'accepted_actual_CA_zero_graph_membership','function':i,'coefficient_nodes':coeff})
    eb=G['CB_f63_zero7_telescoping_coefficients'];bids=G['CB_scalar_odd_generators_normalized_coefficients'][:8]
    ebp=list(map(ref,eb));bp=list(map(ref,bids));fB=cb['actual_zero_odd_I0_generators'][7][1]
    seed[fB]=add(bp[7],*(mul(e,b) for e,b in zip(ebp,bp[:7])))
    identity_records.append({'kind':'accepted_actual_CB_f63_tail','function':fB,'gB':G['CB_g_restriction_node'],'coordinate_quotient_nodes':eb})
    memb=dict(G['all87_actual_odd_membership_lifts']);hs={}
    for r in rows[8:]:
        ids=memb[r][:8];hs[r]=[ref(i) if i else Z for i in ids]
        b=cb['complete_normalized_physical_w_coefficients'][0][r][1]
        seed[b]=add(*(mul(h,source(i)) for h,(_,i) in zip(hs[r],cb['actual_zero_odd_I0_generators'])))
        identity_records.append({'kind':'accepted_actual_CB_zero_graph_membership','physical_row':r,'function':b,'coefficient_nodes':ids})
    aa=[];ab=[];ba=[];bb=[]
    for r in rows:
        x,y=ca['complete_normalized_physical_w_coefficients'][0][r]
        u,v=cb['complete_normalized_physical_w_coefficients'][0][r]
        if x not in seed:ref(x)
        if u not in seed:ref(u)
        for i in [aliases.get(y,y),aliases.get(v,v)]:
            if i not in seed:ref(i)
        aa.append(source(x));ab.append(source(u));ba.append(source(y));bb.append(source(v))
    nuA=ref(ca['actual_tracefree_nu']);nuB=ref(cb['actual_tracefree_nu'])
    assert identity_records==E['accepted_actual_regular_function_identities']
    typed=E['operator_scalar_nodes'];memo={};support=0
    def op(i):
        nonlocal support
        if i in memo:return memo[i]
        n=typed[i];tag=n[0]
        if tag=='c':p={():F(n[1],n[2])} if n[1] else {}
        elif tag=='actual_C0_scalar_function':
            assert n[1] in opaque,('unexpected_opaque_function',n[1]);p=atom(n[1])
        elif tag=='sum':p=add(*(sc(op(j),F(a,b)) for j,a,b in n[1]))
        elif tag=='product':p=sc(mul(*(op(j) for j in n[3])),F(n[1],n[2]))
        else:raise AssertionError(n)
        support+=len(p);assert support<150000
        memo[i]=p;return p
    # Each closure entry must be literal source negation. It is not a relation
    # guessed from vanishing at a point, and it must match the saved operator.
    closures=E['negative_boundary_closure'];seen=set()
    for rec in closures:
        i,j=rec['opaque_function_node'],rec['negative_child_node']
        assert full[i]==['n',j] and (i,j) not in seen;seen.add((i,j))
        val=sc(source(i),-1)
        if j in seed:assert seed[j]==val
        else:seed[j]=val
        assert op(rec['child_bound_operator'])==val
    assert len(closures)==24
    for family in [ca,cb]:
        for r,_ in family['actual_zero_odd_I0_generators']:
            i=family['complete_normalized_physical_w_coefficients'][0][r][1]
            assert any(k==i for k,j in seen)
    def vm(v,w):
        q=[Z]*4
        for i,x in enumerate(v):
            for j,y in enumerate(w):
                q[i^j]=add(q[i^j],mul(x,y,nuA if i&j&1 else O,nuB if i&j&2 else O))
        return q
    def vector(ids):return [op(i) for i in ids]
    F12=[[add(x,sc(y,-1)),b,sc(c,-1),Z] for x,y,b,c in zip(aa,ab,ba,bb)]
    assert list(map(vector,E['actual_F12_C4']))==F12
    assert list(map(op,E['nu_A_B_operator_nodes']))==[nuA,nuB]
    oldmult=[scalar(O)]+[[Z,source(i),Z,Z] for _,i in ca['actual_zero_odd_I0_generators']]
    assert list(map(vector,E['old_actual_CA_Q0_multipliers9_C4']))==oldmult
    assert oldmult[8]!=[Z,source(cg[7]),Z,Z] # f63 is not silently gA.
    grad=[[O if i==j else Z for j in range(12)] for i in range(8)]
    for i in range(7):grad[7][i]=sc(ebp[i],-1)
    assert [[op(v) for v in row] for row in E['normalized_target_gradient8_by12']]==grad
    assert list(map(op,E['CB_g_target_unipotent_coefficients7']))==ebp
    YA=[[add(x,sc(y,-1)),b,Z,Z] for x,y,b in zip(aa,ab,ba)]
    ZA=YA[:8];ZA[7]=va(YA[7],*(vs(YA[i],sc(ebp[i],-1)) for i in range(7)))
    ZB=[[Z,Z,b,Z] for b in bp]
    records=[r for r in G['all99_CB_kernel_and87_graph_generators_evaluated_on_CA'] if r['weight']==0]
    assert len(records)==68 and records==E['CB_polynomial68_source_records']
    values=[];hadamards=[];sem=[]
    for k,r in enumerate(records):
        lab=r.get('label','')
        if lab.startswith('regular_Koszul_') or lab.startswith('cross_'):
            i,j=map(int,lab.split('_')[-2:]);
            H=[scalar(add(mul(bp[i],grad[j][t]),sc(mul(bp[j],grad[i][t]),-1))) for t in range(12)]
            val=va(vs(ZA[j],bp[i]),vs(ZA[i],sc(bp[j],-1)));desc={'kind':'linear_Koszul','i':i,'j':j}
        elif lab.startswith('product_'):
            i,j=map(int,lab.split('_')[-2:]);
            H=[va(vs(ZA[j],grad[i][t]),vs(ZB[i],grad[j][t])) for t in range(12)]
            val=va(vm(ZA[i],ZA[j]),scalar(sc(mul(nuB,bp[i],bp[j]),-1)));desc={'kind':'quadratic_product','i':i,'j':j}
        else:
            r0=r['raw_row'];j=rows.index(r0)
            H=[scalar(hs[r0][t]) if t<8 else scalar(sc(O,-1) if t==j else Z) for t in range(12)]
            val=va(scalar(add(ab[j],sc(aa[j],-1))),[Z,sc(ba[j],-1),Z,Z],*(vs(YA[t],hs[r0][t]) for t in range(8)))
            desc={'kind':'linear_graph','physical_row':r0}
        assert list(map(vector,E['actual_Hadamard68_by12_C4'][k]))==H
        assert E['CB_polynomial68_semantics'][k]==desc
        lhs=[source(r['normalized_CA_trace_pair'][0]),source(r['normalized_CA_trace_pair'][1]),Z,Z]
        assert val==lhs
        rhs=va(*(vm(f,h) for f,h in zip(F12,H)));assert rhs==val
        values.append(val);hadamards.append(H);sem.append(desc)
    # Bind the pruned normalizer graph to the exact SAME original scalar IDs.
    mp=dict(N['original_node_to_pruned_node']);ns={mp[i]:source(i) for i in set(seed)|set(aliases) if i in mp}
    ns[0]=Z;ns[1]=O
    def normal(i):
        if i in ns:return ns[i]
        n=N['DAG_scalar_nodes'][i]
        if n[0]=='c':p={():F(n[1],n[2])} if n[1] else {}
        elif n[0]=='a':p=add(normal(n[1]),normal(n[2]))
        elif n[0]=='n':p=sc(normal(n[1]),-1)
        elif n[0]=='m':p=mul(normal(n[1]),normal(n[2]))
        else:raise AssertionError(('unbound_pruned_source',i,n))
        ns[i]=p;return p
    Q=[]
    for j in range(12):
        for k in range(4):
            e=[Z]*4;e[k]=O;v=vm(F12[j],e)
            assert v==[normal(i) for i in N['normalizer_tensor_relation_columns4_by48'][4*j+k]]
            assert v==vector(E['exact_existing_Q48_C4'][4*j+k]);Q.append(v)
    U=[[normal(i) for i in c] for c in N['actual_C0_free9_to_E4_columns']]
    assert U==list(map(vector,E['exact_existing_U9_C4']))
    P=C['honest_D0_relation_columns'];oldpairs=G['actual_Qn_R0_squared_generator_matrices']['0']
    assert P[:28]==C['all28_C0_Koszul_kernel_columns'] and P[28:]==C['all612_Q0_lifted_columns']
    assert len(P)==640 and len(oldpairs)==612
    B=E['B48_by640_sparse_columns'];assert B[:28]==[[]]*28 and len(B)==640
    for j,c in enumerate(P[:28]):assert va(*(vs(u,source(v)) for u,v in zip(U,c)))==[Z]*4
    count=0
    for k,(val,H) in enumerate(zip(values,hadamards)):
        for m,mult in enumerate(oldmult):
            t=9*k+m;j=28+t;pair=oldpairs[t]
            expect=vm(val,mult);actual=[source(pair[0]),source(pair[1]),Z,Z]
            assert expect==actual
            lifted=va(*(vs(u,source(v)) for u,v in zip(U,P[j])));assert lifted==actual
            rec=E['all612_source_pair_and_honest_lift_bindings'][t]
            assert rec=={'D0_column':j,'CB_polynomial':k,'old_CA_multiplier':m,'actual_Q0_pair':pair,'honest_C0_lift':P[j]}
            predicted={}
            for l,h in enumerate(H):
                for c,p in enumerate(vm(h,mult)):
                    if p:predicted[4*l+c]=p
            actualB={i:op(v) for i,v in B[j]};assert actualB==predicted
            rhs=va(*(vs(Q[i],p) for i,p in actualB.items()));assert rhs==lifted
            count+=len(actualB)
    assert count==2808 and len(typed)==17849
    first=ROOT/'reports/AMBIENT_B18_NORMALIZER_CHAIN_LIFT_640_20261104_v63/report.json'
    assert sha(first)==A['first_failed_attempt_report_sha256'] and not json.loads(first.read_text())['all_passed']
    result={'all68_source_CA_pairs_verified':True,'all68_target_Hadamards_verified':True,
      'all48_exact_Q_columns_verified':True,'all9_exact_U_columns_verified':True,
      'all612_old_f63_multiplier_and_honest_C0_lift_bindings_verified':True,
      'all640_UP_QB_identities_verified':True,'first28_B_columns_zero':True,
      'all26_accepted_regular_identity_source_bindings_verified':True,
      'all16_source_odd_aliases_exact_QQ_linear':True,'negative_boundary_closure_count':24,
      'actual_regular_r_g_nodes_not_redivided_or_set_zero':True,
      'typed_operator_nodes':len(typed),'B_nonzero_entries':count,
      'independent_polynomial_memo_nodes':len(memo),'independent_polynomial_memo_total_terms':support,
      'source_translation_memo_nodes':len(seed),'source_opaque_function_ids':sorted(opaque),
      'no_complete_specialization_kernel_K_torsion_or_D0_depth_result':True}
    ev=out/'independent-source-and-operator-checks.json';ev.write_text(json.dumps(result,indent=2)+'\n')
    note=ROOT/'research/notes'/f'{STEM}.md'
    files=[rp,ep,first,Path(__file__),note,ROOT/'code/ambient_B18_normalizer_chain_lift_640_negative_boundary_20261104_v63.py']
    files += [ROOT/r['path'] for r in A['input_hashes'].values()]
    for n in ['phase0-negative-boundary-map.json','phase1-Hadamard68-operator.json']:
        files.append(folder/n)
    report={'id':STEM,'all_passed':True,'actual_chain_identity_UP640_equals_Q48B_accepted':True,
      'independent_original_source_ID_and_sparse_QQ_operator_arithmetic':True,
      'B_shape':[48,640],'B_nonzero_entries':count,'operator_nodes':len(typed),
      'first_failed_cap_preserved_and_not_used_as_theorem':True,
      'no_author_worker_source_solve_jet_or_scalar_inverse_search':True,
      'K_torsion_specialized_kernel_D0_depth_claimed':False,
      'input_sha256':{str(p.relative_to(ROOT)):sha(p) for p in files},'evidence_sha256':{ev.name:sha(ev)},
      'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:report[k] for k in ['id','all_passed','seconds','max_rss_kib']}))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
