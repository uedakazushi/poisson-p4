#!/usr/bin/env python3
"""ONE rigorous honest9 enlarged degree2 image after actual CA zero7 restriction.

Only a negative truncated membership can imply actual nonmembership.  Unknown
regular-division derivatives are arbitrary formal deviations, never set to0.
"""
import argparse
from fractions import Fraction as F
import hashlib,json
from pathlib import Path
import resource,signal,time
from ambient_B18_CB_g_division_origins_20261104_v63 import evaluate,MON
from ambient_B18_top_graph_telescoping_20261104_v63 import ancestors
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_T3_HONEST9_LOW2_OBSTRUCTION_20261104_v63'
PM=[(),(0,),(1,),(0,0),(0,1),(1,1)];PI={m:i for i,m in enumerate(PM)}
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main(out):
    started=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));out.mkdir(parents=True,exist_ok=False)
    source=ROOT/'reports/AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63/actual-source-coefficients.json'
    fullp=ROOT/'reports/AMBIENT_B18_FULL_RELATIVE_INPUT_20261104_v63/full-relative186-DAG.json'
    jp=ROOT/'reports/AMBIENT_B18_T3_JOINT_CURVE_RESTRICTIONS_20261104_v63/joint-degree3-certificate.json'
    originp=ROOT/'reports/AMBIENT_B18_T3_DOUBLE_ANNIHILATOR_PREFLIGHT_20261104_v63/exact-residues-and-two-ideals.json'
    lp=ROOT/'reports/AMBIENT_B18_T3_DOUBLE_ANNIHILATOR_LENGTH3_20261104_v63/report.json'
    hp=ROOT/'reports/AMBIENT_B18_T3_HONEST10_PRESENTATION_20261104_v63/honest10-presentation-DAG.json'
    up=ROOT/'reports/AMBIENT_B18_T3_HONEST10_UNIT_REDUCTION_20261104_v63/unit-reduced-presentation-DAG.json'
    ur=up.parent/'report.json'
    old=read(source);full=read(fullp);joint=read(jp);orig=read(originp);honest=read(hp);unit=read(up)
    assert read(lp)['all_passed'] and read(ur)['all_passed']
    nodes=unit['DAG_scalar_nodes'];oldnodes=old['DAG_scalar_nodes'];assert nodes[:len(oldnodes)]==oldnodes
    sol={i:tuple(map(F,v)) for i,v in joint['all105_unknown_jets']}
    Q3=unit['reduced_D3_relations_after_dropping_only_literal_zero_columns'];top=unit['all12_transformed_top_columns']
    assert len(Q3)==2069 and len(top)==12 and all(len(c)==9 for c in Q3+top)
    newcoeff={q['node']:(t['original_source_function'],j,q['origin'])
              for t in honest['all_pure_source_Hadamard_membership_telescopes']
              for j,q in enumerate(t['all9_regular_coordinate_divisions'])}
    assert len(newcoeff)==792
    atomic={r['node'] for r in full['all_coordinate_divisions']}|{r['node'] for r in full['all4_known_regular_g_divisions']}
    higher={i:2 for i in full['CB_g_split_rA_plus_sB_nodes']}
    gtail=full['CB_g_restriction_node'];root=[i for pair in Q3+top for i in pair]
    needed=set();terminal=set();stack=root[:]
    while stack:
        i=stack.pop()
        if i in needed:continue
        needed.add(i)
        if i<len(oldnodes):terminal.add(i);continue
        if i in atomic or i in newcoeff or i==gtail:continue
        n=nodes[i]
        if n[0] in ['a','m']:stack.extend(n[1:])
        elif n[0] in ['n','i']:stack.append(n[1])
        elif n[0]=='c':pass
        else:raise AssertionError(('uncontrolled implicit node',i,n))
    usedold=sorted(ancestors(oldnodes,list(terminal)+[r for r,j,o in newcoeff.values()]))
    oldvals=evaluate(oldnodes,usedold,sol)
    atomids={i:i+2 for i in atomic};weights={0:1,1:1,**{i+2:higher.get(i,1) for i in atomic}}
    def degree(m):return sum(weights[i] for i in m)
    def add(a,b):
        c=dict(a)
        for m,v in b.items():
            c[m]=c.get(m,F(0))+v
            if not c[m]:c.pop(m)
        return c
    def scale(a,s):return {m:v*s for m,v in a.items() if v*s}
    def mul(a,b):
        c={}
        for m,v in a.items():
            for n,w in b.items():
                z=tuple(sorted(m+n))
                if degree(z)>2:continue
                c[z]=c.get(z,F(0))+v*w
                if not c[z]:c.pop(z)
        return c
    exact_new_coefficients={}
    for i,(root,j,origin) in newcoeff.items():
        if j<7:p={():F(origin)} if F(origin) else{}
        else:
            p={}
            for m,v in zip(MON,oldvals[root]):
                if not v:continue
                a,b=m
                if j==7 and b==0 and a>=1:dest=tuple([0]*(a-1))
                elif j==8 and b>=1:dest=tuple([0]*a+[1]*(b-1))
                else:continue
                if len(dest)<=2:p[dest]=v
        assert p.get((),F(0))==F(origin)
        exact_new_coefficients[i]=p
    pol={}
    for i in sorted(needed):
        if i<len(oldnodes):
            p={tuple([0]*m[0]+[1]*m[1]):v for m,v in zip(MON,oldvals[i]) if sum(m)<=2 and v}
        elif i in newcoeff:p=exact_new_coefficients[i]
        elif i in atomic:
            c=F(orig['all339669_scalar_origin_values_evaluated'][i]);p={(atomids[i],):F(1)}
            if c:p[()]=c
        elif i==gtail:p={} # actual CB g has order3 in its zero7 coordinate plane
        else:
            n=nodes[i];op=n[0]
            if op=='c':p={():F(n[1],n[2])} if n[1] else{}
            elif op=='n':p=scale(pol[n[1]],-1)
            elif op=='a':p=add(pol[n[1]],pol[n[2]])
            elif op=='m':p=mul(pol[n[1]],pol[n[2]])
            elif op=='i':
                a=pol[n[1]];c=a.get((),F(0));assert c
                h=add(a,{():-c});p=add({():1/c},add(scale(h,-1/c**2),scale(mul(h,h),1/c**3)))
            else:raise ValueError((i,n))
        pol[i]=p
    # Exact known part plus a safe vector-space enlargement of all unknown jets.
    def vectors(pair):
        both={}
        for row,node in enumerate(pair):
            for m,v in pol[node].items():
                both.setdefault(m,[F(0)]*9)[row]+=v
        known=[F(0)]*54;unc=[]
        for m,co in both.items():
            phy=tuple(i for i in m if i<2);unknown=tuple(i for i in m if i>=2)
            if not unknown:
                assert phy in PI
                for row in range(9):known[6*row+PI[phy]]+=co[row]
            else:
                lo=sum(weights[i] for i in unknown)
                for extra in PM:
                    if lo<=len(extra)<=2-len(phy):
                        dest=tuple(sorted(phy+extra));vec=[F(0)]*54
                        for row in range(9):vec[6*row+PI[dest]]=co[row]
                        if any(vec):unc.append(vec)
        return known,unc
    def shift(v,m):
        w=[F(0)]*54
        for row in range(9):
            for k,n in enumerate(PM):
                dest=tuple(sorted(m+n))
                if dest in PI:w[6*row+PI[dest]]+=v[6*row+k]
        return w
    basis=[];provenance=[]
    def insert(v,label):
        v=v[:]
        for pivot,row in basis:
            if v[pivot]:
                c=v[pivot];v=[x-c*y for x,y in zip(v,row)]
        pivot=next((i for i,x in enumerate(v) if x),None)
        if pivot is None:return False
        c=v[pivot];v=[x/c for x in v];basis.append((pivot,v));provenance.append(label);return True
    qrecords=[];allQvectors=[]
    for index,pair in enumerate(Q3):
        known,unc=vectors(pair)
        assert all(known[6*row]==0 for row in range(9))
        qrecords.append({'column':index,'known_part':list(map(str,known)),'uncertainty_vectors':[list(map(str,u)) for u in unc]})
        for kind,v in [('known',known)]+[('unknown',u) for u in unc]:
            for m in PM:
                w=shift(v,m)
                if any(w):insert(w,[index,kind,list(m)]);allQvectors.append(w)
    def residual(v,B):
        v=v[:]
        for p,b in B:
            if v[p]:c=v[p];v=[x-c*y for x,y in zip(v,b)]
        return v
    top_records=[];negative=[]
    for i,pair in enumerate(top):
        known,unc=vectors(pair);oldbasis=[(p,b[:]) for p,b in basis]
        # Include independent unknown variations of this top too: never assume them0.
        backup=list(basis);bp=list(provenance)
        for j,u in enumerate(unc):insert(u,['top_uncertainty',i,j])
        r=residual(known,basis)
        if any(r):negative.append(i)
        top_records.append({'column':i,'known_part':list(map(str,known)),
                            'uncertainty_vectors':[list(map(str,u)) for u in unc],
                            'outside_even_enlarged_Q_and_own_uncertainty':bool(any(r)),
                            'remainder_after_enlarged_exact_linear_elimination':list(map(str,r)),
                            'formal_pair_polynomial_terms':[[[list(m),str(v)] for m,v in sorted(pol[node].items())] for node in pair]})
        basis[:]=backup;provenance[:]=bp
    ev={'specialization':'honest9 presentation on actual CA zero7 plane with original105 source3jets, then degree<=2',
        'physical_monomials':['1','r','s','r^2','rs','s^2'],'vector_row_order':'nine generator rows, each6monomials',
        'all792_pure_CA_Hadamard_coefficients_exact_specialized_low2':{
            str(i):[[list(m),str(v)] for m,v in p.items()] for i,p in exact_new_coefficients.items()},
        'new_coefficient_rules':'z-prefix coefficients are EXACTconstants; h_r=b(r,0)/r,h_s=(b(r,s)-b(r,0))/s,source3=>h2',
        'atomic_unknown_function_deviations':[[i,atomids[i],higher.get(i,1),orig['all339669_scalar_origin_values_evaluated'][i]] for i in sorted(atomic&needed)],
        'unknown_deviations_not_set_to_zero':True,'CB_g_tail_known_order3_node':gtail,
        'all2069_Q_column_known_and_enlarged_uncertain_vectors':qrecords,
        'Q_enlarged_span_basis':[[p,list(map(str,b))] for p,b in basis],'Q_basis_provenance':provenance,
        'all12_top_results':top_records,'negative_nonmembership_columns':negative,
        'original_source_ancestor_count':len(usedold),'post_source_needed_node_count':len(needed-terminal),
        'literal_formal_top2_all_zero':all(not pol[node] for pair in top for node in pair)}
    encoded=json.dumps(ev,separators=(',',':'))+'\n';assert len(encoded.encode())<32*1024*1024
    (out/'honest9-low2-certificate.json').write_text(encoded)
    report={'id':STEM,'all_passed':True,'status':'rigorous_negative_test_only',
            'unknown_regular_division_derivatives_kept_as_independent_deviations':True,
            'degree_tested':2,'variables_after_actual_negative_specialization':2,
            'honest_module_generator_rows':9,'jet_vector_dimension':54,
            'all2069_relations_retained_and_all792_new_coefficient_jets_exact':True,
            'Q_enlarged_jet_span_rank':len(basis),'negative_top_columns':negative,
            'exact_formal_top2_all_zero':ev['literal_formal_top2_all_zero'],
            'T3_nonzero_and_depth19_proved':bool(negative),
            'if_no_negative_no_actual_membership_or_T3_zero_inferred':True,
            'no_Tor_basechange_positive_membership_or_blind_order_raise':True,
            'input_sha256':{str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__),source,fullp,jp,originp,lp,hp,up,ur]},
            'evidence_sha256':{'honest9-low2-certificate.json':sha(out/'honest9-low2-certificate.json')},
            'certificate_bytes':len(encoded.encode()),'seconds':time.monotonic()-started,
            'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:report[k] for k in ['all_passed','Q_enlarged_jet_span_rank','negative_top_columns','exact_formal_top2_all_zero','seconds','max_rss_kib']}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
