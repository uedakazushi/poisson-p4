#!/usr/bin/env python3
"""ONE full nine-variable honest9 first-jet negative test, retaining zero7.

Two pure57 source inverses to degree2 determine coordinate-quotient firstjets.
The four g derivatives remain independent.  No higher honest/source jets.
"""
import argparse
from fractions import Fraction as F
import hashlib,itertools as it,json
from pathlib import Path
import resource,signal,time
from ambient_B18_top_graph_telescoping_20261104_v63 import forward,ancestors,linear_solve
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_T3_FULL9_FIRSTJET_20261104_v63'
M=[()]+[(i,) for i in range(9)]+list(it.combinations_with_replacement(range(9),2));MP={a:i for i,a in enumerate(M)}
CV=[[(i,j) for i,a in enumerate(M) for j,b in enumerate(M) if tuple(sorted(a+b))==m] for m in M];ZERO=(F(0),)*55
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def enc(A):return [[str(v) for v in r] for r in A]
def jm(a,b):return tuple(sum((a[i]*b[j] for i,j in t if a[i] and b[j]),F(0)) for t in CV)
def ev(nodes,needed,bind):
    V={}
    for i in needed:
        n=nodes[i];o=n[0]
        if o=='c':v=(F(n[1],n[2]),)+(F(0),)*54
        elif o=='v':v=bind[i]
        elif o=='a':v=tuple(a+b for a,b in zip(V[n[1]],V[n[2]]))
        elif o=='n':v=tuple(-a for a in V[n[1]])
        elif o=='m':v=jm(V[n[1]],V[n[2]])
        elif o=='i':
            a=V[n[1]];assert a[0];h=(F(0),)+a[1:];hh=jm(h,h);v=(1/a[0],)+tuple(-a[j]/a[0]**2+hh[j]/a[0]**3 for j in range(1,55))
        else:raise ValueError(o)
        V[i]=v
    return V
def main(out):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));out.mkdir(parents=True,exist_ok=False)
    sp=ROOT/'reports/AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63/actual-source-coefficients.json'
    fp=ROOT/'reports/AMBIENT_B18_FULL_RELATIVE_INPUT_20261104_v63/full-relative186-DAG.json'
    hp=ROOT/'reports/AMBIENT_B18_T3_HONEST10_PRESENTATION_20261104_v63/honest10-presentation-DAG.json'
    up=ROOT/'reports/AMBIENT_B18_T3_HONEST10_UNIT_REDUCTION_20261104_v63/unit-reduced-presentation-DAG.json'
    op=ROOT/'reports/AMBIENT_B18_T3_DOUBLE_ANNIHILATOR_PREFLIGHT_20261104_v63/exact-residues-and-two-ideals.json'
    lp=ROOT/'reports/AMBIENT_B18_T3_DOUBLE_ANNIHILATOR_LENGTH3_20261104_v63/report.json'
    old,full,honest,unit,origin=map(read,[sp,fp,hp,up,op]);assert read(lp)['all_passed']
    oldnodes=old['DAG_scalar_nodes'];N=len(oldnodes);nodes=unit['DAG_scalar_nodes'];var={n[1]:i for i,n in enumerate(oldnodes) if n[0]=='v'}
    xids=[var[f'x{i}'] for i in old['physical_zero_base_rows9']];sel=[0,1,2,3,4,5,6,8,9]
    families={};sols={};source2={};stage='source57'
    # Roots of all pure functions used in the old/new Hadamard coefficients.
    for name in ['CA','CB']:
        f=old['families'][name];z=[i for p in f['coefficient_unknown_pairs'].values() for i in p]+[f['p_quadratic_relation']['s'],f['p_quadratic_relation']['t']]
        unknown=z+xids;gens=f['actual_zero_odd_I0_generators']+f['actual_positive_odd_J0_generators'];Y=[gens[i][1] for i in sel];f63=gens[7][1]
        if name=='CA':roots=[t['original_source_function'] for t in honest['all_pure_source_Hadamard_membership_telescopes']]
        else:
            rows=[q['raw_row'] for q in full['all99_CB_kernel_and87_graph_generators_evaluated_on_CA'] if q['kind']=='graph_CB_minus_CA']
            roots=[f63]+[f['complete_normalized_physical_w_coefficients'][next(q['weight'] for q in full['all99_CB_kernel_and87_graph_generators_evaluated_on_CA'] if q.get('raw_row')==r)][r][1] for r in rows]
        roots=list(dict.fromkeys(roots));eq=f['all48_zero_IFT_quadratic_coefficient_equations']+Y;needed=sorted(ancestors(oldnodes,eq+roots))
        vals,g=forward(oldnodes,set(unknown));assert all(vals[i]==0 for i in eq+roots)
        J=[[g[i].get(v,F(0)) for v in unknown] for i in eq];Ji,det=linear_solve(J,[[F(int(i==j)) for j in range(57)] for i in range(57)])
        sol={v:(F(0),)+tuple(Ji[k][48+j] for j in range(9))+(F(0),)*45 for k,v in enumerate(unknown)}
        V=ev(oldnodes,needed,sol);res=[[V[i][j] for j in range(10,55)] for i in eq]
        corr=[[-sum(Ji[i][k]*res[k][j] for k in range(57)) for j in range(45)] for i in range(57)]
        sol={v:sol[v][:10]+tuple(corr[k]) for k,v in enumerate(unknown)}
        V=ev(oldnodes,needed,sol)
        assert all(V[i][j]==int(k>=48 and M[j]==(k-48,)) for k,i in enumerate(eq) for j in range(55))
        source2[name]={'unknown57':unknown,'equations57':eq,'F9_nodes':Y,'f63':f63,'fixed_J57':enc(J),'fixed_J57_inverse':enc(Ji),'determinant':str(det),
                       'all57_source2jets':[[i,list(map(str,sol[i]))] for i in unknown],'pure_function2jets':{str(i):list(map(str,V[i])) for i in roots},'all57_equations_all55_coefficients_verified':True}
        sols[name]=sol;families[name]={'f':f,'z':z,'Y':Y,'f63':f63,'V':V,'roots':roots}
    # First-order change YB(YA), computed through the common physical x9.
    XA=[[sols['CA'][x][1+j] for j in range(9)] for x in xids]
    XB=[[sols['CB'][x][1+j] for j in range(9)] for x in xids]
    XBi,_=linear_solve(XB,[[F(int(i==j)) for j in range(9)] for i in range(9)])
    adapter=[[sum(XBi[i][k]*XA[k][j] for k in range(9)) for j in range(9)] for i in range(9)]
    _,adet=linear_solve(adapter,[[F(int(i==j)) for j in range(9)] for i in range(9)]);assert adet==F(1,4096)
    binding={}
    for i,v in sols['CA'].items():binding[i]=[F(0)]+list(v[1:10])+[F(0)]*4
    for i,v in sols['CB'].items():
        jet=[F(0)]+[sum(v[1+k]*adapter[k][j] for k in range(9)) for j in range(9)]+[F(0)]*4
        if i in binding:assert binding[i]==jet
        binding[i]=jet
    def plus(a,b):return [u+v for u,v in zip(a,b)]
    def times(a,b):return [a[0]*b[0]]+[a[0]*b[i]+b[0]*a[i] for i in range(1,14)]
    def inverse(a):assert a[0];return [1/a[0]]+[-v/a[0]**2 for v in a[1:]]
    original={}
    for i,n in enumerate(oldnodes):
        o=n[0]
        if o=='c':v=[F(n[1],n[2])]+[F(0)]*13
        elif o=='v':v=binding[i]
        elif o=='a':v=plus(original[n[1]],original[n[2]])
        elif o=='n':v=[-x for x in original[n[1]]]
        elif o=='m':v=times(original[n[1]],original[n[2]])
        elif o=='i':v=inverse(original[n[1]])
        else:raise ValueError(o)
        original[i]=v
    def hadamard(poly,j,tail=False):
        h=[poly[1+j]]+[F(0)]*9
        for k in range(9):
            if k<=j or (tail and k>=7):h[1+k]=poly[MP[tuple(sorted((j,k)))]]
        return h
    known={};newh={}
    for rec in honest['all_pure_source_Hadamard_membership_telescopes']:
        p=families['CA']['V'][rec['original_source_function']]
        for j,q in enumerate(rec['all9_regular_coordinate_divisions']):
            h=hadamard(p,j);assert h[0]==F(q['origin']);known[q['node']]=h+[F(0)]*4;newh[q['node']]=h
    oldh={};cb=families['CB'];Braw={q['raw_row']:cb['f']['complete_normalized_physical_w_coefficients'][q['weight']][q['raw_row']][1] for q in full['all99_CB_kernel_and87_graph_generators_evaluated_on_CA'] if q['kind']=='graph_CB_minus_CA'}
    for rec in full['all_coordinate_divisions']:
        i=rec['node'];label=full['DAG_scalar_nodes'][i][3]
        if label.startswith('g_split_'):h=[F(0)]*10
        else:
            head,coord=label.split(':');j=int(coord[1:])
            root=cb['f63'] if head=='f63' else Braw[int(head.replace('zero_graph','').replace('graph',''))]
            h=hadamard(cb['V'][root],j,coord.startswith('z'))
        assert h[0]==F(rec['origin_value'])
        hc=[h[0]]+[sum(h[1+k]*adapter[k][j] for k in range(9)) for j in range(9)]
        known[i]=hc+[F(0)]*4;oldh[i]=hc
    assert len(newh)==792 and len(oldh)==784
    gnodes=[q['node'] for q in full['all4_known_regular_g_divisions']]
    for j,i in enumerate(gnodes):v=[F(0)]*14;v[10+j]=1;known[i]=v
    known[full['CB_g_restriction_node']]=[F(0)]*14
    A=unit['reduced_D3_relations_after_dropping_only_literal_zero_columns'];T=unit['all12_transformed_top_columns']
    needed=set();stack=[i for col in A+T for i in col]
    while stack:
        i=stack.pop()
        if i in needed:continue
        needed.add(i)
        if i<N or i in known:continue
        n=nodes[i]
        if n[0] in ['a','m']:stack.extend(n[1:])
        elif n[0] in ['n','i']:stack.append(n[1])
        elif n[0]!='c':raise AssertionError(('uncontrolled_node',i,n))
    V={}
    for i in sorted(needed):
        if i<N:v=original[i]
        elif i in known:v=known[i]
        else:
            n=nodes[i];o=n[0]
            if o=='c':v=[F(n[1],n[2])]+[F(0)]*13
            elif o=='n':v=[-x for x in V[n[1]]]
            elif o=='a':v=plus(V[n[1]],V[n[2]])
            elif o=='m':v=times(V[n[1]],V[n[2]])
            elif o=='i':v=inverse(V[n[1]])
            else:raise ValueError(o)
        V[i]=v
    def vectors(col):
        x=[a for node in col for a in V[node][:10]];sens=[[V[node][10+j] for node in col] for j in range(4)];u=[]
        for j,s in enumerate(sens):
            if any(s):
                for k in range(9):v=[F(0)]*90;v[1+k::10]=s;u.append(v)
        return x,sens,u
    basis=[];provenance=[]
    def insert(v,label):
        v=v[:]
        for p,b in basis:
            if v[p]:c=v[p];v=[x-c*y for x,y in zip(v,b)]
        p=next((i for i,x in enumerate(v) if x),None)
        if p is not None:
            c=v[p];basis.append((p,[x/c for x in v]));provenance.append(label)
    qvec=[];qunc=[];qsens=[]
    for j,col in enumerate(A):
        x,s,u=vectors(col);assert all(x[10*i]==0 for i in range(9));qvec.append(x);qsens.append(s);qunc+=u;insert(x,['Q',j])
        for k,v in enumerate(u):insert(v,['Qg_unc',j,k])
    rank=len(basis);results=[];negative=[]
    for j,col in enumerate(T):
        x,s,u=vectors(col);save=list(basis);pr=list(provenance)
        for k,v in enumerate(u):insert(v,['top_g_unc',j,k])
        r=x[:]
        for p,b in basis:
            if r[p]:c=r[p];r=[a-c*z for a,z in zip(r,b)]
        witness=None
        if any(r):
            negative.append(j);k=next(i for i,z in enumerate(r) if z);lam=[F(0)]*90;lam[k]=1
            for p,b in reversed(basis):lam[p]=-sum(b[i]*lam[i] for i in range(p+1,90))
            assert all(sum(a*b for a,b in zip(lam,v))==0 for v in qvec+qunc+u)
            value=sum(a*b for a,b in zip(lam,x));assert value
            witness={'left_annihilator90':list(map(str,lam)),'value_on_top':str(value),'all2069relations_and_g_uncertainties_annihilated':True}
        results.append({'top':j,'known_firstjet90':list(map(str,x)),'four_g_sensitivity9_by4':enc(s),'negative':bool(any(r)),'witness':witness})
        basis[:]=save;provenance[:]=pr
    evidence={'full9_CA_coordinate_order':families['CA']['Y'],'source2_monomials':list(map(list,M)),'both57_source2':source2,
              'YB_to_YA_linear':enc(adapter),'adapter_determinant':str(adet),'all792_CA_h_firstjets':{str(i):list(map(str,v)) for i,v in newh.items()},
              'all784_CB_h_firstjets_in_YA':{str(i):list(map(str,v)) for i,v in oldh.items()},'four_g_derivatives_independent_node_ids':gnodes,
              'all2069_known_firstjet90_columns':enc(qvec),'all2069_four_g_sensitivities':[[list(map(str,r)) for r in s] for s in qsens],
              'all_relation_centers_zero_only_constant_multiples_needed':True,'enlarged_relation_span_rank':rank,
              'relation_span_basis':[[p,list(map(str,b))] for p,b in basis],'basis_provenance':provenance,'all12_top_results':results,
              'negative_top_columns':negative,'four_g_derivatives_drop_from_all_relations':all(not any(r) for s in qsens for r in s),
              'four_g_derivatives_drop_from_all_tops':all(not any(map(F,r)) for q in results for r in q['four_g_sensitivity9_by4'])}
    encoded=json.dumps(evidence,separators=(',',':'))+'\n';assert len(encoded.encode())<32*1024*1024
    (out/'full9-firstjet-certificate.json').write_text(encoded)
    rep={'id':STEM,'all_passed':True,'status':'full9_honest_firstjet_negative_test','enlarged_relation_span_rank':rank,
         'negative_top_columns':negative,'actual_T3_nonzero_proved':bool(negative),'cone_depth19_proved_if_independent_gate_accepted':bool(negative),
         'four_g_derivatives_drop_from_relations':evidence['four_g_derivatives_drop_from_all_relations'],
         'four_g_derivatives_drop_from_tops':evidence['four_g_derivatives_drop_from_all_tops'],
         'no_positive_membership_or_zero_inferred_if_inconclusive':True,'no_source_order_above2_Honest_order_above1_GB_or_new_points':True,
         'input_sha256':{str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__),sp,fp,hp,up,op,lp]},
         'evidence_sha256':{'full9-firstjet-certificate.json':sha(out/'full9-firstjet-certificate.json')},
         'certificate_bytes':len(encoded.encode()),'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    (out/'report.json').write_text(json.dumps(rep,indent=2)+'\n');print(json.dumps({k:rep[k] for k in ['all_passed','enlarged_relation_span_rank','negative_top_columns','four_g_derivatives_drop_from_relations','four_g_derivatives_drop_from_tops','seconds','max_rss_kib']}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
