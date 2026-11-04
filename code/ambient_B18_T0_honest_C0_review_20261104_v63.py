#!/usr/bin/env python3
"""Independent fixed certificate review of honest C0, all612 lifts and map."""
import argparse
from fractions import Fraction as F
import hashlib,json
from pathlib import Path
import resource,signal,time

ROOT=Path(__file__).resolve().parents[1];BASE='AMBIENT_B18_T0_HONEST_C0_20261104_v63';ID=BASE+'_REVIEW'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

class Ops:
    def __init__(self,nodes):
        self.nodes=nodes;self.lookup={tuple(n):i for i,n in enumerate(nodes)}
        # Frozen parent contains duplicate display names/constants; node IDs,
        # rather than variable labels, specify its distinct implicit systems.
        self.zero=0;self.one=1
        assert nodes[0]==['c',0,1] and nodes[1]==['c',1,1]
    def c(self,q):
        q=F(q);return self.lookup['c',q.numerator,q.denominator]
    def n(self,a):
        n=self.nodes[a]
        if n[0]=='c':return self.c(-F(n[1],n[2]))
        if n[0]=='n':return n[1]
        return self.lookup['n',a]
    def a(self,a,b):
        if a==self.zero:return b
        if b==self.zero:return a
        if self.nodes[a]==['n',b] or self.nodes[b]==['n',a]:return self.zero
        if self.nodes[a][0]==self.nodes[b][0]=='c':return self.c(F(*self.nodes[a][1:])+F(*self.nodes[b][1:]))
        return self.lookup[('a',*sorted((a,b)))]
    def m(self,a,b):
        if self.zero in [a,b]:return self.zero
        if a==self.one:return b
        if b==self.one:return a
        if self.nodes[a][0]==self.nodes[b][0]=='c':return self.c(F(*self.nodes[a][1:])*F(*self.nodes[b][1:]))
        return self.lookup[('m',*sorted((a,b)))]
    def inv(self,a):
        if self.nodes[a][0]=='c':return self.c(1/F(*self.nodes[a][1:]))
        return self.lookup['i',a]
    def total(self,vs):
        z=self.zero
        for v in vs:z=self.a(z,v)
        return z

def main(out):
    start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
    def timeout(*_):raise TimeoutError('C0 independent fixed review60s cap')
    signal.signal(signal.SIGALRM,timeout);signal.alarm(60);out.mkdir(parents=True,exist_ok=False)
    folder=ROOT/'reports'/BASE;rp=folder/'report.json';primary=read(rp);checked={}
    for key in ['input_sha256','evidence_sha256']:
        for n,h in primary[key].items():
            p=ROOT/n if key=='input_sha256' else folder/n;assert sha(p)==h,n;checked[str(p.relative_to(ROOT))]=h
    source=read(ROOT/'reports/AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63/actual-source-coefficients.json')
    full=read(ROOT/'reports/AMBIENT_B18_FULL_RELATIVE_INPUT_20261104_v63/full-relative186-DAG.json')
    honest=read(ROOT/'reports/AMBIENT_B18_T3_HONEST10_PRESENTATION_20261104_v63/honest10-presentation-DAG.json')
    unit=read(ROOT/'reports/AMBIENT_B18_T3_HONEST10_UNIT_REDUCTION_20261104_v63/unit-reduced-presentation-DAG.json')
    first=read(ROOT/'reports/AMBIENT_B18_T3_FULL9_FIRSTJET_20261104_v63/full9-firstjet-certificate.json')
    joint=read(ROOT/'reports/AMBIENT_B18_T3_JOINT_CURVE_RESTRICTIONS_20261104_v63/joint-degree3-certificate.json')
    origin=read(ROOT/'reports/AMBIENT_B18_T3_DOUBLE_ANNIHILATOR_PREFLIGHT_20261104_v63/exact-residues-and-two-ideals.json')
    certpath=folder/'honest-C0-D0-D4-DAG.json';c=read(certpath);nodes=c['DAG_scalar_nodes'];D=Ops(nodes)
    assert nodes[:len(unit['DAG_scalar_nodes'])]==unit['DAG_scalar_nodes']
    old=source['DAG_scalar_nodes'];ca=source['families']['CA'];gens=ca['actual_zero_odd_I0_generators']+ca['actual_positive_odd_J0_generators']
    Y=[gens[i][1] for i in [0,1,2,3,4,5,6,8,9]];f63=gens[7][1]
    rows=[q['raw_row'] for q in full['all99_CB_kernel_and87_graph_generators_evaluated_on_CA'] if q['kind']=='graph_CB_minus_CA' and q['weight']==0]
    assert rows==[64,68,73,74]
    B=[ca['complete_normalized_physical_w_coefficients'][0][r][1] for r in rows];functions=[f63]+B
    def ancestors(roots):
        todo=list(roots);seen=set()
        while todo:
            i=todo.pop()
            if i in seen:continue
            seen.add(i);n=old[i]
            if n[0] in ['a','m']:todo+=n[1:]
            elif n[0] in ['n','i']:todo.append(n[1])
            else:assert n[0] in ['c','v']
        return sorted(seen)
    eq48=ca['all48_zero_IFT_quadratic_coefficient_equations'];needed=ancestors(eq48+Y+functions)
    prefixes=[]
    for k,rec in enumerate(c['all7_pure_CA_tail_prefix57_systems']):
        assert rec['retained_zero_prefix_length']==k and rec['unit57_determinant']==first['both57_source2']['CA']['determinant']
        m=dict(rec['bindings'])
        for i,j in m.items():assert nodes[j]==['v',f'C0_CA_TAIL_PREFIX{k}_{old[i][1]}']
        for i in needed:
            n=old[i];op=n[0]
            if op=='v':assert i in m
            elif op=='c':m[i]=D.c(F(n[1],n[2]))
            elif op=='a':m[i]=D.a(m[n[1]],m[n[2]])
            elif op=='m':m[i]=D.m(m[n[1]],m[n[2]])
            elif op=='n':m[i]=D.n(m[n[1]])
            else:m[i]=D.inv(m[n[1]])
        eq=[m[i] for i in eq48]+[D.a(m[i],D.n(Y[j] if j<k or j>=7 else D.zero)) for j,i in enumerate(Y)]
        assert eq==rec['all57_unit_equations'];prefixes.append(m)
    assert len(prefixes)==7
    rorig={};coeff=[]
    for root,tel in zip(functions,c['all35_regular_coordinate_divisions']):
        assert tel['original_pure_CA_source_function']==root and tel['restriction_r_s_node']==prefixes[0][root]
        previous=prefixes[0][root];hs=[]
        for i,r in enumerate(tel['all7_regular_coordinate_divisions']):
            current=prefixes[i+1][root] if i<6 else root
            num=D.a(current,D.n(previous));value=F(first['both57_source2']['CA']['pure_function2jets'][str(root)][1+i])
            assert nodes[r['node']]==['r',num,Y[i],f'C0_source{root}:zero{i}']
            assert (r['numerator'],r['coordinate'],r['previous'],r['current'])==(num,Y[i],previous,current)
            assert F(r['origin'])==value;rorig[r['node']]=value;hs.append(r['node']);previous=current
        assert len(hs)==7;coeff.append(hs)
    # Own ten-dimensional truncated polynomial arithmetic, no Newton solve.
    mons=[tuple(m) for m in joint['monomials']];pos={m:i for i,m in enumerate(mons)};zero=(F(0),)*10
    products=[(i,j,pos[(a[0]+b[0],a[1]+b[1])]) for i,a in enumerate(mons) for j,b in enumerate(mons) if sum(a)+sum(b)<=3]
    def mul(a,b):
        z=[F(0)]*10
        for i,j,k in products:
            if a[i] and b[j]:z[k]+=a[i]*b[j]
        return tuple(z)
    jets={i:tuple(map(F,v)) for i,v in joint['all105_unknown_jets']}
    for i in ancestors(functions):
        n=old[i];op=n[0]
        if op=='v':assert i in jets
        elif op=='c':jets[i]=(F(n[1],n[2]),)+zero[1:]
        elif op=='a':jets[i]=tuple(a+b for a,b in zip(jets[n[1]],jets[n[2]]))
        elif op=='n':jets[i]=tuple(-a for a in jets[n[1]])
        elif op=='m':jets[i]=mul(jets[n[1]],jets[n[2]])
        else:
            a=jets[n[1]];assert a[0];h=tuple(-v/a[0] if k else F(0) for k,v in enumerate(a))
            h2=mul(h,h);h3=mul(h2,h);jets[i]=tuple((F(k==0)+h[k]+h2[k]+h3[k])/a[0] for k in range(10))
    gjet=jets[f63];assert list(map(str,gjet))==c['cached_actual_CA_gA_all_degree3_coefficients']
    assert all(not a for a,m in zip(gjet,mons) if sum(m)<3) and any(gjet)
    assert all(jets[b]==zero for b in B)
    gA=prefixes[0][f63];H=[]
    for row,b,rec in zip(rows,B,c['all4_known_regular_gA_divisions']):
        assert nodes[rec['node']]==['g',prefixes[0][b],gA,f'C0_zero_graph{row}:known_regular_zero_image']
        assert F(rec['origin'])==0;H.append(rec['node']);rorig[rec['node']]=F(0)
    basis=Y[:7]+[gA];assert basis==c['actual_C0_basis8_scalar_nodes']
    prim={D.zero:(D.zero,)*8,f63:tuple(coeff[0]+[D.one])}
    for i,z in enumerate(Y[:7]):v=[D.zero]*8;v[i]=D.one;prim[z]=tuple(v)
    for b,h,co in zip(B,H,coeff[1:]):prim[b]=tuple(co+[h])
    rules={int(i):r for i,r in c['C0_odd_membership_propagation_rules'].items()};mem={}
    def lift(i):
        if i in mem:return mem[i]
        rule=rules[i]
        if i in prim:v=prim[i]
        elif rule[0]=='negative_of_certified':
            assert nodes[rule[1]]==['n',i];v=tuple(D.n(a) for a in lift(rule[1]))
        elif rule[0]=='neg':
            assert nodes[i]==['n',rule[1]];v=tuple(D.n(a) for a in lift(rule[1]))
        elif rule[0]=='sum':
            assert nodes[i]==['a',rule[1],rule[2]];v=tuple(D.a(a,b) for a,b in zip(lift(rule[1]),lift(rule[2])))
        else:
            assert rule[0] in ['product_known_left','product_known_right'] and nodes[i][0]=='m'
            assert set(nodes[i][1:])==set(rule[1:])
            known,scalar=(rule[1],rule[2]) if rule[0]=='product_known_left' else (rule[2],rule[1])
            v=tuple(D.m(scalar,a) for a in lift(known))
        mem[i]=v;return v
    q=full['actual_Qn_R0_squared_generator_matrices']['0'];assert len(q)==612
    lifts=[[a]+list(lift(b)) for a,b in q];assert lifts==c['all612_Q0_lifted_columns']
    koszul=[]
    for i in range(8):
        for j in range(i+1,8):v=[D.zero]*9;v[1+i]=D.n(basis[j]);v[1+j]=basis[i];koszul.append(v)
    assert koszul==c['all28_C0_Koszul_kernel_columns'] and koszul+lifts==c['honest_D0_relation_columns']
    H63=next(r['membership_vector9'] for r in honest['all_pure_source_Hadamard_membership_telescopes'] if r['original_source_function']==f63)
    mapping=[[D.one]+[D.zero]*9]
    for i in range(7):v=[D.zero]*10;v[i+1]=D.one;mapping.append(v)
    v=[D.zero]+H63[:]
    for i in range(7):v[i+1]=D.a(v[i+1],D.n(coeff[0][i]))
    mapping.append(v);assert mapping==c['C0_free9_to_H_free10_columns']
    P=unit['exact_projection_matrix_remaining_by_original10']
    image=[[D.total(D.m(P[i][j],v[j]) for j in range(10)) for i in range(9)] for v in mapping]
    assert image==c['C0_free9_to_D4_honest9_columns']
    vals=list(map(F,origin['all339669_scalar_origin_values_evaluated']))
    rorig.update({r['node']:F(r['origin']) for t in honest['all_pure_source_Hadamard_membership_telescopes'] for r in t['all9_regular_coordinate_divisions']})
    for i in range(len(vals),len(nodes)):
        n=nodes[i];op=n[0]
        if op=='c':v=F(n[1],n[2])
        elif op=='v':v=F(0)
        elif op=='a':v=vals[n[1]]+vals[n[2]]
        elif op=='m':v=vals[n[1]]*vals[n[2]]
        elif op=='n':v=-vals[n[1]]
        elif op=='i':assert vals[n[1]];v=1/vals[n[1]]
        else:assert op in ['r','g'] and i in rorig;v=rorig[i]
        vals.append(v)
    assert all(vals[i]==0 for col in koszul+lifts for i in col)
    M=[[vals[i] for i in col] for col in image];rank=0
    for j in range(9):
        k=next((k for k in range(rank,9) if M[k][j]),None)
        if k is None:continue
        M[k],M[rank]=M[rank],M[k];p=M[rank][j];M[rank]=[v/p for v in M[rank]]
        for k in range(rank+1,9):
            q=M[k][j];M[k]=[v-q*w for v,w in zip(M[k],M[rank])]
        rank+=1
    assert rank==7
    ep=out/'independent-C0-membership-and-origin-proof.json'
    ep.write_text(json.dumps({'all7_fixed57_prefix_inputs_rebuilt':True,'all35_regular_coordinate_quotients_checked':True,
        'cached_joint105_own_degree3_arithmetic_gA':list(map(str,gjet)),
        'all4_new_regular_g_quotient_origins':[0]*4,'all612_ideal_memberships_checked':True,
        'complete28_Koszul_checked':True,'D0_minimal_generators':9,'C0_to_D4_origin_rank':rank,
        'T0_kernel_or_socle_inferred_from_origin_rank':False},indent=2)+'\n')
    paths=[Path(__file__),ROOT/'research/notes'/f'{ID}.md',rp,certpath]
    report={'id':ID,'all_passed':True,'all612_actual_Q0_lifts_and_complete28_Koszul_rebuilt':True,
        'all7_prefix57_and35_coordinate_quotients_rebuilt':True,'all4_known_regular_g_origins_independently_zero':True,
        'actual_D0_free9_relations640_and_D4_map_checked':True,'T0_socle_depth_or_kernel_computed':False,
        'input_sha256':{str(p.relative_to(ROOT)):sha(p) for p in paths},'checked_dependency_sha256':checked,
        'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');signal.alarm(0)
    print(json.dumps({k:v for k,v in report.items() if 'sha256' not in k}))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
