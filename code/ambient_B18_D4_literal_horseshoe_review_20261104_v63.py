#!/usr/bin/env python3
"""Independent D4 horseshoe audit; extractor is an explicitly external premise."""
import argparse,itertools,json,hashlib,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
BASE='AMBIENT_B18_D4_LITERAL_HORSESHOE_20261104_v63';ID=BASE+'_REVIEW'
def read(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

class Lookup:
    """Independent syntax, origin arithmetic and lookup-only composition."""
    def __init__(self,nodes):
        self.nodes=[tuple(x) for x in nodes];self.idx={n:i for i,n in enumerate(self.nodes)};assert len(self.idx)==len(nodes);self.values=[]
        for i,n in enumerate(self.nodes):
            o=n[0]
            if o=='c':v=F(n[1],n[2]);assert (v.numerator,v.denominator)==n[1:]
            elif o=='v':v=F(0)
            else:
                assert all(isinstance(j,int) and 0<=j<i for j in n[1:])
                if o=='a':v=self.values[n[1]]+self.values[n[2]]
                elif o=='m':v=self.values[n[1]]*self.values[n[2]]
                elif o=='n':v=-self.values[n[1]]
                else:assert o=='i' and self.values[n[1]];v=1/self.values[n[1]]
            self.values.append(v)
        self.zero=self.const(0);self.one=self.const(1)
    def const(self,q):
        q=F(q);return self.idx['c',q.numerator,q.denominator]
    def var(self,name):return self.idx['v',name]
    def neg(self,a):
        n=self.nodes[a]
        if a==self.zero:return a
        if n[0]=='c':return self.const(-self.values[a])
        if n[0]=='n':return n[1]
        return self.idx['n',a]
    def add(self,a,b):
        if a==self.zero:return b
        if b==self.zero:return a
        if self.nodes[a]==('n',b) or self.nodes[b]==('n',a):return self.zero
        if self.nodes[a][0]==self.nodes[b][0]=='c':return self.const(self.values[a]+self.values[b])
        return self.idx[('a',*sorted([a,b]))]
    def mul(self,a,b):
        if a==self.zero or b==self.zero:return self.zero
        if a==self.one:return b
        if b==self.one:return a
        if self.nodes[a][0]==self.nodes[b][0]=='c':return self.const(self.values[a]*self.values[b])
        return self.idx[('m',*sorted([a,b]))]
    def inv(self,a):
        assert self.values[a]
        return self.const(1/self.values[a]) if self.nodes[a][0]=='c' else self.idx['i',a]
    def total(self,roots):
        out=self.zero
        for i in roots:out=self.add(out,i)
        return out
    def grad(self,root,variables):
        # Reverse-mode exact first derivative at the pointed origin. Zero
        # sensitivities are discarded; no source or higher jet is solved.
        adj={root:F(1)};answer={}
        def put(i,q):
            if q:
                q=adj.get(i,F(0))+q
                if q:adj[i]=q
                else:adj.pop(i,None)
        for i in range(root,-1,-1):
            v=adj.pop(i,F(0))
            if not v:continue
            n=self.nodes[i];o=n[0]
            if o=='v':
                if i in variables:answer[i]=answer.get(i,F(0))+v
            elif o=='a':put(n[1],v);put(n[2],v)
            elif o=='m':put(n[1],v*self.values[n[2]]);put(n[2],v*self.values[n[1]])
            elif o=='n':put(n[1],-v)
            elif o=='i':put(n[1],-v/self.values[n[1]]**2)
            else:assert o=='c'
        return [answer.get(i,F(0)) for i in variables]

def det(A):
    A=[row[:] for row in A];v=F(1)
    for j in range(len(A)):
        i=next(i for i in range(j,len(A)) if A[i][j])
        if i!=j:A[i],A[j]=A[j],A[i];v=-v
        pivot=A[j][j];v*=pivot
        for i in range(j+1,len(A)):
            if not A[i][j]:continue
            co=A[i][j]/pivot
            for k in range(j+1,len(A)):A[i][k]-=co*A[j][k]
            A[i][j]=F(0)
    return v
def add(a,b):
    c=dict(a)
    for m,v in b.items():c[m]=c.get(m,0)+v
    return {m:v for m,v in c.items() if v}
def scale(a,q):return {m:v*q for m,v in a.items() if v*q}
def mul(a,b):
    c={}
    for m,v in a.items():
        for n,w in b.items():e=tuple(sorted(m+n));c[e]=c.get(e,0)+v*w
    return {m:v for m,v in c.items() if v}
def var(i):return {(i,):1}
def total(polys):
    q={}
    for p in polys:q=add(q,p)
    return q
def column(col):return {int(i):{tuple(mon):int(co) for co,mon in entries if co} for i,entries in col}
def compose(cols,col):
    out={}
    for j,p in col.items():
        for i,q in cols[j].items():out[i]=add(out.get(i,{}),mul(p,q))
    return {i:p for i,p in out.items() if p}
def same(a,b):assert {int(k):v for k,v in a.items()}==b

def main(out):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False)
    folder=ROOT/'reports'/BASE;report=read(folder/'report.json');assert report['all_passed'];hashes={}
    for name,h in report['input_sha256'].items():assert sha(ROOT/name)==h;hashes[name]=h
    for name,h in report['evidence_sha256'].items():assert sha(folder/name)==h;hashes[str((folder/name).relative_to(ROOT))]=h
    data=read(folder/'actual-D4-horseshoe-DAG.json');np=ROOT/'reports/AMBIENT_B18_N_LITERAL_MODULE_CONE_20261104_v63/actual-N-literal-cone-DAG.json';nd=read(np)
    dp=ROOT/'reports/AMBIENT_B18_NIL_DERIVATION_EXTRACTOR_20261104_v63/nil-delta-extractor-DAG.json';dd=read(dp)
    gate=ROOT/'reports/AMBIENT_B18_N_LITERAL_MODULE_CONE_20261104_v63_REVIEW/report.json';assert read(gate)['all_passed'];hashes[str(gate.relative_to(ROOT))]=sha(gate)
    extractor_gate=ROOT/'reports/AMBIENT_B18_NIL_DERIVATION_EXTRACTOR_20261104_v63_REVIEW/report.json'
    eg=read(extractor_gate);assert eg['all_passed'] and eg['input_sha256'][str(dp.relative_to(ROOT))]==sha(dp)
    for name,h in eg['input_sha256'].items():assert sha(ROOT/name)==h
    for name,h in eg['evidence_sha256'].items():assert sha(extractor_gate.parent/name)==h
    hashes[str(extractor_gate.relative_to(ROOT))]=sha(extractor_gate)
    assert data['DAG_scalar_nodes'][:len(nd['DAG_scalar_nodes'])]==nd['DAG_scalar_nodes'];D=Lookup(data['DAG_scalar_nodes'])
    mapping={};variables={};Z=dd['curve_parameter_Z_node'];z=nd['z1_through_z8_actual_odd_coefficient_nodes'];vv=nd['actual_v_node']
    for i,n in enumerate(dd['DAG_scalar_nodes']):
        o=n[0]
        if i==Z:assert o=='v';q=z[0]
        elif o=='v':q=D.var('D4_EXTRACTOR_'+n[1]);variables[i]=q
        elif o=='c':q=D.const(F(n[1],n[2]))
        elif o=='a':q=D.add(mapping[n[1]],mapping[n[2]])
        elif o=='m':q=D.mul(mapping[n[1]],mapping[n[2]])
        elif o=='n':q=D.neg(mapping[n[1]])
        else:assert o=='i';q=D.inv(mapping[n[1]])
        mapping[i]=q
    assert data['extractor_variable_node_map']==[[i,j] for i,j in sorted(variables.items())]
    cu=[mapping[i] for i in dd['curve77_unknowns_group17_free4_N56']];du=[mapping[i] for i in dd['linear_dual77_centered_unknowns']]
    assert len(cu)==len(du)==77 and not set(cu)&set(du) and not (set(cu)|set(du))&set(nd['all56_coefficient_unknown_nodes'])
    ce=[mapping[i] for i in dd['curve77_equations']];de=[mapping[i] for i in dd['linear_dual77_equations']]
    assert cu==data['cloned_curve77_unknowns'] and du==data['cloned_linear_dual77_centered_unknowns'];assert ce==data['cloned_curve77_equations'] and de==data['cloned_linear_dual77_equations']
    J=[[F(x) for x in r] for r in dd['curve77_origin_Jacobian']]
    JC=[D.grad(i,cu) for i in ce];JD=[D.grad(i,du) for i in de];assert JC==JD==J
    assert det(JC)==F(-47775744000)==F(data['curve77_unit_determinant'])
    assert all(D.values[i]==0 for i in ce+de)
    for k in ['same_actual_R0_coordinate_nodes_by_row'] :assert data[k]==nd[k]
    assert data['N56_unit_unknowns']==nd['all56_coefficient_unknown_nodes'] and data['N56_unit_equations']==nd['all56_actual_coefficient_Hensel_equations']
    assert data['actual_z1_through_z8_nodes']==z and data['actual_v_node']==vv and data['actual_nu_node']==nd['nu_actual_node']
    c=[mapping[i] for i in dd['positive_nil_residue_coefficients_c36']];psi=[mapping[i] for i in dd['nilcurve_graph_ideal_P8_right_sides_psi']]
    assert c==data['positive_c36_actual_nodes'] and psi==data['nilcurve_psi8_actual_nodes']
    assert len(c)==36 and len(psi)==8 and all(D.values[i]==0 for i in c+psi)
    P=[D.add(a,D.neg(b)) for a,b in zip(z[1:]+[vv],psi)];assert P==data['nil_annihilator8_nodes']
    n1=[{int(i):v for i,v in col.items()} for col in nd['R0_sparse_d1_8_by36_columns']];n2=[{int(i):v for i,v in col.items()} for col in nd['R0_sparse_d2_36_by84_columns']]
    U1=list(map(column,nd['d1_universal_sparse_polynomial_columns']));U2=list(map(column,nd['d2_universal_sparse_polynomial_columns']));linear=[]
    for col in U2:
        outcol={}
        for alpha,p in col.items():
            assert all(len(mon)==1 and mon[0] in [*range(8),9] for mon in p)
            outcol[alpha]={mon[0]:co for mon,co in p.items()}
        linear.append(outcol)
    d1=[{0:p} for p in P]+[{0:D.neg(c[a]),**{i+1:q for i,q in col.items()}} for a,col in enumerate(n1)]
    for actual,expected in zip(data['D4_d1_9_by44_columns'],d1):same(actual,expected)
    d2=[{j:P[i],i:D.neg(P[j])} for i,j in itertools.combinations(range(8),2)];obstructions=[];off=[];curve={0:z[0],**{j:psi[j-1] for j in range(1,8)},9:psi[7]}
    for beta,col in enumerate(n2):
        h=[D.total(D.mul(D.const(co.get(j,0)),c[a]) for a,co in linear[beta].items() if co.get(j,0)) for j in [*range(1,8),9]]
        obs=D.total(D.mul(c[a],D.total(D.mul(D.const(coef),curve[j]) for j,coef in co.items())) for a,co in linear[beta].items());obstructions.append(obs);off.append(h)
        d2.append({**{i:q for i,q in enumerate(h)},**{8+a:q for a,q in col.items()}})
    assert len(d1)==44 and len(d2)==112
    for actual,expected in zip(data['D4_d2_44_by112_columns'],d2):same(actual,expected)
    assert off==data['D2_offdiagonal8_by84_columns'] and obstructions==data['all84_curve_cocycle_obstruction_nodes']
    # Independent full polynomial matrix identity, not an origin substitution.
    zp=[var(i) for i in range(8)];vp=var(9);cp=[var(10+i) for i in range(36)];psip=[var(46+i) for i in range(8)];pp=[add(a,scale(b,-1)) for a,b in zip(zp[1:]+[vp],psip)]
    Q1=[{0:p} for p in pp]+[{0:scale(cp[a],-1),**{i+1:p for i,p in col.items()}} for a,col in enumerate(U1)]
    Q2=[{j:pp[i],i:scale(pp[j],-1)} for i,j in itertools.combinations(range(8),2)];curvep={0:zp[0],**{j:psip[j-1] for j in range(1,8)},9:psip[7]};obsps=[]
    for beta,col in enumerate(U2):
        h=[total(scale(cp[a],co.get(j,0)) for a,co in linear[beta].items()) for j in [*range(1,8),9]]
        obs=total(mul(cp[a],total(scale(curvep[j],coef) for j,coef in co.items())) for a,co in linear[beta].items());obsps.append(obs)
        Q2.append({**{i:p for i,p in enumerate(h) if p},**{8+a:p for a,p in col.items()}})
    for col in Q2[:28]:assert not compose(Q1,col)
    for b,col in enumerate(Q2[28:]):assert compose(Q1,col)==({0:scale(obsps[b],-1)} if obsps[b] else {})
    assert data['D4_generators_F_then_N8'][0]['actual_target_function']=='F_old=2*aa+v0+z25'
    assert data['D4_generators_F_then_N8'][1:]==[{**q,'index':q['index']+1} for q in nd['eight_N_module_generators']]
    proof=dict(all139490_scalar_nodes_origin_syntax_checked=len(D.nodes)==139490,all83218_extractor_nodes_composed=len(mapping)==83218,cloned77_and_dual77_disjoint_from_ambient_N56=True,both_cloned77_QQ_origin_Jacobians_exact=77,unit_determinant='-47775744000',same_actual_R0_and_even_subtracted_N8=True,all44_and112_columns_source_bound=True,all_N84_entries_linear_z_v_no_nu=True,all84_offdiagonals_constant_linear_combinations_c36=True,full_universal_polynomial_chain_identities=True,actual_obstruction_vanishing_uses_external_extractor_and_free_nil_premise=True,no_unexpanded_c36_or_psi_series_inferred_zero=True,extractor_independent_audit_owned_by_root=True,extractor_primary_not_independently_reaudited_here=True,extractor_independent_gate_passed=True,scope='Actual horseshoe9/44/112 exact through these terms using the separately independently accepted extractor; no honest9/2069 lifts, minimal beta2, T0..2 or depth conclusion')
    (out/'independent-proof.json').write_text(json.dumps(proof,indent=2)+'\n')
    for p in [Path(__file__),ROOT/f'research/notes/{ID}.md',folder/'report.json']:hashes[str(p.relative_to(ROOT))]=sha(p)
    result=dict(id=ID,all_passed=True,horseshoe_independent_construction_and_conditional_exactness_passed=True,extractor_is_external_primary_premise_not_independently_reviewed_here=True,input_sha256=hashes,evidence_sha256={'independent-proof.json':sha(out/'independent-proof.json')},seconds=time.monotonic()-start,max_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(dict(all_passed=True,seconds=result['seconds'],rss=result['max_rss_kib'],report_sha256=sha(out/'report.json'))))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
