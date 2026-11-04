#!/usr/bin/env python3
"""Independent pure CA96/N56 composition and regular-coordinate B_N audit."""
import argparse,itertools,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
from ambient_B18_D4_literal_horseshoe_review_20261104_v63 import Lookup,det,read,sha
from ambient_B18_actual_top_coefficient_preflight_review_20261104_v63 import QuadraticReplay
ROOT=Path(__file__).resolve().parents[1]
BASE='AMBIENT_B18_HONEST9_TO_N8_20261104_v63';ID=BASE+'_REVIEW'
def ancestors(nodes,roots):
    need=set();todo=list(roots)
    while todo:
        i=todo.pop()
        if i in need:continue
        need.add(i);n=nodes[i]
        if n[0] not in ['c','v']:
            assert n[0]!='r';todo.extend(n[1:])
    return sorted(need)

class RegularLookup(Lookup):
    def __init__(self,nodes,origins):
        self.nodes=[tuple(n) for n in nodes];self.idx={n:i for i,n in enumerate(self.nodes)};assert len(self.idx)==len(nodes);self.values=[];rdep=[]
        for i,n in enumerate(self.nodes):
            o=n[0];has=False
            if o=='c':v=F(n[1],n[2])
            elif o=='v':v=F(0)
            elif o=='r':
                assert i in origins and n[1]<i and n[2]<i;v=origins[i];has=True
            else:
                assert all(isinstance(j,int) and 0<=j<i for j in n[1:]);has=any(rdep[j] for j in n[1:])
                if o=='a':v=self.values[n[1]]+self.values[n[2]]
                elif o=='m':v=self.values[n[1]]*self.values[n[2]]
                elif o=='n':v=-self.values[n[1]]
                else:assert o=='i' and self.values[n[1]] and not has;v=1/self.values[n[1]]
            self.values.append(v);rdep.append(has)
        self.zero=self.const(0);self.one=self.const(1);self.rdep=rdep
    inverse=Lookup.inv

def grad_add(a,b,ca,cb):
    out={i:v*ca for i,v in a.items() if v*ca}
    for i,v in b.items():
        q=out.get(i,F(0))+v*cb
        if q:out[i]=q
        else:out.pop(i,None)
    return out
def gradients(D,roots,variables):
    G={};variables=set(variables)
    for i in ancestors(D.nodes,roots):
        n=D.nodes[i];o=n[0]
        if o=='c':g={}
        elif o=='v':g={i:F(1)} if i in variables else {}
        elif o=='a':g=grad_add(G[n[1]],G[n[2]],1,1)
        elif o=='m':g=grad_add(G[n[1]],G[n[2]],D.values[n[2]],D.values[n[1]])
        elif o=='n':g={j:-v for j,v in G[n[1]].items()}
        else:assert o=='i';g={j:-v/D.values[n[1]]**2 for j,v in G[n[1]].items()}
        G[i]=g
    return G
def implicit_witness(J):
    """Only a constant QQ linearized witness; no source/Newton calculation."""
    A=[row[:]+[F(int(56<=i<64 and i-56==j)) for j in range(8)] for i,row in enumerate(J)];n=len(A);determinant=F(1)
    for j in range(n):
        k=next(k for k in range(j,n) if A[k][j])
        if k!=j:A[k],A[j]=A[j],A[k];determinant=-determinant
        pivot=A[j][j];determinant*=pivot
        for i in range(j+1,n):
            if not A[i][j]:continue
            co=A[i][j]/pivot;A[i][j]=F(0)
            for k in range(j+1,n+8):
                if A[j][k]:A[i][k]-=co*A[j][k]
    W=[[F(0)]*8 for _ in range(n)]
    for i in range(n-1,-1,-1):
        for j in range(8):W[i][j]=(A[i][n+j]-sum((A[i][k]*W[k][j] for k in range(i+1,n) if A[i][k]),F(0)))/A[i][i]
    for i,row in enumerate(J):
        for j in range(8):assert sum((v*W[k][j] for k,v in enumerate(row) if v),F(0))==int(56<=i<64 and i-56==j)
    return W,determinant

def main(out):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False)
    folder=ROOT/'reports'/BASE;r=read(folder/'report.json');assert r['all_passed'];hashes={}
    for name,h in r['input_sha256'].items():assert sha(ROOT/name)==h;hashes[name]=h
    for name,h in r['evidence_sha256'].items():assert sha(folder/name)==h;hashes[str((folder/name).relative_to(ROOT))]=h
    data=read(folder/'honest9-to-N8-DAG.json');cone=read(ROOT/'reports/AMBIENT_B18_N_LITERAL_MODULE_CONE_20261104_v63/actual-N-literal-cone-DAG.json');N=read(ROOT/'reports/AMBIENT_B18_N_SAME_R0_BINDING_20261104_v63/actual-N-R0-binding-DAG.json')
    old=read(ROOT/'reports/AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63/actual-source-coefficients.json');unit=read(ROOT/'reports/AMBIENT_B18_T3_HONEST10_UNIT_REDUCTION_20261104_v63/unit-reduced-presentation-DAG.json')
    gates=['AMBIENT_B18_N_LITERAL_MODULE_CONE_20261104_v63_REVIEW/report.json','AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63_REVIEW/report.json','AMBIENT_B18_T3_HONEST10_UNIT_REDUCTION_20261104_v63_REVIEW/report.json','AMBIENT_B18_W1_TAIL_BINDING_20261104_v63_REVIEW/report.json','AMBIENT_B18_W1_POLYNOMIAL_LAURENT_SCOPE_20261104_v63.json']
    for name in gates:
        p=ROOT/'reports'/name;g=read(p);assert g.get('all_passed',g.get('passed',False));hashes[str(p.relative_to(ROOT))]=sha(p)
    assert data['DAG_scalar_nodes'][:len(cone['DAG_scalar_nodes'])]==cone['DAG_scalar_nodes']
    origins={q['node']:F(q['origin']) for telescope in data['all64_regular_coordinate_divisions'] for q in telescope['all8_coordinate_divisions']};assert len(origins)==64
    D=RegularLookup(data['DAG_scalar_nodes'],origins);ca=old['families']['CA'];oldnodes=old['DAG_scalar_nodes'];oldvars={n[1]:i for i,n in enumerate(oldnodes) if n[0]=='v'}
    az=[i for pair in ca['coefficient_unknown_pairs'].values() for i in pair]+[ca['p_quadratic_relation']['s'],ca['p_quadratic_relation']['t']];assert len(az)==48 and data['CA_scalar48_unknown_nodes']==az
    AA=[D.var(f'HN_CA_A_{i}') for i in az];BB=[D.var(f'HN_CA_B_{i}') for i in az];unknown96=AA+BB;assert unknown96==data['CA96_quadratic_N_basechange_unknown_nodes']
    assert data['CA96_coefficient_pairs_of_original48']==[[i,[a,b]] for i,a,b in zip(az,AA,BB)]
    sourcegens=ca['actual_zero_odd_I0_generators']+ca['actual_positive_odd_J0_generators'];Yrows=[sourcegens[i][0] for i in [0,1,2,3,4,5,6,8,9]]
    honest=[None if i==0 else Yrows[i-1] for i in unit['remaining_original_free_generator_indices']];assert honest==data['honest9_generator_physical_rows']==[None,18,19,28,29,67,72,46,95]
    even={row:[ca['complete_normalized_physical_w_coefficients'][degree][row][0] for degree in range(4)] for row in honest[1:]}
    xids=[oldvars[f'x{row}'] for row in old['physical_zero_base_rows9']];binding=dict(zip(xids,map(tuple,N['all9_R0_quadratic_p_readback_pairs'])));binding.update({i:(a,b) for i,a,b in zip(az,AA,BB)})
    S,T=N['quadratic_trace_S_T_nodes'];Q=QuadraticReplay(D,S,T);eq48=ca['all48_zero_IFT_quadratic_coefficient_equations'];need=ancestors(oldnodes,eq48+[v for roots in even.values() for v in roots]);mp={}
    for i in need:
        n=oldnodes[i];o=n[0]
        if o=='v':assert i in binding;q=binding[i]
        elif o=='c':q=(D.const(F(n[1],n[2])),D.zero)
        elif o=='a':q=Q.add(mp[n[1]],mp[n[2]])
        elif o=='m':q=Q.mul(mp[n[1]],mp[n[2]])
        elif o=='n':q=Q.neg(mp[n[1]])
        else:assert o=='i';q=Q.inverse(mp[n[1]])
        mp[i]=q
    eq96=[q for i in eq48 for q in mp[i]];assert eq96==data['all96_unit_basechange_equations'];images=[]
    for row,record in zip(honest[1:],data['all8_honest_nonunit_generators_on_N']):
        ce=(D.zero,D.zero)
        for i in even[row]:ce=Q.add(ce,mp[i])
        p=Q.add(tuple(cone['all150_actual_quadratic_p_pairs'][row]),Q.neg(ce));a=D.add(p[0],D.mul(D.const(F(1,2)),D.mul(S,p[1])));b=p[1]
        assert record==dict(honest_physical_row=row,CA_even_scalar_basechanged_to_N_p_pair=list(ce),N_tracefree_even_odd=[a,b]);assert not D.values[a] and not D.values[b];images.append((a,b))
    EN=cone['all56_actual_coefficient_Hensel_equations'];z=cone['z1_through_z8_actual_odd_coefficient_nodes'];base=cone['actual_R8_parameter_nodes8'];variables=cone['all56_coefficient_unknown_nodes']+base+unknown96;eq160=EN+z+eq96
    assert len(variables)==len(eq160)==160 and variables==data['pure_unit160_variable_order'] and eq160==data['pure_unit160_equation_order'];assert data['actual_N_z8_coordinate_nodes']==z and data['actual_R0_v_node']==cone['actual_v_node']
    assert all(not D.values[i] for i in eq160)
    roots=eq160+[b for a,b in images];G=gradients(D,roots,variables);J=[[G[q].get(i,F(0)) for i in variables] for q in eq160];J96=[[G[q].get(i,F(0)) for i in unknown96] for q in eq96]
    assert J==[[F(v) for v in row] for row in data['unit160_Jacobian']] and J96==[[F(v) for v in row] for row in data['unit96_Jacobian']]
    d96=det(J96);assert d96==F(data['unit96_determinant'])==F(ca['unit48_block_determinant_even24_then_odd24_vars_a23_t_b23_s'])**2
    W,d160=implicit_witness(J);assert d160==F(data['unit160_determinant'])==F(cone['unit56_minor'])*F(cone['unit8_odd_minor'])*d96
    derivatives=[[sum((v*W[variables.index(i)][j] for i,v in G[b].items()),F(0)) for j in range(8)] for a,b in images]
    # All7 prefix inverse systems are cloned as functions, with only coordinate
    # targets changed. No implicit solve or finite jet is recomputed here.
    need2=ancestors(D.nodes,roots);assert {i for i in need2 if D.nodes[i][0]=='v'}<=set(variables);prefix=[]
    for level,system in enumerate(data['all7_prefix160_systems'],1):
        assert system['level']==level and F(system['unit160_determinant'])==d160;bind={int(i):j for i,j in system['all160_variable_bindings']};assert set(bind)==set(variables)
        assert len(set(bind.values()))==160 and not set(bind.values())&set(variables)
        for i,j in bind.items():assert D.nodes[j]==('v',f'HN_PREFIX{level}_{D.nodes[i][1]}')
        m={}
        for i in need2:
            n=D.nodes[i];o=n[0]
            if o=='v':q=bind[i]
            elif o=='c':q=D.const(F(n[1],n[2]))
            elif o=='a':q=D.add(m[n[1]],m[n[2]])
            elif o=='m':q=D.mul(m[n[1]],m[n[2]])
            elif o=='n':q=D.neg(m[n[1]])
            else:assert o=='i';q=D.inv(m[n[1]])
            m[i]=q
        expected=[m[i] for i in EN]+[D.add(m[i],D.neg(z[j] if j<level else D.zero)) for j,i in enumerate(z)]+[m[i] for i in eq96]
        assert expected==system['all160_unit_equations'];prefix.append({b:m[b] for a,b in images})
    v=cone['actual_v_node'];columns=[[D.one]+[D.zero]*7]
    for j,((a,b),telescope) in enumerate(zip(images,data['all64_regular_coordinate_divisions'])):
        assert telescope['honest_index']==j+1 and telescope['odd_source_function']==b;previous=D.zero;coeffs=[]
        for i,item in enumerate(telescope['all8_coordinate_divisions']):
            current=prefix[i][b] if i<7 else b;numerator=D.add(current,D.neg(previous));node=item['node']
            assert D.nodes[node]==('r',numerator,z[i],f'HN_honest{j+1}:z{i+1}')
            assert item['numerator']==numerator and item['coordinate']==z[i] and item['current']==current and item['previous']==previous
            assert F(item['origin'])==derivatives[j][i]==D.values[node];coeffs.append(node);previous=current
        assert coeffs==telescope['membership_coefficients8'] and previous==b
        columns.append([D.add(a,D.mul(coeffs[0],v))]+coeffs[1:])
    assert columns==data['B_N_8_by9_columns'];constant=[[D.values[q] for q in col] for col in columns];assert constant==[[F(v) for v in col] for col in data['B_N_origin_matrix_columns']]
    minor=[[constant[j][i] for j in range(8)] for i in range(8)];unit=det(minor);assert unit in [1,-1] and all(constant[7][i]+constant[8][i]==0 for i in range(8));assert data['B_N_origin_rank']==8
    assert data['N8_module_generator_order']==cone['eight_N_module_generators']
    proof=dict(DAG_nodes=len(D.nodes),all675735_scalar_nodes_origin_and_unit_syntax_verified=len(D.nodes)==675735,pure_CA48_basechange_all96_equations_bound=True,all8_actual_even_and_odd_pairs_bound=True,CA_tau_not_identified_with_N_tau=True,all160_QQ_Jacobian_entries_independently_AD_verified=True,unit96_determinant=str(d96),unit160_determinant=str(d160),all1280_constant_implicit_derivative_witness_entries_verified=True,all7_prefix160_systems_rebuilt_without_source_solve=True,all64_coordinate_quotients_regular_by_unit_uniqueness=True,all64_quotient_origins_exactly_derived=True,all8by9_map_entries_bound=True,origin_first8_columns_minor=str(unit),origin_kernel='span(e7+e8)',old_CB_r_g_never_composed=True,finite_w_Laurent_scope_gate_reused=True,no_nilrow_full_transition_or_T0_conclusion=True)
    (out/'independent-proof.json').write_text(json.dumps(proof,indent=2)+'\n')
    for p in [Path(__file__),ROOT/'code/ambient_B18_D4_literal_horseshoe_review_20261104_v63.py',ROOT/'code/ambient_B18_actual_top_coefficient_preflight_review_20261104_v63.py',ROOT/f'research/notes/{ID}.md',folder/'report.json']:hashes[str(p.relative_to(ROOT))]=sha(p)
    result=dict(id=ID,all_passed=True,actual_honest9_to_same_R0_N8_map_independently_verified=True,no_author_source_solver_newjets_GB_or_nil_transition=True,input_sha256=hashes,evidence_sha256={'independent-proof.json':sha(out/'independent-proof.json')},seconds=time.monotonic()-start,max_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(dict(all_passed=True,seconds=result['seconds'],rss=result['max_rss_kib'],report_sha256=sha(out/'report.json'))))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
