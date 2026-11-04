#!/usr/bin/env python3
"""Independent read-only audit of the actual finite cokernel Z, not T0.

Cached source2 is already independently accepted. No producer, IFT solve,
source jet, GB, or point is run here. Sparse dual numbers reconstruct the
actual C0 map; descending-pivot QQ elimination audits the finite cokernel.
"""
import argparse,hashlib,json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
AUTHOR='AMBIENT_B18_T0_COKERNEL_Z_20261104_v63'
ID=AUTHOR+'_REVIEW'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

# Dual numbers with a sparse differential, square-zero in all 13 directions.
# The final four directions are the UNASSIGNED old regular-g derivatives.
def add(a,b):
    p=dict(a[1])
    for k,v in b[1].items():
        p[k]=p.get(k,Q(0))+v
        if not p[k]:del p[k]
    return a[0]+b[0],p
def scale(a,c):return a[0]*c,{k:v*c for k,v in a[1].items() if v*c}
def mul(a,b):return add(scale((Q(0),a[1]),b[0]),(a[0]*b[0],{k:v*a[0] for k,v in b[1].items() if v*a[0]}))
def inv(a):
    assert a[0]
    return 1/a[0],{k:-v/a[0]**2 for k,v in a[1].items()}
def dual(vec):return Q(vec[0]),{i:Q(v) for i,v in enumerate(vec[1:]) if Q(v)}
def dense(a):return [a[0]]+[a[1].get(i,Q(0)) for i in range(13)]

def quotient_first(poly,monomials,j,tail):
    """Generic prefix restriction difference divided by its coordinate.

    Work directly with the saved degree<=2 monomials, rather than using a
    closed Hessian-index formula from the producer.
    """
    value=Q(0); differential={}
    for mon,c in zip(monomials,poly):
        c=Q(c)
        if not c or j not in mon:continue
        if not all(k<=j or (tail and k>=7) for k in mon):continue
        rem=list(mon);rem.remove(j)
        if not rem:value+=c
        else:
            assert len(rem)==1
            differential[rem[0]]=differential.get(rem[0],Q(0))+c
    return value,{k:v for k,v in differential.items() if v}

class DescendingSpan:
    """Sparse exact row reduction with the OPPOSITE pivot order to author."""
    def __init__(self):self.rows={}
    def reduce(self,vec):
        v={i:x for i,x in (vec.items() if isinstance(vec,dict) else enumerate(vec)) if x}
        for k in sorted(self.rows,reverse=True):
            if k in v:
                c=v[k]
                for j,x in self.rows[k].items():
                    v[j]=v.get(j,Q(0))-c*x
                    if not v[j]:del v[j]
        return v
    def insert(self,vec):
        v=self.reduce(vec)
        if not v:return False
        k=max(v);c=v[k];self.rows[k]={j:x/c for j,x in v.items()}
        return True

def main(out):
    started=time.monotonic();signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
    out.mkdir(parents=True,exist_ok=False)
    folder=ROOT/'reports'/AUTHOR;report=read(folder/'report.json')
    assert report['all_passed']
    hashes={}
    for p,h in report['input_sha256'].items():assert sha(ROOT/p)==h;hashes[p]=h
    for p,h in report['evidence_sha256'].items():assert sha(folder/p)==h;hashes[str((folder/p).relative_to(ROOT))]=h
    hashes[str((folder/'report.json').relative_to(ROOT))]=sha(folder/'report.json')
    # Numerical cache is certified by a DIFFERENT full firstjet review, and
    # Q1=Jinf by its independently accepted 44-dimensional quadratic gate.
    gates=[
      'reports/AMBIENT_B18_T0_HONEST_C0_20261104_v63_REVIEW/report.json',
      'reports/AMBIENT_B18_T1_T2_QUADRATIC_PHI_20261104_v63_REVIEW/report.json',
      'reports/AMBIENT_B18_T3_ZERO_STRUCTURAL_20261104_v63/report.json',
      'reports/AMBIENT_B18_T3_DOUBLE_ANNIHILATOR_LENGTH3_20261104_v63_REVIEW.json',
      'reports/AMBIENT_B18_W1_TAIL_BINDING_20261104_v63_REVIEW/report.json',
    ]
    for p in gates:
        assert read(ROOT/p)['all_passed'];hashes[p]=sha(ROOT/p)
    data=read(folder/'finite-cokernel-Z-and-actions.json')
    c=read(ROOT/'reports/AMBIENT_B18_T0_HONEST_C0_20261104_v63/honest-C0-D0-D4-DAG.json')
    source=read(ROOT/'reports/AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63/actual-source-coefficients.json')
    full=read(ROOT/'reports/AMBIENT_B18_FULL_RELATIVE_INPUT_20261104_v63/full-relative186-DAG.json')
    cache=read(ROOT/'reports/AMBIENT_B18_T3_FULL9_FIRSTJET_20261104_v63/full9-firstjet-certificate.json')
    q1=read(ROOT/'reports/AMBIENT_B18_T1_T2_FIRSTJET_20261104_v63/Q1-Q2-input-and-firstjets.json')
    nodes=c['DAG_scalar_nodes'];original=source['DAG_scalar_nodes'];N=len(original)
    assert nodes[:N]==original
    adapter=[[Q(v) for v in row] for row in cache['YB_to_YA_linear']]
    monomials=list(map(tuple,cache['source2_monomials']))
    assert len(monomials)==55
    functions={name:{int(k):list(map(Q,v)) for k,v in rec['pure_function2jets'].items()} for name,rec in cache['both57_source2'].items()}
    # Source variable firstjets, with actual separate CA/CB trace bases.
    values={}
    for name in ('CA','CB'):
        for i,v in cache['both57_source2'][name]['all57_source2jets']:
            lin=list(map(Q,v[1:10]))
            if name=='CB':lin=[sum((lin[k]*adapter[k][j] for k in range(9)),Q(0)) for j in range(9)]
            p=(Q(0),{j:x for j,x in enumerate(lin) if x})
            if i in values:assert values[i]==p
            values[i]=p
    for i,n in enumerate(original):
        if n[0]=='v':assert i in values
        elif n[0]=='c':values[i]=(Q(n[1],n[2]),{})
        elif n[0]=='a':values[i]=add(values[n[1]],values[n[2]])
        elif n[0]=='n':values[i]=scale(values[n[1]],Q(-1))
        elif n[0]=='m':values[i]=mul(values[n[1]],values[n[2]])
        else:assert n[0]=='i';values[i]=inv(values[n[1]])
    special={}
    # Reconstruct each cached CA prefix cofactor from source2, including the
    # coefficients actually used by the 10/3 unit Schur map.
    # The membership metadata is needed only to bind function and coordinate.
    hp=ROOT/'reports/AMBIENT_B18_T3_HONEST10_PRESENTATION_20261104_v63/honest10-presentation-DAG.json'
    honest=read(hp);hashes[str(hp.relative_to(ROOT))]=sha(hp)
    for rec in honest['all_pure_source_Hadamard_membership_telescopes']:
        root=rec['original_source_function']
        for j,q in enumerate(rec['all9_regular_coordinate_divisions']):
            h=quotient_first(functions['CA'][root],monomials,j,False)
            assert h[0]==Q(q['origin'])
            special[q['node']]=h
            assert dense(h)[:10]==list(map(Q,cache['all792_CA_h_firstjets'][str(q['node'])]))
    cb=source['families']['CB']
    cb_f63=functions['CB'][cache['both57_source2']['CB']['f63']]
    assert all(not a for m,a in zip(monomials,cb_f63) if all(k>=7 for k in m))
    raw={q['raw_row']:cb['complete_normalized_physical_w_coefficients'][q['weight']][q['raw_row']][1]
         for q in full['all99_CB_kernel_and87_graph_generators_evaluated_on_CA'] if q['kind']=='graph_CB_minus_CA'}
    for rec in full['all_coordinate_divisions']:
        node=rec['node'];label=full['DAG_scalar_nodes'][node][3]
        if label.startswith('g_split_'):
            # gCB has order3, so each of its coordinate-split cofactors has
            # order>=2. This is a known order statement, not an unknown jet.
            h=(Q(0),{})
        else:
            head,coordinate=label.split(':');j=int(coordinate[1:])
            root=cache['both57_source2']['CB']['f63'] if head=='f63' else raw[int(head.replace('zero_graph','').replace('graph',''))]
            h=quotient_first(functions['CB'][root],monomials,j,coordinate.startswith('z'))
            h=(h[0],{j:sum((h[1].get(k,Q(0))*adapter[k][j] for k in range(9)),Q(0)) for j in range(9)})
            h=(h[0],{k:v for k,v in h[1].items() if v})
        assert h[0]==Q(rec['origin_value'])
        assert dense(h)[:10]==list(map(Q,cache['all784_CB_h_firstjets_in_YA'][str(node)]))
        special[node]=h
    gnodes=cache['four_g_derivatives_independent_node_ids']
    for j,node in enumerate(gnodes):special[node]=(Q(0),{9+j:Q(1)})
    special[full['CB_g_restriction_node']]=(Q(0),{})
    for rec in c['all35_regular_coordinate_divisions']:
        root=rec['original_pure_CA_source_function']
        for j,q in enumerate(rec['all7_regular_coordinate_divisions']):
            h=quotient_first(functions['CA'][root],monomials,j,True)
            assert h[0]==Q(q['origin']);special[q['node']]=h
    roots=[i for col in c['C0_free9_to_D4_honest9_columns'] for i in col]
    need=set();todo=roots[:]
    while todo:
        i=todo.pop()
        if i in need:continue
        need.add(i)
        if i<N or i in special:continue
        op,*a=nodes[i]
        if op in ('a','m'):todo.extend(a)
        elif op in ('n','i'):todo.append(a[0])
        else:assert op=='c',(i,nodes[i])
    assert not (need & {q['node'] for q in c['all4_known_regular_gA_divisions']})
    for i in sorted(need):
        if i<N:continue
        if i in special:values[i]=special[i];continue
        op,*a=nodes[i]
        if op=='c':p=(Q(a[0],a[1]),{})
        elif op=='a':p=add(values[a[0]],values[a[1]])
        elif op=='m':p=mul(values[a[0]],values[a[1]])
        elif op=='n':p=scale(values[a[0]],Q(-1))
        else:assert op=='i';p=inv(values[a[0]])
        values[i]=p
    M=[]
    for j,col in enumerate(c['C0_free9_to_D4_honest9_columns']):
        this=[]
        for k,i in enumerate(col):
            h=values[i]
            assert not any(h[1].get(t,Q(0)) for t in range(9,13))
            a=dense(h)[:10];this.append(a)
            assert a==list(map(Q,data['all9_C0_map_exact_full_firstjet_columns'][j][k]))
        M.append(this)
    positions=q1['families']['1']['honest9_column_source_positions_in_Q3_reduced']
    J=[]
    for i in positions:
        a=list(map(Q,cache['all2069_known_firstjet90_columns'][i]))
        assert len(a)==90
        J.append([[a[10*k+j] for j in (0,8,9)] for k in range(9)])
    C=[[[a[j] for j in (0,8,9)] for a in col] for col in M]
    assert len(J)==1548 and len(C)==9
    assert J==[[list(map(Q,a)) for a in col] for col in data['all1548_D4_Q1_relation_columns_over_A3']]
    assert C==[[list(map(Q,a)) for a in col] for col in data['all9_C0_map_columns_over_A3']]
    def flat(col):return [x for a in col for x in a]
    def act(v,k):
        out=[Q(0)]*27
        for i in range(9):out[3*i+k]=v[3*i]
        return out
    space=DescendingSpan()
    for col in J+C:
        v=flat(col)
        for a in (v,act(v,1),act(v,2)):space.insert(a)
    rank=len(space.rows);length=27-rank
    assert rank==data['QQ_relation_rank']==report['QQ_relation_rank_in27']
    assert length==data['Z_length']==report['Z_length']
    for b in space.rows.values():
        v=[b.get(i,Q(0)) for i in range(27)]
        assert not space.reduce(act(v,1)) and not space.reduce(act(v,2))
    # Independently verify the author's basis, instead of copying its pivot
    # algorithm: equality of subspaces follows from both dimensions and spans.
    recorded=DescendingSpan()
    for pivot,pairs in data['QQ_relation_basis']:
        v={i:Q(a) for i,a in pairs}
        assert not space.reduce(v)
        assert recorded.insert(v)
    assert len(recorded.rows)==rank
    free=data['quotient_basis_free_coordinate_indices'];assert len(free)==length
    for j in free:
        e=[Q(0)]*27;e[j]=1
        assert recorded.insert(e)
    assert len(recorded.rows)==27
    actions=[[[Q(x) for x in col] for col in data[key]] for key in ('r_action_columns','s_action_columns')]
    assert all(len(A)==length and all(len(c)==length for c in A) for A in actions)
    for k,A in enumerate(actions,1):
        for j,index in enumerate(free):
            e=[Q(0)]*27;e[index]=1;v=act(e,k)
            for q,x in zip(free,A[j]):v[q]-=x
            assert not space.reduce(v)
    def compose(A,B):
        return [[sum((A[k][i]*B[j][k] for k in range(length)),Q(0)) for i in range(length)] for j in range(length)]
    assert all(not any(x for c in compose(A,B) for x in c) for A in actions for B in actions)
    radical=DescendingSpan()
    for A in actions:
        for col in A:radical.insert(col)
    mu=length-len(radical.rows)
    assert mu==2==data['Z_minimal_R0_generators']==report['Z_minimal_generators']
    assert length<=6
    assert length==3 and free==[21,24,26]
    # Explicit comparison with the dual of A3, not merely the Hilbert function.
    Ar,As=actions
    pair=[[Ar[0][2],Ar[1][2]],[As[0][2],As[1][2]]]
    determinant=pair[0][0]*pair[1][1]-pair[0][1]*pair[1][0]
    assert determinant==Q(36,25)==Q(data['r_s_pairing_determinant'])
    dualmap=[[Q(0),pair[0][0],pair[1][0]],
             [Q(0),pair[0][1],pair[1][1]],[Q(1),Q(0),Q(0)]]
    assert dualmap==[[Q(x) for x in col] for col in data['omega_isomorphism_columns_basis_eps0_epsr_epss']]
    canonical_r=[[Q(0)]*3,[Q(1),Q(0),Q(0)],[Q(0)]*3]
    canonical_s=[[Q(0)]*3,[Q(0)]*3,[Q(1),Q(0),Q(0)]]
    assert compose(dualmap,Ar)==compose(canonical_r,dualmap)
    assert compose(dualmap,As)==compose(canonical_s,dualmap)
    # (2+3t+t^2)(1+t)^7 is the full minimal resolution of omega_A3 over R0.
    from math import comb
    betti=[sum(c*comb(7,k-j) for j,c in enumerate((2,3,1)) if 0<=k-j<=7) for k in range(10)]
    assert betti==data['omega_A3_R0_minimal_resolution_ranks']
    # The exact free lift of D4 -> Z is the last-two-coordinate projection.
    # The seven killed basis vectors are literal columns, not just their jets.
    for k,j in zip((0,1,2,3,4,6,7),range(7)):
        assert c['C0_free9_to_D4_honest9_columns'][k]==[int(i==j) for i in range(9)]
    projected=DescendingSpan();projected_Q1=DescendingSpan()
    allfirst=cache['all2069_known_firstjet90_columns']
    for i,col in enumerate(allfirst):
        p=[Q(col[10*j+k]) for j in (7,8) for k in range(1,10)]
        projected.insert(p)
        if i in positions:projected_Q1.insert(p)
    assert len(projected.rows)==len(projected_Q1.rows)==16==data['projected_J_full9_linear_rank']
    for v in projected.rows.values():assert not projected_Q1.reduce(v)
    recorded_linear=DescendingSpan()
    for pivot,pairs in data['projected_J_full9_basis']:
        v={i:Q(x) for i,x in pairs}
        assert not projected.reduce(v) and recorded_linear.insert(v)
    assert len(recorded_linear.rows)==16
    assert data['Tor1_D4_to_Z_rank']==16
    mu_B0=betti[1]-16+9-2
    assert mu_B0==8==data['image_B0_minimal_generators']==report['actual_image_B0_minimal_generators']
    # Exact Tor sequence: this proves nonzero T0 from mu(D0)=9, not from the
    # origin rank7 of its map to D4.  Distinguish leading class from a function.
    qg=[Q(0)]*8+[Q(1)]
    qt=[Q(0)]*5+[Q(1),Q(1,2),Q(-1),Q(0)]
    originspace=DescendingSpan()
    for col in M:originspace.insert([a[0] for a in col])
    assert len(originspace.rows)==7
    for name,q,saved in zip(('gA_tau','f60_plus_half_f67_minus_f72'),(qg,qt),data['unique_visible_T0_class_in_D0_fiber']):
        assert saved['label']==name and q==list(map(Q,saved['constant_C0_vector']))
        assert all(sum((q[i]*M[i][j][0] for i in range(9)),Q(0))==0 for j in range(9))
        p=[sum((q[i]*M[i][j][k] for i in range(9)),Q(0)) for j in (7,8) for k in range(1,10)]
        assert p==list(map(Q,saved['projected_firstjet18']))
        outside=bool(projected.reduce(p))
        assert outside==saved['nonzero_in_Tor1_Z_mod_Tor1_D4']==(name!='gA_tau')
    vars={n[1]:i for i,n in enumerate(original) if n[0]=='v'}
    sol={i:list(map(Q,v)) for i,v in cache['both57_source2']['CA']['all57_source2jets']}
    physical=data['all9_actual_physical_R0_actions']
    for row,record in zip(source['physical_zero_base_rows9'],physical):
        a,b=sol[vars[f'x{row}']][8:10]
        assert record['row']==row and [a,b]==list(map(Q,record['linear_r_s']))
        expect=[[a*x+b*y for x,y in zip(c,d)] for c,d in zip(*actions)]
        assert expect==[[Q(x) for x in c] for c in record['action_columns']]
    proof=dict(both_actual_zero_ideals_annihilate_Z=True,
      right_exact_complete_cokernel_basechange_no_kernel_or_Tor_shortcut=True,
      actual_A3_length=3,independent_C0_full_firstjets=81,
      all35_CA_tail_quotient_firstjets_and792_plus784_prior_prefix_jets_checked=True,
      old_four_g_first_derivatives_retained_and_drop_exactly=True,
      new_CA_four_g_quotients_not_in_map_dependency_closure=True,
      all1548_Q1_columns_plus9_C0_columns_exact=True,
      independent_descending_relation_rank=rank,Z_exact_length=length,Z_minimal_generators=mu,
      Z_is_canonical_module_of_A3_with_pairing_det=str(determinant),
      omega_A3_minimal_resolution_ranks=betti,
      exact_free_projection_and_independent_Tor1_rank=16,
      exact_image_B0_minimal_generators=mu_B0,T0_nonzero_by_mu_9_vs8=True,
      unique_visible_fiber_class_gA_not_a_literal_kernel_function=True,
      author_basis_equal_same_subspace=True,all9_physical_R0_actions_verified=True,
      r_s_actions_commute_and_maxideal_squared_zero=True,
      no_author_producer_or_source_solver_run=True,no_new_source_IFT_jets_GB_points=True,
      T0_finite_socle_or_depth_not_claimed=True)
    (out/'independent-Z-proof.json').write_text(json.dumps(proof,indent=2)+'\n')
    for p in (Path(__file__),ROOT/'research/notes'/f'{ID}.md'):
        hashes[str(p.relative_to(ROOT))]=sha(p)
    review=dict(id=ID,all_passed=True,input_sha256=hashes,
      evidence_sha256={'independent-Z-proof.json':sha(out/'independent-Z-proof.json')},
      Z_length=length,Z_minimal_generators=mu,QQ_relation_rank_in27=rank,
      Z_is_omega_A3=True,Tor1_D4_to_Z_rank=16,image_B0_minimal_generators=mu_B0,
      actual_T0_nonzero_independently_accepted=True,T0_finite_socle_or_depth_not_claimed=True,
      seconds=time.monotonic()-started,max_rss_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    (out/'report.json').write_text(json.dumps(review,indent=2)+'\n')
    assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
    print(json.dumps(review))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True)
    main(p.parse_args().output)
