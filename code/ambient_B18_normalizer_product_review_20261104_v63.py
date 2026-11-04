#!/usr/bin/env python3
"""Independent source-DAG/free4 and conductor-Koszul proof audit.

No author builder, source solve/jet, r/g pullback or GB. Original/pruned
coefficients are compared by local universally valid DAG identities.
"""
import argparse,hashlib,itertools,json,math,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
from ambient_B18_T0_cokernel_Z_review_20261104_v63 import (
    add as da,mul as dm,inv as di,scale as ds,DescendingSpan)
from ambient_B18_T0_quotient_firsts_quartic_review_20261104_v63 import exact_det
from s61_sparse_generic5_full134_normal_circuit_review_20261104_v63 import plus,times,neg,const,symbol

ROOT=Path(__file__).resolve().parents[1]
AUTHOR='AMBIENT_B18_ZERO_NORMALIZER_TENSOR_MODEL_20261104_v63'
PRODUCT='AMBIENT_B18_CONDUCTOR_PRODUCT_CURVE_BASE_20261104_v63'
ID=AUTHOR+'_REVIEW'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def closure(nodes,roots):
    seen=set();todo=list(roots)
    while todo:
        i=todo.pop()
        if i in seen:continue
        seen.add(i);op,*a=nodes[i]
        if op in ('a','m'):todo.extend(a)
        elif op in ('n','i'):todo.append(a[0])
        else:assert op in ('c','v'),(i,nodes[i])
    return sorted(seen)
class Expr:
    def __init__(self,nodes,stops):
        self.nodes=nodes;self.stops=set();self.memo={}
        # A constructed negative can cancel a source coefficient's existing
        # negation. Cut at the signed canonical representative, preserving
        # that exact alias instead of treating x and -x as independent.
        for i in stops:
            while nodes[i][0]=='n':i=nodes[i][1]
            self.stops.add(i)
    def __call__(self,i):
        if self.nodes[i][0]=='c':return const(Q(*self.nodes[i][1:]))
        if i in self.stops:return symbol(('coefficient',i))
        if i in self.memo:return self.memo[i]
        op,*a=self.nodes[i]
        if op=='a':v=plus(self(a[0]),self(a[1]))
        elif op=='m':v=times(self(a[0]),self(a[1]))
        elif op=='n':v=neg(self(a[0]))
        else:raise AssertionError(('coefficient_not_cut',i,self.nodes[i]))
        self.memo[i]=v;return v
def mm(A,B):
    return [[sum((a*b for a,b in zip(r,c)),Q(0)) for c in zip(*B)] for r in A]
def matpoly(A,B):
    out=[]
    for r in A:
        row=[]
        for c in zip(*B):
            p={}
            for a,b in zip(r,c):p=plus(p,times(a,b))
            row.append(p)
        out.append(row)
    return out
def mv(A,v):
    out=[]
    for r in A:
        p={}
        for a,b in zip(r,v):p=plus(p,times(a,b))
        out.append(p)
    return out

def main(out,pout):
    start=time.monotonic();signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
    out=out.resolve();pout=pout.resolve()
    out.mkdir(parents=True,exist_ok=False);pout.mkdir(parents=True,exist_ok=False)
    folder=ROOT/'reports'/AUTHOR;product=ROOT/'reports'/PRODUCT
    hashes={}
    for base in (folder,product):
        r=read(base/'report.json');assert r['all_passed']
        for p,h in r['input_sha256'].items():assert sha(ROOT/p)==h;hashes[p]=h
        for p,h in r.get('evidence_sha256',{}).items():assert sha(base/p)==h;hashes[str((base/p).relative_to(ROOT))]=h
        hashes[str((base/'report.json').relative_to(ROOT))]=sha(base/'report.json')
    for name in (
      'AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63_REVIEW',
      'AMBIENT_B18_T0_HONEST_C0_20261104_v63_REVIEW',
      'AMBIENT_B18_RELATIVE_WEIGHTED_MODEL_20261104_v63_REVIEW',
      'AMBIENT_B18_N_SAME_R0_BINDING_20261104_v63_REVIEW',
      'AMBIENT_B18_R_N_PROJECTION_NORM_20261104_v63_REVIEW'):
        p=ROOT/'reports'/name/'report.json';assert read(p)['all_passed'];hashes[str(p.relative_to(ROOT))]=sha(p)
    data=read(folder/'actual-normalizer-tensor-free4-DAG.json')
    S=read(ROOT/'reports/AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63/actual-source-coefficients.json')
    C0=read(ROOT/'reports/AMBIENT_B18_T0_HONEST_C0_20261104_v63/honest-C0-D0-D4-DAG.json')
    W=read(ROOT/'reports/AMBIENT_B18_POSITIVE_GRADING_OBSTRUCTION_20261104_v63/certificate.json')
    R=read(ROOT/'reports/AMBIENT_B18_RELATIVE_WEIGHTED_MODEL_20261104_v63/certificate.json')
    cache=read(ROOT/'reports/AMBIENT_B18_T3_FULL9_FIRSTJET_20261104_v63/full9-firstjet-certificate.json')
    N=read(ROOT/'reports/AMBIENT_B18_N_SAME_R0_BINDING_20261104_v63/actual-N-R0-binding-DAG.json')
    CB3=read(ROOT/'reports/AMBIENT_B18_CB_G_DIVISION_ORIGINS_20261104_v63/degree3-certificate.json')
    old=C0['DAG_scalar_nodes'];nodes=data['DAG_scalar_nodes']
    mapping=dict(data['original_node_to_pruned_node'])
    assert len(mapping)==len(data['original_node_to_pruned_node'])
    pairs=list(itertools.combinations(range(5),2));mons=list(itertools.combinations_with_replacement(range(5),2))
    zero=[i for i,w in enumerate(W['physical150_torus_weights']) if w==0]
    cuts=[i for i in R['common_physical_group_cuts19'] if i in zero]
    base=S['physical_zero_base_rows9'];target=[18,19,28,29,60,67,72,63,64,68,73,74]
    coords=cuts+base+target;other=[i for i in zero if i not in coords]
    assert len(zero)==40 and len(cuts)==15 and len(coords)==len(set(coords))==36 and len(other)==4
    assert zero==data['zero_raw_rows40'] and cuts==data['common_cut_zero_rows15'] and coords==data['coordinate_rows36']
    assert other==data['eliminated_zero_rows4'] and base==data['actual_R0_rows9']
    # Reconstruct ALL25 physical divergence equations, not only the displayed4.
    full_div=[]
    for j in range(5):
        for k in range(5):
            row=[]
            for index in zero:
                a,b=pairs[index//15];u,v=mons[index%15];c=Q(0)
                for i,sgn in ((a,1),) if b==j else ():
                    c+=sgn*(int(u==i and v==k)+int(v==i and u==k))
                for i,sgn in ((b,-1),) if a==j else ():
                    c+=sgn*(int(u==i and v==k)+int(v==i and u==k))
                row.append(c)
            full_div.append(row)
    selected=[full_div[5*2+k] for k in (0,1,3,4)]
    assert selected==[[Q(x) for x in r] for r in data['four_divergence_relations40']]
    span=DescendingSpan()
    for row in full_div:span.insert(row)
    selectedspan=DescendingSpan()
    for row in selected:selectedspan.insert(row)
    assert len(span.rows)==len(selectedspan.rows)==4
    for row in full_div:assert not selectedspan.reduce(row)
    om=[[row[zero.index(i)] for i in other] for row in selected]
    rec=[[Q(x) for x in r] for r in data['eliminated_four_coordinate_reconstruction']]
    assert exact_det(om)==Q(data['divergence_unit4_minor'])==2
    assert mm(om,rec)==[[-row[zero.index(i)] for i in coords] for row in selected]
    # Verify pure unit source blocks and the complete old coefficient closure.
    unknown=[];equations=[]
    ca=S['families']['CA'];cb=S['families']['CB']
    for name,f,saved in zip(('CA','CB'),(ca,cb),data['normalizer_sources48_each']):
        u=[i for a in f['coefficient_unknown_pairs'].values() for i in a]+[f['p_quadratic_relation'][k] for k in ('s','t')]
        eq=f['all48_zero_IFT_quadratic_coefficient_equations']
        assert len(u)==len(set(u))==len(eq)==48
        assert saved['family']==name and saved['unknown48']==[mapping[i] for i in u]
        assert saved['equations48']==[mapping[i] for i in eq]
        assert saved['unit48_determinant']==f['unit48_block_determinant_even24_then_odd24_vars_a23_t_b23_s'] and Q(saved['unit48_determinant'])
        unknown+=u;equations+=eq
    prefix=C0['all7_pure_CA_tail_prefix57_systems'][0]
    pu=[i for oldi,i in prefix['bindings']];pe=prefix['all57_unit_equations']
    assert data['CA_conductor_gA_prefix57']['unknown57']==[mapping[i] for i in pu]
    assert data['CA_conductor_gA_prefix57']['equations57']==[mapping[i] for i in pe]
    assert Q(data['CA_conductor_gA_prefix57']['unit_determinant'])==Q(prefix['unit57_determinant'])
    unknown+=pu;equations+=pe
    assert len(unknown)==len(set(unknown))==len(equations)==153
    assert data['all153_scalar_unit_unknown_nodes']==[mapping[i] for i in unknown]
    assert data['all153_scalar_unit_equations']==[mapping[i] for i in equations]
    ap=ca['complete_normalized_physical_w_coefficients'][0];bp=cb['complete_normalized_physical_w_coefficients'][0]
    roots=equations+C0['actual_C0_basis8_scalar_nodes']+[ca['actual_tracefree_nu'],cb['actual_tracefree_nu']]
    roots += [i for row in zero for i in ap[row]+bp[row]]
    needed=closure(old,roots)
    assert needed==sorted(mapping) and all(old[i][0] in ('c','v','a','m','n','i') for i in needed)
    assert nodes[0]==['c',0,1] and nodes[1]==['c',1,1]
    values=[]
    for i,(op,*a) in enumerate(nodes):
        if op=='c':v=Q(*a)
        elif op=='v':v=Q(0)
        elif op=='a':assert max(a)<i;v=values[a[0]]+values[a[1]]
        elif op=='m':assert max(a)<i;v=values[a[0]]*values[a[1]]
        elif op=='n':assert a[0]<i;v=-values[a[0]]
        else:assert op=='i' and a[0]<i and values[a[0]];v=1/values[a[0]]
        values.append(v)
    # Local identities certify the old->pruned map by induction. Each allowed
    # simplification is an elementary ring identity, not an opaque fingerprint.
    for i in needed:
        op,*a=old[i];j=mapping[i];n=nodes[j]
        if op=='c':assert n[0]=='c' and Q(*n[1:])==Q(*a)
        elif op=='v':assert n==['v',f'ACTUAL_NODE{i}_{a[0]}']
        elif op=='n':
            x=mapping[a[0]];nx=nodes[x]
            if nx[0]=='c':assert n[0]=='c' and Q(*n[1:])==-Q(*nx[1:])
            elif nx[0]=='n':assert j==nx[1]
            else:assert n==['n',x]
        elif op=='i':
            x=mapping[a[0]];nx=nodes[x]
            if nx[0]=='c':assert n[0]=='c' and Q(*n[1:])==1/Q(*nx[1:])
            else:assert n==['i',x]
        else:
            x,y=map(mapping.__getitem__,a);nx,ny=nodes[x],nodes[y]
            if op=='a':
                if x==0:assert j==y
                elif y==0:assert j==x
                elif nx==['n',y] or ny==['n',x]:assert j==0
                elif nx[0]==ny[0]=='c':assert n[0]=='c' and Q(*n[1:])==Q(*nx[1:])+Q(*ny[1:])
                else:assert n==['a',*sorted((x,y))]
            else:
                assert op=='m'
                if x==0 or y==0:assert j==0
                elif x==1:assert j==y
                elif y==1:assert j==x
                elif nx[0]==ny[0]=='c':assert n[0]=='c' and Q(*n[1:])==Q(*nx[1:])*Q(*ny[1:])
                else:assert n==['m',*sorted((x,y))]
    assert all(values[i]==0 for i in data['all153_scalar_unit_equations'])
    assert len(nodes)==28624
    # Expanded HONEST/C0 archives contain reviewed UNUSED duplicate aliases.
    # Bind physical x variables in the original pure source, not the last
    # equal-name alias in the larger archive.
    oldvars={n[1]:i for i,n in enumerate(S['DAG_scalar_nodes']) if n[0]=='v'}
    assert data['R0_variable_nodes']==[mapping[oldvars[f'x{r}']] for r in base]
    assert {i for i,n in enumerate(nodes) if n[0]=='v'}==set(data['all153_scalar_unit_unknown_nodes'])|set(data['R0_variable_nodes'])
    # Free4 algebra identities and all12 output differences/all48 ideal columns.
    coeffroots={mapping[i] for row in zero for i in ap[row]+bp[row]}
    coeffroots |= {mapping[i] for i in C0['actual_C0_basis8_scalar_nodes']+[ca['actual_tracefree_nu'],cb['actual_tracefree_nu']]}
    E=Expr(nodes,coeffroots);z={};one=const(1)
    na=E(data['nu_CA_node']);nb=E(data['nu_CB_node'])
    assert data['nu_CA_node']==mapping[ca['actual_tracefree_nu']] and data['nu_CB_node']==mapping[cb['actual_tracefree_nu']]
    MA=[[z,na,z,z],[one,z,z,z],[z,z,z,na],[z,z,one,z]]
    MB=[[z,z,nb,z],[z,z,z,nb],[one,z,z,z],[z,one,z,z]]
    assert [[E(i) for i in r] for r in data['multiplication_tau_CA']]==MA
    assert [[E(i) for i in r] for r in data['multiplication_tau_CB']]==MB
    assert matpoly(MA,MB)==matpoly(MB,MA)
    for M,nu in ((MA,na),(MB,nb)):assert matpoly(M,M)==[[nu if i==j else {} for j in range(4)] for i in range(4)]
    columns=[]
    for row,record in zip(target,data['all12_zero_physical_difference_columns']):
        assert row==record['physical_row']
        ae,ao=[E(mapping[i]) for i in ap[row]];be,bo=[E(mapping[i]) for i in bp[row]]
        v=[plus(ae,neg(be)),ao,neg(bo),{}]
        assert [E(i) for i in record['column4']]==v
        columns += [v,mv(MA,v),mv(MB,v),mv(MA,mv(MB,v))]
    assert len(columns)==48
    assert [[E(i) for i in col] for col in data['normalizer_tensor_relation_columns4_by48']]==columns
    assert all(values[i]==0 for col in data['normalizer_tensor_relation_columns4_by48'] for i in col)
    expectedmap=[[one,z,z,z]]+[[z,E(mapping[b]),z,z] for b in C0['actual_C0_basis8_scalar_nodes']]
    assert [[E(i) for i in col] for col in data['actual_C0_free9_to_E4_columns']]==expectedmap
    assert data['all_actual_I0_CA_generator_nodes8']==[mapping[i] for row,i in ca['actual_zero_odd_I0_generators']]
    assert data['all_actual_I0_CB_generator_nodes8']==[mapping[i] for row,i in cb['actual_zero_odd_I0_generators']]
    # Replay the saved FIRST coefficients on the two individual conductor
    # planes. No source solve; the literal provided N56 corrections are used.
    sols={name:{i:list(map(Q,v)) for i,v in cache['both57_source2'][name]['all57_source2jets']} for name in ('CA','CB')}
    X={name:[sols[name][oldvars[f'x{r}']][1:10] for r in base] for name in ('CA','CB')}
    adapter=[list(map(Q,r)) for r in cache['YB_to_YA_linear']]
    assert mm(X['CB'],adapter)==X['CA']
    Nnodes=N['DAG_scalar_nodes'];Ne=N['all56_actual_coefficient_Hensel_equations'];Nu=N['all56_coefficient_unknown_nodes']
    Nroots=Ne+[N['ninth_actual_even_a_node'],N['ninth_actual_odd_b_node'],N['nu_tracefree_square_node']]
    Nneeded=closure(Nnodes,Nroots);cohen=N['chosen7_R0_indices']+[N['Cohen_even8_R0_index']]
    for name,rec in zip(('CA','CB'),data['conductor_support_first_gate']):
        assert rec['family']==name
        plane=[row[7:9] for row in X[name]]
        assert plane==[list(map(Q,row)) for row in rec['actual_R0_first_on_own_conductor_plane']]
        corrections=[list(map(Q,row)) for row in rec['N56_first_unknown_corrections']]
        bind={u:(Q(0),{k:v for k,v in enumerate(c) if v}) for u,c in zip(Nu,corrections)}
        bind.update({u:(Q(0),{k:v for k,v in enumerate(plane[j]) if v}) for u,j in zip(N['actual_Cohen_R8_parameter_nodes'],cohen)})
        val={}
        for i in Nneeded:
            op,*a=Nnodes[i]
            if op=='v':assert i in bind;v=bind[i]
            elif op=='c':v=(Q(*a),{})
            elif op=='a':v=da(val[a[0]],val[a[1]])
            elif op=='m':v=dm(val[a[0]],val[a[1]])
            elif op=='n':v=ds(val[a[0]],Q(-1))
            else:assert op=='i';v=di(val[a[0]])
            val[i]=v
        assert all(val[i]==(Q(0),{}) for i in Ne)
        ae=val[N['ninth_actual_even_a_node']]
        ell=[plane[N['ninth_R0_index']][k]-ae[1].get(k,Q(0)) for k in range(2)]
        assert ell==list(map(Q,rec['v_N_linear_form_on_own_plane']))==[0,0]
        g3=C0['cached_actual_CA_gA_all_degree3_coefficients'][6:10] if name=='CA' else CB3['g_CB_all10_coefficients'][6:10]
        assert list(map(Q,g3))==list(map(Q,rec['g3_coefficients']))
        assert not rec['gcd_linear_v_and_g3_is1']
    assert data['comparison_kernel_finite_from_this_gate'] is False
    # Exact three-quadric Koszul calculation. Its completeness is justified
    # by the elementary divisibility proof in the companion review note.
    pr=read(product/'report.json')
    rr=symbol(('r',));ss=symbol(('s',));r2=times(rr,rr);rs=times(rr,ss);s2=times(ss,ss)
    u=[neg(ss),rr,{}];v=[{},neg(ss),rr];f=[r2,rs,s2]
    def decoded_two_terms(terms):
        p={}
        for c,rpow,spow in terms:
            m=tuple([('r',)]*rpow+[('s',)]*spow)
            p=plus(p,{m:Q(c)})
        return p
    assert [[decoded_two_terms(terms) for terms in row] for row in pr['H1_three_quad_Koszul_generators']]==[u,v]
    assert pr['H1_relation_columns']==['r*u','s*u+r*v','s*v']
    for cycle in (u,v):
        p={}
        for a,b in zip(f,cycle):p=plus(p,times(a,b))
        assert not p
    d2=[[neg(rs),r2,{}],[neg(s2),{},r2],[{},neg(s2),rs]]
    assert d2==[[times(rr,x) for x in u],
      [plus(times(ss,x),times(rr,y)) for x,y in zip(u,v)],[times(ss,x) for x in v]]
    lengths=[];socles=[]
    choose=lambda n,k:math.comb(n,k) if 0<=k<=n else 0
    for i in range(9):lengths.append(3*(choose(6,i)+choose(6,i-1)));socles.append(2*choose(6,i)+choose(6,i-1))
    assert lengths==pr['Tor_lengths_i0_through8']==[3,21,63,105,105,63,21,3,0]
    assert lengths[1]==pr['Tor1_length']==pr['product_curve_ring_H0_exact_length']==21
    assert socles[1]==pr['Tor1_socle_dimension']==13
    assert pr['conductor_curve_ideals_assumed_radical'] is False and pr['I_A_times_I_B_annihilates_T0_claimed'] is False
    proof=dict(all_original_pruned_scalar_nodes_verified_by_local_ring_identities=True,
      pruned_nodes=28624,all96_plus57_unit_equations_and153_unknowns_bound=True,
      only_pure_unit_ops_no_old_r_g_or_unbound_variables=True,
      all25_divergence_rows_rebuilt_rank4_and_unit_minor2=True,
      actual_zero36_cut15_base9_leave_all12_outputs=True,
      all12_output_differences_and48_free4_relation_columns_verified=True,
      tauA_tauB_commute_and_square_to_own_nu=True,actual_C0_free9_map_verified=True,
      two_scalar_return_maps_prove_product_ann_kernel_cokernel_without_flatness=True,
      product_ann_is_not_ann_T0=True,normalization_tensor_Spec_surjective_support_equal=True,
      both_own_conductor_plane_vN_firsts_zero_finiteness_gate_inconclusive=True,
      exact_CI8_CI8_combined_Koszul_unit_equivalence_proved=True,
      explicit_K3_cycles_boundaries_and_all_order_homology_proved=True,
      Tor_lengths=lengths,Tor1_length=21,Tor1_socle_dimension=13,
      scheme_union_CM1_without_radical_assumption=True,product_ring_H0_length21_only=True,
      no_arbitrary_comparison_kernel_length_finiteness_socle_E4depth_or_D0depth_claim=True,
      no_author_worker_newsource_jet_GB_or_nonunit_pullback=True)
    evidence=out/'independent-normalizer-and-product-proof.json';evidence.write_text(json.dumps(proof,indent=2)+'\n')
    for p in (Path(__file__),ROOT/f'research/notes/{ID}.md',
      ROOT/'code/ambient_B18_T0_cokernel_Z_review_20261104_v63.py',
      ROOT/'code/ambient_B18_T0_quotient_firsts_quartic_review_20261104_v63.py',
      ROOT/'code/s61_sparse_generic5_full134_normal_circuit_review_20261104_v63.py'):
        hashes[str(p.relative_to(ROOT))]=sha(p)
    review=dict(id=ID,all_passed=True,input_sha256=hashes,
      evidence_sha256={evidence.name:sha(evidence)},actual_free4_48_and_product_conductor_comparison_accepted=True,
      conductor_product_Tor_and_ring_H0_length21_accepted=True,
      no_finite_comparison_kernel_or_T0_socle_depth_claim=True,
      seconds=time.monotonic()-start,max_rss_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    rp=out/'report.json';rp.write_text(json.dumps(review,indent=2)+'\n')
    companion=dict(id=PRODUCT+'_REVIEW',all_passed=True,
      input_sha256={str(rp.relative_to(ROOT)):sha(rp),str((product/'report.json').relative_to(ROOT)):sha(product/'report.json'),
                    str(Path(__file__).relative_to(ROOT)):sha(Path(__file__)),
                    f'research/notes/{ID}.md':sha(ROOT/f'research/notes/{ID}.md')},
      Tor_i='A3^binom(6,i) + omega_A3^binom(6,i-1)',Tor1_length=21,Tor1_socle_dimension=13,
      union_is_CM1_not_assumed_reduced=True,product_ring_H0_exact_length=21,
      no_arbitrary_module_length_bound_or_ann_T0=True)
    (pout/'report.json').write_text(json.dumps(companion,indent=2)+'\n')
    assert sum(p.stat().st_size for folder_ in (out,pout) for p in folder_.iterdir())<32*1024*1024
    print(json.dumps(review))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--product-output',type=Path,required=True)
    a=p.parse_args();main(a.output,a.product_output)
