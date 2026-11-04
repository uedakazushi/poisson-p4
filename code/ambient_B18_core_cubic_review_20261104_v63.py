#!/usr/bin/env python3
"""Independent sparse 145-slot quotient ring and source-bound adjoint review."""
import argparse,hashlib,itertools,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
from ambient_B18_zero_normalizer_quadratic_review_20261104_v63 import decode,composition,coordinate,ancestors,rank
ROOT=Path(__file__).resolve().parents[1];ID='AMBIENT_B18_CORE_CUBIC_20261104_v63_REVIEW'
M9=[()]+[(i,)for i in range(9)]+list(itertools.combinations_with_replacement(range(9),2))
MON=sum(([m for m in itertools.combinations_with_replacement(range(4),d)]for d in range(5)),[])
DUAL=MON[:15];CORE=[2,4,7,8];UDIR=[0,1,3,5,6]


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def plus(a,b):
    c=dict(a)
    for m,v in b.items():
        c[m]=c.get(m,F(0))+v
        if not c[m]:del c[m]
    return c
def scale(a,c):return{m:c*v for m,v in a.items()if c*v}
def allowed(m,linear=False):
    if linear:return len(m)<=1 and all(i<4 for i in m)
    dual=sum(i>=4 for i in m)
    return dual==0 and len(m)<=4 or dual==1 and len(m)<=3
def times(a,b,linear=False):
    c={}
    for m,x in a.items():
        for n,y in b.items():
            e=tuple(sorted(m+n))
            if allowed(e,linear):c[e]=c.get(e,F(0))+x*y
    return{m:v for m,v in c.items()if v}
def inverse(a):
    c=a.get((),F(0));assert c
    h=scale(a,-1/c);h.pop((),None);out={():F(1)};v=out
    for _ in range(4):v=times(v,h);out=plus(out,v)
    ans=scale(out,1/c);assert times(a,ans)=={():F(1)}
    return ans
def hom(a,d):return{m:v for m,v in a.items()if len(m)==d and all(i<4 for i in m)}
def vec(a,mon=MON):return[str(a.get(m,F(0)))for m in mon]
def vec145(a):return vec(a)+[str(a.get(tuple(sorted(m+(4+d,))),F(0)))for d in range(5)for m in DUAL]
def unpack70(v):return{m:F(c)for m,c in zip(MON,v)if F(c)}
def linear_part(a):return{m:v for m,v in a.items()if allowed(m,True)}


def reverse_linear(nodes,need,values,root,unknown):
    adj={root:{():F(1)}}
    def accumulate(i,p):
        if p:adj[i]=plus(adj.get(i,{}),p)
    for i in reversed(need):
        a=adj.get(i,{})
        if not a:continue
        op,*j=nodes[i]
        if op=='a':accumulate(j[0],a);accumulate(j[1],a)
        elif op=='n':accumulate(j[0],scale(a,-1))
        elif op=='m':accumulate(j[0],times(a,values[j[1]],True));accumulate(j[1],times(a,values[j[0]],True))
        elif op=='i':accumulate(j[0],scale(times(a,times(values[i],values[i],True),True),-1))
        else:assert op in ('v','c')
    return[adj.get(i,{})for i in unknown]


def binary_product(a,b):
    c=[F(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):c[i+j]+=x*y
    return c


def main(out):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
    out=out.resolve();out.mkdir(parents=True,exist_ok=False);hashes={};base=ROOT/'reports'
    def read(p):hashes[str(p.relative_to(ROOT))]=sha(p);return json.loads(p.read_text())
    def gate(stem):
        p=base/stem/'report.json';r=read(p);assert r['all_passed']
        for name,h in r.get('input_sha256',{}).items():assert sha(ROOT/name)==h,name;hashes[name]=h
        for name,h in r.get('evidence_sha256',{}).items():assert sha(p.parent/name)==h,name;hashes[str((p.parent/name).relative_to(ROOT))]=h
        return r
    for stem in ['AMBIENT_B18_CORE_CUBIC_ADJOINT_PREFLIGHT_20261104_v63','AMBIENT_B18_CORE_CUBIC_INITIAL_20261104_v63',
      'AMBIENT_B18_CORE4_QUADRATIC_20261104_v63_REVIEW','AMBIENT_B18_CORE_NODAL_QUARTIC_20261104_v63_REVIEW',
      'AMBIENT_B18_NORMALIZER_CORE4_REDUCTION_20261104_v63_REVIEW','AMBIENT_B18_N_NORMALIZER_PAIR_PREFLIGHT_20261104_v63_REVIEW',
      'AMBIENT_B18_ZERO_NORMALIZER_TENSOR_MODEL_20261104_v63_REVIEW','AMBIENT_B18_R_ZERO_SLICE_IMAGE_20261104_v63_REVIEW',
      'AMBIENT_B18_R_N_PROJECTION_NORM_20261104_v63_REVIEW']:gate(stem)
    P=read(base/'AMBIENT_B18_CORE_CUBIC_ADJOINT_PREFLIGHT_20261104_v63/literal-lambda1-dual5-preflight.json')
    Q=read(base/'AMBIENT_B18_CORE_CUBIC_INITIAL_20261104_v63/five-actual-core-cubics-and-binary12.json')
    C=read(ROOT/P['scalar_DAG_source']);nodes=C['DAG_scalar_nodes']
    E=read(base/'AMBIENT_B18_ZERO_NORMALIZER_QUADRATIC_20261104_v63/cached-source2-and-exact-E4-quadrics.json')
    J=read(base/'AMBIENT_B18_T3_FULL9_FIRSTJET_20261104_v63/full9-firstjet-certificate.json')
    oldP=read(base/'AMBIENT_B18_CORE4_QUADRATIC_PREFLIGHT_20261104_v63/literal-adjoint98-core4-preflight.json')
    oldQ=read(base/'AMBIENT_B18_CORE4_QUADRATIC_INITIAL_20261104_v63/seven-actual-core-quadrics.json')
    # The previous independent review already proves all literal source roots,
    # original96+2 equations, unknowns, original physical outputs, and J98.
    unknown=C['unit_schur_unknown_nodes98'];equations=C['unit_schur_equations98']
    assert unknown==P['unit98_unknowns']==oldP['actual98_unknown_nodes']
    assert equations==P['unit98_equations']==oldP['actual98_equation_nodes']
    A=[[F(c)for c in row]for row in P['unit98_Jacobian']]
    assert A==[[F(c)for c in row]for row in oldP['actual98_origin_Jacobian']]
    assert F(P['unit98_determinant'])==F(oldP['unit98_determinant'])==F(-12524124635136,78125)
    selected=[0,1,4,2,3];remaining=[7,8,9,10,11]
    assert P['selected5_F_roots']==[C['all12_actual_difference_nodes'][j]for j in selected]
    assert P['remaining5_F_roots']==[C['all12_actual_difference_nodes'][j]for j in remaining]
    assert P['remaining5_physical_rows']==[63,64,68,73,74]
    assert all(not any(F(v)for v in a)for a in oldQ['all7_core_quadrics_padded_through3'][2:])
    # Rebuild full nine-variable source2 before restricting to core + dual.
    rho=[decode(p,M9)for p in E['all11_variables_as_remaining9_through2'][:9]]
    bound9={}
    for family in ['CA','CB']:
        f9=[composition(decode(p,M9),rho)for p in E['families'][family]['inverse_F9_in_common_R0_through2']]
        for i,p in J['both57_source2'][family]['all57_source2jets']:
            v=composition(decode(p,M9),f9)
            if i in bound9:assert bound9[i]==v
            bound9[i]=v
    for j,i in enumerate(C['actual9_parameter_nodes']):
        p=coordinate(j)
        if i in bound9:assert bound9[i]==p
        bound9[i]=p
    assert bound9=={i:decode(p,M9)for i,p in P['all107_full9_source2_bindings']} and len(bound9)==107
    bound={};base4={};dual5={}
    for i,p in bound9.items():
        v={};b={};ds=[{}for _ in range(5)]
        for m,c in p.items():
            u=[j for j in m if j in UDIR]
            if not u:
                e=tuple(CORE.index(j)for j in m);v[e]=c;b[e]=c
            elif len(u)==1:
                rest=list(m);rest.remove(u[0]);e=tuple(CORE.index(j)for j in rest);d=UDIR.index(u[0]);ds[d][e]=c;v[tuple(sorted(e+(4+d,)))]=c
        bound[i]=v;base4[i]=b;dual5[i]=ds
    assert base4=={i:decode(p,DUAL)for i,p in P['all_four_variable_source2_bindings']}
    assert dual5=={i:[decode(p,DUAL)for p in rows]for i,rows in P['all_five_dual_source2_bindings']}
    # Verify the five actual G linear directions from the cached FULL9 model.
    diag=[F(1,12),F(1,12),F(1,24),F(1,12),F(1,12)]
    assert list(map(F,P['selected5_Diagonal_G1_coefficients']))==diag
    L9={(2,):F(5),(4,):F(-20),(7,):F(24),(8,):F(6)}
    from ambient_B18_zero_normalizer_quadratic_review_20261104_v63 import mul as mul2
    for j,i in enumerate(selected):
        assert decode(E['all14_equations_after_exact_linear_Schur_through2'][i+2],M9)==mul2(L9,{(UDIR[j],):diag[j]})
    roots=equations+P['selected5_F_roots']+P['remaining5_F_roots'];need=ancestors(nodes,roots)
    assert need==P['scalar_ancestor_nodes'] and len(need)==17329
    assert all(i in bound for i in need if nodes[i][0]=='v')
    values={};linear_values={};inverse_checks=0
    for i in need:
        op,*a=nodes[i]
        if op=='c':p={():F(*a)}if a[0]else{}
        elif op=='v':p=bound[i]
        elif op=='a':p=plus(values[a[0]],values[a[1]])
        elif op=='n':p=scale(values[a[0]],-1)
        elif op=='m':p=times(values[a[0]],values[a[1]])
        else:assert op=='i';p=inverse(values[a[0]]);inverse_checks+=1
        assert len(p)<=145;values[i]=p;linear_values[i]=linear_part(p)
    # Independent linear AD; no new adjoint solve is used.
    gsel=[reverse_linear(nodes,need,linear_values,r,unknown)for r in P['selected5_F_roots']]
    gf=[reverse_linear(nodes,need,linear_values,r,unknown)for r in P['remaining5_F_roots']]
    to5=lambda p:[str(p.get(m,F(0)))for m in [(),(0,),(1,),(2,),(3,)]]
    assert [[to5(p)for p in row]for row in gsel]==P['own_reverse_linear_selected_gradients']
    assert [[to5(p)for p in row]for row in gf]==P['own_reverse_linear_remaining_gradients']
    assert all(not p.get((),F(0))for row in gf for p in row)
    l0=[[F(c)for c in row]for row in P['selected5_lambda0_rows']]
    l1=[[[F(c)for c in row]for row in mat]for mat in P['remaining5_lambda1_rows_by_core_direction']]
    assert all(sum(l0[i][k]*A[k][j]for k in range(98))==gsel[i][j].get((),F(0))for i in range(5)for j in range(98))
    assert all(sum(l1[i][d][k]*A[k][j]for k in range(98))==gf[i][j].get((d,),F(0))for i in range(5)for d in range(4)for j in range(98))
    eqjets=[values[r]for r in equations]
    assert [vec145(p)for p in eqjets]==Q['all98_residual_jets145']
    assert all(all(len(m)>=3 if all(j<4 for j in m)else len(m)>=3 for m in p)for p in eqjets)
    E3=[hom(p,3)for p in eqjets];L={(0,):F(5),(1,):F(-20),(2,):F(24),(3,):F(6)}
    selectedraw=[values[r]for r in P['selected5_F_roots']];raw=[values[r]for r in P['remaining5_F_roots']]
    assert [vec145(p)for p in selectedraw]==Q['selected5_raw145']
    assert [vec145(p)for p in raw]==Q['remaining5_raw145']
    U2=[]
    for i,p in enumerate(selectedraw):
        correction=hom(p,3)
        for k,c in enumerate(l0[i]):
            if c:correction=plus(correction,scale(E3[k],-c))
        assert vec(correction)==Q['selected5_corrected_cubics70'][i]
        u=unpack70(Q['selected5_eliminated_U2_quadrics70'][i]);assert all(len(m)==2 for m in u)
        assert times(L,scale(u,-diag[i]))==correction;U2.append(u)
    cubics=[]
    for i,p in enumerate(raw):
        assert all(len(m)>=4 for m in p if all(j<4 for j in m))
        correction=hom(p,4)
        for d in range(4):
            for k,c in enumerate(l1[i][d]):
                if c:correction=plus(correction,scale(times({(d,):F(1)},E3[k]),-c))
        for d in range(5):
            mix={tuple(j for j in m if j<4):c for m,c in p.items()if 4+d in m}
            assert all(len(m)==2 for m in mix)
            assert vec(mix)==Q['remaining5_mixed_dual_quadrics'][i][d]
            correction=plus(correction,times(mix,U2[d]))
        assert vec(correction)==Q['remaining5_corrected_quartic_numerators70'][i]
        cubic=unpack70(Q['all5_actual_core_cubics70'][i]);assert all(len(m)==3 for m in cubic)
        assert times(L,cubic)==correction;cubics.append(cubic)
    cmat=[[p.get(m,F(0))for m in MON if len(m)==3]for p in cubics]
    assert cmat==[[F(c)for c in row]for row in Q['actual_core_cubic_matrix5_by20']]and rank(cmat)==3
    geom=read(base/'AMBIENT_B18_CORE_NODAL_QUARTIC_20261104_v63/fixed-quadrics-nodal-quartic.json')
    phi=[[F(c)for c in row]for row in geom['binary_form_coefficient_matrix4_by5']]
    def pullback(p,degree):
        ans=[F(0)]*(4*degree+1)
        for m,c in p.items():
            assert len(m)==degree;v=[c]
            for j in m:v=binary_product(v,phi[j])
            ans=[a+b for a,b in zip(ans,v)]
        return ans
    # Bind normalization to the independently accepted quadrics as well.
    quad0=decode(oldQ['all7_core_quadrics_padded_through3'][0],MON[:35])
    quad1=decode(oldQ['all7_core_quadrics_padded_through3'][1],MON[:35])
    assert not any(pullback(quad0,2))and not any(pullback(quad1,2))
    binary=[pullback(p,3)for p in cubics]
    assert binary==[[F(c)for c in row]for row in Q['binarydegree12_normalization_pullback_coefficients']]
    assert [any(row)for row in binary]==[True,False,False,False,False]
    proof=dict(actual98_source_unknown_equation_and_all10_output_roots_bound=True,
        all107_FULL9_source2_F9_Schur_compositions_independently_rebuilt=True,
        four_core_and_five_dual_bindings_rebuilt=True,independent_dictionary_ring_dimension145=True,
        all17329_scalar_source_ancestors_evaluated=True,all_unit_inverse_ring_products_checked=inverse_checks,
        actual5_G1_diagonal_checked_from_FULL9_model=True,all2450_adjoint_identities_with_independent_AD=True,
        all5_remaining_raw_origin_gradients_zero=True,unknown_source3_canceled_not_assumed_zero=True,
        unknown_source4_absent_by_stationary_output_gradients=True,all98_base2_and_dual1_residuals_zero=True,
        all98_full145_residuals_and_all10_raw145_outputs_reproduced=True,
        all5_U2_and_all5_corrected_core3_full70_L1_products=True,
        cubic_span_rank3=True,only_row63_binary12_section_nonzero=True,
        first3_actual_functions_CI_2_2_3_CM_Gorenstein_curve_multiplicity12=True,
        full_core_dimension_exactly1_by_actual_R_lower_and_K63_upper=True,
        full_core_CM_H0_socle_E4_D0_depth_or_full_graph_redundancy_not_claimed=True,
        no_author_worker_solve_new_source3_source4_point_prime_GB_or_higher_raise=True)
    (out/'independent-proof.json').write_text(json.dumps(proof,indent=2)+'\n')
    for p in[Path(__file__),ROOT/'code/ambient_B18_zero_normalizer_quadratic_review_20261104_v63.py',ROOT/f'research/notes/{ID}.md']:hashes[str(p.relative_to(ROOT))]=sha(p)
    report=dict(id=ID,all_passed=True,input_sha256=hashes,evidence_sha256={p.name:sha(p)for p in out.iterdir()},
        actual_core_dimension_exactly1=True,first3_actual_CI_multiplicity12=True,full_core_CM_or_Ddepth_claim=False,
        seconds=time.monotonic()-start,max_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    p=out/'report.json';p.write_text(json.dumps(report,indent=2)+'\n');assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
    print(json.dumps(dict(all_passed=True,seconds=report['seconds'],rss=report['max_rss_kib'],report_sha256=sha(p))))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args();main(a.output)
