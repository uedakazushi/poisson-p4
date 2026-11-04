#!/usr/bin/env python3
"""Independent Fraction AD, physical interface and finite Res4/Koszul audit."""
import argparse
from collections import defaultdict
from fractions import Fraction as F
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import resource
import signal
import time

ROOT=Path(__file__).resolve().parents[1]
ID='AMBIENT_B18_COMPLETE_NORMALIZER_SYZYGY_OPERATOR_20261104_v63_REVIEW'
FI='AMBIENT_B18_PROPER_PHYSICAL_INTERFACE_FIRSTJET_20261104_v63_REVIEW'


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def M(rows):return [[F(c) for c in row] for row in rows]
def Z(n,m):return [[F(0) for _ in range(m)] for _ in range(n)]
def tr(A):return [list(r) for r in zip(*A)]
def plus(A,B):return [[a+b for a,b in zip(r,s)] for r,s in zip(A,B)]
def neg(A):return [[-a for a in r] for r in A]
def mm(A,B):
    C=Z(len(A),len(B[0]))
    for i,row in enumerate(A):
        for k,a in enumerate(row):
            if a:
                for j,b in enumerate(B[k]):
                    if b:C[i][j]+=a*b
    return C
def take(A,rows,cols=None):return [[A[i][j] for j in (range(len(A[0])) if cols is None else cols)] for i in rows]
def em(A):return [[str(c) for c in r] for r in A]
def determinant(A):
    A=[r[:] for r in A];ans=F(1)
    for k in range(len(A)):
        j=next((j for j in range(k,len(A)) if A[j][k]),None)
        if j is None:return F(0)
        if j!=k:A[k],A[j]=A[j],A[k];ans=-ans
        p=A[k][k];ans*=p
        for i in range(k+1,len(A)):
            if not A[i][k]:continue
            c=A[i][k]/p
            for j in range(k+1,len(A)):A[i][j]-=c*A[k][j]
            A[i][k]=0
    return ans
def parity(p):return (-1)**sum(p[i]>p[j] for i in range(len(p)) for j in range(i+1,len(p)))


def main(out):
    start=time.monotonic();signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
    resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024))
    out=out.resolve();out.mkdir(parents=True,exist_ok=False)
    interface_out=ROOT/'reports'/FI;interface_out.mkdir(parents=True,exist_ok=False)
    hashes={}
    def read(p):
        hashes[str(p.relative_to(ROOT))]=sha(p)
        return json.loads(gzip.decompress(p.read_bytes()) if p.suffix=='.gz' else p.read_text())
    def gate(stem):
        p=ROOT/'reports'/stem/'report.json';r=read(p);assert r['all_passed']
        for name,h in r.get('input_sha256',{}).items():assert sha(ROOT/name)==h,name;hashes[name]=h
        for name,h in r.get('evidence_sha256',{}).items():
            path=p.parent/name;assert sha(path)==h,name;hashes[str(path.relative_to(ROOT))]=h
        return r
    # The accepted combined model review covers its preflight and builder.
    for stem in ['AMBIENT_B18_PROPER_POLYNOMIAL_MODEL',
        'AMBIENT_B18_PROPER_RATIONAL_ELIMINATION','AMBIENT_B18_PROPER_CORE_ELIMINATION',
        'AMBIENT_B18_PROPER_CORE_CI_MEMBERSHIP','AMBIENT_B18_ZERO_NORMALIZER_TENSOR_MODEL',
        'AMBIENT_B18_ZERO_NORMALIZER_QUADRATIC']:
        gate(stem+'_20261104_v63_REVIEW')
    for stem in ['AMBIENT_B18_PROPER_PHYSICAL_INTERFACE_FIRSTJET','AMBIENT_B18_COMPLETE_NORMALIZER_SYZYGY_OPERATOR']:
        gate(stem+'_20261104_v63')
    b=ROOT/'reports'
    I=read(b/'AMBIENT_B18_PROPER_PHYSICAL_INTERFACE_FIRSTJET_20261104_v63/actual-source-physical-firstjet-interface.json')
    L=read(b/'AMBIENT_B18_COMPLETE_NORMALIZER_SYZYGY_OPERATOR_20261104_v63/actual-complete-Res4-syzygy-operator.json')
    E=read(b/'AMBIENT_B18_ZERO_NORMALIZER_TENSOR_MODEL_20261104_v63/actual-normalizer-tensor-free4-DAG.json')
    P=read(b/'AMBIENT_B18_PROPER_POLYNOMIAL_MODEL_PREFLIGHT_20261104_v63/literal-proper-source14-preflight.json')
    A17=read(b/'AMBIENT_B18_PROPER_POLYNOMIAL_MODEL_20261104_v63/actual-proper-source14-equations17.json')
    H=read(b/'AMBIENT_B18_PROPER_RATIONAL_ELIMINATION_20261104_v63/actual12-flag-divided-polynomials.json')
    K=read(b/'AMBIENT_B18_PROPER_CORE_ELIMINATION_20261104_v63/actual4variable7core-polynomials.json')
    C=read(b/'AMBIENT_B18_PROPER_CORE_CI_MEMBERSHIP_20261104_v63/membership-and-source-cofactors.json')
    prep=read(b/'AMBIENT_B18_PROPER_CORE_CI_MEMBERSHIP_20261104_v63/preparation.json')
    R=read(b/'AMBIENT_B18_QUADRATIC_IMAGE_REDUCTION_20261103_v63/report.json')
    O=read(b/'AMBIENT_B18_ACTUAL_INTERSECTION_ORBIT_PRODUCT_20261103_v63/certificate.json.gz')
    S=read(b/'AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63/actual-source-coefficients.json')
    Q=read(b/'AMBIENT_B18_ZERO_NORMALIZER_QUADRATIC_20261104_v63/cached-source2-and-exact-E4-quadrics.json')
    nodes=E['DAG_scalar_nodes'];rnodes=E['R0_variable_nodes']
    oldmap=dict(E['original_node_to_pruned_node'])
    preflight=dict(source_DAG_nodes=len(nodes),two_source_unit_systems48=True,
        fixed_constant_determinants=[48,48,15,9,12],largest_new_polynomial_space='8 independent G scalars',
        source_interface_bytes=242121,typed_operator_bytes=291886,
        no_new_Hensel_solution_source_jet_GB_or_scalar_coefficient_expansion=True,
        limits=dict(seconds=60,address_bytes=2_000_000_000,output_bytes=32*1024*1024,expanded_terms=150000))
    (out/'preflight.json').write_text(json.dumps(preflight,indent=2)+'\n')

    # Fresh differentiation of the actual finite source prefix. No producer AD,
    # derivative matrix, inverse or source solver is imported.
    def differentiate(roots,inputs):
        pos={i:j for j,i in enumerate(inputs)};needed=set();todo=list(roots)
        while todo:
            i=todo.pop()
            if i in needed:continue
            needed.add(i);op=nodes[i][0]
            if op not in ['c','v']:todo.extend(nodes[i][1:])
        values={};grads={}
        for i in sorted(needed):
            op,*a=nodes[i]
            if op=='c':val=F(a[0],a[1]);g={}
            elif op=='v':val=F(0);g={pos[i]:F(1)} if i in pos else {}
            elif op=='n':val=-values[a[0]];g={j:-c for j,c in grads[a[0]].items()}
            elif op=='i':
                val=1/values[a[0]];g={j:-c*val*val for j,c in grads[a[0]].items()}
            elif op in ['a','m']:
                l,r=a;gl,gr=grads[l],grads[r]
                if op=='a':val=values[l]+values[r];g={j:gl.get(j,F(0))+gr.get(j,F(0)) for j in gl.keys()|gr.keys()}
                else:val=values[l]*values[r];g={j:gl.get(j,F(0))*values[r]+values[l]*gr.get(j,F(0)) for j in gl.keys()|gr.keys()}
                g={j:c for j,c in g.items() if c}
            else:raise AssertionError((i,op))
            values[i]=val;grads[i]=g
        return values,grads,len(needed)

    base=Z(9,9)
    for j in range(7):base[j+1][j]=1
    assert M(I['physical_R0_first9_by9'])==base
    orbit=M(O['physical150_by19_scalar_then_GL'])
    allgl=[tuple(p) for p in P['GL4_generator_positions16']]
    chosen=[tuple(p) for p in P['actual_G15_generator_positions']]
    orbits=M(P['all150_GL4_orbit16']);orb15=take(orbits,range(150),[allgl.index(p) for p in chosen])
    bad=[(1,0),(1,3),(1,4),(0,3),(0,4)]
    bad_orbit=take(orbits,range(150),[allgl.index(p) for p in bad])
    left=M(P['CA_left_inverse19']);selected19=[P['zero_coordinate_rows36'][j] for j in P['CA_constant_left_inverse_rows19']]
    assert mm(left,take(M(P['all150_CA_fixed_flag19']),selected19))==[[F(i==j) for j in range(19)] for i in range(19)]
    Jg=M(P['source_gauge15_flag5_then_fiber10_Jacobian']);assert determinant(Jg)==F(-2,3)
    physical=[];frames=[];flags=[];source_results=[];monic=[]
    for family_i,family in enumerate(['CA','CB']):
        sys=next(s for s in E['normalizer_sources48_each'] if s['family']==family)
        rec=next(s for s in I['families'] if s['family']==family)
        unknown=sys['unknown48'];eq=sys['equations48'];nu=E['nu_CA_node' if family=='CA' else 'nu_CB_node']
        vals,grad,ncount=differentiate(eq+[nu],unknown+rnodes)
        assert all(vals[e]==0 for e in eq)
        J=[[grad[e].get(j,F(0)) for j in range(48)] for e in eq]
        B=[[grad[e].get(48+j,F(0)) for j in range(9)] for e in eq]
        assert J==M(rec['actual_J48']) and B==M(rec['actual_input_rhs48_by9'])
        det=determinant(J);assert det==-F(sys['unit48_determinant'])==F(rec['actual_literal_order_J48_determinant'])
        cp=list(range(0,46,2))+[47]+list(range(1,46,2))+[46]
        rp=list(range(0,48,2))+list(range(1,48,2))
        assert cp==rec['block_column_permutation'] and rp==rec['block_row_permutation']
        assert parity(cp)*parity(rp)*det==F(sys['unit48_determinant'])
        sol=M(rec['solution48_by9']);assert plus(mm(J,sol),B)==Z(48,9)
        nd=[grad[nu].get(48+j,F(0))+sum(grad[nu].get(k,F(0))*sol[k][j] for k in range(48)) for j in range(9)]
        assert [-v for v in nd]+[F(0),F(0)]==list(map(F,Q['linear_Jacobian14_by11'][family_i]))
        monic.append([-v for v in nd])
        Dsol=mm(sol,base);up={n:j for j,n in enumerate(unknown)}
        source=S['families'][family];names=source['zero_hensel_unknowns23_order']
        first={n:Dsol[up[oldmap[source['coefficient_unknown_pairs'][n][0]]]] for n in names}
        mean=Dsol[up[oldmap[source['p_quadratic_relation']['s']]]]
        pfirst=[c/2+F(j==7+family_i) for j,c in enumerate(mean)]
        assert [pfirst]==M(rec['p_first'])
        gam=Z(19,9);k=Z(9,9)
        for name,values in first.items():
            if name.startswith('gamma'):gam[int(name[5:])]=values
            else:k[int(name[1:])]=values
        frame=[Z(5,5) for _ in range(9)]
        for n in range(1,19):
            i,j=divmod(O['GL18_generator_indices'][n-1],5)
            for z in range(9):frame[z][i][j]+=gam[n][z]
        core=M(R['families'][family]['source_kernel_Lie_matrix'])
        for z in range(9):frame[z]=plus(frame[z],[[pfirst[z]*c for c in row] for row in core])
        assert [em(m) for m in frame]==rec['source_frame5_by5_by9_first']
        pp=plus(mm(orbit,gam),mm(take(M(R['families'][family]['linear_J']),range(150),range(19,28)),k))
        assert pp==M(rec['physical150_by9_first'])
        assert take(pp,E['common_cut_zero_rows15'])==Z(15,9)
        assert take(pp,E['actual_R0_rows9'])==base
        flag=[[-frame[z][i][j] for z in range(9)] for i,j in bad]
        assert flag==M(rec['normalizer_flag5_by9_first'])
        physical.append(pp);frames.append(frame);flags.append(flag)
        source_results.append(dict(family=family,DAG_ancestor_count=ncount,
            all48_J_B_and_solution_equations_verified=True,literal_order_determinant=str(det),
            source_mean_shift_and_physical150_by9_reconstructed=True))
    assert determinant(take(monic,[0,1],[0,8]))==F(-20,9)
    assert mm(monic,base)==Z(2,9)
    before=mm(left,take(plus(physical[0],mm(bad_orbit,flags[0])),selected19))
    fixed=P['CA_fiber10_fixed_coefficient_indices'];free=P['CA_free_fiber9_coefficient_indices']
    constraints=flags[0]+take(before,fixed)
    gauge=M(I['group15_first9']);assert plus(mm(Jg,gauge),constraints)==Z(15,9)
    PA=plus(physical[0],mm(orb15,gauge));proper=mm(left,take(PA,selected19))
    assert take(proper,fixed)==Z(10,9)
    G=[Z(5,5) for _ in range(9)]
    for j,(i,k) in enumerate(chosen):
        for z in range(9):G[z][i][k]=gauge[j][z]
    relative=[[-frames[1][z][i][j]-G[z][i][j] for z in range(9)] for i,j in bad]
    actual14=take(proper,free)+relative
    assert actual14==M(I['source14_first_by_physical9'])
    retain=[9,12,2,3,4,5,6,7,8];J9=take(actual14,retain)
    assert retain==I['retained_source9_indices'] and J9==M(I['actual_source9_by_physical9_Jacobian'])
    assert determinant(J9)==F(5,48)==F(I['actual_source9_unit_determinant'])
    ann=M(P['CB_annihilator17_by36']);rows=P['zero_coordinate_rows36']
    targets=[q['physical_row'] for q in E['all12_zero_physical_difference_columns']]
    assert targets==[18,19,28,29,60,67,72,63,64,68,73,74]
    assert mm(ann,take(M(P['all150_CB_fixed_flag19']),rows))==Z(17,19)
    R17=take(ann,range(17),[rows.index(q) for q in targets])
    assert R17==M(I['source17_by_physicalF12_origin'])
    J17=M(A17['Jacobian17_by14']);piv=[0,1,10,11,13];sel=list(range(5));rest=list(range(5,17))
    JS=take(J17,sel,piv);subtract=M(I['source_remaining12_subtract_first5_origin'])
    assert JS==M(I['source_first5_pivot_Jacobian']) and determinant(JS)==F(8,3)
    assert mm(subtract,JS)==take(J17,rest,piv)
    RR=plus(take(R17,rest),neg(mm(subtract,take(R17,sel))))
    norms=[F(q['constant_denominator'],q['primitive_content_unit']) for q in H['all12_actual_remaining_and_divided_records']]
    R12=[[norms[i]*c for c in row] for i,row in enumerate(RR)]
    assert R12==M(I['source_dH12_by_physicalF12_origin'])
    assert determinant(R12)==F(687194767360000)==F(I['source_dH12_row_unit_determinant'])
    assert mm(J17,actual14)==Z(17,9)

    # Bind every actual source argument of the typed all-order recipe.
    typed=I['source_bound_finite_matrix_typed_operator'];ops={q['id']:q for q in typed['operations']}
    for path in typed['input_reference_files'].values():assert ROOT.joinpath(path).is_file();hashes[path]=sha(ROOT/path)
    assert ops['parameter_pairs']['unknown_and_equation_records']==E['normalizer_sources48_each']
    assert ops['G']['positions']==P['actual_G15_generator_positions']
    assert M(ops['gauge15']['origin_Jacobian15'])==Jg
    assert ops['actual14']['free9']==free and ops['kappa9']['indices']==retain
    assert ops['R17']['Ann_CB']==P['CB_annihilator17_by36'] and ops['R17']['F12_physical_rows']==targets
    assert ops['root5']['roots']==H['five_explicit_rational_roots']
    assert M(ops['A5_B12']['origin_A5'])==JS and M(ops['R12']['origin_matrix'])==R12
    assert [F(q['QQ_factor']) for q in ops['R12']['N_i']]==norms
    assert L['source_bound_interface']==typed
    T=L['source_adapter_T']
    assert T['source_H12_polynomials']==[q['divided_primitive_polynomial'] for q in H['all12_actual_remaining_and_divided_records']]
    assert T['second5_core_roots']==K['five_actual_rational_roots']
    assert T['graph_cofactors3_by4']==C['cofactor3_by4'] and T['CI_generators3']==prep['generators3']
    assert T['permutation']==[0,4,1,2,3,5,6,7,8,9,10,11]
    lam=T['residual_lambda']
    assert lam==[dict(QQ=str(F(q['QQ_unit_content_numerator'],q['QQ_unit_content_denominator'])),
        D_denominator_power=q['D_denominator_power'],P_denominator_power=q['P_denominator_power']) for q in K['remaining7_actual_core_polynomials']]
    lam0=[F(q['QQ'])/5**q['P_denominator_power'] for q in lam]
    detT=F(2,5)/F(10*8*80*16*192)
    for c in lam0:detT/=c
    assert detT==F(5**10,2359296)
    assert T['T']=='E4*E3*E2*E1*Pi' and T['full_N']=='T(kappa)*R12'
    assert T['full_N_inverse']=='R12^-1*T(kappa)^-1' and T['d_never_inverted']
    assert L['actual_d2']['shape']==[48,128] and L['actual_d3']['shape']==[128,224]
    assert L['actual_d2']['exact_block_recipe']=='block(i,j)=mu(sum_k N[k,i]*S[k,j]); native12rowblock then4basis columns'
    odd=[row[7] for row in gauge]
    assert odd==list(map(F,[-1,0,0,2,0,0,-2,0,-2,0,0,-1,0,-4,0]))
    gr=L['group30_to120']
    assert gr['auxiliary_etaA_distinct_from_physical_tauA']
    assert list(map(F,gr['group30_odd_origin']))==odd and gr['group30_even_origin']==[0]*15
    assert M(gr['group120_origin'])==Z(15,4)+[[q,0,0,0] for q in odd]
    assert F(gr['original_group30_det_origin'])==F(-2,3)**2
    assert F(gr['restricted120_det_origin'])==F(-2,3)**8==F(256,6561)

    # Independent monomial algebra: basis multiplication, actual Q columns,
    # and both Koszul products. All dictionaries are small free QQ-polynomials.
    def add(a,b):
        c=dict(a)
        for m,x in b.items():
            c[m]=c.get(m,F(0))+x
            if not c[m]:del c[m]
        return c
    def mul(a,b):
        c=defaultdict(F)
        for p,x in a.items():
            for q,y in b.items():c[tuple(i+j for i,j in zip(p,q))]+=x*y
        return {m:x for m,x in c.items() if x}
    def bp(i,j):
        A=(i%2)+(j%2);B=(i//2)+(j//2)
        return (A%2)+2*(B%2),(A//2,B//2)
    expected=[dict(left=i,right=j,output=bp(i,j)[0],nuA_power=bp(i,j)[1][0],nuB_power=bp(i,j)[1][1]) for i,j in itertools.product(range(4),repeat=2)]
    assert expected==L['structure16']
    for i,j,k in itertools.product(range(4),repeat=3):
        a,p=bp(i,j);a,q=bp(a,k);b,r=bp(j,k);b,s=bp(i,b)
        assert a==b and tuple(x+y for x,y in zip(p,q))==tuple(x+y for x,y in zip(r,s))
    mus=[]
    for i in range(4):
        A=[[{} for _ in range(4)] for _ in range(4)]
        for j in range(4):k,p=bp(i,j);A[k][j]={p:F(1)}
        mus.append(A)
    for a,b in itertools.product(range(4),repeat=2):
        c,p=bp(a,b)
        for i,j in itertools.product(range(4),repeat=2):
            val={}
            for k in range(4):val=add(val,mul(mus[a][i][k],mus[b][k][j]))
            assert val==mul({p:F(1)},mus[c][i][j])
    atoms=[]
    for q in E['all12_zero_physical_difference_columns']:
        for a in q['column4']:
            if a not in [0,1] and a not in atoms:atoms.append(a)
    nus=[E['nu_CA_node'],E['nu_CB_node']]
    for a in nus:
        if a not in atoms:atoms.append(a)
    apos={a:j for j,a in enumerate(atoms)};zero=(0,)*len(atoms);cache={}
    def atom(i):
        if i in cache:return cache[i]
        if i in apos:
            m=list(zero);m[apos[i]]=1;v={tuple(m):F(1)}
        else:
            op,*q=nodes[i]
            if op=='c':v={zero:F(q[0],q[1])} if q[0] else {}
            elif op=='a':v=add(atom(q[0]),atom(q[1]))
            elif op=='m':v=mul(atom(q[0]),atom(q[1]))
            elif op=='n':v={m:-c for m,c in atom(q[0]).items()}
            else:raise AssertionError((i,op))
        cache[i]=v;return v
    assert L['actual_Q4_by48_column_node_ids']==E['normalizer_tensor_relation_columns4_by48']
    for i,rec in enumerate(E['all12_zero_physical_difference_columns']):
        for j in range(4):
            target=[{} for _ in range(4)]
            for k in range(4):
                row,p=bp(k,j);c=atom(rec['column4'][k])
                for nu,count in zip(nus,p):
                    for _ in range(count):c=mul(c,atom(nu))
                target[row]=add(target[row],c)
            assert [atom(k) for k in E['normalizer_tensor_relation_columns4_by48'][4*i+j]]==target
    pairs=list(itertools.combinations(range(8),2));triples=list(itertools.combinations(range(8),3))
    k2=[]
    for c,(i,j) in enumerate(pairs):k2.extend([[i,c,j,1],[j,c,i,-1]])
    k3=[]
    for c,(i,j,k) in enumerate(triples):
        k3.extend([[pairs.index((j,k)),c,i,1],[pairs.index((i,k)),c,j,-1],[pairs.index((i,j)),c,k,1]])
    assert k2==L['native_S12_by32']['koszul28_in_first8_rows']
    assert L['native_S12_by32']['free4_zero_source_columns']==[[8+i,28+i,1] for i in range(4)]
    assert k3==L['native_V32_by56']['koszul3_in_first28_rows']
    zeros=(0,)*8
    def g(i,c=1):m=list(zeros);m[i]=1;return {tuple(m):F(c)}
    A=[[{} for _ in range(28)] for _ in range(8)]
    B=[[{} for _ in range(56)] for _ in range(28)]
    for i,j,k,c in k2:A[i][j]=g(k,c)
    for i,j,k,c in k3:B[i][j]=g(k,c)
    for j in range(28):
        val={}
        for i in range(8):val=add(val,mul(g(i),A[i][j]))
        assert not val
    for i,j in itertools.product(range(8),range(56)):
        val={}
        for k in range(28):val=add(val,mul(A[i][k],B[k][j]))
        assert not val
    ranks=[4,48,128,224,280,224,112,32,4]
    assert L['restricted_resolution_ranks']==ranks

    physical_proof=dict(all_passed=True,arithmetic='Independent Fraction AD and constant matrix multiplication/determinants; saved derivative solutions verified, never re-solved.',
        both_source48_details=source_results,monic2_minor='-20/9',
        all_source_frame_means_and_actual150_derivatives_checked=True,
        actual_gauge15_determinant='-2/3',all15_times9_gauge_equations=True,
        all14_source_readbacks_and_kappa9_reconstructed=True,actual_kappa9_determinant='5/48',
        source_R17_annihilator_and_first5_Schur_exact=True,actual_R12_determinant='687194767360000',
        allorder_source_row_identity='dH12(kappa)=R12 F12',
        typed_source_operator_arguments_and_exact_root_signs_bound=True,
        allorder_proof='The own CA flags/fiber readbacks define the actual unit group gauge. The relative CB frame annihilates its own fixed source plane; its action on Delta is linear in the twelve retained physical differences. Finite polynomial telescoping at the explicit five roots gives R12; its checked constant determinant is a unit. The checked kappa9 minor gives an actual Cohen-coordinate isomorphism, not only an abstract quotient-type comparison.',
        scalar_series_coefficients_expanded=False,new_source_solve_jet_or_GB=False)
    (interface_out/'independent-proof.json').write_text(json.dumps(physical_proof,indent=2)+'\n')
    proof=dict(all_passed=True,actual_interface_independently_verified=True,
        all_H12_core_roots_and_CI_cofactor_arrays_bound=True,detT_origin=str(detT),
        N_orientation='N=T(kappa)R12, NF=(dG,0)^t; old syzygies=N^t S',
        sixteen_basis_products_and64_associativity_and256_matrix_entries=True,
        all48_actual_Qcolumns192entries_bound=True,Koszul56_and168_sparse_entries_bound=True,
        GK2_zero_columns=28,K2K3_zero_entries=448,
        group_auxiliary_eta_distinct_and120_origins_checked=True,group120_determinant='256/6561',
        completeness='C is regular9 by the actual monic and kappa unit minors. G8 is regular by the accepted actual CI source presentation. d is nonzero in the regular domain C, hence ker(dG)=ker(G), without inverting d. Koszul exactness plus four zero source rows gives the entire relation module S=K2 directsum I4 and next kernel [K3;0]. N is unit, so N^t S is complete in the original physical frame. Restriction from R0-free rank4 C is exact and gives the stated full 48x128 and128x224 operators.',
        restricted_resolution_ranks=ranks,not_claimed_minimal=True,
        no_evaluated_R0_scalar_functions_comparison_kernel_torsion_socle_D0_depth_or_old_Hadamard_lifts=True)
    (out/'independent-proof.json').write_text(json.dumps(proof,indent=2)+'\n')
    for p in [Path(__file__),ROOT/f'research/notes/{ID}.md',ROOT/f'research/notes/{FI}.md']:
        hashes[str(p.relative_to(ROOT))]=sha(p)
    common=dict(all_passed=True,input_sha256=hashes,
        no_author_producer_solver_or_new_source_arithmetic=True,
        seconds=time.monotonic()-start,max_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    first=dict(common,id=FI,actual_source9_unit_minor='5/48',actual_relation12_unit_minor='687194767360000',
        evidence_sha256={p.name:sha(p) for p in interface_out.iterdir()})
    (interface_out/'report.json').write_text(json.dumps(first,indent=2)+'\n')
    final=dict(common,id=ID,complete_actual_SyzQ_48_by128=True,complete_next128_by224=True,
        physical_interface_review_sha256=sha(interface_out/'report.json'),
        physical_interface_review_path=str((interface_out/'report.json').relative_to(ROOT)),
        evidence_sha256={p.name:sha(p) for p in out.iterdir()})
    (out/'report.json').write_text(json.dumps(final,indent=2)+'\n')
    assert sum(p.stat().st_size for p in out.iterdir())+sum(p.stat().st_size for p in interface_out.iterdir())<32*1024*1024
    print(json.dumps(dict(all_passed=True,seconds=common['seconds'],interface_review_sha256=sha(interface_out/'report.json'),operator_review_sha256=sha(out/'report.json'))))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    main(parser.parse_args().output)
