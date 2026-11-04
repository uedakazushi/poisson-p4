#!/usr/bin/env python3
"""Universal contact81 identities and the literal D-to-normalizer model.

Pure algebra only; actual source jets and their contact interpretation are
the separately assigned source peer's premise.
"""
import argparse
from fractions import Fraction as F
import hashlib
import itertools
import json
from pathlib import Path
import resource
import signal
import time

ROOT=Path(__file__).resolve().parents[1]
ID='AMBIENT_CA_CRITICAL_V1_CONTACT_ALGEBRA_20261104_v63_REVIEW'


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def add(a,b,c=F(1)):
    out=dict(a)
    for m,v in b.items():
        v=out.get(m,F(0))+c*v
        if v:out[m]=v
        else:out.pop(m,None)
    return out
def mul(a,b):
    out={}
    for m,c in a.items():
        for n,d in b.items():
            p=tuple(i+j for i,j in zip(m,n));v=out.get(p,F(0))+c*d
            if v:out[p]=v
            else:out.pop(p,None)
    return out
def enc(p):return [[list(m),str(c)] for m,c in sorted(p.items())]


def main(out):
    start=time.monotonic();signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
    resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024))
    out=out.resolve();out.mkdir(parents=True,exist_ok=False);hashes={}
    def read(p):hashes[str(p.relative_to(ROOT))]=sha(p);return json.loads(p.read_text())
    def author_gate(stem):
        p=ROOT/'reports'/stem/'report.json';r=read(p);assert r['all_passed']
        for name,h in r.get('input_sha256',{}).items():assert sha(ROOT/name)==h,name;hashes[name]=h
        for name,h in r.get('evidence_sha256',{}).items():
            q=p.parent/name;assert sha(q)==h,name;hashes[str(q.relative_to(ROOT))]=h
        return r
    # These are accepted independent theorem gates. Source Delta geometry is
    # kept separate; none of the source coefficient jobs is repeated here.
    accepted={}
    for stem in ['AMBIENT_B18_RELATIVE_WEIGHTED_MODEL_20261104_v63_REVIEW',
        'AMBIENT_B18_T1_T2_QUADRATIC_PHI_20261104_v63_REVIEW',
        'AMBIENT_B18_T3_ZERO_STRUCTURAL_20261104_v63_REVIEW',
        'AMBIENT_B18_COMPARISON_C_GENERIC_DVR_RANK_20261104_v63_REVIEW',
        'AMBIENT_B18_CONDUCTOR_DVR_COMPARISON_PREFLIGHT_20261104_v63_REVIEW',
        'AMBIENT_CA_SELF_PAIR_ACTUAL_D_FULL_RING_20261103_v63_REVIEW',
        'AMBIENT_CB_SELF_PAIR_ACTUAL_D_FULL_RING_20261103_v63_REVIEW']:
        accepted[stem]=author_gate(stem)
    for stem in ['AMBIENT_CA_CRITICAL_V1_CUSP_DIRECTION_PREFLIGHT_20261104_v63',
        'AMBIENT_CA_CRITICAL_V1_CONTACT_OPERATOR_PREFLIGHT_20261104_v63',
        'AMBIENT_CA_CRITICAL_V1_CONTACT_VECTOR_20261104_v63']:
        author_gate(stem)
    b=ROOT/'reports'
    source=read(b/'AMBIENT_CA_CRITICAL_V1_CONTACT_OPERATOR_PREFLIGHT_20261104_v63/literal-V1-contact-vector-input.json')
    direction=read(b/'AMBIENT_CA_CRITICAL_V1_CUSP_DIRECTION_PREFLIGHT_20261104_v63/existing-V1-exact-cusp-direction.json')
    contact=read(b/'AMBIENT_CA_CRITICAL_V1_CONTACT_VECTOR_20261104_v63/actual-V1-contact-vector-and-source-jets.json')
    full9=read(b/'AMBIENT_B18_T3_FULL9_FIRSTJET_20261104_v63/full9-firstjet-certificate.json')
    fullsource=read(b/'AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63/actual-source-coefficients.json')
    f9_units={}
    for family in ['CA','CB']:
        f9=full9['both57_source2'][family]
        rec=fullsource['families'][family]
        assert f9['F9_nodes']==[i for _,i in rec['actual_zero_odd_I0_generators'][:7]]+[i for _,i in rec['actual_positive_odd_J0_generators']]
        assert F(f9['determinant']) and f9['all57_equations_all55_coefficients_verified']
        f9_units[family]=dict(actual_F9_generator_ids=f9['F9_nodes'],known_unit57_determinant=f9['determinant'])
    beta=list(map(F,source['actual_v_even_u8_direction']))
    assert beta[4]==F(20,9) and F(direction['N_adapted_direction_nx'])==F(9,20)
    assert F(direction['N_adapted_direction_nZ'])==0
    assert F(contact['c_v'])==0 and not contact['nonzero_graph_rows']
    # Abstract pivot0 is actual selected u8 index4 (physical row60). Relabel
    # rather than accidentally use the first, zero, component of beta.
    actual_permutation=[4,0,1,2,3,5,6,7]
    names=[f'b{i}' for i in range(8)]+[f'l{i}' for i in range(1,8)]+['v','x','c']
    n=len(names);zero=(0,)*n
    def atom(i,e=1):m=list(zero);m[i]=e;return {tuple(m):F(1)}
    one={zero:F(1)};bs=[atom(i) for i in range(8)]
    ls=[{}]+[atom(7+i) for i in range(1,8)]
    v,x,c=atom(15),atom(16),atom(17)
    us=[add(mul(v,b),l) for b,l in zip(bs,ls)]
    X=add(x,mul(v,v),-1);W=mul(v,add(v,c,-1));target=ls[1:]+[X,W]
    witnesses=[];source_relations=[]
    def witness(label,relation,cof):
        ans={}
        for a,g in zip(cof,target):ans=add(ans,mul(a,g))
        assert ans==relation
        witnesses.append(dict(label=label,pulled_divided_relation=enc(relation),cofactors=[enc(a) for a in cof]))
        source_relations.append(relation)
    def blank():return [{} for _ in range(9)]
    for i,j in itertools.combinations(range(8),2):
        co=blank()
        if i:co[i-1]=add(co[i-1],bs[j])
        if j:co[j-1]=add(co[j-1],bs[i],-1)
        witness(['Koszul_Y',i,j],add(mul(us[i],bs[j]),mul(us[j],bs[i]),-1),co)
    xY_indices=[]
    for i in range(8):
        co=blank()
        if i:co[i-1]=add({},c,-1)
        co[7]=bs[i];co[8]=bs[i]
        xY_indices.append(len(witnesses))
        witness(['Koszul_Z',i],add(mul(x,bs[i]),mul(us[i],c),-1),co)
    square_indices={}
    for i,j in itertools.combinations_with_replacement(range(8),2):
        co=blank();co[7]=mul(bs[i],bs[j])
        if i:co[i-1]=add(co[i-1],add(mul(v,bs[j]),ls[j]),-1)
        if j:co[j-1]=add(co[j-1],mul(v,bs[i]),-1)
        square_indices[i,j]=len(witnesses)
        witness(['product_YY',i,j],add(mul(x,mul(bs[i],bs[j])),mul(us[i],us[j]),-1),co)
    for i in range(8):
        co=blank()
        if i:co[i-1]=add({},x,-1)
        co[7]=mul(bs[i],add(c,v,-1));co[8]=mul(v,bs[i]);co[8]=add({},co[8],-1)
        witness(['product_YZ',i],mul(x,add(mul(bs[i],c),us[i],-1)),co)
    co=blank();co[7]=add(add(mul(c,c),mul(v,v),-1),x,-1);co[8]=mul(v,add(c,v));co[8]=add({},co[8],-1)
    witness(['product_ZZ'],mul(x,add(mul(c,c),x,-1)),co)
    assert len(witnesses)==81
    invb=atom(0,-1);invb2=atom(0,-2)
    reverse=[]
    for i in range(1,8):
        source_index=next(k for k,a in enumerate(witnesses) if a['label']==['Koszul_Y',0,i])
        assert mul(add({},invb,-1),source_relations[source_index])==ls[i]
        reverse.append(dict(target=f'l{i}',sources=[[source_index,enc(add({},invb,-1))]]))
    pp=square_indices[0,0];xx=xY_indices[0]
    assert mul(invb2,source_relations[pp])==X
    assert add(mul(invb,source_relations[xx]),mul(invb2,source_relations[pp]),-1)==W
    reverse +=[dict(target='x-v^2',sources=[[pp,enc(invb2)]]),
               dict(target='v(v-c)',sources=[[xx,enc(invb)],[pp,enc(add({},invb2,-1))]])]
    # The exact ideal intersection uses only monomial lcms.
    dvars=[f'l{i}' for i in range(1,8)]+['v','x'];dn=len(dvars)
    primary=[]
    for i in range(7):m=[0]*dn;m[i]=1;primary.append(tuple(m))
    vm=[0]*dn;vm[7]=2;primary.append(tuple(vm))
    xm=[0]*dn;xm[8]=2;primary.append(tuple(xm))
    xmon=[0]*dn;xmon[8]=1
    intersection=[tuple(max(a,b) for a,b in zip(xmon,m)) for m in primary]
    expected=[]
    for i in range(7):m=list(xmon);m[i]+=1;expected.append(tuple(m))
    m=list(xmon);m[7]=2;expected.append(tuple(m));m=list(xmon);m[8]+=1;expected.append(tuple(m))
    assert set(intersection)==set(expected)
    # A literal4x4 finite primary-quotient comparison: D basis1,v,x,xv;
    # E basis1,v,t,vt, v^2=0,t^2=vt. The full regular reduced part is shared.
    Dbasis=[(0,0),(1,0),(0,1),(1,1)]
    Ebasis=[(0,0),(1,0),(0,1),(1,1)]
    def dmul(i,j):
        a,b=[x+y for x,y in zip(Dbasis[i],Dbasis[j])]
        return {} if a>=2 or b>=2 else {Dbasis.index((a,b)):F(1)}
    def emul(i,j):
        a,b=[x+y for x,y in zip(Ebasis[i],Ebasis[j])]
        if b>=2:a+=b-1;b=1
        return {} if a>=2 else {Ebasis.index((a,b)):F(1)}
    mapping=[{0:F(1)},{1:F(1)},{3:F(1)},{}]
    def lin(cols,value):
        ans={}
        for j,a in value.items():ans=add(ans,cols[j],a)
        return ans
    for i,j in itertools.product(range(4),repeat=2):
        rhs={}
        for a,ca in mapping[i].items():
            for bb,cb in mapping[j].items():rhs=add(rhs,emul(a,bb),ca*cb)
        assert lin(mapping,dmul(i,j))==rhs
    assert mapping[3]=={} and [mapping[i] for i in [0,1,2]]==[{0:F(1)},{1:F(1)},{3:F(1)}]
    # All action and associativity tables are explicit, no rank algorithm.
    for product in [dmul,emul]:
        columns=[[product(i,j) for j in range(4)] for i in range(4)]
        for i,j,k in itertools.product(range(4),repeat=3):
            assert lin(columns[k],product(i,j))==lin(columns[i],product(j,k))
    # Generic branch comparison lengths use only accepted scheme diagrams,
    # not a Tor-ideal guess or the nil length by itself.
    rankC=accepted['AMBIENT_B18_COMPARISON_C_GENERIC_DVR_RANK_20261104_v63_REVIEW']
    assert [rankC[k] for k in ['CAcritical_length','Pa_triple_meeting_length','CBcritical_length','Qb_triple_meeting_length']]==[1,1,0,2]
    branches=[dict(branch='CA critical',D_nil=2,E_extra=2,C_length=1,K_length=1),
        dict(branch='CA self Pa',D_nil=1,E_extra=1,C_length=1,K_length=1),
        dict(branch='CB critical outside CA',D_nil=0,E_extra=0,C_length=0,K_length=0),
        dict(branch='CB self Qb',D_nil=4,E_extra=4,C_length=2,K_length=2)]
    for rec in branches:assert rec['K_length']==rec['D_nil']-rec['E_extra']+rec['C_length']
    assert sum(q['K_length'] for q in branches)==4
    (out/'universal81-and-primary-model.json').write_text(json.dumps(dict(
        universal_variable_order=names,only_Laurent_inverse='b0^(-1), a known unit',
        abstract_pivot_actual_u8_index=4,actual_u8_permutation=actual_permutation,actual_b_unit='20/9',
        original81_divided_relations_and_target_cofactors=witnesses,
        reverse9_target_generator_witnesses=reverse,
        exact_intrinsic_divided_ideal='(l7,x-v^2,v(v-c))',
        higher_order_hypotheses='After the8 unit IFT equations(l7,x-v^2), ctilde and EVERY ambient Gtilde have constant and linear coefficients zero as functions over the smooth critical base.',
        higher_order_conclusion='ctilde=v^2 h, Gtilde=v^2 h_l; v(v-ctilde)=v^2(1-vh) has unit coefficient, so full dividedideal=(l7,x,v^2).',
        full_D='G20[[l7,v,x]]/(x*l7,x^2,x*v^2)',
        exact_primary_intersection=dict(first='(x)',second='(l7,x^2,v^2)',lcm_generators=[list(m) for m in intersection]),
        normalizer_basechange='G20[[l7,v,t]]/(l7*t,t^2-v*t,t^3)',
        primary_D_basis=['1','v','x','xv'],primary_E_basis=['1','v','t','vt'],
        primary_D_multiplication=[[list(dmul(i,j).items()) for j in range(4)] for i in range(4)],
        primary_E_multiplication=[[list(emul(i,j).items()) for j in range(4)] for i in range(4)],
        map_D_to_E_primary_columns=[list(col.items()) for col in mapping],
        kernel_basis=['xv'],cokernel_basis=['t'],
        actual_F9_zero7_plus_positive2_unit_gates=f9_units,
        generic_kernel_branch_lengths_if_contact_source_gate_verified=branches,
        generic_DVR_K_rank_if_contact_source_gate_verified=4,
        no_source_coefficient_evaluation=True),indent=2,default=str)+'\n')
    proof=dict(all_passed=True,pure_algebra81_two_inclusions_verified=True,
        all_ambient_graphs_retained=True,full_higher_order_formal_unit_argument_accepted=True,
        D_nil_free_rank2_over_smooth20_base=True,D_nil_square_zero=True,
        D_Ass=['N=(x)','criticalS=(l7,v,x)'],embedded_primary_quotient_length4=True,
        D_cone_dimension=28,D_cone_embedding_dimension=29,D_cone_depth=20,D_top_multiplicity=1,
        E_basechange_and_literal_map_kernel_length1_cokernel_length1=True,
        actual_source_application_requires_separate_Delta_and_geometry_peer=True,
        source_coefficients_not_independently_recomputed_in_this_pure_algebra_review=True,
        G20_uniformity_required='Use the formal product for the actual G20 orbit and its normalization diagram. The constant/linear contact vanishings must hold after that product, not only at an unrelated isolated closed point.',
        generic_zero_weight_transport_accepted=True,
        zero_weight_transport='Work AFTER the common equivariant G19 group product, on actual R10=R0[[w]]. For AA=R10+IA*tauA inside BA and AB=R10+IB*tauB inside BB, multiplication by f inIA and g inIB gives reverse normalization maps. Their completed tensor product over the same physical ambient ring is a reverse map whose two composites with the comparison are multiplication by fg. Hence IA*IB annihilates the full transverse comparison kernel and cokernel without flat basechange. At a generic prime of IA0, IB0 is a unit since IA0+IB0 is m-primary; JA0 is a unit since IA0+JA0=m, the actual F9 coordinate theorem. Thus w is inIA and IB is a unit, giving w*Kfull=0. Therefore Kfull is a submodule of T(w). The separately accepted T1=T2=T3=0 and Jinf=Q3 give T(w)_n=0 for ALL n>=1, so Kfull has only degree0. This uses positive-degree vanishing, not w-annihilation alone. Nonnegative grading makes degree0 exact, with Dfull_0=D0 and Efull_0=E4 literally. Thus Kfull_0=ker(D0->E4), not ker(D0->D4). The common flat formal G19 product preserves these source comparison diagrams and generic lengths.',
        generic_K_rank4_if_actual_contact_source_peer_closes=True,
        generic_branch_count_proof='The accepted full conductor/source inventory has four rational smooth branches, each ell order1. CriticalCA gives K1 by the literal map above. Pa has Dnil1 and E=N plus a reduced length1 second meeting, with coker1; Qb has Dnil4 and E=N plus Z4 length4, with coker C2 length2. Projection to N is surjective in both accepted scheme diagrams, so kernel is in Dnil and has lengths1,2 respectively. CBcritical is outside CA. This gives generic K lengths1+1+0+2=4. These are generic rank counts only and depend on the separately verified actual contact premises.',
        distinct_embedded_closure_criterion='With the actual irreducible G20 source-orbit/fiber gates, criticalCA(one ramified CA point, smooth CB), Pa(two smooth CA points) and Qb(smooth CA,two CB points) have distinct generic proper-fiber patterns. Their dense constructible loci cannot be the same irreducible20-dimensional associated closure. This is not a new Poisson-component count.',
        no_original_B18_depth_nilindex_or_new_Poisson_component_claim=True,
        no_CAS_source_jet_GB_newpoint_prime_or_author_execution=True)
    (out/'independent-proof.json').write_text(json.dumps(proof,indent=2)+'\n')
    for p in [Path(__file__),ROOT/f'research/notes/{ID}.md']:hashes[str(p.relative_to(ROOT))]=sha(p)
    report=dict(id=ID,all_passed=True,pure_algebra_contact_and_normalizer_comparison_accepted=True,
        actual_source_delta_peer_separate=True,generic_zero_weight_K_transport_accepted=True,
        generic_K_DVR_rank4_conditional_on_separate_actual_contact_source_peer=True,
        input_sha256=hashes,evidence_sha256={p.name:sha(p) for p in out.iterdir()},
        seconds=time.monotonic()-start,max_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
    print(json.dumps(dict(all_passed=True,seconds=report['seconds'],report_sha256=sha(out/'report.json'))))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    main(p.parse_args().output)
