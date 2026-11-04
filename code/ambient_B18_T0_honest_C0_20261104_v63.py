#!/usr/bin/env python3
"""ONE honest C0=R0+(z7,gA)tau free9/Koszul28, plus612 Q0 lifts."""
import argparse
from fractions import Fraction as F
import hashlib,json
from pathlib import Path
import resource,signal,time
from ambient_B18_actual_top_coefficient_preflight_20261104_v63 import DAG
from ambient_B18_top_graph_telescoping_20261104_v63 import ancestors,independent_rows
from ambient_B18_CB_g_division_origins_20261104_v63 import evaluate,MON
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_T0_HONEST_C0_20261104_v63'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main(out):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));out.mkdir(parents=True,exist_ok=False)
    sp=ROOT/'reports/AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63/actual-source-coefficients.json'
    fp=ROOT/'reports/AMBIENT_B18_FULL_RELATIVE_INPUT_20261104_v63/full-relative186-DAG.json'
    hp=ROOT/'reports/AMBIENT_B18_T3_HONEST10_PRESENTATION_20261104_v63/honest10-presentation-DAG.json'
    up=ROOT/'reports/AMBIENT_B18_T3_HONEST10_UNIT_REDUCTION_20261104_v63/unit-reduced-presentation-DAG.json'
    op=ROOT/'reports/AMBIENT_B18_T3_DOUBLE_ANNIHILATOR_PREFLIGHT_20261104_v63/exact-residues-and-two-ideals.json'
    jp=ROOT/'reports/AMBIENT_B18_T3_JOINT_CURVE_RESTRICTIONS_20261104_v63/joint-degree3-certificate.json'
    firstp=ROOT/'reports/AMBIENT_B18_T3_FULL9_FIRSTJET_20261104_v63/full9-firstjet-certificate.json'
    old,full,honest,unit,origin,joint,first=map(read,[sp,fp,hp,up,op,jp,firstp])
    oldnodes=old['DAG_scalar_nodes'];ca=old['families']['CA'];var={n[1]:i for i,n in enumerate(oldnodes) if n[0]=='v'}
    x=[var[f'x{r}'] for r in old['physical_zero_base_rows9']]
    az=[i for p in ca['coefficient_unknown_pairs'].values() for i in p]+[ca['p_quadratic_relation']['s'],ca['p_quadratic_relation']['t']]
    fca=[i for r,i in ca['actual_zero_odd_I0_generators']];gens=ca['actual_zero_odd_I0_generators']+ca['actual_positive_odd_J0_generators']
    sel=[0,1,2,3,4,5,6,8,9];Y=[gens[i][1] for i in sel];f63=fca[7]
    rows=[q['raw_row'] for q in full['all99_CB_kernel_and87_graph_generators_evaluated_on_CA'] if q['kind']=='graph_CB_minus_CA' and q['weight']==0]
    assert rows==[64,68,73,74]
    Bnodes=[ca['complete_normalized_physical_w_coefficients'][0][r][1] for r in rows]
    # Reuse the already-computed combined105 degree3 solution; no new IFT/jet.
    sol={i:tuple(map(F,v)) for i,v in joint['all105_unknown_jets']}
    needjets=sorted(ancestors(oldnodes,[f63]+Bnodes));jj=evaluate(oldnodes,needjets,sol)
    gj=jj[f63];assert list(map(str,gj))==joint['g_CA_all10_coefficients']
    assert all(v==0 for m,v in zip(MON,gj) if sum(m)<3);lead=next(i for i,(m,v) in enumerate(zip(MON,gj)) if sum(m)==3 and v)
    H0=[]
    for b in Bnodes:
        bj=jj[b];assert all(v==0 for m,v in zip(MON,bj) if sum(m)<3)
        h=bj[lead]/gj[lead];assert bj==tuple(h*v for v in gj);H0.append(h)
    D=DAG();D.nodes=[tuple(n) for n in unit['DAG_scalar_nodes']];D.values=list(map(F,origin['all339669_scalar_origin_values_evaluated']))
    rvals={q['node']:F(q['origin']) for a in honest['all_pure_source_Hadamard_membership_telescopes'] for q in a['all9_regular_coordinate_divisions']}
    for i in range(len(D.values),len(D.nodes)):
        n=D.nodes[i];o=n[0]
        if o=='c':v=F(n[1],n[2])
        elif o=='v':v=F(0)
        elif o=='a':v=D.values[n[1]]+D.values[n[2]]
        elif o=='m':v=D.values[n[1]]*D.values[n[2]]
        elif o=='n':v=-D.values[n[1]]
        elif o=='i':assert D.values[n[1]];v=1/D.values[n[1]]
        elif o=='r':v=rvals[i]
        else:raise ValueError((i,n))
        D.values.append(v)
    D.cache={n:i for i,n in enumerate(D.nodes)}
    eq48=ca['all48_zero_IFT_quadratic_coefficient_equations'];needed=sorted(ancestors(oldnodes,eq48+Y+[f63]+Bnodes))
    assert {i for i in needed if oldnodes[i][0]=='v'}<=set(az+x)
    def substitute(binding):
        m={}
        for i in needed:
            n=oldnodes[i];o=n[0]
            if o=='c':v=D.const(F(n[1],n[2]))
            elif o=='v':v=binding[i]
            elif o=='a':v=D.add(m[n[1]],m[n[2]])
            elif o=='m':v=D.mul(m[n[1]],m[n[2]])
            elif o=='n':v=D.neg(m[n[1]])
            elif o=='i':v=D.inverse(m[n[1]])
            else:raise ValueError(o)
            m[i]=v
        return m
    systems=[];prefix=[];det57=first['both57_source2']['CA']['determinant']
    for k in range(7):
        binding={i:D.var(f'C0_CA_TAIL_PREFIX{k}_{oldnodes[i][1]}') for i in az+x};m=substitute(binding);prefix.append(m)
        equations=[m[i] for i in eq48]+[D.add(m[i],D.neg(Y[j] if j<k or j>=7 else D.zero)) for j,i in enumerate(Y)]
        assert len(equations)==57 and all(D.values[i]==0 for i in equations)
        systems.append({'retained_zero_prefix_length':k,'bindings':[[i,j] for i,j in binding.items()],'all57_unit_equations':equations,'unit57_determinant':det57})
    functions=[f63]+Bnodes;tel=[];coeffs=[]
    pureCA=first['both57_source2']['CA']['pure_function2jets']
    for root in functions:
        previous=prefix[0][root];hs=[];record=[]
        for i in range(7):
            current=prefix[i+1][root] if i<6 else root;num=D.add(current,D.neg(previous));value=F(pureCA[str(root)][1+i])
            h=D.put(('r',num,Y[i],f'C0_source{root}:zero{i}'),value);hs.append(h)
            record.append({'node':h,'numerator':num,'coordinate':Y[i],'current':current,'previous':previous,'origin':str(value)});previous=current
        assert previous==root
        coeffs.append(hs);tel.append({'original_pure_CA_source_function':root,'restriction_r_s_node':prefix[0][root],'all7_regular_coordinate_divisions':record})
    gA=prefix[0][f63];e=coeffs[0];H=[];gdiv=[]
    for row,b,h0 in zip(rows,Bnodes,H0):
        h=D.put(('g',prefix[0][b],gA,f'C0_zero_graph{row}:known_regular_zero_image'),h0);H.append(h)
        gdiv.append({'node':h,'numerator':prefix[0][b],'denominator':gA,'origin':str(h0),
                     'proof':'actualzeroimageR0+I0tau givesb in(z7,f63); restrictiongivesbbar in(gA); gA nonzero'})
    basis=Y[:7]+[gA];membership={D.zero:(D.zero,)*8};proof={D.zero:['zero']}
    for i,z in enumerate(Y[:7]):v=[D.zero]*8;v[i]=D.one;membership[z]=tuple(v);proof[z]=['coordinate',i]
    membership[f63]=tuple(e+[D.one]);proof[f63]=['source_f63_tail_identity',0]
    for j,b in enumerate(Bnodes):membership[b]=tuple(coeffs[1+j]+[H[j]]);proof[b]=['zero_graph_tail_plus_known_g_division',j]
    for root in list(membership):
        n=D.nodes[root]
        if n[0]=='n' and n[1] not in membership:membership[n[1]]=tuple(D.neg(v) for v in membership[root]);proof[n[1]]=['negative_of_certified',root]
    def get(i):
        if i in membership:return membership[i]
        n=D.nodes[i];o=n[0];a=None;why=None
        if o=='n':
            v=get(n[1])
            if v is not None:a=tuple(D.neg(x) for x in v);why=['neg',n[1]]
        elif o=='a':
            v,w=get(n[1]),get(n[2])
            if v is not None and w is not None:a=tuple(D.add(x,y) for x,y in zip(v,w));why=['sum',n[1],n[2]]
        elif o=='m':
            v=get(n[1])
            if v is not None:a=tuple(D.mul(n[2],x) for x in v);why=['product_known_left',n[1],n[2]]
            else:
                w=get(n[2])
                if w is not None:a=tuple(D.mul(n[1],x) for x in w);why=['product_known_right',n[1],n[2]]
        membership[i]=a
        if why:proof[i]=why
        return a
    Q0=full['actual_Qn_R0_squared_generator_matrices']['0'];assert len(Q0)==612;lift=[]
    for a,b in Q0:
        v=get(b);assert v is not None,('missing_C0_membership',b,D.nodes[b]);lift.append([a]+list(v))
    koszul=[]
    for i in range(8):
        for j in range(i+1,8):v=[D.zero]*9;v[1+i]=D.neg(basis[j]);v[1+j]=basis[i];koszul.append(v)
    assert len(koszul)==28
    H63=next(r['membership_vector9'] for r in honest['all_pure_source_Hadamard_membership_telescopes'] if r['original_source_function']==f63)
    c_to_h10=[[D.one]+[D.zero]*9]
    for i in range(7):v=[D.zero]*10;v[1+i]=D.one;c_to_h10.append(v)
    v=[D.zero]+H63[:]
    for i in range(7):v[1+i]=D.add(v[1+i],D.neg(e[i]))
    c_to_h10.append(v)
    P=unit['exact_projection_matrix_remaining_by_original10']
    toD4=[[D.total(D.mul(P[i][j],v[j]) for j in range(10)) for i in range(9)] for v in c_to_h10]
    originrank=len(independent_rows([[D.values[i] for i in v] for v in koszul+lift]));maprank=len(independent_rows([[D.values[i] for i in v] for v in toD4]))
    ev={'DAG_scalar_nodes':D.nodes,'actual_C0_basis8_scalar_nodes':basis,'actual_C0_free9_images_in_CA_trace_pairs':[[D.one,D.zero]]+[[D.zero,z] for z in basis],
        'all7_pure_CA_tail_prefix57_systems':systems,'all35_regular_coordinate_divisions':tel,'all4_known_regular_gA_divisions':gdiv,
        'cached_actual_CA_gA_all_degree3_coefficients':list(map(str,gj)),
        'cached_actual_CA_four_bbar_all_degree3_coefficients':[list(map(str,jj[b])) for b in Bnodes],
        'all4_gA_quotient_origin_values':list(map(str,H0)),'regular_sequence_proof':'z7arecoordinates,gA nonzeroorder3modz7,so8regularsequence',
        'all28_C0_Koszul_kernel_columns':koszul,'all612_Q0_lifted_columns':lift,'honest_D0_relation_columns':koszul+lift,
        'C0_odd_membership_propagation_rules':{str(i):r for i,r in proof.items()},
        'C0_free9_to_H_free10_columns':c_to_h10,'C0_free9_to_D4_honest9_columns':toD4,
        'D0_relation_origin_rank':originrank,'C0_to_D4_origin_rank':maprank,
        'all_old_CB_r_g_are_untouched_original_scalar_cofactors':True,
        'T0':'kernel(coker(honestD0relations)->D4) viaexplicit9by9map; kernelnotcomputed',
        'no_T0_nonzero_depth_or_support_only_nilcurve_claim':True}
    encoded=json.dumps(ev,separators=(',',':'))+'\n';assert len(encoded.encode())<=32*1024*1024
    (out/'honest-C0-D0-D4-DAG.json').write_text(encoded)
    paths=[Path(__file__),sp,fp,hp,up,op,jp,firstp]
    rep={'id':STEM,'all_passed':True,'status':'actual_honest_C0_D0_presentation_and_D4_map_not_T0_solution',
         'C0_free_generators':9,'complete_Koszul_kernel_columns':28,'Q0_actual_lifted_columns':612,'D0_relation_columns':640,
         'new_pure_CA_prefix57_systems':7,'new_regular_coordinate_quotients':35,'new_known_regular_gA_quotients':4,
         'four_new_gA_quotient_origin_values':list(map(str,H0)),'D0_relation_origin_rank':originrank,'C0_to_D4_origin_rank':maprank,
         'T0_kernel_socle_depth_or_support_classification_solved':False,'no_new_source_IFT_solve_jets_GB_points_or_old_nonunit_pullback':True,
         'input_sha256':{str(p.relative_to(ROOT)):sha(p) for p in paths},
         'evidence_sha256':{'honest-C0-D0-D4-DAG.json':sha(out/'honest-C0-D0-D4-DAG.json')},
         'DAG_nodes':len(D.nodes),'certificate_bytes':len(encoded.encode()),'seconds':time.monotonic()-start,
         'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    (out/'report.json').write_text(json.dumps(rep,indent=2)+'\n');print(json.dumps({k:rep[k] for k in ['all_passed','four_new_gA_quotient_origin_values','D0_relation_origin_rank','C0_to_D4_origin_rank','DAG_nodes','certificate_bytes','seconds','max_rss_kib']}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
