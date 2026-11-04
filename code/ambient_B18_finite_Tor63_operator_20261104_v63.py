#!/usr/bin/env python3
"""ONE all-order typed operator construction of the27 finite Tor2 image maps."""
import argparse,hashlib,itertools,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
from ambient_B18_actual_top_coefficient_preflight_20261104_v63 import DAG
from ambient_B18_top_graph_telescoping_review_20261104_v63 import closure
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_FINITE_TOR63_OPERATOR_20261104_v63'
def read(p):return json.loads(p.read_text())
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()
def enc(a):
    if isinstance(a,F):return str(a)
    if isinstance(a,(list,tuple)):return [enc(x) for x in a]
    if isinstance(a,dict):return {str(k):enc(v) for k,v in a.items()}
    return a
def addp(a,b,sign=1):
    c=dict(a)
    for m,v in b.items():
        x=c.get(m,0)+sign*v
        if x:c[m]=x
        else:c.pop(m,None)
    return c
def mulp(a,b):
    c={}
    for m,x in a.items():
        for n,y in b.items():
            k=tuple(sorted(m+n));c[k]=c.get(k,0)+x*y
    return {m:v for m,v in c.items() if v}
def universal_transform_check(selected,redundant):
    z={};one={():1};counter=0
    def var():
        nonlocal counter
        p={(counter,):1};counter+=1;return p
    U=[[var() for j in range(3)] for i in range(3)]
    H=[[var() for j in range(7)] for i in range(3)]
    Ar=[[var() for j in range(10)] for i in range(6)]
    def det3(M):
        p={}
        for perm in itertools.permutations(range(3)):
            q=one
            for i,j in enumerate(perm):q=mulp(q,M[i][j])
            p=addp(p,q,(-1)**sum(perm[i]>perm[j] for i in range(3) for j in range(i+1,3)))
        return p
    det=det3(U);adj=[]
    for i in range(3):
        row=[]
        for j in range(3):
            rr=[r for r in range(3) if r!=j];cc=[c for c in range(3) if c!=i]
            q=addp(mulp(U[rr[0]][cc[0]],U[rr[1]][cc[1]]),mulp(U[rr[0]][cc[1]],U[rr[1]][cc[0]]),-1)
            row.append({m:(-1)**(i+j)*v for m,v in q.items()})
        adj.append(row)
    A=[[z for j in range(10)] for i in range(16)]
    for i in range(7):A[i][i]=one
    for k,i in enumerate(selected[7:]):A[i]=H[k]+U[k]
    for k,i in enumerate(redundant):A[i]=Ar[k]
    T=[[z for j in range(16)] for i in range(10)]
    for i in range(7):T[i][i]=det
    for i in range(3):
        for j in range(7):
            q={}
            for k in range(3):q=addp(q,mulp(adj[i][k],H[k][j]),-1)
            T[7+i][j]=q
        for k,j in enumerate(selected[7:]):T[7+i][j]=adj[i][k]
    B=T[:]
    for k,original in enumerate(redundant):
        row=[]
        for j in range(16):
            q=det if j==original else z
            for i in range(10):q=addp(q,mulp(Ar[k][i],T[i][j]),-1)
            row.append(q)
        B.append(row)
    Inv=[A[i]+[one if i==j else z for j in redundant] for i in range(16)]
    checked=0
    for L,R in [(B,Inv),(Inv,B)]:
        for i in range(16):
            for j in range(16):
                q={}
                for k in range(16):q=addp(q,mulp(L[i][k],R[k][j]))
                assert q==(det if i==j else z)
                checked+=1
    return {'free_coefficient_symbols':counter,'left_and_right_inverse_entries_verified':checked,
      'identity':'B_numerator*A_inverse=det(U)*I and A_inverse*B_numerator=det(U)*I',
      'no_source_coefficient_expansion':True}
def main(out):
    start=time.monotonic();signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False)
    sp=ROOT/'reports/AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63/actual-source-coefficients.json'
    jp=ROOT/'reports/AMBIENT_B18_T3_FULL9_FIRSTJET_20261104_v63/full9-firstjet-certificate.json'
    cp=ROOT/'reports/AMBIENT_B18_T3_JOINT_CURVE_RESTRICTIONS_20261104_v63/joint-degree3-certificate.json'
    tp=ROOT/'reports/AMBIENT_B18_FINITE_TOR63_IMAGE_PREFLIGHT_20261104_v63/finite63-proof-and-input-preflight.json'
    gates=[tp.parent/'report.json',ROOT/'reports/AMBIENT_B18_ZERO_NORMALIZER_TENSOR_MODEL_20261104_v63_REVIEW/report.json',
      ROOT/'reports/AMBIENT_B18_T0_HONEST_C0_20261104_v63_REVIEW/report.json',
      ROOT/'reports/AMBIENT_B18_T3_DOUBLE_ANNIHILATOR_LENGTH3_20261104_v63_REVIEW.json']
    assert all(read(p)['all_passed'] for p in gates)
    S,J,C,T=map(read,[sp,jp,cp,tp]);old=S['DAG_scalar_nodes'];ca,cb=[S['families'][f] for f in ['CA','CB']]
    f=[p[1] for p in ca['actual_zero_odd_I0_generators']];g=[p[1] for p in cb['actual_zero_odd_I0_generators']]
    rows=[p[0] for p in ca['actual_zero_odd_I0_generators']]
    Y=J['both57_source2']['CA']['F9_nodes'];assert len(Y)==9
    eqCA=ca['all48_zero_IFT_quadratic_coefficient_equations'];eqCB=cb['all48_zero_IFT_quadratic_coefficient_equations']
    proto=C['all105_equation_original_node_order'];us=C['all105_unknown_original_node_order']
    assert proto==eqCA+Y+eqCB and len(us)==105
    roots=proto+[f[7]]+g;needed=sorted(closure(old,roots))
    assert {i for i in needed if old[i][0]=='v'}<=set(us)
    D=DAG();D.nodes=[tuple(n) for n in old];D.cache={n:i for i,n in enumerate(D.nodes)}
    assert len(D.cache)==len(D.nodes)
    vals=[]
    for n in D.nodes:
        op=n[0]
        if op=='c':v=F(n[1],n[2])
        elif op=='v':v=F(0)
        elif op=='a':v=vals[n[1]]+vals[n[2]]
        elif op=='m':v=vals[n[1]]*vals[n[2]]
        elif op=='n':v=-vals[n[1]]
        elif op=='i':assert vals[n[1]];v=1/vals[n[1]]
        else:raise ValueError(op)
        vals.append(v)
    D.values=vals;D.zero=0;D.one=1
    systems=[];maps=[];functions=[f[7]]+g
    specs=[('prefix'+str(k),[Y[j] if j<k or j>=7 else D.zero for j in range(9)]) for k in range(7)]
    specs += [('axis_r',[D.zero]*7+[Y[7],D.zero]),('axis_s',[D.zero]*8+[Y[8]])]
    for name,targets in specs:
        bind={i:D.var(f'TOR63_{name}_source_node{i}') for i in us};m={}
        for i in needed:
            n=old[i];o=n[0]
            if o=='v':v=bind[i]
            elif o=='c':v=D.const(F(n[1],n[2]))
            elif o=='a':v=D.add(m[n[1]],m[n[2]])
            elif o=='m':v=D.mul(m[n[1]],m[n[2]])
            elif o=='n':v=D.neg(m[n[1]])
            elif o=='i':v=D.inverse(m[n[1]])
            else:raise ValueError(o)
            m[i]=v
        equations=[m[i] for i in eqCA]+[D.add(m[i],D.neg(t)) for i,t in zip(Y,targets)]+[m[i] for i in eqCB]
        assert len(equations)==105 and all(D.values[i]==0 for i in equations)
        systems.append({'name':name,'source105_bindings':sorted(bind.items()),'all105_unit_equations':equations,
          'F9_target_nodes':targets,'unit_Jacobian_minor':C['det_joint105'],
          'old_to_clone_node_map':[[i,m[i]] for i in needed]})
        maps.append(m)
        checkpoint={'DAG_scalar_nodes':D.nodes,'completed_unit_systems':systems,'pure_source_needed_node_ids':needed}
        (out/'prefix-operator-checkpoint.json').write_text(json.dumps(checkpoint,separators=(',',':'))+'\n')
        assert (out/'prefix-operator-checkpoint.json').stat().st_size<32*1024*1024
    regular=[]
    def quotient(num,den,origin,why):
        assert D.values[num]==0 and D.values[den]==0
        i=D.put(('q',num,den,f'TOR63_regular_{len(regular)}'),F(origin))
        regular.append({'node':i,'numerator':num,'denominator':den,'origin':str(F(origin)),'all_order_regularity_proof':why})
        return i
    M=[[F(v) for v in row] for row in J['YB_to_YA_linear']]
    fa=[F(v) for v in J['both57_source2']['CA']['pure_function2jets'][str(f[7])][1:10]]
    gb8=[F(v) for v in J['both57_source2']['CB']['pure_function2jets'][str(g[7])][1:10]]
    Gfirst=M[:7]+[[sum((gb8[j]*M[j][i] for j in range(9)),F(0)) for i in range(9)]]
    grad=[fa]+Gfirst
    H=[];telescopes=[]
    for n,root in enumerate(functions):
        h=[];record=[]
        for j in range(7):
            before=maps[j][root];after=maps[j+1][root] if j<6 else root
            num=D.add(after,D.neg(before))
            q=quotient(num,Y[j],grad[n][j],{'type':'coordinate_prefix_difference','function_source_node':root,
              'before_prefix':j,'after_prefix':j+1,'coordinate_index':j,
              'identity_by_pointed_unit_uniqueness':'after|Yj=0 equals before; all other targets coincide'})
            h.append(q);record.append([before,after,q,Y[j]])
        H.append(h);telescopes.append({'function_source_node':root,'plane_node':maps[0][root],
          'full_node':root,'all7_difference_steps':record,'identity':'full=plane+sum z_j*coordinate_quotient_j'})
    planeU=[];planeidentities=[];r,s=Y[7],Y[8]
    source_quad=[[F(0)]*3]+[[F(v) for v in p['all_degree_le3_coefficients'][3:6]] for p in C['all8_CB_zero_restrictions']]
    assert not any(map(F,C['g_CA_all10_coefficients'][:6]))
    for n,root in enumerate(functions):
        p,pr,ps=maps[0][root],maps[7][root],maps[8][root]
        nums=[pr,D.add(D.add(p,D.neg(pr)),D.neg(ps)),ps]
        dens=[D.mul(r,r),D.mul(r,s),D.mul(s,s)]
        proofs=[{'type':'axis_order2','axis':'s=0','source_constant_and_linear_zero_from_saved105':True},
          {'type':'mixed_axis_difference','vanishes_on_r_axis_and_s_axis':True,'regular_ring_intersection':'(r) intersect (s)=(r*s)'},
          {'type':'axis_order2','axis':'r=0','source_constant_and_linear_zero_from_saved105':True}]
        u=[quotient(a,b,c,proof) for a,b,c,proof in zip(nums,dens,source_quad[n],proofs)]
        planeU.append(u);planeidentities.append({'function_source_node':root,'plane_node':p,'axis_r_node':pr,'axis_s_node':ps,
          'quadratic_Hadamard_coefficients':u,'identity':'plane=r^2*U0+r*s*U1+s^2*U2',
          'numerator_identity':'axis_r+(plane-axis_r-axis_s)+axis_s=plane'})
    selected=T['selected10_combined_indices'];redundant=T['redundant6_combined_indices'];assert selected==[0,1,2,3,4,5,6,8,12,13]
    A=[]
    for i in range(7):A.append([D.one if i==j else D.zero for j in range(7)]+[D.zero]*3)
    A.append(H[0]+planeU[0]);A += [H[i+1]+planeU[i+1] for i in range(8)]
    U=[A[i][7:] for i in selected[7:]]
    det=D.zero
    for perm in itertools.permutations(range(3)):
        q=D.one
        for i,j in enumerate(perm):q=D.mul(q,U[i][j])
        if sum(perm[i]>perm[j] for i in range(3) for j in range(i+1,3))%2:q=D.neg(q)
        det=D.add(det,q)
    assert D.values[det]==F(625,768)
    invdet=D.inverse(det);Ui=[]
    for i in range(3):
        row=[]
        for j in range(3):
            rr=[k for k in range(3) if k!=j];cc=[k for k in range(3) if k!=i]
            q=D.add(D.mul(U[rr[0]][cc[0]],U[rr[1]][cc[1]]),D.neg(D.mul(U[rr[0]][cc[1]],U[rr[1]][cc[0]])))
            if (i+j)%2:q=D.neg(q)
            row.append(D.mul(invdet,q))
        Ui.append(row)
    Z=D.zero;O=D.one
    transform10=[[Z]*16 for i in range(10)]
    for i in range(7):transform10[i][i]=O
    for i in range(3):
        for j in range(7):transform10[7+i][j]=D.neg(D.total(D.mul(Ui[i][k],A[selected[7+k]][j]) for k in range(3)))
        for k,gindex in enumerate(selected[7:]):transform10[7+i][gindex]=Ui[i][k]
    qcycles=[]
    for k in redundant:
        qcycles.append([D.add(O if j==k else Z,D.neg(D.total(D.mul(A[k][i],transform10[i][j]) for i in range(10)))) for j in range(16)])
    assert [[D.values[i] for i in q] for q in qcycles]==[[F(x) for x in row] for row in T['all6_zero_column_cycle_origins16']]
    eta_u=[D.add(D.neg(D.mul(s,transform10[7][j])),D.mul(r,transform10[8][j])) for j in range(16)]
    eta_v=[D.add(D.neg(D.mul(s,transform10[8][j])),D.mul(r,transform10[9][j])) for j in range(16)]
    assert all(D.values[i]==0 for i in eta_u+eta_v)
    transform=transform10+qcycles
    inverse=[A[i]+[O if i==k else Z for k in redundant] for i in range(16)]
    generic_proof=universal_transform_check(selected,redundant)
    Xi=[]
    for i in range(7):Xi.append([O if j==i+1 else Z for j in range(9)])
    Xi.append([Z]+H[0]+[O])
    rawca=ca['complete_normalized_physical_w_coefficients'][0];rawcb=cb['complete_normalized_physical_w_coefficients'][0]
    delta=[D.add(rawca[row][0],D.neg(rawcb[row][0])) for row in rows]
    nu=ca['actual_tracefree_nu'];images=[]
    specs=[('q_wedge_q',a,b,qcycles[a],qcycles[b]) for a,b in itertools.combinations(range(6),2)]
    specs += [(name,k,None,eta,qcycles[k]) for name,eta in [('eta_u_wedge_q',eta_u),('eta_v_wedge_q',eta_v)] for k in range(6)]
    for kind,a,b,q,p in specs:
        h=[[D.add(D.mul(q[i],p[8+j]),D.neg(D.mul(p[i],q[8+j]))) for j in range(8)] for i in range(8)]
        even=D.mul(nu,D.total(D.mul(h[i][j],D.mul(f[i],f[j])) for i in range(8) for j in range(8)))
        vector=[even]+[D.total(D.mul(h[i][j],D.mul(delta[j],Xi[i][k])) for i in range(8) for j in range(8)) for k in range(1,9)]
        assert all(D.values[i]==0 for i in vector)
        images.append({'type':kind,'first_cycle_index':a,'second_cycle_index':b,'h8_by8':h,
          'actual_C0_free9_image_nodes':vector,
          'formula_identity':'product of the two partial-syzygy target forms, expanded in CA trace via tauA^2=nuA'})
    assert len(images)==27 and len(regular)==90
    data={'DAG_scalar_nodes':D.nodes,'scalar_node_semantics_extension_q':'known regular quotient with exact numerator=denominator*q and saved all-order proof',
      'actual_common_R0_rows9':S['physical_zero_base_rows9'],'source_original_node_prefix_length':len(old),
      'all96_base_normalizer_unit_unknowns':[[i for p in rec['coefficient_unknown_pairs'].values() for i in p]+[rec['p_quadratic_relation'][k] for k in ['s','t']] for rec in [ca,cb]],
      'all96_base_normalizer_unit_equations':eqCA+eqCB,'all9_new_unit105_systems':systems,
      'all90_known_regular_Hadamard_quotients':regular,'all9_coordinate_telescopes':telescopes,
      'all9_second_order_Hadamard_identities':planeidentities,
      'canonical10_functions':Y[:7]+[D.mul(r,r),D.mul(r,s),D.mul(s,s)],
      'combined_actual16_generator_nodes':f+g,'actual16_by10_generator_expression_matrix':A,
      'selected10_indices':selected,'redundant6_indices':redundant,
      'actual_quadratic3_unit_matrix':U,'actual_unit3_determinant_node':det,'actual_unit3_inverse':Ui,
      'canonical10_from_actual16_matrix':transform10,'actual16_by16_unit_transform':transform,
      'actual16_by16_inverse_transform':inverse,'universal_exact_matrix_inverse_proof':generic_proof,
      'six_actual_zero_column_cycles16':qcycles,'eta_u_actual_cycle16':eta_u,'eta_v_actual_cycle16':eta_v,
      'whole_syzygy_proof':'G=A*h and h=T*G; q=e_red-A_red*T and eta=T^t*(-s,r,0)/(0,-s,r)',
      'C0_basis_odd8_nodes':Y[:7]+[maps[0][f[7]]],'C0_actual_CA8_odd_target_free9_vectors':Xi,
      'C0_basis_identification':'CA block of prefix105 is the existing pointed pureCA57 prefix; uniqueness identifies old gA and all seven coordinate quotients',
      'all8_actual_even_difference_nodes':delta,'all27_finite_Tor2_image_maps_in_C0_free9':images,
      'exact_image_is_m_R0_squared_annihilated_by_Tor2_theorem':True,
      'all_operator_coefficients_preserved_as_regular_functions_not_evaluated_or_expanded':True,
      'old_r_g_functions_used_or_evaluated':False,'image_nonzero_zero_socle_or_depth_solved':False}
    ep=out/'literal-Tor63-transform-and27-image-operators.json';txt=json.dumps(enc(data),separators=(',',':'))+'\n'
    assert len(txt.encode())<32*1024*1024;ep.write_text(txt)
    assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
    paths=[Path(__file__),ROOT/'research/notes'/f'{STEM}.md',sp,jp,cp,tp,*gates,
      ROOT/'code/ambient_B18_actual_top_coefficient_preflight_20261104_v63.py',ROOT/'code/ambient_B18_top_graph_telescoping_review_20261104_v63.py']
    report={'id':STEM,'all_passed':True,'status':'all_order_typed_Tor2_transform_and27_actual_C0_image_operators_only',
      'unit105_systems_constructed':9,'known_regular_Hadamard_quotients':90,
      'exact_16_generator_unit_transform_and_inverse_constructed':True,'generic_inverse_identity_entries_checked':512,
      'unit3_determinant_at_origin':'625/768','actual_C0_free9_Tor2_image_maps':27,
      'source_binding_by_unit_clone_and_regular_identity_not_coefficient_expansion':True,
      'source_jets_old_r_g_evaluation_Gaussian_or_nonmembership_used':False,
      'actual_finite_image_nonzero_or_zero_or_full_D_depth_proved':False,
      'input_sha256':{str(p.relative_to(ROOT)):sha(p) for p in paths},
      'evidence_sha256':{p.name:sha(p) for p in out.iterdir()},
      'DAG_nodes':len(D.nodes),'total_bytes':sum(p.stat().st_size for p in out.iterdir()),
      'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:report[k] for k in ['all_passed','DAG_nodes','total_bytes','seconds','max_rss_kib']}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
