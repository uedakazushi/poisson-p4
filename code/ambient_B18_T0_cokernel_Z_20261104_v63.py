#!/usr/bin/env python3
"""ONE exact finite cokernel Z=coker(D0->D4) over the proven A3 quotient.

This is not T0, and does not compute a kernel by naive base change.
"""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import resource
import signal
import time
from ambient_B18_T3_zero_structural_20261104_v63 import insert_basis

ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_T0_COKERNEL_Z_20261104_v63'
def read(p):return json.loads(p.read_text())
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for q in iter(lambda:f.read(1048576),b''):h.update(q)
    return h.hexdigest()
def add(a,b):return [x+y for x,y in zip(a,b)]
def mul(a,b):return [a[0]*b[0]]+[a[0]*b[j]+b[0]*a[j] for j in range(1,14)]
def inv(a):
    assert a[0]
    return [1/a[0]]+[-x/a[0]**2 for x in a[1:]]

def main(out):
    start=time.monotonic();signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
    out.mkdir(parents=True,exist_ok=False)
    cp=ROOT/'reports/AMBIENT_B18_T0_HONEST_C0_20261104_v63/honest-C0-D0-D4-DAG.json'
    cr=ROOT/'reports/AMBIENT_B18_T0_HONEST_C0_20261104_v63_REVIEW/report.json'
    sp=ROOT/'reports/AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63/actual-source-coefficients.json'
    fp=ROOT/'reports/AMBIENT_B18_FULL_RELATIVE_INPUT_20261104_v63/full-relative186-DAG.json'
    jp=ROOT/'reports/AMBIENT_B18_T3_FULL9_FIRSTJET_20261104_v63/full9-firstjet-certificate.json'
    qp=ROOT/'reports/AMBIENT_B18_T1_T2_FIRSTJET_20261104_v63/Q1-Q2-input-and-firstjets.json'
    ep=ROOT/'reports/AMBIENT_B18_T3_ZERO_STRUCTURAL_20261104_v63/independent-firstjet-and-structural-proof.json'
    gates=[cr,ROOT/'reports/AMBIENT_B18_T1_T2_QUADRATIC_PHI_20261104_v63_REVIEW/report.json',
           ROOT/'reports/AMBIENT_B18_T3_DOUBLE_ANNIHILATOR_LENGTH3_20261104_v63_REVIEW.json',
           ROOT/'reports/AMBIENT_B18_W1_TAIL_BINDING_20261104_v63_REVIEW/report.json']
    assert all(read(p)['all_passed'] for p in gates)
    c,source,full,jet,qdata,own=map(read,[cp,sp,fp,jp,qp,ep])
    nodes=c['DAG_scalar_nodes'];oldnodes=source['DAG_scalar_nodes'];nold=len(oldnodes)
    # Actual source first derivatives in the CA F9 coordinates, previously
    # independently verified through all saved quadratic coefficients.
    adapter=[list(map(F,a)) for a in jet['YB_to_YA_linear']]
    binding={}
    for name in ['CA','CB']:
        for i,a in jet['both57_source2'][name]['all57_source2jets']:
            lin=list(map(F,a[1:10]))
            if name=='CB':lin=[sum((lin[k]*adapter[k][j] for k in range(9)),F(0)) for j in range(9)]
            value=[F(0)]+lin+[F(0)]*4
            if i in binding:assert binding[i]==value
            binding[i]=value
    old=[]
    for i,n in enumerate(oldnodes):
        op=n[0]
        if op=='c':value=[F(n[1],n[2])]+[F(0)]*13
        elif op=='v':value=binding[i]
        elif op=='n':value=[-x for x in old[n[1]]]
        elif op=='a':value=add(old[n[1]],old[n[2]])
        elif op=='m':value=mul(old[n[1]],old[n[2]])
        elif op=='i':value=inv(old[n[1]])
        else:raise ValueError(op)
        old.append(value)
    known={int(i):list(map(F,a))+[F(0)]*4 for i,a in jet['all784_CB_h_firstjets_in_YA'].items()}
    known.update({int(i):list(map(F,a))+[F(0)]*4 for i,a in jet['all792_CA_h_firstjets'].items()})
    for j,node in enumerate(jet['four_g_derivatives_independent_node_ids']):
        a=[F(0)]*14;a[10+j]=1;known[node]=a
    known[full['CB_g_restriction_node']]=[F(0)]*14
    MON=list(map(tuple,jet['source2_monomials']));POS={m:i for i,m in enumerate(MON)}
    ca=jet['both57_source2']['CA']['pure_function2jets']
    # Tail prefixes keep r,s as well as earlier z-coordinates.  Their first
    # jets are fixed by the already saved source2 coefficients.
    for rec in c['all35_regular_coordinate_divisions']:
        p=list(map(F,ca[str(rec['original_pure_CA_source_function'])]))
        for j,q in enumerate(rec['all7_regular_coordinate_divisions']):
            a=[p[1+j]]+[p[POS[tuple(sorted((j,k)))]] if k<=j or k>=7 else F(0) for k in range(9)]
            assert a[0]==F(q['origin']);known[q['node']]=a+[F(0)]*4
    roots=[i for col in c['C0_free9_to_D4_honest9_columns'] for i in col]
    needed=set();stack=roots[:]
    while stack:
        i=stack.pop()
        if i in needed:continue
        needed.add(i)
        if i<nold or i in known:continue
        n=nodes[i];op=n[0]
        if op in ['a','m']:stack+=n[1:]
        elif op in ['n','i']:stack.append(n[1])
        elif op!='c':raise ValueError((i,n))
    assert not ({q['node'] for q in c['all4_known_regular_gA_divisions']} & needed)
    V={}
    for i in sorted(needed):
        if i<nold:a=old[i]
        elif i in known:a=known[i]
        else:
            n=nodes[i];op=n[0]
            if op=='c':a=[F(n[1],n[2])]+[F(0)]*13
            elif op=='n':a=[-x for x in V[n[1]]]
            elif op=='a':a=add(V[n[1]],V[n[2]])
            elif op=='m':a=mul(V[n[1]],V[n[2]])
            elif op=='i':a=inv(V[n[1]])
            else:raise ValueError(op)
        V[i]=a
    M=[]
    for col in c['C0_free9_to_D4_honest9_columns']:
        assert all(not any(V[i][10:]) for i in col)
        M.append([V[i][:10] for i in col])
    origin_basis=[]
    for col in M:insert_basis(origin_basis,[a[0] for a in col])
    assert len(origin_basis)==7
    # Legitimate RIGHT-exact base change of a cokernel, after the independent
    # all-order annihilator proof in the note.
    ids=qdata['families']['1']['honest9_column_source_positions_in_Q3_reduced']
    J=[[[F(own['all2069_rebuilt_exact_firstjet90_columns'][i][10*k+j]) for j in [0,8,9]] for k in range(9)] for i in ids]
    C=[[[a[j] for j in [0,8,9]] for a in col] for col in M]
    cols=J+C;assert len(cols)==1548+9
    def flat(col):return [x for a in col for x in a]
    def act(v,which):
        ans=[F(0)]*27
        for row in range(9):ans[3*row+which]=v[3*row]
        return ans
    basis=[];provenance=[]
    for j,col in enumerate(cols):
        a=flat(col)
        for label,v in [('1',a),('r',act(a,1)),('s',act(a,2))]:
            if insert_basis(basis,v) is not None:provenance.append([j,label])
    def reduce(v):
        a={i:x for i,x in enumerate(v) if x}
        for pivot,b in basis:
            if pivot in a:
                coef=a[pivot]
                for i,x in b.items():
                    q=a.get(i,F(0))-coef*x
                    if q:a[i]=q
                    else:a.pop(i,None)
        return a
    pivot={i for i,_ in basis};free=[i for i in range(27) if i not in pivot]
    dim=len(free);assert 2<=dim<=6
    # Verify action-stability of the whole relation subspace.
    for _,b in basis:
        a=[b.get(i,F(0)) for i in range(27)]
        assert not reduce(act(a,1)) and not reduce(act(a,2))
    actions=[]
    for j in [1,2]:
        matrix=[]
        for coord in free:
            a=[F(0)]*27;a[coord]=1;q=reduce(act(a,j))
            matrix.append([q.get(k,F(0)) for k in free])
        actions.append(matrix)
    def mm(A,B):
        # matrices stored as columns
        return [[sum((A[k][i]*B[j][k] for k in range(dim)),F(0)) for i in range(dim)] for j in range(dim)]
    ZER=[[F(0)]*dim for _ in range(dim)]
    assert all(mm(A,B)==ZER for A in actions for B in actions)
    image=[]
    for A in actions:
        for col in A:insert_basis(image,col)
    mu=dim-len(image);assert mu==2
    assert dim==3 and free==[21,24,26]
    # Identify the whole module, not just its Hilbert function.
    a,b=actions[0][0][2],actions[0][1][2]
    cc,d=actions[1][0][2],actions[1][1][2]
    pairdet=a*d-b*cc;assert pairdet==F(36,25)
    omega_map=[[F(0),a,cc],[F(0),b,d],[F(1),F(0),F(0)]]
    omega_r=[[F(0)]*3,[F(1),F(0),F(0)],[F(0)]*3]
    omega_s=[[F(0)]*3,[F(0)]*3,[F(1),F(0),F(0)]]
    assert mm(omega_map,actions[0])==mm(omega_r,omega_map)
    assert mm(omega_map,actions[1])==mm(omega_s,omega_map)
    # The actual quotient map from honest D4 free9 to Z is projection to E7,E8:
    # the first seven free generators are killed by literal C0 unit columns.
    for k,row in zip([0,1,2,3,4,6,7],range(7)):
        assert c['C0_free9_to_D4_honest9_columns'][k]==[int(i==row) for i in range(9)]
    projected=[];projected_pivots=[]
    for j,col in enumerate(own['all2069_rebuilt_exact_firstjet90_columns']):
        p=[F(col[10*row+k]) for row in [7,8] for k in range(1,10)]
        if insert_basis(projected,p) is not None:projected_pivots.append(j)
    assert len(projected)==16
    muB=17-16+9-2;assert muB==8
    # The two-dimensional kernel at the D4 origin loses exactly one direction
    # in B0's fiber.  Detect it in Tor1(Z)/imageTor1(D4).
    qg=[F(0)]*8+[F(1)]
    qzero=[F(0)]*5+[F(1),F(1,2),F(-1),F(0)]
    visible=[]
    for name,q in [('gA_tau',qg),('f60_plus_half_f67_minus_f72',qzero)]:
        assert all(sum((q[i]*M[i][row][0] for i in range(9)),F(0))==0 for row in range(9))
        p=[sum((q[i]*M[i][row][k] for i in range(9)),F(0)) for row in [7,8] for k in range(1,10)]
        temp=[(i,dict(v)) for i,v in projected]
        outside=insert_basis(temp,p) is not None
        assert outside==(name!='gA_tau')
        visible.append({'constant_C0_vector':list(map(str,q)),'label':name,
                        'projected_firstjet18':list(map(str,p)),
                        'nonzero_in_Tor1_Z_mod_Tor1_D4':outside})
    ca_source={i:list(map(F,a)) for i,a in jet['both57_source2']['CA']['all57_source2jets']}
    oldvars={n[1]:i for i,n in enumerate(oldnodes) if n[0]=='v'}
    physical_actions=[]
    physical_F9=[]
    for row in source['physical_zero_base_rows9']:
        p=ca_source[oldvars[f'x{row}']];a,b=p[8],p[9]
        physical_F9.append([str(x) for x in p[1:10]])
        A=[[a*x+b*y for x,y in zip(ca,cb)] for ca,cb in zip(actions[0],actions[1])]
        physical_actions.append({'row':row,'linear_r_s':[str(a),str(b)],'action_columns':[[str(x) for x in col] for col in A]})
    ev={'A3_coordinates':'CA F9: z1..z7=0, r=Y7,s=Y8; r^2=rs=s^2=0',
      'all9_C0_map_exact_full_firstjet_columns':[[list(map(str,a)) for a in col] for col in M],
      'all1548_D4_Q1_relation_columns_over_A3':[[list(map(str,a)) for a in col] for col in J],
      'all9_C0_map_columns_over_A3':[[list(map(str,a)) for a in col] for col in C],
      'QQ_relation_rank':len(basis),'QQ_relation_basis':[[i,[[j,str(a)] for j,a in sorted(b.items())]] for i,b in basis],
      'basis_provenance':provenance,'quotient_basis_free_coordinate_indices':free,
      'Z_length':dim,'Z_minimal_R0_generators':mu,
      'Z_is_canonical_module_of_A3':True,'omega_isomorphism_columns_basis_eps0_epsr_epss':[[str(x) for x in col] for col in omega_map],
      'r_s_pairing_determinant':str(pairdet),
      'omega_A3_R0_minimal_resolution_ranks':[2,17,64,140,196,182,112,44,10,1],
      'projected_J_full9_linear_rank':16,
      'all1548_Q1_projected_full9_linear_columns':
          [[str(own['all2069_rebuilt_exact_firstjet90_columns'][i][10*row+k]) for row in [7,8] for k in range(1,10)] for i in ids],
      'projected_J_full9_basis':[[i,[[j,str(a)] for j,a in sorted(v.items())]] for i,v in projected],
      'projected_J_full9_pivot_columns':projected_pivots,
      'Tor1_D4_to_Z_rank':16,'image_B0_minimal_generators':muB,
      'T0_nonzero_by_mu_D0_9_vs_mu_B0_8':True,
      'unique_visible_T0_class_in_D0_fiber':visible,
      'visible_gA_tau_is_a_leading_class_not_claimed_literal_kernel_function':True,
      'r_action_columns':[[str(x) for x in a] for a in actions[0]],'s_action_columns':[[str(x) for x in a] for a in actions[1]],
      'all9_actual_physical_R0_actions':physical_actions,
      'physical_R0_rows9':source['physical_zero_base_rows9'],
      'actual_physical_R0_firstjets_on_CA_F9_9by9':physical_F9,
      'all_actions_commute_and_all_maxideal_products_zero':True,
      'all4_unknown_old_g_first_derivative_sensitivities_drop_from_map':True,
      'new_CA_g_quotients_not_in_map_closure':True,
      'finite_module_is_cokernel_Z_not_T0':True,
      'T0_socle_or_depth_not_computed':True}
    evidence=out/'finite-cokernel-Z-and-actions.json';evidence.write_text(json.dumps(ev,separators=(',',':'))+'\n')
    paths=[Path(__file__),ROOT/'research/notes'/f'{STEM}.md',cp,sp,fp,jp,qp,ep,*gates,
      ROOT/'code/ambient_B18_T3_zero_structural_20261104_v63.py']
    report={'id':STEM,'all_passed':True,'status':'exact_finite_cokernel_Z_actual_A3_module_not_T0_solution',
      'Z_length':dim,'Z_minimal_generators':mu,'QQ_relation_rank_in27':len(basis),
      'Z_is_omega_A3_proved':True,'Tor1_D4_to_Z_rank':16,'actual_image_B0_minimal_generators':muB,
      'actual_T0_nonzero_proved_by_mu_comparison':True,
      'actual_both_zero_ideals_annihilate_Z_by_homogeneous_module_identity':True,
      'right_exact_cokernel_basechange_only_no_kernel_Tor_dropped':True,
      'all_actual_R0_commuting_actions_saved':True,
      'T0_finite_socle_or_full_depth_claimed':False,
      'no_new_source_IFT_jets_GB_points_primes_or_g_derivative_assumptions':True,
      'input_sha256':{str(p.relative_to(ROOT)):sha(p) for p in paths},
      'evidence_sha256':{evidence.name:sha(evidence)},'seconds':time.monotonic()-start,
      'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'certificate_bytes':evidence.stat().st_size}
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:report[k] for k in ['all_passed','Z_length','Z_minimal_generators','QQ_relation_rank_in27','seconds','max_rss_kib','certificate_bytes']}))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
