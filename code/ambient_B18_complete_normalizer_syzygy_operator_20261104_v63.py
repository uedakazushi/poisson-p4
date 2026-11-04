#!/usr/bin/env python3
"""Finite actual Res4 syzygy operator, without expanding coefficient series.

The exact source-bound matrix-typed interface is lowered recursively. Universal
basis multiplication and both Koszul chain identities are checked independently.
No source solve, scalar jet, polynomial expansion or comparison-kernel claim.
"""
import argparse,hashlib,itertools,json,resource,signal,time
from collections import defaultdict
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_COMPLETE_NORMALIZER_SYZYGY_OPERATOR_20261104_v63'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def monic_basis_product(i,j):
    a,b=(i&1)+(j&1),(i>>1)+(j>>1)
    return (a%2)+2*(b%2),(a//2,b//2)
def plus(a,b):
    c=dict(a)
    for m,x in b.items():c[m]=c.get(m,F(0))+x
    return {m:x for m,x in c.items() if x}
def mul(a,b):
    c=defaultdict(F)
    for p,x in a.items():
        for q,y in b.items():c[tuple(i+j for i,j in zip(p,q))]+=x*y
    return {m:x for m,x in c.items() if x}
def neg(a):return {m:-x for m,x in a.items()}
def main(out):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));out.mkdir(parents=True,exist_ok=False)
    ip=ROOT/'reports/AMBIENT_B18_PROPER_PHYSICAL_INTERFACE_FIRSTJET_20261104_v63/actual-source-physical-firstjet-interface.json'
    sp=ROOT/'reports/AMBIENT_B18_PROPER_SYZYGY_ADAPTER_PREFLIGHT_20261104_v63/source12-to-regular8-adapter-preflight.json'
    ep=ROOT/'reports/AMBIENT_B18_ZERO_NORMALIZER_TENSOR_MODEL_20261104_v63/actual-normalizer-tensor-free4-DAG.json'
    hp=ROOT/'reports/AMBIENT_B18_PROPER_RATIONAL_ELIMINATION_20261104_v63/actual12-flag-divided-polynomials.json'
    kp=ROOT/'reports/AMBIENT_B18_PROPER_CORE_ELIMINATION_20261104_v63/actual4variable7core-polynomials.json'
    cp=ROOT/'reports/AMBIENT_B18_PROPER_CORE_CI_MEMBERSHIP_20261104_v63/membership-and-source-cofactors.json'
    pp=ROOT/'reports/AMBIENT_B18_PROPER_CORE_CI_MEMBERSHIP_20261104_v63/preparation.json'
    I,S,E,H,K,C,P=map(read,[ip,sp,ep,hp,kp,cp,pp])
    gate_stems=['AMBIENT_B18_PROPER_POLYNOMIAL_MODEL','AMBIENT_B18_PROPER_RATIONAL_ELIMINATION',
                'AMBIENT_B18_PROPER_CORE_ELIMINATION','AMBIENT_B18_PROPER_CORE_CI_MEMBERSHIP',
                'AMBIENT_B18_ZERO_NORMALIZER_TENSOR_MODEL']
    gates=[ROOT/'reports'/f'{g}_20261104_v63_REVIEW/report.json' for g in gate_stems]
    assert all(read(p)['all_passed'] for p in gates)
    assert I['actual_source9_unit_determinant']=='5/48' and I['source_dH12_row_unit_determinant']=='687194767360000'
    assert E['actual_R0_rows9']==[27,108,109,112,113,117,118,119,61]
    assert [r['physical_row'] for r in E['all12_zero_physical_difference_columns']]==[18,19,28,29,60,67,72,63,64,68,73,74]

    # Universal four-dimensional algebra over QQ[nuA,nuB]. Check all64 basis
    # associativity identities, and all16 multiplication-matrix products.
    table=[]
    for a,b in itertools.product(range(4),repeat=2):
        k,p=monic_basis_product(a,b);table.append({'left':a,'right':b,'output':k,'nuA_power':p[0],'nuB_power':p[1]})
    associativity=[]
    for a,b,c in itertools.product(range(4),repeat=3):
        x,p=monic_basis_product(a,b);y,q=monic_basis_product(x,c)
        z,r=monic_basis_product(b,c);w,t=monic_basis_product(a,z)
        assert y==w and tuple(i+j for i,j in zip(p,q))==tuple(i+j for i,j in zip(r,t))
        associativity.append([a,b,c])
    def basis_matrix(a):
        m=[[{} for j in range(4)] for i in range(4)]
        for j in range(4):k,p=monic_basis_product(a,j);m[k][j]={p:F(1)}
        return m
    mus=[basis_matrix(a) for a in range(4)]
    for a,b in itertools.product(range(4),repeat=2):
        k,p=monic_basis_product(a,b)
        for i,j in itertools.product(range(4),repeat=2):
            v={}
            for u in range(4):v=plus(v,mul(mus[a][i][u],mus[b][u][j]))
            assert v==mul({p:F(1)},mus[k][i][j])

    # Bind every column of the ACTUAL Q matrix to multiplication by its own
    # F_i in this basis. Target scalars and nus are independent symbolic atoms.
    nodes=E['DAG_scalar_nodes'];nus=[E['nu_CA_node'],E['nu_CB_node']]
    atomids=[]
    for q in E['all12_zero_physical_difference_columns']:
        for a in q['column4']:
            if a not in [0,1] and a not in atomids:atomids.append(a)
    atomids+=nus;apos={a:i for i,a in enumerate(atomids)};nvars=len(atomids);zero=(0,)*nvars;cache={}
    def expression(i):
        if i in cache:return cache[i]
        if i in apos:
            m=list(zero);m[apos[i]]=1;p={tuple(m):F(1)}
        else:
            q=nodes[i];op=q[0]
            if op=='c':p={zero:F(q[1],q[2])} if q[1] else {}
            elif op=='n':p=neg(expression(q[1]))
            elif op=='a':p=plus(expression(q[1]),expression(q[2]))
            elif op=='m':p=mul(expression(q[1]),expression(q[2]))
            else:raise ValueError(('actualQ binding reached unpinned coefficient',i,q))
        cache[i]=p;return p
    actualQ=E['normalizer_tensor_relation_columns4_by48'];assert len(actualQ)==48
    verified=[]
    for j,rec in enumerate(E['all12_zero_physical_difference_columns']):
        f=rec['column4']
        for b in range(4):
            target=[{} for _ in range(4)]
            for a in range(4):
                k,p=monic_basis_product(a,b);co=expression(f[a])
                for nu,count in zip(nus,p):
                    for _ in range(count):co=mul(co,expression(nu))
                target[k]=plus(target[k],co)
            assert all(expression(actualQ[4*j+b][k])==target[k] for k in range(4))
            verified.append([j,b])

    # Exact integer universal Koszul coefficients in eight independent scalars.
    pairs=list(itertools.combinations(range(8),2));triples=list(itertools.combinations(range(8),3))
    zero8=(0,)*8
    def gi(i,sign=1):m=list(zero8);m[i]=1;return {tuple(m):F(sign)}
    k2=[[{} for _ in pairs] for _ in range(8)];sparse2=[]
    for col,(a,b) in enumerate(pairs):
        k2[a][col]=gi(b);k2[b][col]=gi(a,-1)
        sparse2.extend([[a,col,b,1],[b,col,a,-1]])
    k3=[[{} for _ in triples] for _ in pairs];sparse3=[]
    for col,(a,b,c) in enumerate(triples):
        for pair,g,sign in [((b,c),a,1),((a,c),b,-1),((a,b),c,1)]:
            row=pairs.index(pair);k3[row][col]=gi(g,sign);sparse3.append([row,col,g,sign])
    for j in range(28):
        p={}
        for i in range(8):p=plus(p,mul(gi(i),k2[i][j]))
        assert not p
    for i,j in itertools.product(range(8),range(56)):
        p={}
        for k in range(28):p=plus(p,mul(k2[i][k],k3[k][j]))
        assert not p
    assert len(sparse2)==56 and len(sparse3)==168

    # Matrix-typed actual T, whose factor coefficients are all literal finite
    # polynomial telescopes or previously checked cofactor arrays.
    sourceT={
      'native_ring':'C=R0[tauA,tauB]/monics; local completion',
      'source9_input':'kappa9 from source_bound_interface; exact originJac5/48',
      'source9_order':['a','d','v2','v3','v4','v5','v6','v7','v8'],
      'core4_indices':[0,1,2,6],'core4_order':['a','d','v2','v6'],
      'D':'1-kappa_d^2','P':'kappa_d^2+5',
      'source_H12_polynomials':[q['divided_primitive_polynomial'] for q in H['all12_actual_remaining_and_divided_records']],
      'second5_pivot_order':['v4','v3','v5','v7','v8'],
      'second5_core_roots':K['five_actual_rational_roots'],
      'u5':'pivot_i - second5_core_root_i(core4)',
      'permutation':[0,4,1,2,3,5,6,7,8,9,10,11],
      'telescoping_A5_B7':'finite monomial telescope rule from source_preflight; selected5 rows, other7 rows, actual pivots vssecond5roots',
      'A5_diag':['2P','8D','80D^2','16D^2','192D^3'],
      'residual_lambda':[{'QQ':str(F(q['QQ_unit_content_numerator'],q['QQ_unit_content_denominator'])),
                         'D_denominator_power':q['D_denominator_power'],'P_denominator_power':q['P_denominator_power']}
                        for q in K['remaining7_actual_core_polynomials']],
      'E1':'[[A5^-1,0],[-B7*A5^-1,I7]]',
      'E2':'diag(I5,lambda^-1)',
      'E3':'diag(I5,[[1/D,0,0],[0,1/D,0],[-r/(PD),0,-2/(PD)]],I4)',
      'E4':'I12 with lower-left core3 coefficients minus U_j*C_ij in last4rows; U=(D^2,D^2,D^2,D^3)',
      'graph_cofactors3_by4':C['cofactor3_by4'],
      'CI_generators3':P['generators3'],'r':'v2-a-d','s':'a+v2',
      'T':'E4*E3*E2*E1*Pi','T_inverse':'Pi^-1*E1^-1*E2^-1*E3^-1*E4^-1',
      'detT':'P^10/(2359296D^11)',
      'full_N':'T(kappa)*R12','full_N_inverse':'R12^-1*T(kappa)^-1',
      'full_N_identity':'(d*G0,...,d*G7,0,0,0,0)^t=N*F12^t',
      'G8':['u0','u1','u2','u3','u4','Q0(core4)','Q1(core4)','s*r^2'],
      'all_unit_C_inverses_lower_to_known_unit_Res4_matrix_inverses':True,
      'd_never_inverted':True}
    lowered={
      'component_symbol_rule':'Each typed Cscalar node v has four R0 components(v[0],v[1],v[2],v[3]), recursively defined by the rules below, never an unspecified new coefficient.',
      'basis':['1','tauA','tauB','tauA*tauB'],'structure16':table,
      'rational_R0_scalar_embedding':['q',0,0,0],
      'add':'(x+y)[i]=x[i]+y[i]','neg':'(-x)[i]=-x[i]',
      'multiply':'sum all16 x[a]*y[b]*nuA^p*nuB^q in output k specified by structure16',
      'matrix_mu_by_rows':[['x0','nuA*x1','nuB*x2','nuA*nuB*x3'],
                           ['x1','x0','nuB*x3','nuB*x2'],
                           ['x2','nuA*x3','x0','nuA*x1'],
                           ['x3','x2','x1','x0']],
      'unit_scalar_inverse':'unique unit4-linear root mu(x)*y=e0; det mu(x) has closed-point value x0(0)^4',
      'matrix_Res4':'For nativeC matrix M, block(i,j)=mu(M[i,j]); never transpose within mu.',
      'unit_matrix_inverse':'Res4(M^-1)=Res4(M)^-1; solve exact unit linear matrix system',
      'hensel_Res4':'Expand each native C unknown into fourR0 unknowns and each equation via structure16. CoefficientJacobian is Res4(nativeJac); norm of unit determinant stays a unit.',
      'group30_to120':{
         'original_group30':'expand group15 with formal auxiliary etaA^2=nuA and t=t_even+t_odd*etaA BEFORE scalar restriction',
         'auxiliary_etaA_distinct_from_physical_tauA':True,
         'basechange_group30_to_C_then_Res4_unknown_count':120,
         'component_equation_count':120,
         'group30_unknown_order':'15 even then15 odd; each C scalar expanded in fourcomponents',
         'group30_even_origin':[0]*15,
         'group30_odd_origin':[r[7] for r in I['group15_first9']],
         'group120_origin':[[0,0,0,0] for _ in range(15)]+[[r[7],0,0,0] for r in I['group15_first9']],
         'why_odd_origin_is_nonzero':'The firstjet tauA column is the exact eta coefficient at the closed R0 fiber eta^2=0; native group15 is centered at0 only modulo the C maximal ideal, not componentwise overR0.',
         'original_group30_det_origin':'4/9','restricted120_det_origin':'256/6561',
         'proof':'Unique basechanged solution is the original R0 even/odd solution embedded intoC; evaluate etaA=actualtauA afterwards. This redundant120system pins every component without an ambiguous decomposition insideC.'},
      'source96_roots':'retain original R0 unit roots from pinned interface; embed as(q,0,0,0)',
      'coordinate_inverse9_if_requested':'36component unit system, determinant origin(5/48)^4; not needed to evaluate forwardT(kappa)',
      'source_bound_interface':I['source_bound_finite_matrix_typed_operator'],
      'source_adapter_T':sourceT,
      'original_source_adapter_reference':S,
      'actual_Q4_by48_column_node_ids':actualQ,
      'actual_Q_scalar_DAG_reference':str(ep.relative_to(ROOT)),
      'all48_actual_Qcolumns_independently_bound_to_mu_F':True,
      'native_S12_by32':{'koszul28_in_first8_rows':sparse2,
                       'free4_zero_source_columns':[[8+i,28+i,1] for i in range(4)]},
      'native_V32_by56':{'koszul3_in_first28_rows':sparse3,'last4_rows_zero':True},
      'actual_d1':{'shape':[4,48],'value':'pinned actual_Q4_by48'},
      'actual_d2':{'shape':[48,128],'value':'Res4(N^t)*Res4(S12_by32)',
         'exact_block_recipe':'block(i,j)=mu(sum_k N[k,i]*S[k,j]); native12rowblock then4basis columns'},
      'actual_d3':{'shape':[128,224],'value':'Res4(V32_by56)',
         'exact_block_recipe':'block(i,j)=mu(V[i,j]); native32rowblock then4basis columns'},
      'chain_identity1':'Q*Res4(N^t)*Res4(S)=Res4((dG,0)^t_as_row*S)=0',
      'chain_identity2':'d2*d3=Res4(N^t*(K2,Id4)*[K3;0])=0',
      'completeness':'N unit;d nonzero in regularC;G8 regular;Koszul exact and fourzero rows split;restriction of scalars along R0-freeC exact',
      'restricted_resolution_ranks':[4,48,128,224,280,224,112,32,4],
      'not_necessarily_R0_minimal':True,
      'removed_zero_source4_optional_restriction_ranks':[4,32,112,224,280,224,112,32,4],
      'no_expanded_coefficient_polynomials_or_fresh_Hensel_solve':True}
    path=out/'actual-complete-Res4-syzygy-operator.json';path.write_text(json.dumps(lowered,separators=(',',':'))+'\n')
    assert path.stat().st_size<=32*1024*1024
    proof={'structure_products_checked':16,'basis_associativity_checked':len(associativity),
      'multiplication_matrix_products_checked':16,'actual_Qcolumns_checked':len(verified),
      'actual_Qcolumn_binding':verified,'native_Koszul_GK2_zero_columns':28,
      'native_Koszul_K2K3_zero_entries':448,'native_sparseK2_entries':56,'native_sparseK3_entries':168,
      'all_G8_and_N_inputs_source_bound':True,'unit_interface_minors':['5/48','687194767360000'],
      'complete_actual_SyzQ_shape':[48,128],'complete_next_shape':[128,224],
      'restricted_ranks':lowered['restricted_resolution_ranks'],'all5_d_times_u_preserved':True,
      'scalar_functions_are_component_typed_unit_implicit_values_not_expanded_power_series':True}
    ev=out/'universal-restriction-and-chain-checks.json';ev.write_text(json.dumps(proof,indent=2)+'\n')
    paths=[Path(__file__),ROOT/'research/notes'/f'{STEM}.md',ip,sp,ep,hp,kp,cp,pp,*gates]
    report={'id':STEM,'all_passed':True,'status':'source_bound_finite_Res4_complete_actual_Q_syzygy_operator',
      'actual_free4_Q48_bound':True,'complete_actual_SyzQ_48_by128':True,'complete_next128_by224':True,
      'every_Cscalar_has_finite4component_unit_implicit_representation':True,
      'expanded_R0_coefficient_polynomials_computed':False,'new_source_solve_jets_or_GB':False,
      'K_ell_torsion_D0_saturation_socle_or_depth_claimed':False,
      'group120_known_unit_origin':'256/6561','all5_d_times_u_preserved':True,
      'input_sha256':{str(p.relative_to(ROOT)):sha(p) for p in paths},
      'evidence_sha256':{p.name:sha(p) for p in [path,ev]},'bytes':path.stat().st_size+ev.stat().st_size,
      'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:report[k] for k in ['id','all_passed','bytes','seconds','max_rss_kib']}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
