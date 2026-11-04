#!/usr/bin/env python3
"""ONE actual H_A-valued W1 normalization and finite missing-gauge preflight.

No old49 identification is inferred from the central point or an abstract type.
"""
import argparse
from fractions import Fraction as F
import gzip,hashlib,json
from pathlib import Path
import resource,signal,time
import sympy as sy
from ambient_B18_actual_top_coefficient_preflight_20261104_v63 import DAG,Quadratic,action,PAIRS,MONS,PI,MI
from ambient_B18_top_graph_telescoping_20261104_v63 import independent_rows,linear_solve
from s61_model import S61,vertical_rref
from weight1_structure import sp_basis
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_W1_TAIL_COORDINATE_PREFLIGHT_20261104_v63'
def read(p):return json.loads(gzip.decompress(p.read_bytes()) if p.suffix=='.gz' else p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main(out):
    started=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));out.mkdir(parents=True,exist_ok=False)
    sp=ROOT/'reports/AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63/actual-source-coefficients.json'
    op=ROOT/'reports/AMBIENT_B18_T3_DOUBLE_ANNIHILATOR_PREFLIGHT_20261104_v63/exact-residues-and-two-ideals.json'
    wp=ROOT/'reports/AMBIENT_W1_S61_INTERSECTION_RECOGNITION_20261103_v63/recognition.json.gz'
    wr=wp.parent/'report.json'
    weights=ROOT/'reports/AMBIENT_B18_RELATIVE_WEIGHTED_MODEL_20261104_v63/certificate.json'
    bundlep=ROOT/'reports/AMBIENT_FIBER_BUNDLE_20261103_v63/report.json'
    bp=ROOT/'reports/AMBIENT_BASE_PINCHED_SURFACE_20261103_v63/report.json'
    src,origin,wrec,weight,bundle,base=map(read,[sp,op,wp,weights,bundlep,bp]);assert read(wr)['all_passed']
    ca=src['families']['CA'];D=DAG();D.nodes=[tuple(n) for n in src['DAG_scalar_nodes']]
    D.values=list(map(F,origin['all339669_scalar_origin_values_evaluated'][:len(D.nodes)]));D.cache={n:i for i,n in enumerate(D.nodes)}
    Q=object.__new__(Quadratic);Q.D=D;Q.S=D.zero;Q.T=ca['actual_tracefree_nu'];Q.zero=(D.zero,D.zero);Q.one=(D.one,D.zero);Q.p=(D.zero,D.one)
    def C(v):return Q.scalar(D.const(str(v)))
    def matrix(M):return [[C(M[i,j]) for j in range(M.cols)] for i in range(M.rows)]
    U=[Q.total(ca['complete_normalized_physical_w_coefficients'][d][r] for d in range(4)) for r in range(150)]
    U0=[Q.center(p)[0] for p in U];assert all(Q.center(p)[1]==0 for p in U)
    frame=sy.Matrix(wrec['exact_homogeneous_frame']);fixed=action(Q,U,matrix(frame),matrix(frame.inv()),Q.one)
    value=[fixed[15*PI[pair]+MI[(0,0)]] for pair in [(1,2),(1,3),(1,4),(2,3),(2,4),(3,4)]]
    a,b,c,d,e,f=value;ai=Q.inv(a);delta=Q.add(Q.add(Q.mul(a,f),Q.neg(Q.mul(b,e))),Q.mul(c,d));di=Q.inv(delta)
    assert Q.center(a)==(F(1),F(0)) and Q.center(delta)==(F(1),F(0))
    FF=[[Q.zero for j in range(5)] for i in range(5)];FI=[[Q.zero for j in range(5)] for i in range(5)]
    for i in range(5):FF[i][i]=FI[i][i]=Q.one
    FF[2][2]=ai;FF[3][1]=Q.mul(d,ai);FF[3][2]=Q.neg(Q.mul(b,ai))
    FF[4][1]=Q.mul(e,di);FF[4][2]=Q.neg(Q.mul(c,di));FF[4][4]=Q.mul(a,di)
    FI[2][2]=a;FI[3][1]=Q.neg(Q.mul(d,ai));FI[3][2]=b
    FI[4][1]=Q.neg(Q.mul(e,ai));FI[4][2]=c;FI[4][4]=Q.mul(delta,ai)
    based=action(Q,fixed,FF,FI,Q.one)
    model=S61();rr,piv=vertical_rref()
    def quotient(V):return [Q.add(v,Q.neg(Q.total(Q.scale(V[p],rr[i,k]) for i,p in enumerate(piv) if rr[i,k]))) for k,v in enumerate(V)]
    initial=quotient(based);response=model.perspective_matrix();ri=response.inv()
    alpha=[Q.neg(Q.total(Q.scale(initial[p],ri[i,j]) for j,p in enumerate(model.certificate['perspective_orbit_pivots_raw']) if ri[i,j])) for i in range(4)]
    P=[[Q.one if i==j else Q.zero for j in range(5)] for i in range(5)];Pinv=[row[:] for row in P]
    for j in range(1,5):P[0][j]=alpha[j-1];Pinv[0][j]=Q.neg(alpha[j-1])
    normal=quotient(action(Q,based,P,Pinv,Q.one));u=[normal[i] for i in model.certificate['raw_free_indices']]
    assert [str(Q.center(v)[0]) for v in u]==wrec['exact_normalized61'] and all(Q.center(v)[1]==0 for v in u)
    bi=bundle['beta_indices_zero_based'];qi=bundle['q_indices_zero_based'];adapter=sy.Matrix(bundle['beta16_to_intrinsic24_adapter'])
    beta=[Q.total(Q.scale(u[bi[j]],adapter[i,j]) for j in range(len(bi)) if adapter[i,j]) for i in range(24)]
    assert [str(Q.center(v)[0]) for v in beta]==wrec['full_beta24']
    # Genuine original-point group response, not a new source/jet point.
    def lie(U,i,j):
        V=[F(0)]*150
        def add(k,v):V[k]+=v
        for pair,(a,b) in enumerate(PAIRS):
            for mon,(r,s) in enumerate(MONS):
                co=U[15*pair+mon]
                if not co:continue
                if a==j and i!=b:add(15*PI[tuple(sorted((i,b)))]+mon,co*(1 if i<b else-1))
                if b==j and a!=i:add(15*PI[tuple(sorted((a,i)))]+mon,co*(1 if a<i else-1))
                for q,other in [(r,s),(s,r)]:
                    if q==i:add(15*pair+MI[tuple(sorted((j,other)))],-co)
        return V
    section=[];labels=[]
    for i,j in [(2,2),(4,4),(3,2),(4,2),(3,1),(4,1)]+[(0,j) for j in range(1,5)]:
        M=sy.zeros(5);M[i,j]=1;section.append(M);labels.append(f'value/perspective E{i}{j}')
    chosen=base['orbit_constant_block_witness']['columns_zero_based']
    for j in chosen:
        M=sy.zeros(5);M[1:5,1:5]=sp_basis()[j];section.append(M);labels.append(f'Sp4 orbit basis{j}')
    section=[frame.inv()*M*frame for M in section]
    L={(i,j):lie(U0,i,j) for i in range(5) for j in range(5)}
    columns=[[sum(F(str(M[i,j]))*L[i,j][k] for i in range(5) for j in range(5) if M[i,j]) for k in range(150)] for M in section]
    cuts=weight['common_physical_group_cuts19']+[140];R=[[col[k] for col in columns] for k in cuts]
    selected=independent_rows(R);rank=len(selected);det=None
    if rank==17:_,det=linear_solve([R[i] for i in selected],[[F(int(i==j)) for j in range(17)] for i in range(17)])
    units=[x for x in D.divisions if x[0]>=len(src['DAG_scalar_nodes'])]
    ev={'DAG_scalar_nodes':D.nodes,'CA_trace_square_nu_node':Q.T,
        'actual_physical150_tail_w1_trace_pairs':U,'fixed_GL5':[[str(v) for v in row] for row in frame.tolist()],
        'actual_variable_value6_pairs':value,'value_pivot_a_pair':a,'value_Pfaffian_delta_pair':delta,
        'exact_value_Darboux_F5_pairs':FF,'exact_value_Darboux_F5_inverse_pairs':FI,
        'perspective_parameters4_pairs':alpha,'all150_normalized_pairs':normal,'all61_normalized_pairs':u,
        'all24_primitive_beta_pairs':beta,'all45_q_pairs':[u[i] for i in qi],
        'all_new_unit_inverse_norms':units,'all61_center_matches_frozen_W1':True,
        'original_physical_E1_coefficients':list(map(str,U0)),
        'chosen_group17_labels':labels,'chosen_group17_matrices_original_frame':[[[str(v) for v in row] for row in M.tolist()] for M in section],
        'physical150_by17_group_response':[[str(col[k]) for col in columns] for k in range(150)],
        'common19cuts_plus_w140':cuts,'all20_by17_response':[[str(v) for v in row] for row in R],
        'response_rank':rank,'response_selected17_rows':selected,'response_minor':str(det) if det else None,
        'missing_Sp4_target_gauge7_equations':'beta[0,4,5,6,7,8,16]=frozen beta0 constants, selectedSp4 orbitminor -80',
        'old49_base_extraction_AFTER_gauge':{'ss':'beta11/2','aa':'beta2','bb':'-beta14/2','cc':'-beta3/2'},
        'old49_q_extraction_AFTER_gauge':'commonqselected10 v, remaining35 z, peroldbuilder',
        'additional_slice_equations_required':3,
        'whole_tail_binding_nil_pullback_or_Jinf_syzygies_not_claimed':True}
    encoded=json.dumps(ev,separators=(',',':'))+'\n';assert len(encoded.encode())<32*1024*1024
    (out/'tail-coordinate-preflight-DAG.json').write_text(encoded)
    report={'id':STEM,'all_passed':True,'status':'actual_W1_to_S61_finite_map_plus_remaining49_gauge_gates',
            'actual_H_valued_w1_tensor_all150_reconstructed':True,
            'generic_value_pivot_and_Pfaffian_both_units_at_E1':True,
            'fixed_GL_then_variable_Darboux_then_perspective_all61_constructed':True,
            'all61_and_beta24_exactly_match_frozen_E1_center':True,
            'group17_response_on_common19cuts_and_w_rank':rank,'selected17_minor':str(det) if det else None,
            'old49_gauge7_and_residual3_slice_not_constructed':True,
            'no_whole_tail_type_nil_pullback_or_syzygy_claim':True,
            'no_new_source_jets_GB_points_or_p_inverse':True,
            'input_sha256':{str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__),sp,op,wp,wr,weights,bundlep,bp,ROOT/'code/s61_model.py']},
            'evidence_sha256':{'tail-coordinate-preflight-DAG.json':sha(out/'tail-coordinate-preflight-DAG.json')},
            'DAG_nodes':len(D.nodes),'certificate_bytes':len(encoded.encode()),'seconds':time.monotonic()-started,
            'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:report[k] for k in ['all_passed','group17_response_on_common19cuts_and_w_rank','selected17_minor','DAG_nodes','certificate_bytes','seconds','max_rss_kib']}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
