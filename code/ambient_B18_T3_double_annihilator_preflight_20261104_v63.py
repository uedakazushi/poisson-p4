#!/usr/bin/env python3
"""ONE bounded exact-origin override and two actual zero-odd tangent ideals.

No new jet or module-membership inference from the residue matrix.
"""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import resource
import signal
import time
from ambient_B18_top_graph_telescoping_20261104_v63 import independent_rows,linear_solve
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_T3_DOUBLE_ANNIHILATOR_PREFLIGHT_20261104_v63'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def enc(M):return [[str(v) for v in row] for row in M]
def nullspace(M):
    A=[r[:] for r in M]; n=len(A[0]); row=0;piv=[]
    for col in range(n):
        i=next((i for i in range(row,len(A)) if A[i][col]),None)
        if i is None:continue
        A[row],A[i]=A[i],A[row];q=A[row][col];A[row]=[v/q for v in A[row]]
        for j in range(len(A)):
            if j!=row and A[j][col]:
                q=A[j][col];A[j]=[v-q*w for v,w in zip(A[j],A[row])]
        piv.append(col);row+=1
        if row==len(A):break
    free=[j for j in range(n) if j not in piv];answer=[]
    for f in free:
        v=[F(0)]*n;v[f]=1
        for i,p in enumerate(piv):v[p]=-A[i][f]
        assert all(sum(a*b for a,b in zip(r,v))==0 for r in M)
        answer.append(v)
    return answer
def main(out):
    started=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000))
    out.mkdir(parents=True,exist_ok=False)
    fullp=ROOT/'reports/AMBIENT_B18_FULL_RELATIVE_INPUT_20261104_v63/full-relative186-DAG.json'
    cap=ROOT/'reports/AMBIENT_B18_T3_CURVE_ANNIHILATOR_20261104_v63/annihilator-and-parameter.json'
    cbp=ROOT/'reports/AMBIENT_B18_TOP_GRAPH_TELESCOPING_20261104_v63/rank-gate.json'
    origins=ROOT/'reports/AMBIENT_B18_CB_G_DIVISION_ORIGINS_20261104_v63/degree3-certificate.json'
    originrep=origins.parent/'report.json'
    d,ca,cb,ori=map(read,[fullp,cap,cbp,origins]);assert read(originrep)['all_passed']
    overrides={q['known_regular_full_DAG_g_quotient_node']:F(q['quotient_origin_value']) for q in ori['all4_regular_quotient_origin_values']}
    assert len(overrides)==4 and all(v==0 for v in overrides.values())
    divs={q['node']:F(q['origin_value']) for q in d['all_coordinate_divisions']}
    vals=[];inverse_denominators=[]
    for i,n in enumerate(d['DAG_scalar_nodes']):
        op=n[0]
        if op=='c':v=F(n[1],n[2])
        elif op=='v':v=F(0)
        elif op=='a':v=vals[n[1]]+vals[n[2]]
        elif op=='n':v=-vals[n[1]]
        elif op=='m':v=vals[n[1]]*vals[n[2]]
        elif op=='i':
            assert vals[n[1]];v=1/vals[n[1]];inverse_denominators.append([i,n[1],str(vals[n[1]])])
        elif op=='r':v=divs[i]
        elif op=='g':v=overrides[i]
        else:raise ValueError(op)
        vals.append(v)
    residue_Q=[[vals[a],vals[b]] for a,b in d['actual_Qn_R0_squared_generator_matrices']['3']]
    residue_top=[[vals[a],vals[b]] for a,b in d['actual12_top_T3_additional_columns']]
    caM=[[F(v) for v in r] for r in ca['actual_CA10_by9_Jacobian']]
    cbM=[[F(v) for v in r] for r in cb['actual_CB_10_by9_odd_Jacobian']]
    selected=[0,1,2,3,4,5,6,8,9]
    invca,detca=linear_solve([caM[i] for i in selected],[[F(int(i==j)) for j in range(9)] for i in range(9)])
    assert detca==F(128,3)
    cb_in_ca=[[sum(cbM[i][k]*invca[k][j] for k in range(9)) for j in range(9)] for i in range(8)]
    restriction=[r[7:] for r in cb_in_ca]
    union=caM[:8]+cbM[:8];unionindices=independent_rows(union)
    data={'four_origin_overrides':{str(i):str(v) for i,v in overrides.items()},
          'all339669_scalar_origin_values_evaluated':list(map(str,vals)),
          'all_unit_inverse_denominators_checked':inverse_denominators,
          'Q3_all2034_residue_columns':enc(residue_Q),'top12_all_residue_columns':enc(residue_top),
          'Q3_residue_column_rank':len(independent_rows(residue_Q)),
          'top12_residue_column_rank':len(independent_rows(residue_top)),
          'actual_CA_zero_odd8_firstjets':enc(caM[:8]),'actual_CB_zero_odd8_firstjets':enc(cbM[:8]),
          'union_firstjet_rank':len(unionindices),'union_independent_row_indices':unionindices,
          'union_common_tangent_kernel_basis':enc(nullspace(union)),
          'CA_F9_inverse_jacobian':enc(invca),'CA_F9_determinant':str(detca),
          'actual_CB_zero8_in_CA_F9_firstjets':enc(cb_in_ca),
          'CB_zero8_restriction_to_CA_r_s_firstjets':enc(restriction),
          'restricted_rank':len(independent_rows(restriction)),
          'restricted_common_tangent_kernel':enc(nullspace(restriction)),
          'no_T3_origin_fiber_inferred_from_ambient2row_residues':True,
          'no_colon_depth_or_nil_answer_computed':True}
    encoded=json.dumps(data,separators=(',',':'))+'\n';assert len(encoded.encode())<32*1024*1024
    (out/'exact-residues-and-two-ideals.json').write_text(encoded)
    report={'id':STEM,'all_passed':True,'status':'actual_origin_override_and_double_zero_ideal_preflight',
            'all4_proved_zero_origins_applied_without_mutating_frozen_DAG':True,
            'DAG_origin_values_count':len(vals),'Q3_center_rank':data['Q3_residue_column_rank'],
            'top12_center_rank':data['top12_residue_column_rank'],
            'actual_I0_CA_plus_I0_CB_tangent_rank':len(unionindices),
            'CB_zero8_restriction_rank_on_actual_CA_curve_tangent_plane':data['restricted_rank'],
            'no_Tor_unsafe_residue_module_or_depth_claim':True,
            'no_new_jets_GB_new_point_source_regeneration_or_p_inverse':True,
            'input_sha256':{str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__),fullp,cap,cbp,origins,originrep]},
            'evidence_sha256':{'exact-residues-and-two-ideals.json':sha(out/'exact-residues-and-two-ideals.json')},
            'certificate_bytes':len(encoded.encode()),'seconds':time.monotonic()-started,
            'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:report[k] for k in ['all_passed','Q3_center_rank','top12_center_rank','actual_I0_CA_plus_I0_CB_tangent_rank','CB_zero8_restriction_rank_on_actual_CA_curve_tangent_plane','seconds','max_rss_kib']}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
