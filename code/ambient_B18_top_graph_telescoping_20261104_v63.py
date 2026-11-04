#!/usr/bin/env python3
"""ONE exact derivative gate, then actual regular-algebraic graph lifts.

No source/core regeneration, GB, point or formal-series truncation.
"""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import resource
import signal
import time
from ambient_B18_actual_top_coefficient_preflight_20261104_v63 import DAG
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_TOP_GRAPH_TELESCOPING_20261104_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def merge(a,b,scale=F(1)):
    z=dict(a)
    for k,v in b.items():
        z[k]=z.get(k,F(0))+scale*v
        if not z[k]:z.pop(k)
    return z

def forward(nodes,relevant):
    vals=[];grads=[]
    for n,node in enumerate(nodes):
        op=node[0]
        if op=='c':v=F(node[1],node[2]);g={}
        elif op=='v':v=F(0);g={n:F(1)}if n in relevant else{}
        elif op=='n':v=-vals[node[1]];g={k:-x for k,x in grads[node[1]].items()}
        elif op=='a':
            i,j=node[1:];v=vals[i]+vals[j];g=merge(grads[i],grads[j])
        elif op=='m':
            i,j=node[1:];v=vals[i]*vals[j]
            g=merge({k:vals[j]*x for k,x in grads[i].items()if vals[j]},grads[j],vals[i])
        elif op=='i':
            i=node[1];assert vals[i];v=1/vals[i];g={k:-x/(vals[i]*vals[i])for k,x in grads[i].items()}
        else:raise ValueError(op)
        vals.append(v);grads.append(g)
    return vals,grads

def linear_solve(M,B):
    n=len(M);m=len(B[0]);a=[list(r)+list(b)for r,b in zip(M,B)];det=F(1)
    for j in range(n):
        i=next(i for i in range(j,n)if a[i][j])
        if i!=j:a[i],a[j]=a[j],a[i];det=-det
        pivot=a[j][j];det*=pivot
        a[j]=[x/pivot for x in a[j]]
        for i in range(n):
            if i==j or not a[i][j]:continue
            c=a[i][j];a[i]=[x-c*y for x,y in zip(a[i],a[j])]
    return[row[n:]for row in a],det

def independent_rows(M):
    a=[list(row)for row in zip(*M)];nr=len(a);nc=len(a[0]);row=0;piv=[]
    for j in range(nc):
        i=next((i for i in range(row,nr)if a[i][j]),None)
        if i is None:continue
        a[row],a[i]=a[i],a[row];p=a[row][j];a[row]=[v/p for v in a[row]]
        for i in range(nr):
            if i!=row and a[i][j]:
                c=a[i][j];a[i]=[x-c*y for x,y in zip(a[i],a[row])]
        piv.append(j);row+=1
        if row==nr:break
    return piv

def encode(M):return[[str(x)for x in row]for row in M]

def replay(D,nodes,bindings,needed=None):
    mapping={}
    for old,node in enumerate(nodes):
        if needed is not None and old not in needed:continue
        op=node[0]
        if op=='c':new=D.const(F(node[1],node[2]))
        elif op=='v':new=bindings.get(old,D.var(node[1]))
        elif op=='n':new=D.neg(mapping[node[1]])
        elif op=='a':new=D.add(mapping[node[1]],mapping[node[2]])
        elif op=='m':new=D.mul(mapping[node[1]],mapping[node[2]])
        elif op=='i':new=D.inverse(mapping[node[1]])
        else:raise ValueError(op)
        mapping[old]=new
    return mapping

def ancestors(nodes,roots):
    wanted=set();stack=list(roots)
    while stack:
        i=stack.pop()
        if i in wanted:continue
        wanted.add(i);node=nodes[i]
        if node[0]in['n','i']:stack.append(node[1])
        elif node[0]in['a','m']:stack.extend(node[1:])
    return wanted

def main(out):
    started=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));out.mkdir(parents=True,exist_ok=False)
    folder=ROOT/'reports/AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63'
    path=folder/'actual-source-coefficients.json';data=json.loads(path.read_text());nodes=data['DAG_scalar_nodes']
    ca,cb=data['families']['CA'],data['families']['CB'];base=data['physical_zero_base_rows9']
    varids={node[1]:i for i,node in enumerate(nodes)if node[0]=='v'}
    xids=[varids[f'x{i}']for i in base]
    zids=[i for pair in cb['coefficient_unknown_pairs'].values()for i in pair]
    zids +=[cb['p_quadratic_relation']['s'],cb['p_quadratic_relation']['t']]
    assert len(zids)==48 and len(set(zids))==48
    vals,grads=forward(nodes,set(xids+zids))
    eq=cb['all48_zero_IFT_quadratic_coefficient_equations']
    J=[[grads[i].get(v,F(0))for v in zids]for i in eq]
    B=[[-grads[i].get(v,F(0))for v in xids]for i in eq]
    zx,detJ=linear_solve(J,B)
    assert detJ and abs(detJ)==F(cb['unit48_block_determinant_even24_then_odd24_vars_a23_t_b23_s'])
    def derivative(i):
        answer=[grads[i].get(v,F(0))for v in xids]
        for k,v in enumerate(zids):
            c=grads[i].get(v,F(0))
            if c:answer=[a+c*b for a,b in zip(answer,zx[k])]
        return answer
    generators=cb['actual_zero_odd_I0_generators']+cb['actual_positive_odd_J0_generators']
    gids=[i for row,i in generators];assert all(vals[i]==0 for i in gids)
    matrix=[derivative(i)for i in gids];selected=independent_rows(matrix)
    gate={'actual_CB_10_by9_odd_Jacobian':encode(matrix),'generator_physical_rows':[row for row,node in generators],
          'source48_actual_AD_Jacobian':encode(J),'source48_dz_dx':encode(zx),
          'source48_Jacobian_determinant':str(detJ),'selected_generator_indices':selected,'rank':len(selected)}
    if len(selected)!=9:
        (out/'rank-gate.json').write_text(json.dumps(gate,indent=2)+'\n')
        (out/'report.json').write_text(json.dumps({'id':STEM,'all_passed':False,'status':'actual_rank_not9_STOP','rank':len(selected)},indent=2)+'\n')
        print(json.dumps({'rank':len(selected),'status':'STOP'}));return
    F9=[matrix[i]for i in selected]
    inverse,detF=linear_solve(F9,[[F(int(i==j))for j in range(9)]for i in range(9)])
    gate['selected9_Jacobian_determinant']=str(detF)
    D=DAG();original=replay(D,nodes,{})
    assert all(original[i]==i for i in range(len(nodes)))
    Y=[original[gids[i]]for i in selected]
    cubic_B=cb['actual_individual_physical_cubic_trace_pairs']
    brows=[row for row,pair in cubic_B];bids=[pair[1]for row,pair in cubic_B]
    assert all(vals[i]==0 for i in bids)
    hcenters=[]
    for b in bids:
        dx=derivative(b);hcenters.append([sum(dx[k]*inverse[k][j]for k in range(9))for j in range(9)])
    roots=eq+[gids[i]for i in selected]+bids
    needed=ancestors(nodes,roots)
    assert all(nodes[i][1]in {nodes[v][1]for v in xids+zids}for i in needed if nodes[i][0]=='v')
    inverse_instances=[];b_partial=[]
    for level in range(1,10):
        bindings={old:D.var(f'INV{level}_{nodes[old][1]}')for old in xids+zids}
        instance=replay(D,nodes,bindings,needed)
        equations=[instance[i]for i in eq]
        equations +=[D.add(instance[gids[j]],D.neg(Y[k]if k<level else D.zero))for k,j in enumerate(selected)]
        assert len(equations)==57 and all(D.values[i]==0 for i in equations)
        inverse_instances.append({'level':level,'variable_bindings':[[old,new]for old,new in bindings.items()],
                                  'all57_equations':equations,
                                  'unit57_Jacobian_determinant_in_z_then_x_order':str(detJ*detF)})
        b_partial.append([instance[i]for i in bids])
    h_records=[];relative=[];top=[]
    ca_zero=ca['actual_zero_odd_I0_generators'];ca_pos=ca['actual_positive_odd_J0_generators']
    zero_rows=[row for row,node in ca_zero];pos_rows=[row for row,node in ca_pos]
    assert zero_rows==[row for row,node in cb['actual_zero_odd_I0_generators']]
    for alpha,row in enumerate(brows):
        hs=[D.zero]*10;division_records=[];previous=D.zero
        for i,index in enumerate(selected):
            current=b_partial[i][alpha];numerator=D.add(current,D.neg(previous))
            # Not a unit inverse: exact coordinate-restriction divisibility.
            h=D.put(('r',numerator,Y[i],alpha,i),hcenters[alpha][i])
            division_records.append({'h_node':h,'numerator':numerator,'coordinate':Y[i],
                                     'proof':'two consecutive unit57 inverse instances agree after coordinate Yi=0',
                                     'origin_value':str(hcenters[alpha][i])})
            hs[index]=h;previous=current
        h_records.append({'raw_cubic_row':row,'membership_coefficients10':hs,
                          'regular_coordinate_divisions9':division_records,
                          'exact_identity':'sum selected Fi*hi = b9-b0 = actual cubic odd coefficient, by telescoping and unit inverse uniqueness'})
        Aca,Bca=ca['actual_individual_physical_cubic_trace_pairs'][alpha][1]
        Acb,Bcb=cb['actual_individual_physical_cubic_trace_pairs'][alpha][1]
        even=D.add(Acb,D.neg(Aca));odd=D.neg(Bca)
        for j,(zrow,fca)in enumerate(ca_zero):
            aA=ca['complete_normalized_physical_w_coefficients'][0][zrow][0]
            aB=cb['complete_normalized_physical_w_coefficients'][0][zrow][0]
            even=D.add(even,D.mul(hs[j],D.add(aA,D.neg(aB))))
            odd=D.add(odd,D.mul(hs[j],fca))
        for e,(prow,jca)in enumerate(ca_pos):
            aA=ca['complete_normalized_physical_w_coefficients'][1][prow][0]
            aB=cb['complete_normalized_physical_w_coefficients'][1][prow][0]
            even=D.add(even,D.mul(hs[8+e],D.add(aA,D.neg(aB))))
            odd=D.add(odd,D.mul(hs[8+e],jca))
        relative.append({'raw_cubic_row':row,'same_CA_tau_basis_CB_graph_minus_CA_pair':[even,odd]})
        for e,(prow,jca)in enumerate(ca_pos):
            left=D.mul(ca['actual_tracefree_nu'],D.mul(jca,odd));right=D.mul(jca,even)
            top.append({'positive_odd_raw_row':prow,'cubic_raw_row':row,'pair_nodes':[left,right]})
    gate['b_cubic_origin_values']=[str(vals[i])for i in bids]
    certificate={'DAG_scalar_nodes':D.nodes,
                 'operation_r':'regular coordinate quotient numerator/Yi, algebraic with exact restriction-zero certificate; NOT arbitrary nonunit inverse',
                 'original_unit48_equations_CA':ca['all48_zero_IFT_quadratic_coefficient_equations'],
                 'original_unit48_equations_CB':eq,
                 'coordinate_gate':gate,'selected_Y_coordinate_nodes':Y,
                 'all9_pointed_unit57_inverse_instances':inverse_instances,
                 'all6_actual_CB_odd_membership_lifts':h_records,
                 'all6_actual_relative_cubic_pairs_in_same_CA_trace_basis':relative,
                 'all12_actual_top_T3_pairs':top,
                 'CA_tau_square_nu_node':ca['actual_tracefree_nu'],
                 'old_DAG_input_sha256':sha(path),
                 'all_lower_CB_kernel_and_graph_relations_must_still_be_kept':True}
    encoded=json.dumps(certificate,separators=(',',':'))+'\n';assert len(encoded.encode())<=32*1024*1024
    (out/'actual-top12-DAG.json').write_text(encoded)
    (out/'rank-gate.json').write_text(json.dumps(gate,indent=2)+'\n')
    result={'id':STEM,'all_passed':True,
            'actual_CB_Hodd_10_generators_Jacobian_rank':9,
            'actual_selected9_unit_determinant':str(detF),
            'all6_memberships_exact_by_finite_telescoping':True,
            'all6_CB_graph_minus_CA_pairs_in_same_CA_trace_basis_constructed':True,
            'all12_actual_top_T3_pair_finite_algebraic_circuits_constructed':True,
            'nonunit_regular_divisions_only_with_exact_coordinate_restriction_zero_proof':54,
            'no_p_inverse_arbitrary_unknown_membership_coeff_or_series_truncation':True,
            'depth_nil_T3_membership_or_socle_answer_claimed':False,
            'no_source_core_GB_new_point_or_ungauged_R9_rerun':True,
            'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in [Path(__file__),path,folder/'report.json']},
            'evidence_sha256':{name:sha(out/name)for name in ['actual-top12-DAG.json','rank-gate.json']},
            'DAG_nodes':len(D.nodes),'certificate_bytes':len(encoded.encode()),
            'seconds':time.monotonic()-started,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k]for k in ['all_passed','actual_selected9_unit_determinant','DAG_nodes','certificate_bytes','seconds','max_rss_kib']}))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True)
    main(p.parse_args().output)
