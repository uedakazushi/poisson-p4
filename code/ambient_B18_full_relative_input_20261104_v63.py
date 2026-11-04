#!/usr/bin/env python3
"""ONE finite algebraic-DAG completion of all99+87 actual relative inputs."""
import argparse
from collections import Counter
from fractions import Fraction as F
import gzip
import hashlib
import itertools as it
import json
from pathlib import Path
import resource
import signal
import time
from ambient_B18_actual_top_coefficient_preflight_20261104_v63 import DAG
from ambient_B18_top_graph_telescoping_20261104_v63 import forward,replay,ancestors,linear_solve
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_FULL_RELATIVE_INPUT_20261104_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(gzip.decompress(p.read_bytes())if p.suffix=='.gz'else p.read_text())

class RegularDAG(DAG):
    """A known regular quotient may have an unevaluated, but unique, origin value."""
    def neg(self,a):
        if a==self.zero:return a
        n=self.nodes[a]
        if n[0]=='c':return self.const(-self.values[a])
        if n[0]=='n':return n[1]
        return self.put(('n',a),None if self.values[a]is None else-self.values[a])
    def add(self,a,b):
        if a==self.zero:return b
        if b==self.zero:return a
        if self.nodes[a]==('n',b)or self.nodes[b]==('n',a):return self.zero
        if self.nodes[a][0]==self.nodes[b][0]=='c':return self.const(self.values[a]+self.values[b])
        a,b=sorted((a,b));v=None if self.values[a]is None or self.values[b]is None else self.values[a]+self.values[b]
        return self.put(('a',a,b),v)
    def mul(self,a,b):
        if a==self.zero or b==self.zero:return self.zero
        if a==self.one:return b
        if b==self.one:return a
        if self.nodes[a][0]==self.nodes[b][0]=='c':return self.const(self.values[a]*self.values[b])
        a,b=sorted((a,b));v=None if self.values[a]is None or self.values[b]is None else self.values[a]*self.values[b]
        return self.put(('m',a,b),v)

def main(out):
    started=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));out.mkdir(parents=True,exist_ok=False)
    source=ROOT/'reports/AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63/actual-source-coefficients.json'
    topgate=ROOT/'reports/AMBIENT_B18_TOP_GRAPH_TELESCOPING_20261104_v63/rank-gate.json'
    graphpath=ROOT/'reports/AMBIENT_B18_NOETHER_FIBER_EVEN_QUADRATIC_20261103_v63/certificate.json.gz'
    weightpath=ROOT/'reports/AMBIENT_B18_POSITIVE_GRADING_OBSTRUCTION_20261104_v63/certificate.json'
    mupath=ROOT/'reports/AMBIENT_CB_ORIGIN_CONDUCTOR_MULTIPLICITY4_20261103_v63_REVIEW/report.json'
    data,gate,graph,weights,mureview=map(read,[source,topgate,graphpath,weightpath,mupath])
    assert mureview['all_passed']and mureview['CB_actual_I_minimal10_and_image_depth21_edim39_e2_accepted']
    nodes=data['DAG_scalar_nodes'];ca,cb=data['families']['CA'],data['families']['CB'];varids={n[1]:i for i,n in enumerate(nodes)if n[0]=='v'}
    xids=[varids[f'x{i}']for i in data['physical_zero_base_rows9']]
    zids=[i for p in cb['coefficient_unknown_pairs'].values()for i in p]+[cb['p_quadratic_relation']['s'],cb['p_quadratic_relation']['t']]
    vals,grads=forward(nodes,set(xids+zids));zx=[[F(x)for x in r]for r in gate['source48_dz_dx']]
    def deriv(node):
        answer=[grads[node].get(i,F(0))for i in xids]
        for k,i in enumerate(zids):
            q=grads[node].get(i,F(0))
            if q:answer=[a+q*b for a,b in zip(answer,zx[k])]
        return answer
    gens=cb['actual_zero_odd_I0_generators']+cb['actual_positive_odd_J0_generators'];gids=[i for row,i in gens]
    selected=gate['selected_generator_indices'];assert selected==[0,1,2,3,4,5,6,8,9]
    inverse,detF=linear_solve([deriv(gids[i])for i in selected],[[F(int(i==j))for j in range(9)]for i in range(9)])
    assert detF==F(1,96)
    def derivY(node):
        dx=deriv(node);return[sum(dx[k]*inverse[k][j]for k in range(9))for j in range(9)]
    rows=graph['independent87_physical_graph_rows'];ws=weights['physical150_torus_weights']
    assert len(rows)==87 and Counter(ws[i]for i in rows)=={0:4,1:53,2:24,3:6}
    assert graph['common_odd10_physical_rows']==[r for r,i in ca['actual_zero_odd_I0_generators'][:7]]+[46,95,63]
    Bnodes={row:cb['complete_normalized_physical_w_coefficients'][ws[row]][row][1]for row in rows}
    assert all(vals[i]==0 for i in Bnodes.values())
    f63=gids[7];assert gens[7][0]==63 and vals[f63]==0
    D=RegularDAG();original=replay(D,nodes,{});assert all(original[i]==i for i in range(len(nodes)))
    Y=[gids[i]for i in selected];r,s=Y[7:]
    allset=tuple(range(9));empty=tuple()
    subsets=set(tuple(range(i))for i in range(1,9))
    subsets.update(tuple(list(range(i))+[7,8])for i in range(7));subsets.add((8,))
    roots=cb['all48_zero_IFT_quadratic_coefficient_equations']+[gids[i]for i in selected]+list(Bnodes.values())+[f63]
    needed=ancestors(nodes,roots)
    instances={allset:{'mapping':original},empty:{'mapping':None}};equation_records=[]
    for number,subset in enumerate(sorted(subsets),1):
        bindings={i:D.var(f'FULLINV{number}_{nodes[i][1]}')for i in xids+zids}
        mapping=replay(D,nodes,bindings,needed)
        equations=[mapping[i]for i in cb['all48_zero_IFT_quadratic_coefficient_equations']]
        equations +=[D.add(mapping[gids[index]],D.neg(Y[j]if j in subset else D.zero))for j,index in enumerate(selected)]
        assert len(equations)==57 and all(D.values[i]==0 for i in equations)
        instances[subset]={'mapping':mapping}
        equation_records.append({'retained_Y_indices':list(subset),'bindings':[[i,j]for i,j in bindings.items()],
                                 'all57_pointed_equations':equations,'unit_Jacobian':str(F(gate['source48_Jacobian_determinant'])*detF)})
    def at(node,subset):
        subset=tuple(sorted(subset))
        if not subset:assert vals[node]==0;return D.zero
        return instances[subset]['mapping'][node]
    divisions=[]
    def coordinate_div(num,coordinate,label,value):
        node=D.put(('r',num,coordinate,label),value)
        divisions.append({'node':node,'numerator':num,'coordinate':coordinate,'proof':'restriction-zero by paired unique unit57 inverse instances','origin_value':str(value)})
        return node
    tail=(7,8);g=at(f63,tail)
    # Zero7 telescoping with the final two coordinates retained throughout.
    def telescope_zero7(node,label):
        previous=at(node,tail);coeff=[];dy=derivY(node)
        for i in range(7):
            current=at(node,tuple(list(range(i+1))+[7,8]))
            coeff.append(coordinate_div(D.add(current,D.neg(previous)),Y[i],f'{label}:z{i}',dy[i]));previous=current
        return coeff
    e=telescope_zero7(f63,'f63')
    g_s=at(f63,(8,))
    AA=coordinate_div(D.add(g,D.neg(g_s)),r,'g_split_r',F(0))
    BB=coordinate_div(g_s,s,'g_split_s',F(0))
    memberships={};graph_pairs=[];graph_by_weight=Counter();regular_g_divs=[]
    for row in rows:
        degree=ws[row];node=Bnodes[row];hs=[D.zero]*10
        if degree:
            previous=D.zero;dy=derivY(node)
            for i,index in enumerate(selected):
                current=at(node,tuple(range(i+1)))
                hs[index]=coordinate_div(D.add(current,D.neg(previous)),Y[i],f'graph{row}:Y{i}',dy[i]);previous=current
        else:
            bbar=at(node,tail)
            h8=D.put(('g',bbar,g,f'graph{row}:known_I0_membership'),None)
            regular_g_divs.append({'node':h8,'numerator':bbar,'denominator_g':g,
                                  'proof':'actual zero-weight output odd coefficient belongs I0=(z7,f63); restriction gives (g)',
                                  'origin_value':'not evaluated; no unit use'})
            hi=telescope_zero7(node,f'zero_graph{row}')
            for i in range(7):hs[i]=D.add(hi[i],D.neg(D.mul(h8,e[i])))
            hs[7]=h8
        memberships[row]=hs
        Acb,Bcb=cb['complete_normalized_physical_w_coefficients'][degree][row]
        Aca,Bca=ca['complete_normalized_physical_w_coefficients'][degree][row]
        even=D.add(Acb,D.neg(Aca));odd=D.neg(Bca)
        for i,(raw,fi)in enumerate(ca['actual_zero_odd_I0_generators']):
            aA=ca['complete_normalized_physical_w_coefficients'][0][raw][0]
            aB=cb['complete_normalized_physical_w_coefficients'][0][raw][0]
            even=D.add(even,D.mul(hs[i],D.add(aA,D.neg(aB))));odd=D.add(odd,D.mul(hs[i],fi))
        for j,(raw,ji)in enumerate(ca['actual_positive_odd_J0_generators']):
            aA=ca['complete_normalized_physical_w_coefficients'][1][raw][0]
            aB=cb['complete_normalized_physical_w_coefficients'][1][raw][0]
            even=D.add(even,D.mul(hs[8+j],D.add(aA,D.neg(aB))));odd=D.add(odd,D.mul(hs[8+j],ji))
        graph_pairs.append({'kind':'graph_CB_minus_CA','raw_row':row,'weight':degree,'normalized_CA_trace_pair':[even,odd]});graph_by_weight[degree]+=1
    assert graph_by_weight=={0:4,1:53,2:24,3:6}
    # Actual unipotent odd generator change, evaluated on CA physical coordinates.
    oddraw=[raw for raw,node in ca['actual_zero_odd_I0_generators']]+[raw for raw,node in ca['actual_positive_odd_J0_generators']]
    oo=[]
    for i,raw in enumerate(oddraw):
        degree=int(i>=8)
        aA,bA=ca['complete_normalized_physical_w_coefficients'][degree][raw]
        aB,bB=cb['complete_normalized_physical_w_coefficients'][degree][raw]
        oo.append((D.add(aA,D.neg(aB)),bA))
    def pa(a,b):return(D.add(a[0],b[0]),D.add(a[1],b[1]))
    def pn(a):return(D.neg(a[0]),D.neg(a[1]))
    def ps(a,z):return(D.mul(a[0],z),D.mul(a[1],z))
    nuA=ca['actual_tracefree_nu'];nuB=cb['actual_tracefree_nu']
    def pm(a,b):return(D.add(D.mul(a[0],b[0]),D.mul(nuA,D.mul(a[1],b[1]))),D.add(D.mul(a[0],b[1]),D.mul(a[1],b[0])))
    Z=oo.copy()
    for i in range(7):Z[7]=pa(Z[7],pn(ps(oo[i],e[i])))
    bn=Y[:7]+[g,r,s];kernel=[]
    def save(label,degree,pair):kernel.append({'kind':'CB_kernel','label':label,'weight':degree,'normalized_CA_trace_pair':list(pair)})
    for i,j in it.combinations(range(7),2):save(f'regular_Koszul_{i}_{j}',0,pa(ps(Z[j],bn[i]),pn(ps(Z[i],bn[j]))))
    for i in range(7):
        for j in [7,8,9]:save(f'cross_{i}_{j}',int(j>=8),pa(ps(Z[j],bn[i]),pn(ps(Z[i],bn[j]))))
    save('HB_wg_minus_A_wr_minus_B_ws',1,pa(pa(Z[7],pn(ps(Z[8],AA))),pn(ps(Z[9],BB))))
    save('HB_s_wr_minus_r_ws',1,pa(ps(Z[8],s),pn(ps(Z[9],r))))
    assert len(kernel)==44
    for i in range(10):
        for j in range(i,10):
            pair=pm(Z[i],Z[j]);pair=(D.add(pair[0],D.neg(D.mul(nuB,D.mul(bn[i],bn[j])))),pair[1])
            save(f'product_{i}_{j}',int(i>=8)+int(j>=8),pair)
    assert len(kernel)==99 and Counter(x['weight']for x in kernel)=={0:64,1:32,2:3}
    Qgens=kernel+graph_pairs;assert len(Qgens)==186 and Counter(x['weight']for x in Qgens)=={0:68,1:85,2:27,3:6}
    # Fully literal finite R0^2 matrices for the four-module criterion.
    H=[(D.one,D.zero)]+[(D.zero,i)for raw,i in ca['actual_zero_odd_I0_generators']+ca['actual_positive_odd_J0_generators']]
    C0=H[:9]
    def span(indices,multipliers):return[list(pm(tuple(Qgens[i]['normalized_CA_trace_pair']),h))for i in indices for h in multipliers]
    idx={d:[i for i,q in enumerate(Qgens)if q['weight']==d]for d in range(4)}
    Jinf=span(list(range(186)),H)
    Qn={str(n):span([i for d in range(n)for i in idx[d]],H)+span(idx[n],C0)for n in range(4)}
    assert [len(Qn[str(n)])for n in range(4)]==[612,1513,1926,2034] and len(Jinf)==2046
    top=[list(pm(tuple(q['normalized_CA_trace_pair']),h))for q in graph_pairs if q['weight']==3 for h in H[-2:]]
    assert len(top)==12
    certificate={'DAG_scalar_nodes':D.nodes,
                 'r_operator':'known regular coordinate quotient, restriction-zero certified',
                 'g_operator':'known regular quotient by nonzero g; regularity follows actual zero-output membership in I0; centre intentionally not evaluated',
                 'original_unit48_equations_CA':ca['all48_zero_IFT_quadratic_coefficient_equations'],
                 'original_unit48_equations_CB':cb['all48_zero_IFT_quadratic_coefficient_equations'],
                 'all16_pointed_unit57_inverse_instances':equation_records,
                 'CB_coordinate_Y9_nodes':Y,'CB_f63_node':f63,'CB_g_restriction_node':g,
                 'CB_f63_zero7_telescoping_coefficients':e,
                 'CB_g_split_rA_plus_sB_nodes':[AA,BB],
                 'nonzero_g_proof':'g=0 would reduce I0+wJ to9 generators, contradict independently accepted actual muI=10',
                 'all_coordinate_divisions':divisions,'all4_known_regular_g_divisions':regular_g_divs,
                 'all87_actual_odd_membership_lifts':[[row,hs]for row,hs in memberships.items()],
                 'CB_new_odd_generators10_evaluated_on_CA_normalized_pairs':[list(z)for z in Z],
                 'CB_scalar_odd_generators_normalized_coefficients':bn,
                 'CB_scalar_odd_generator_weights':[0]*8+[1,1],
                 'CB_HB_matrix':'[[w,0],[-A,s],[-B,-r]], minors(g,-wr,ws)',
                 'all99_CB_kernel_and87_graph_generators_evaluated_on_CA':Qgens,
                 'actual_CA_C0_generators':[list(h)for h in C0],'actual_CA_H_generators':[list(h)for h in H],
                 'actual_Qn_R0_squared_generator_matrices':Qn,
                 'actual_Jinf_R0_squared_generator_matrix':Jinf,
                 'actual12_top_T3_additional_columns':top,
                 'no_module_membership_colon_socle_or_depth_computed':True}
    encoded=json.dumps(certificate,separators=(',',':'))+'\n';assert len(encoded.encode())<=32*1024*1024
    (out/'full-relative186-DAG.json').write_text(encoded)
    paths=[source,topgate,graphpath,weightpath,mupath,Path(__file__)]
    report={'id':STEM,'all_passed':True,'status':'complete_actual_relative_input_not_depth_solution',
            'all99_CB_kernel_inputs_complete_by_regular7_and_HB_plus55_products':True,
            'all87_actual_same_target_graph_differences_constructed':True,
            'actual_Q_generator_count':186,'Q_weight_distribution':dict(Counter(x['weight']for x in Qgens)),
            'four_finite_R0_squared_Qn_column_counts':[612,1513,1926,2034],
            'Jinf_column_count':2046,'additional_top_T3_column_count':12,
            'noncoordinate_regular_g_divisions4_have_proved_graphmembership_premise':True,
            'g_division_origin_constants_evaluated':False,
            'full_nil_depth_socle_T3_membership_or_primary_answer_claimed':False,
            'no_source_core_GB_new_point_series_truncation_or_p_inverse':True,
            'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},
            'evidence_sha256':{'full-relative186-DAG.json':sha(out/'full-relative186-DAG.json')},
            'DAG_nodes':len(D.nodes),'certificate_bytes':len(encoded.encode()),
            'seconds':time.monotonic()-started,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:report[k]for k in ['all_passed','actual_Q_generator_count','DAG_nodes','certificate_bytes','seconds','max_rss_kib']}))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True)
    main(p.parse_args().output)
