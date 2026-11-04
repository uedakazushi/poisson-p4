#!/usr/bin/env python3
"""Independent unit-clone, regular quotient and actual512 inverse audit."""
import argparse,itertools,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
from ambient_B18_D4_literal_horseshoe_review_20261104_v63 import Lookup,det,read,sha
from ambient_B18_honest9_to_N8_review_20261104_v63 import gradients,ancestors
from s61_sparse_generic5_full134_normal_circuit_review_20261104_v63 import plus,times,neg,const,symbol
ROOT=Path(__file__).resolve().parents[1];AUTHOR='AMBIENT_B18_FINITE_TOR63_OPERATOR_20261104_v63';ID=AUTHOR+'_REVIEW'
def total(ps):
    out={}
    for p in ps:out=plus(out,p)
    return out
class QLookup(Lookup):
    def __init__(self,nodes,quotients):
        self.nodes=[tuple(n)for n in nodes];self.idx={n:i for i,n in enumerate(self.nodes)};assert len(self.idx)==len(nodes)
        self.values=[]
        for i,n in enumerate(self.nodes):
            op=n[0]
            if op=='c':v=F(n[1],n[2]);assert n[1:]==(v.numerator,v.denominator)
            elif op=='v':v=F(0)
            elif op=='q':
                assert i in quotients and 0<=n[1]<i and 0<=n[2]<i and isinstance(n[3],str)
                assert self.values[n[1]]==self.values[n[2]]==0;v=quotients[i]
            else:
                assert all(isinstance(j,int)and 0<=j<i for j in n[1:])
                if op=='a':v=self.values[n[1]]+self.values[n[2]]
                elif op=='m':v=self.values[n[1]]*self.values[n[2]]
                elif op=='n':v=-self.values[n[1]]
                else:assert op=='i'and self.values[n[1]];v=1/self.values[n[1]]
            self.values.append(v)
        assert {i for i,n in enumerate(self.nodes)if n[0]=='q'}==set(quotients)
        self.zero=self.const(0);self.one=self.const(1)
def main(out):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False)
    b=ROOT/'reports';folder=b/AUTHOR;hashes={}
    def bind(p):hashes[str(p.relative_to(ROOT))]=sha(p)
    def audit(p):
        d=read(p);assert d['all_passed'];bind(p)
        for n,v in d.get('input_sha256',{}).items():assert sha(ROOT/n)==v,n;hashes[n]=v
        for n,v in d.get('evidence_sha256',{}).items():assert sha(p.parent/n)==v,n;bind(p.parent/n)
        return d
    audit(folder/'report.json')
    for name in ['AMBIENT_B18_FINITE_TOR63_IMAGE_20261104_v63_REVIEW/report.json','AMBIENT_B18_ZERO_NORMALIZER_TENSOR_MODEL_20261104_v63_REVIEW/report.json',
                 'AMBIENT_B18_T0_HONEST_C0_20261104_v63_REVIEW/report.json','AMBIENT_B18_T3_DOUBLE_ANNIHILATOR_LENGTH3_20261104_v63_REVIEW.json']:
        audit(b/name)
    data=read(folder/'literal-Tor63-transform-and27-image-operators.json');checkpoint=read(folder/'prefix-operator-checkpoint.json')
    S=read(b/'AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63/actual-source-coefficients.json')
    J=read(b/'AMBIENT_B18_T3_FULL9_FIRSTJET_20261104_v63/full9-firstjet-certificate.json')
    C=read(b/'AMBIENT_B18_T3_JOINT_CURVE_RESTRICTIONS_20261104_v63/joint-degree3-certificate.json')
    T=read(b/'AMBIENT_B18_FINITE_TOR63_IMAGE_PREFLIGHT_20261104_v63/finite63-proof-and-input-preflight.json')
    nodes=data['DAG_scalar_nodes'];old=S['DAG_scalar_nodes'];assert len(nodes)==209337 and nodes[:len(old)]==old
    assert nodes[:len(checkpoint['DAG_scalar_nodes'])]==checkpoint['DAG_scalar_nodes']
    systems=data['all9_new_unit105_systems'];assert len(systems)==9 and checkpoint['completed_unit_systems']==systems
    regular=data['all90_known_regular_Hadamard_quotients'];assert len(regular)==90
    origins={r['node']:F(r['origin'])for r in regular};D=QLookup(nodes,origins)
    ca,cb=[S['families'][f]for f in ('CA','CB')];rows=[18,19,28,29,60,67,72,63]
    f=[p[1]for p in ca['actual_zero_odd_I0_generators']];g=[p[1]for p in cb['actual_zero_odd_I0_generators']]
    Y=J['both57_source2']['CA']['F9_nodes'];proto=ca['all48_zero_IFT_quadratic_coefficient_equations']+Y+cb['all48_zero_IFT_quadratic_coefficient_equations']
    us=C['all105_unknown_original_node_order'];assert proto==C['all105_equation_original_node_order'] and len(us)==105
    funcs=[f[7]]+g;need=ancestors(old,proto+funcs);assert need==checkpoint['pure_source_needed_node_ids']
    G=gradients(D,proto,us);JJ=[[G[r].get(v,F(0))for v in us]for r in proto]
    assert JJ==[list(map(F,row))for row in C['fixed105_Jacobian']]
    inverse=[list(map(F,row))for row in C['fixed105_Jacobian_inverse']]
    for A,B in ((JJ,inverse),(inverse,JJ)):
        assert all(sum((x*B[k][j]for k,x in enumerate(A[i])if x),F(0))==int(i==j)for i in range(105)for j in range(105))
    determinant=det(JJ);assert determinant==F(C['det_joint105'])
    specs=[('prefix'+str(k),[Y[j]if j<k or j>=7 else D.zero for j in range(9)])for k in range(7)]
    specs += [('axis_r',[D.zero]*7+[Y[7],D.zero]),('axis_s',[D.zero]*8+[Y[8]])]
    maps=[];used=set(us)
    for system,(name,targets)in zip(systems,specs):
        assert system['name']==name and system['F9_target_nodes']==targets and F(system['unit_Jacobian_minor'])==determinant
        bindings=dict(system['source105_bindings']);assert set(bindings)==set(us)and len(set(bindings.values()))==105 and not used&set(bindings.values())
        used.update(bindings.values())
        for i,j in bindings.items():assert D.nodes[j]==('v',f'TOR63_{name}_source_node{i}')
        m={}
        for i in need:
            op,*a=old[i]
            if op=='v':v=bindings[i]
            elif op=='c':v=D.const(F(*a))
            elif op=='a':v=D.add(m[a[0]],m[a[1]])
            elif op=='m':v=D.mul(m[a[0]],m[a[1]])
            elif op=='n':v=D.neg(m[a[0]])
            else:assert op=='i';v=D.inv(m[a[0]])
            m[i]=v
        assert system['old_to_clone_node_map']==[[i,m[i]]for i in need]
        eq=[m[i]for i in ca['all48_zero_IFT_quadratic_coefficient_equations']]+[D.add(m[i],D.neg(t))for i,t in zip(Y,targets)]+[m[i]for i in cb['all48_zero_IFT_quadratic_coefficient_equations']]
        assert eq==system['all105_unit_equations']and all(D.values[i]==0 for i in eq);maps.append(m)
    A0=[list(map(F,r))for r in T['actual_constant16_by10_generator_matrix']]
    grad=[A0[7][:7]]+[r[:7]for r in A0[8:]];quad=[A0[7][7:]]+[r[7:]for r in A0[8:]]
    H=[];plane=[]
    assert len(data['all9_coordinate_telescopes'])==len(data['all9_second_order_Hadamard_identities'])==9
    for n,root in enumerate(funcs):
        tel=data['all9_coordinate_telescopes'][n];assert tel['function_source_node']==root and tel['plane_node']==maps[0][root]and tel['full_node']==root
        assert len(tel['all7_difference_steps'])==7
        h=[]
        for j,record in enumerate(tel['all7_difference_steps']):
            before=maps[j][root];after=maps[j+1][root]if j<6 else root;reg=regular[7*n+j]
            numerator=D.add(after,D.neg(before));ref=reg['node']
            assert record==[before,after,ref,Y[j]]and D.nodes[ref]==('q',numerator,Y[j],f'TOR63_regular_{7*n+j}')
            assert reg['numerator']==numerator and reg['denominator']==Y[j]
            proof=reg['all_order_regularity_proof'];assert proof['type']=='coordinate_prefix_difference'and proof['function_source_node']==root
            assert [proof['before_prefix'],proof['after_prefix'],proof['coordinate_index']]==[j,j+1,j]
            assert F(reg['origin'])==grad[n][j];h.append(ref)
        H.append(h)
        p,pr,ps=maps[0][root],maps[7][root],maps[8][root];nums=[pr,D.add(D.add(p,D.neg(pr)),D.neg(ps)),ps]
        dens=[D.mul(Y[7],Y[7]),D.mul(Y[7],Y[8]),D.mul(Y[8],Y[8])];u=[]
        for j in range(3):
            reg=regular[63+3*n+j];ref=reg['node']
            assert D.nodes[ref]==('q',nums[j],dens[j],f'TOR63_regular_{63+3*n+j}')and F(reg['origin'])==quad[n][j]
            assert reg['numerator']==nums[j]and reg['denominator']==dens[j]
            typ=reg['all_order_regularity_proof']['type'];assert typ==('mixed_axis_difference'if j==1 else'axis_order2');u.append(ref)
        plane.append(u);ident=data['all9_second_order_Hadamard_identities'][n]
        assert [ident['function_source_node'],ident['plane_node'],ident['axis_r_node'],ident['axis_s_node']]==[root,p,pr,ps]and ident['quadratic_Hadamard_coefficients']==u
    assert not any(map(F,C['g_CA_all10_coefficients'][:6]))
    assert all(not any(map(F,p['all_degree_le3_coefficients'][:3]))for p in C['all8_CB_zero_restrictions'])
    # Quotient relation/telescope identities use no inverse of Y,r,s.
    t=[symbol(('level',i))for i in range(8)];assert total([t[0],*[plus(t[i+1],neg(t[i]))for i in range(7)]])==t[7]
    pp,pr,ps=[symbol(('plane',i))for i in range(3)];assert total([pr,pp,neg(pr),neg(ps),ps])==pp
    A=[[D.one if i==j else D.zero for j in range(7)]+[D.zero]*3 for i in range(7)]+[H[0]+plane[0]]+[H[i+1]+plane[i+1]for i in range(8)]
    assert A==data['actual16_by10_generator_expression_matrix']and data['combined_actual16_generator_nodes']==f+g
    selected=T['selected10_combined_indices'];redundant=T['redundant6_combined_indices'];assert selected==data['selected10_indices']and redundant==data['redundant6_indices']
    U=[A[i][7:]for i in selected[7:]];assert U==data['actual_quadratic3_unit_matrix']
    detroot=data['actual_unit3_determinant_node'];assert D.values[detroot]==F(625,768)
    invdet=D.inv(detroot);assert [i for i,n in enumerate(D.nodes)if i>=len(checkpoint['DAG_scalar_nodes'])and n[0]=='i']==[invdet]
    # A embeds 90 actual regular q-functions into a free coefficient ring.
    qrefs=[i for row in A for i in row if D.nodes[i][0]=='q'];assert len(qrefs)==len(set(qrefs))==90 and set(qrefs)==set(origins)
    memo={};tt=('inverse_det',)
    def free(i):
        if i in memo:return memo[i]
        op,*a=D.nodes[i]
        if op=='q':v=symbol(('q',i))
        elif op=='c':v=const(F(*a))
        elif op=='a':v=plus(free(a[0]),free(a[1]))
        elif op=='m':v=times(free(a[0]),free(a[1]))
        elif op=='n':v=neg(free(a[0]))
        else:assert op=='i'and i==invdet;v=symbol(tt)
        memo[i]=v;return v
    uu=[[free(i)for i in row]for row in U]
    def cofactor(i,j):
        rr=[k for k in range(3)if k!=i];cc=[k for k in range(3)if k!=j]
        p=plus(times(uu[rr[0]][cc[0]],uu[rr[1]][cc[1]]),neg(times(uu[rr[0]][cc[1]],uu[rr[1]][cc[0]])))
        return p if(i+j)%2==0 else neg(p)
    delta=total(times(uu[0][j],cofactor(0,j))for j in range(3));assert free(detroot)==delta
    for i in range(3):
        for j in range(3):assert free(data['actual_unit3_inverse'][i][j])==times(symbol(tt),cofactor(j,i))
    transform=data['actual16_by16_unit_transform'];inverse=data['actual16_by16_inverse_transform'];T10=data['canonical10_from_actual16_matrix'];qcycles=data['six_actual_zero_column_cycles16']
    assert transform==T10+qcycles and inverse==[row+[D.one if i==j else D.zero for j in redundant]for i,row in enumerate(A)]
    def clear(p):
        out={}
        for mon,c in p.items():
            power=mon.count(tt);assert power<=1;base=list(mon)
            if power:base.remove(tt)
            term={tuple(base):c};out=plus(out,term if power else times(delta,term))
        return out
    BB=[[clear(free(i))for i in row]for row in transform];CC=[[free(i)for i in row]for row in inverse]
    checks=0
    for X,Y0 in ((BB,CC),(CC,BB)):
        for i in range(16):
            for j in range(16):
                assert total(times(X[i][k],Y0[k][j])for k in range(16))==(delta if i==j else{});checks+=1
    assert checks==512
    assert [[D.values[i]for i in row]for row in qcycles]==[list(map(F,row))for row in T['all6_zero_column_cycle_origins16']]
    r,s=Y[7:];eta_u=[D.add(D.neg(D.mul(s,T10[7][j])),D.mul(r,T10[8][j]))for j in range(16)]
    eta_v=[D.add(D.neg(D.mul(s,T10[8][j])),D.mul(r,T10[9][j]))for j in range(16)]
    assert eta_u==data['eta_u_actual_cycle16']and eta_v==data['eta_v_actual_cycle16']
    assert data['canonical10_functions']==Y[:7]+[D.mul(r,r),D.mul(r,s),D.mul(s,s)]
    Xi=[[D.one if j==i+1 else D.zero for j in range(9)]for i in range(7)]+[[D.zero]+H[0]+[D.one]]
    assert Xi==data['C0_actual_CA8_odd_target_free9_vectors']and data['C0_basis_odd8_nodes']==Y[:7]+[maps[0][f[7]]]
    delta8=[D.add(ca['complete_normalized_physical_w_coefficients'][0][row][0],D.neg(cb['complete_normalized_physical_w_coefficients'][0][row][0]))for row in rows]
    assert delta8==data['all8_actual_even_difference_nodes']
    specs=[('q_wedge_q',a,b,qcycles[a],qcycles[b])for a,b in itertools.combinations(range(6),2)]
    specs +=[(name,k,None,eta,qcycles[k])for name,eta in [('eta_u_wedge_q',eta_u),('eta_v_wedge_q',eta_v)]for k in range(6)]
    images=data['all27_finite_Tor2_image_maps_in_C0_free9'];assert len(images)==27
    for record,(kind,a,b0,q,p)in zip(images,specs):
        assert [record['type'],record['first_cycle_index'],record['second_cycle_index']]==[kind,a,b0]
        h=[[D.add(D.mul(q[i],p[8+j]),D.neg(D.mul(p[i],q[8+j])))for j in range(8)]for i in range(8)]
        assert h==record['h8_by8']
        even=D.mul(ca['actual_tracefree_nu'],D.total(D.mul(h[i][j],D.mul(f[i],f[j]))for i in range(8)for j in range(8)))
        col=[even]+[D.total(D.mul(h[i][j],D.mul(delta8[j],Xi[i][k]))for i in range(8)for j in range(8))for k in range(1,9)]
        assert col==record['actual_C0_free9_image_nodes']
    proof=dict(all209337_scalar_syntax_and_unit_origins_checked=True,all9_source105_homomorphic_clones_and_targets_bound=True,
        independent_original105_AD_and_both_inverse_products=True,unit105_determinant=str(determinant),
        all63_coordinate_and27_second_Hadamard_regular_interfaces_bound=True,
        regularity_from_unit_uniqueness_and_axis_order2_not_zero_denominator_inversion=True,
        actual16_generator_decomposition_bound=True,actual90_coefficients_free_ring_embedding=True,
        actual512_denominator_cleared_inverse_entries_independently_verified=True,
        unit3_origin='625/768',all8_actual_cycles_and27_C0_image_maps_source_bound=True,
        C0_basis_identified_by_CA57_prefix_uniqueness=True,finite_m_R0_squared_image_theorem_reused=True,
        no_source_coefficient_expansion_newjets_author_worker_membership_nonzero_socle_or_depth=True)
    (out/'independent-proof.json').write_text(json.dumps(proof,indent=2)+'\n')
    for p in (Path(__file__),ROOT/'code/ambient_B18_D4_literal_horseshoe_review_20261104_v63.py',ROOT/'code/ambient_B18_honest9_to_N8_review_20261104_v63.py',ROOT/'code/s61_sparse_generic5_full134_normal_circuit_review_20261104_v63.py',ROOT/f'research/notes/{ID}.md'):bind(p)
    report=dict(id=ID,all_passed=True,input_sha256=hashes,evidence_sha256={p.name:sha(p)for p in out.iterdir()},seconds=time.monotonic()-start,
        max_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,scope='All-order typed Tor63 transform and27 actual C0 image operators only; no coefficient expansion or image zero/nonzero/depth decision.')
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(dict(all_passed=True,seconds=report['seconds'],rss=report['max_rss_kib'],report_sha256=sha(out/'report.json'))))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args();main(a.output)
