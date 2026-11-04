#!/usr/bin/env python3
"""Independent literal R source, finite group, unit curve and weight audit."""
import argparse,ast,gzip,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
from ambient_B18_actual_top_coefficient_preflight_review_20261104_v63 import DAGReview,lie
from ambient_B18_W1_tail_binding_review_20261104_v63 import push,PAIRS,MONS,PI,MI
from ambient_W1_S61_intersection_recognition_review_20261103_v63 import inverse
ROOT=Path(__file__).resolve().parents[1]
BASE='AMBIENT_B18_R_ZERO_SLICE_IMAGE_20261104_v63';ID=BASE+'_REVIEW'
def read(p):return json.loads(gzip.decompress(Path(p).read_bytes()) if str(p).endswith('.gz') else Path(p).read_text())
def sha(p):import hashlib;return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def mat(A):return [[F(v) for v in r] for r in A]
def mm(A,B):return [[sum((a*b for a,b in zip(r,c)),F(0)) for c in zip(*B)] for r in A]
def mv(A,b):return [sum((x*y for x,y in zip(r,b)),F(0)) for r in A]
class Scalar:
    def __init__(self,D):self.D=D
    def add(self,a,b):return self.D.add(a,b)
    def neg(self,a):return self.D.neg(a)
    def mul(self,a,b):return self.D.mul(a,b)
    def total(self,items):return self.D.total(items)
def padd(a,b):
    q=dict(a)
    for e,c in b.items():q[e]=q.get(e,F(0))+c
    return {e:c for e,c in q.items() if c}
def pmul(a,b):
    q={}
    for e,c in a.items():
        for f,d in b.items():q[e+f]=q.get(e+f,F(0))+c*d
    return {e:c for e,c in q.items() if c}
def main(out):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False)
    folder=ROOT/'reports'/BASE;pr=read(folder/'report.json');assert pr['all_passed'];hashes={}
    for p,h in pr['input_sha256'].items():assert sha(ROOT/p)==h;hashes[p]=h
    for p,h in pr['evidence_sha256'].items():assert sha(folder/p)==h;hashes[str((folder/p).relative_to(ROOT))]=h
    data=read(folder/'actual-R-zero-slice-image-DAG.json');curve=read(ROOT/'reports/AMBIENT_B18_R_CURVE_TANGENT_20261103_v63/certificate.json.gz');orbit=read(ROOT/'reports/AMBIENT_B18_ACTUAL_INTERSECTION_ORBIT_PRODUCT_20261103_v63/certificate.json.gz');torus=read(ROOT/'reports/AMBIENT_B18_POSITIVE_GRADING_OBSTRUCTION_20261104_v63/certificate.json');weight=read(ROOT/'reports/AMBIENT_B18_RELATIVE_WEIGHTED_MODEL_20261104_v63/certificate.json')
    D=DAGReview(data['DAG_scalar_nodes']);T=data['frozen_actual_R_curve_parameter_node'];groups=data['group19_unknown_nodes'];assert len(groups)==19 and D.nodes[T]==('v','R_ACTUAL_CURVE_T')
    for i,g in enumerate(groups):assert D.nodes[g]==('v',f'R_SLICE_GROUP_{i}')
    # Independent full physical construction from e_y wedge V, not a derivative.
    source=[{} for _ in range(150)]
    def insert(vector_index,mon,degree,coefficient):
        pair=tuple(sorted((2,vector_index)));row=15*PI[pair]+MI[tuple(sorted(mon))];sign=1 if 2<vector_index else -1;source[row][degree]=source[row].get(degree,F(0))+sign*F(coefficient)
    for r,m,e,c in [(0,(1,1),0,1),(0,(0,1),1,-13),(0,(1,3),2,25),(1,(1,1),1,12),(3,(0,1),0,1),(3,(1,3),1,-13),(4,(0,0),0,1),(4,(3,3),0,1),(4,(1,4),1,2)]:insert(r,m,e,c)
    def parse(node):
        if isinstance(node,ast.Constant):return ({0:F(node.value)} if node.value else {},D.const(node.value))
        if isinstance(node,ast.Name):assert node.id=='t';return {1:F(1)},T
        if isinstance(node,ast.UnaryOp):assert isinstance(node.op,ast.USub);p,q=parse(node.operand);return {e:-c for e,c in p.items()},D.neg(q)
        assert isinstance(node,ast.BinOp);p,q=parse(node.left)
        if isinstance(node.op,ast.Pow):
            assert isinstance(node.right,ast.Constant) and 0<=node.right.value<=2;z={0:F(1)};v=D.one
            for _ in range(node.right.value):z=pmul(z,p);v=D.mul(v,q)
            return z,v
        p2,q2=parse(node.right)
        if isinstance(node.op,ast.Add):return padd(p,p2),D.add(q,q2)
        if isinstance(node.op,ast.Sub):return padd(p,{e:-c for e,c in p2.items()}),D.add(q,D.neg(q2))
        if isinstance(node.op,ast.Mult):return pmul(p,p2),D.mul(q,q2)
        assert isinstance(node.op,ast.Div) and set(p2)=={0};return {e:c/p2[0] for e,c in p.items()},D.mul(q,D.inverse(q2))
    U=[]
    for i,text in enumerate(curve['physical_curve150']):p,q=parse(ast.parse(text,mode='eval').body);assert p==source[i];U.append(q)
    point=[p.get(0,F(0)) for p in source];velocity0=[p.get(1,F(0)) for p in source];assert source[140]=={}
    GL=orbit['GL18_generator_indices'];assert GL==data['actual_G19_GL18_indices'];G=[[D.one if i==j else D.zero for j in range(5)] for i in range(5)]
    for g,k in zip(groups[1:],GL):i,j=divmod(k,5);G[i][j]=D.add(G[i][j],g)
    assert G==data['group_matrix5'];A=[row[:]+[D.one if i==j else D.zero for j in range(5)] for i,row in enumerate(G)];trace=[]
    for j in range(5):
        pivot=A[j][j];assert D.values[pivot];trace.append([j,pivot,str(D.values[pivot])]);inv=D.inverse(pivot);A[j]=[D.mul(inv,q) for q in A[j]]
        for i in range(5):
            if i==j:continue
            c=A[i][j];A[i]=[D.add(x,D.neg(D.mul(c,y))) for x,y in zip(A[i],A[j])]
    Gi=[row[5:] for row in A];assert Gi==data['group_inverse5'] and trace==data['group_inverse_unit_Gaussian_trace']
    physical=[D.mul(D.add(D.one,groups[0]),q) for q in push(Scalar(D),U,G,Gi)];assert physical==data['all150_transformed_physical_curve_nodes']
    centered=[D.add(q,D.const(-p)) for q,p in zip(physical,point)];assert centered==data['all150_centered_curve_nodes']
    O=[point]+[lie(point,*divmod(k,5)) for k in GL];O=[list(row) for row in zip(*O)]
    assert O==mat(orbit['physical150_by19_scalar_then_GL']);assert [[D.grad[q].get(g,F(0)) for g in groups] for q in physical]==O
    assert [D.grad[q].get(T,F(0)) for q in physical]==velocity0==list(map(F,curve['physical_R_curve_derivative150']))
    cuts=weight['common_physical_group_cuts19'];assert cuts==data['same_common_cut_rows19'];equations=[centered[i] for i in cuts];assert equations==data['all19_cut_equations']
    J=[O[i] for i in cuts];Ji,d19=inverse(J);assert J==mat(data['all19_cut_Jacobian']) and d19==F(-32)==F(data['unit19_minor'])
    correction=mv(Ji,[-velocity0[i] for i in cuts]);assert correction==list(map(F,data['first_group_correction']))
    rows=[27,108,109,112,113,117,118,119,61];assert rows==data['same_actual_R0_rows9'];rho=[centered[i] for i in rows];assert rho==data['actual_R0_readback_nodes9']
    velocity=[velocity0[i]+sum((a*b for a,b in zip(O[i],correction)),F(0)) for i in rows]
    assert velocity==list(map(F,data['actual_R0_slice_curve_velocity9']))==[0,0,0,-2,0,-12,0,0,0]
    parameter=3;assert data['DVR_parameter_actual_R0_row']==112 and F(data['DVR_parameter_derivative'])==-2
    unknown=groups+[T];assert data['all20_unknowns']==unknown
    J20=[[D.grad[q].get(g,F(0)) for g in unknown] for q in equations+[rho[parameter]]];_,d20=inverse(J20);assert d20==F(64)==d19*velocity[parameter]==F(data['unit20_minor'])
    external=data['external_R0_coordinate_nodes9'];assert data['all20_pointed_unit_equations']==equations+[D.add(rho[parameter],D.neg(external[parameter]))]
    assert data['all8_regular_graph_ideal_equations_in_actual_R0']==[D.add(external[i],D.neg(rho[i])) for i in range(9) if i!=parameter]
    # Exact torus weights of the entire finite rational map, not just tangent.
    diagonal=[0,0,-1,0,0];weights=[1+diagonal[i]+diagonal[j]-diagonal[a]-diagonal[b] for i,j in PAIRS for a,b in MONS]
    assert weights==torus['physical150_torus_weights'];assert all(weights[i]==0 for i,p in enumerate(source) if p)
    gw=[0]+[diagonal[i]-diagonal[j] for i,j in map(lambda k:divmod(k,5),GL)];assert gw==data['group19_weights']==torus['orbit19_column_weights'];assert Counter(gw)=={0:15,1:4}
    variableweights={T:0,**dict(zip(groups,gw)),**{x:0 for x in external}};nodeweights=[]
    for i,node in enumerate(D.nodes):
        o=node[0]
        if o=='c':w=0 if node[1] else None
        elif o=='v':w=variableweights[i]
        elif o=='a':
            x,y=nodeweights[node[1]],nodeweights[node[2]];assert x is None or y is None or x==y;w=y if x is None else x
        elif o=='m':
            x,y=nodeweights[node[1]],nodeweights[node[2]];w=None if x is None or y is None else x+y
        elif o=='n':w=nodeweights[node[1]]
        else:assert o=='i' and nodeweights[node[1]]==0;w=0
        nodeweights.append(w)
    assert all(nodeweights[q] in [None,weights[i]] for i,q in enumerate(physical));assert all(nodeweights[e] in [None,weights[row]] for e,row in zip(equations,cuts))
    assert all(not correction[i] for i,w in enumerate(gw) if w)
    # Cached actual unit57 first-coordinate adapters, only linear algebra.
    jet=read(ROOT/'reports/AMBIENT_B18_T3_FULL9_FIRSTJET_20261104_v63/full9-firstjet-certificate.json');src=read(ROOT/'reports/AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63/actual-source-coefficients.json');xids=[next(i for i,n in enumerate(src['DAG_scalar_nodes']) if n==['v',f'x{row}']) for row in rows]
    L=[]
    for name in ['CA','CB']:
        sol={int(i):list(map(F,v)) for i,v in jet['both57_source2'][name]['all57_source2jets']};L.append([[sol[i][1+j] for j in range(9)] for i in xids])
    LAi,_=inverse(L[0]);adapter=mat(jet['YB_to_YA_linear']);assert mm(L[1],adapter)==L[0]
    ya=mv(LAi,velocity);yb=mv(adapter,ya);assert ya[:7]==yb[:7]==[0]*7
    c0=read(ROOT/'reports/AMBIENT_B18_T0_HONEST_C0_20261104_v63/honest-C0-D0-D4-DAG.json');ga=list(map(F,c0['cached_actual_CA_gA_all_degree3_coefficients']));gb=list(map(F,read(ROOT/'reports/AMBIENT_B18_CB_G_DIVISION_ORIGINS_20261104_v63/degree3-certificate.json')['g_CB_all10_coefficients']));assert ga[:6]==gb[:6]==[0]*6
    def cubic(coeff,v):r,s=v[-2:];return coeff[6]*r**3+coeff[7]*r*r*s+coeff[8]*r*s*s+coeff[9]*s**3
    gA=cubic(ga,ya);gB=cubic(gb,yb);assert gA==7500 and gB==F(-1280000,9)
    tests=data['first_conductor_tests'];assert list(map(F,tests['CA_zero7_first_on_actual_R']))==ya[:7] and list(map(F,tests['CB_zero7_first_on_actual_R']))==yb[:7] and F(tests['CA_g3_on_positive_first'])==gA and F(tests['CB_g3_on_positive_first'])==gB
    gates=['AMBIENT_B18_ACTUAL_INTERSECTION_ORBIT_PRODUCT_20261103_v63_REVIEW.json','AMBIENT_B18_N_R_DVR_CUBIC_GRAPH_20261103_v63_REVIEW/report.json','AMBIENT_B18_W1_POLYNOMIAL_LAURENT_SCOPE_20261104_v63.json','AMBIENT_B18_T0_COKERNEL_Z_20261104_v63_REVIEW/report.json']
    for name in gates:
        p=ROOT/'reports'/name;assert read(p)['all_passed'];hashes[str(p.relative_to(ROOT))]=sha(p)
    proof=dict(actual_e_y_wedge_V_all150_source_coefficients_rebuilt=True,all2194_finite_DAG_nodes_and_group_routing=len(D.nodes)==2194,all150_transformed_outputs_and19_Lie_columns_bound=True,cut19_det='-32',augmented20_det='64',actual_R0_velocity=list(map(str,velocity)),actual_uniformizer_row=112,full_finite_torus_weights_not_only_firstjets=True,all_order_positive_coordinates_zero_by_unique_equivariant_cut_solution=True,regular_graph_ideal8_over_same_R0=True,CA_g3_along_R=str(gA),CB_g3_along_R=str(gB),R_generic_T0_support_not_generic_length2=True,D0_R_prime_not_automatically_R0_projected_prime=True,no_missed_N_projection_or_T0_length_depth_socle_claim=True)
    (out/'independent-proof.json').write_text(json.dumps(proof,indent=2)+'\n')
    paths=[Path(__file__),ROOT/'code/ambient_B18_W1_tail_binding_review_20261104_v63.py',ROOT/'code/ambient_B18_actual_top_coefficient_preflight_review_20261104_v63.py',ROOT/f'research/notes/{ID}.md',folder/'report.json']
    for p in paths:hashes[str(p.relative_to(ROOT))]=sha(p)
    report=dict(id=ID,all_passed=True,actual_R_zero_slice_DVR_and_weighted_normalization_independently_verified=True,R_generic_support_in_T0_but_no_length2_promotion=True,input_sha256=hashes,evidence_sha256={'independent-proof.json':sha(out/'independent-proof.json')},seconds=time.monotonic()-start,max_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(dict(all_passed=True,seconds=report['seconds'],rss=report['max_rss_kib'],report_sha256=sha(out/'report.json'))))
if __name__=='__main__':
    from collections import Counter
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
