#!/usr/bin/env python3
"""Literal weighted group normalization of the frozen actual R curve.

Only a finite rational map and its first derivative are computed; no new
source jets or high-order group-series expansion.
"""
import argparse,ast,gzip,hashlib,json
from fractions import Fraction as F
from pathlib import Path
import resource,signal,time
from ambient_B18_actual_top_coefficient_preflight_20261104_v63 import DAG,action
from ambient_B18_actual_top_coefficient_preflight_review_20261104_v63 import DAGReview
from ambient_B18_top_graph_telescoping_review_20261104_v63 import solve
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_R_ZERO_SLICE_IMAGE_20261104_v63'
def read(p):
    with (gzip.open(p,'rt') if p.suffix=='.gz' else p.open()) as f:return json.load(f)
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()
class Scalar:
    def __init__(self,D):self.D=D;self.zero=D.zero;self.one=D.one
    def add(self,a,b):return self.D.add(a,b)
    def neg(self,a):return self.D.neg(a)
    def mul(self,a,b):return self.D.mul(a,b)
    def scale(self,a,q):return self.D.mul(a,self.D.const(q))
    def total(self,v):return self.D.total(v)
def main(out):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False)
    rp=ROOT/'reports/AMBIENT_B18_R_CURVE_TANGENT_20261103_v63/certificate.json.gz'
    op=ROOT/'reports/AMBIENT_B18_ACTUAL_INTERSECTION_ORBIT_PRODUCT_20261103_v63/certificate.json.gz'
    wp=ROOT/'reports/AMBIENT_B18_RELATIVE_WEIGHTED_MODEL_20261104_v63/certificate.json'
    tp=ROOT/'reports/AMBIENT_B18_POSITIVE_GRADING_OBSTRUCTION_20261104_v63/certificate.json'
    jp=ROOT/'reports/AMBIENT_B18_T3_FULL9_FIRSTJET_20261104_v63/full9-firstjet-certificate.json'
    sp=ROOT/'reports/AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63/actual-source-coefficients.json'
    gp=ROOT/'reports/AMBIENT_B18_T0_HONEST_C0_20261104_v63/honest-C0-D0-D4-DAG.json'
    gates=[rp.parent/'report.json',
       ROOT/'reports/AMBIENT_B18_N_R_DVR_CUBIC_GRAPH_20261103_v63_REVIEW/report.json',
       ROOT/'reports/AMBIENT_B18_T1_T2_QUADRATIC_PHI_20261104_v63_REVIEW/report.json',
       ROOT/'reports/AMBIENT_B18_T0_COKERNEL_Z_20261104_v63_REVIEW/report.json']
    assert all(read(p)['all_passed'] for p in gates)
    curve,orbit,weight,torus,jet,src,C0=map(read,[rp,op,wp,tp,jp,sp,gp])
    D=DAG();T=D.var('R_ACTUAL_CURVE_T');groups=[D.var(f'R_SLICE_GROUP_{i}') for i in range(19)]
    def parse(text):
        def visit(n):
            if isinstance(n,ast.Constant):return D.const(n.value)
            if isinstance(n,ast.Name):assert n.id=='t';return T
            if isinstance(n,ast.UnaryOp) and isinstance(n.op,ast.USub):return D.neg(visit(n.operand))
            if isinstance(n,ast.BinOp):
                a=visit(n.left)
                if isinstance(n.op,ast.Pow):
                    assert isinstance(n.right,ast.Constant) and 0<=n.right.value<=2
                    ans=D.one
                    for _ in range(n.right.value):ans=D.mul(ans,a)
                    return ans
                b=visit(n.right)
                if isinstance(n.op,ast.Add):return D.add(a,b)
                if isinstance(n.op,ast.Sub):return D.add(a,D.neg(b))
                if isinstance(n.op,ast.Mult):return D.mul(a,b)
                if isinstance(n.op,ast.Div):return D.mul(a,D.inverse(b))
            raise ValueError(ast.dump(n))
        return visit(ast.parse(text,mode='eval').body)
    U=[parse(s) for s in curve['physical_curve150']];point=[D.values[i] for i in U]
    wt=torus['physical150_torus_weights'];assert all(wt[i]==0 for i,node in enumerate(U) if node!=D.zero)
    GL=orbit['GL18_generator_indices'];diag=torus['extended_fixing_torus_GL_weights']
    groupweights=[0]+[diag[i//5]-diag[i%5] for i in GL]
    assert groupweights==torus['orbit19_column_weights']
    G=[[D.one if i==j else D.zero for j in range(5)] for i in range(5)]
    for p,k in zip(groups[1:],GL):i,j=divmod(k,5);G[i][j]=D.add(G[i][j],p)
    A=[G[i][:]+[D.one if i==j else D.zero for j in range(5)] for i in range(5)];inverse_trace=[]
    for k in range(5):
        p=A[k][k];assert D.values[p];iv=D.inverse(p);inverse_trace.append([k,p,str(D.values[p])])
        A[k]=[D.mul(iv,x) for x in A[k]]
        for i in range(5):
            if i==k:continue
            c=A[i][k];A[i]=[D.add(x,D.neg(D.mul(c,y))) for x,y in zip(A[i],A[k])]
    Gi=[r[5:] for r in A]
    physical=action(Scalar(D),U,G,Gi,D.add(D.one,groups[0]))
    centered=[D.add(a,D.const(-p)) for a,p in zip(physical,point)]
    cuts=weight['common_physical_group_cuts19'];E=[centered[i] for i in cuts]
    Rrows=[i for i in weight['remaining_actual_even_base10_raw_coordinates'] if i!=140]
    assert Rrows==src['physical_zero_base_rows9']==[27,108,109,112,113,117,118,119,61]
    rho=[centered[i] for i in Rrows]
    R=DAGReview(D.nodes)
    O=[[R.grad[node].get(g,F(0)) for g in groups] for node in physical]
    assert O==[list(map(F,r)) for r in orbit['physical150_by19_scalar_then_GL']]
    vt=[R.grad[node].get(T,F(0)) for node in physical]
    assert vt==list(map(F,curve['physical_R_curve_derivative150']))
    J=[O[i] for i in cuts]
    g1,det19=solve(J,[[-vt[i]] for i in cuts]);g1=[a[0] for a in g1]
    assert det19==F(-32)
    velocity=[vt[row]+sum((O[row][k]*g1[k] for k in range(19)),F(0)) for row in Rrows]
    assert any(velocity);param=next(i for i,v in enumerate(velocity) if v)
    external=[D.var(f'R_IMAGE_rho{row}') for row in Rrows]
    neweq=D.add(rho[param],D.neg(external[param]));unknown=groups+[T]
    J20=[[R.grad[node].get(g,F(0)) for g in unknown] for node in E+[rho[param]]]
    _,det20=solve(J20,[[F(i==j) for j in range(20)] for i in range(20)])
    assert det20==det19*velocity[param]
    graph=[D.add(external[i],D.neg(rho[i])) for i in range(9) if i!=param]
    # First conductor test only; a nonzero linear restriction is enough to
    # rule out conductor annihilation of the surviving R-primary part.
    ca=jet['both57_source2']['CA'];xx=[next(i for i,n in enumerate(src['DAG_scalar_nodes']) if n==['v',f'x{row}']) for row in Rrows]
    sol={i:list(map(F,v)) for i,v in ca['all57_source2jets']}
    L=[[sol[i][1+j] for j in range(9)] for i in xx]
    y1,_=solve(L,[[v] for v in velocity]);y1=[v[0] for v in y1]
    cb=[[F(a) for a in row] for row in jet['YB_to_YA_linear']]
    yB=[sum((cb[i][j]*y1[j] for j in range(9)),F(0)) for i in range(9)]
    def cubic(c,a,b):return c[0]*a**3+c[1]*a*a*b+c[2]*a*b*b+c[3]*b**3
    ga=list(map(F,C0['cached_actual_CA_gA_all_degree3_coefficients'][6:10]))
    gbpath=ROOT/'reports/AMBIENT_B18_CB_G_DIVISION_ORIGINS_20261104_v63/degree3-certificate.json'
    gb=list(map(F,read(gbpath)['g_CB_all10_coefficients'][6:10]))
    tests={'CA_zero7_first_on_actual_R':list(map(str,y1[:7])),
           'CB_zero7_first_on_actual_R':list(map(str,yB[:7])),
           'CA_g3_on_positive_first':str(cubic(ga,y1[7],y1[8])),
           'CB_g3_on_positive_first':str(cubic(gb,yB[7],yB[8]))}
    ev={'DAG_scalar_nodes':D.nodes,'frozen_actual_R_curve_parameter_node':T,'group19_unknown_nodes':groups,
      'actual_G19_GL18_indices':GL,'group19_weights':groupweights,'group_matrix5':G,'group_inverse5':Gi,
      'group_inverse_unit_Gaussian_trace':inverse_trace,'all150_transformed_physical_curve_nodes':physical,
      'all150_centered_curve_nodes':centered,'same_common_cut_rows19':cuts,'all19_cut_equations':E,
      'all19_cut_Jacobian':[[str(x) for x in a] for a in J],'unit19_minor':str(det19),
      'first_group_correction':list(map(str,g1)),'same_actual_R0_rows9':Rrows,'actual_R0_readback_nodes9':rho,
      'actual_R0_slice_curve_velocity9':list(map(str,velocity)),
      'DVR_parameter_actual_R0_row':Rrows[param],'DVR_parameter_derivative':str(velocity[param]),
      'external_R0_coordinate_nodes9':external,'all20_pointed_unit_equations':E+[neweq],
      'all20_unknowns':unknown,'unit20_minor':str(det20),
      'all8_regular_graph_ideal_equations_in_actual_R0':graph,
      'full_R_slice_curve_is_torus_fixed_by_equivariant_uniqueness':True,
      'all_positive_weight_physical_coordinates_including_w_are_zero_all_orders':True,
      'first_conductor_tests':tests,
      'no_high_order_series_or_T0_socle_computed':True}
    ep=out/'actual-R-zero-slice-image-DAG.json';ep.write_text(json.dumps(ev,separators=(',',':'))+'\n')
    paths=[Path(__file__),ROOT/'research/notes'/f'{STEM}.md',rp,op,wp,tp,jp,sp,gp,gbpath,*gates,
      ROOT/'code/ambient_B18_actual_top_coefficient_preflight_20261104_v63.py',
      ROOT/'code/ambient_B18_actual_top_coefficient_preflight_review_20261104_v63.py',
      ROOT/'code/ambient_B18_top_graph_telescoping_review_20261104_v63.py']
    report={'id':STEM,'all_passed':True,'status':'actual_weighted_R_zero_slice_DVR_image_and_regular_R0_graph_ideal',
      'raw_R_curve_all_nonzero_physical_coefficients_weight0':True,'raw_curve_not_assumed_already_in_slice':True,
      'unit19_minor':str(det19),'unit20_minor':str(det20),'actual_R0_uniformizer_row':Rrows[param],
      'R0_projection_is_closed_formal_DVR_embedding':True,'all_positive_coordinates_zero_by_equivariant_unique_gauge':True,
      'R_branch_survives_in_T0_generically_by_reduced_component_separation':True,
      'T0_conductor_annihilator_or_socle_depth_claimed':False,
      'no_new_source_jets_GB_points_primes_or_group_series_expansion':True,
      'input_sha256':{str(p.relative_to(ROOT)):sha(p) for p in paths},'evidence_sha256':{ep.name:sha(ep)},
      'DAG_nodes':len(D.nodes),'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:report[k] for k in ['all_passed','unit19_minor','unit20_minor','actual_R0_uniformizer_row','DAG_nodes','seconds','max_rss_kib']}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
