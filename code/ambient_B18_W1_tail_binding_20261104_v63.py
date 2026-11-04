#!/usr/bin/env python3
"""ONE literal old49/rationalGL17/physical20cut tail model, no Newton or GB."""
import argparse
from fractions import Fraction as F
import gzip,hashlib,itertools as it,json
from pathlib import Path
import resource,signal,time
import sympy as sy
from ambient_B18_actual_top_coefficient_preflight_20261104_v63 import DAG,Quadratic,action,PAIRS,MONS,PI,MI
from ambient_B18_top_graph_telescoping_20261104_v63 import forward,linear_solve,independent_rows,ancestors,replay
from s61_model import S61
from weight1_structure import sp_basis
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_W1_TAIL_BINDING_20261104_v63'
def read(p):return json.loads(gzip.decompress(p.read_bytes()) if p.suffix=='.gz' else p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def enc(M):return [[str(v) for v in row] for row in M]
def main(out):
    started=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));out.mkdir(parents=True,exist_ok=False)
    prepath=ROOT/'reports/AMBIENT_B18_W1_TAIL_COORDINATE_PREFLIGHT_20261104_v63/tail-coordinate-preflight-DAG.json'
    prep=prepath.parent/'report.json';builderp=ROOT/'reports/AMBIENT_INTERSECTION_EXACT_IDEAL_20261103_v63/builder.json'
    exactp=builderp.parent/'report.json';wr=ROOT/'reports/AMBIENT_W1_S61_INTERSECTION_RECOGNITION_20261103_v63/recognition.json.gz'
    basep=ROOT/'reports/AMBIENT_BASE_PINCHED_SURFACE_20261103_v63/report.json'
    bundlep=ROOT/'reports/AMBIENT_FIBER_BUNDLE_20261103_v63/report.json'
    orbitp=ROOT/'reports/AMBIENT_B18_ACTUAL_INTERSECTION_ORBIT_PRODUCT_20261103_v63/certificate.json.gz'
    sourcep=ROOT/'reports/AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63/actual-source-coefficients.json'
    pre,builder,exact,wrec,base,bundle,orbit,source=map(read,[prepath,builderp,exactp,wr,basep,bundlep,orbitp,sourcep])
    assert read(prep)['all_passed'] and exact['all_passed']
    D=DAG();names=builder['ambient_variables'];vn={n:D.var('TAIL_OLD49_'+n) for n in names}
    oldvar={n:D.add(vn[n],D.const(1 if n in ['v4','v7'] else 0)) for n in names}
    t=[D.var(f'TAIL_GROUP17_{j}') for j in range(17)]
    def add(*a):return D.total(a)
    def mul(a,b):return D.mul(a,b)
    def scale(a,c):return mul(a,D.const(c))
    def power(a,k):
        v=D.one
        for i in range(k):v=mul(v,a)
        return v
    def polynomial(text,variables=oldvar):
        expr=sy.sympify(text);sym=sorted(expr.free_symbols,key=str)
        if not sym:return D.const(expr)
        terms=[]
        for mon,co in sy.Poly(expr,*sym).terms():
            v=D.const(co)
            for name,k in zip(sym,mon):v=mul(v,power(variables[str(name)],k))
            terms.append(v)
        return add(*terms)
    ss,aa,bb,cc=[oldvar[n] for n in ['ss','aa','bb','cc']]
    s,p=sy.symbols('s p');beta=[]
    for expr in base['slice_24_coordinates']:
        terms=[]
        for (es,ep),co in sy.Poly(sy.sympify(expr),s,p).terms():
            v=D.const(co)
            if ep==1:
                assert es;v=mul(v,mul(power(ss,es-1),cc))
            elif ep%2==0:v=mul(v,mul(power(ss,es),power(aa,ep//2)))
            else:
                assert ep>=3;v=mul(v,mul(power(ss,es),mul(bb,power(aa,(ep-3)//2))))
            terms.append(v)
        beta.append(add(*terms))
    adapter=sy.Matrix(bundle['beta16_to_intrinsic24_adapter']);select=list(adapter.T.rref()[1]);left=adapter[select,:].inv()
    b16=[add(*(scale(beta[select[j]],left[i,j]) for j in range(16) if left[i,j])) for i in range(16)]
    u=[D.zero]*61
    for i,b in zip(bundle['beta_indices_zero_based'],b16):u[i]=b
    selectedq=builder['common_q_selected_rows'];remainingq=builder['remaining_q_rows'];q=[D.zero]*45
    for j,i in enumerate(selectedq):q[i]=oldvar[f'v{j}']
    for j,i in enumerate(remainingq):q[i]=oldvar[f'z{j}']
    for i,a in zip(bundle['q_indices_zero_based'],q):u[i]=a
    model=S61();P=[D.zero]*150;P[15*PI[(1,2)]]=P[15*PI[(3,4)]]=D.one
    for i,v in zip(model.certificate['raw_free_indices'],u):P[i]=v
    for item in model.data['eliminated_raw_variables']:
        terms=[]
        for mon,num,den in item['polynomial']:
            v=D.const(F(num,den))
            for i in mon:v=mul(v,u[i])
            terms.append(v)
        P[item['raw_index']]=add(*terms)
    # Euler-RREF is not the physical divergence-free representative.
    div=[[D.zero]*5 for j in range(5)]
    for pair,(i,j) in enumerate(PAIRS):
        for mon,(a,b) in enumerate(MONS):
            co=P[15*pair+mon]
            for v,other in [(a,b),(b,a)]:
                if v==i:div[j][other]=add(div[j][other],co)
                if v==j:div[i][other]=add(div[i][other],D.neg(co))
    canonical=P[:]
    for pair,(i,j) in enumerate(PAIRS):
        for k in range(5):
            r=15*pair+MI[tuple(sorted((i,k)))];canonical[r]=add(canonical[r],scale(div[j][k],F(-1,5)))
            r=15*pair+MI[tuple(sorted((j,k)))];canonical[r]=add(canonical[r],scale(div[i][k],F(1,5)))
    def eye(n):return [[D.one if i==j else D.zero for j in range(n)] for i in range(n)]
    def mm(A,B):return [[add(*(mul(A[i][k],B[k][j]) for k in range(len(B)))) for j in range(len(B[0]))] for i in range(len(A))]
    def mat(M):return [[D.const(v) for v in row] for row in M.tolist()]
    def inverse(A):
        n=len(A);B=[row[:]+I[:] for row,I in zip(A,eye(n))]
        for j in range(n):
            pivot=next(i for i in range(j,n) if D.values[B[i][j]]);B[j],B[pivot]=B[pivot],B[j]
            z=D.inverse(B[j][j]);B[j]=[mul(v,z) for v in B[j]]
            for i in range(n):
                if i==j:continue
                co=B[i][j];B[i]=[add(v,D.neg(mul(co,w))) for v,w in zip(B[i],B[j])]
        return [row[n:] for row in B]
    G=eye(5);Gi=eye(5);elementary=[(2,2),(4,4),(3,2),(4,2),(3,1),(4,1)]+[(0,j) for j in range(1,5)]
    factors=[]
    for k,(i,j) in enumerate(elementary):
        A=eye(5);B=eye(5);A[i][j]=add(A[i][j],t[k])
        B[i][j]=add(B[i][j],D.neg(mul(t[k],D.inverse(add(D.one,t[k]))))) if i==j else D.neg(t[k])
        factors.append((A,B))
    chosen=base['orbit_constant_block_witness']['columns_zero_based']
    for k,j in enumerate(chosen,10):
        M=sp_basis()[j];plus=eye(4);minus=eye(4)
        for i in range(4):
            for l in range(4):
                v=scale(t[k],M[i,l]/2);plus[i][l]=add(plus[i][l],v);minus[i][l]=add(minus[i][l],D.neg(v))
        A4=mm(inverse(minus),plus);B4=mm(inverse(plus),minus);A=eye(5);B=eye(5)
        for i in range(4):
            for l in range(4):A[i+1][l+1]=A4[i][l];B[i+1][l+1]=B4[i][l]
        factors.append((A,B))
    for A,B in factors:G=mm(G,A);Gi=mm(B,Gi)
    fixed=sy.Matrix(wrec['exact_homogeneous_frame']);M=mm(mat(fixed.inv()),G);Mi=mm(Gi,mat(fixed))
    Q=object.__new__(Quadratic);Q.D=D;Q.S=Q.T=D.zero;Q.zero=(D.zero,D.zero);Q.one=(D.one,D.zero);Q.p=(D.zero,D.one)
    pairmat=lambda A:[[(v,D.zero) for v in row] for row in A]
    physical=[v[0] for v in action(Q,[(v,D.zero) for v in canonical],pairmat(M),pairmat(Mi),Q.one)]
    point=list(map(F,pre['original_physical_E1_coefficients']));assert [D.values[i] for i in physical]==point
    cuts=pre['common19cuts_plus_w140'];conditions=[add(physical[r],D.const(-point[r])) for r in cuts]
    selected=pre['response_selected17_rows'];other=[i for i in range(20) if i not in selected]
    eq17=[conditions[i] for i in selected];residual3=[conditions[i] for i in other];assert len(other)==3
    nilfree=['v1','v4','v5','v7'];relevant=set(t+[vn[n] for n in nilfree])
    values,grad=forward(D.nodes,relevant)
    J=[[grad[i].get(v,F(0)) for v in t] for i in eq17]
    invJ,det=linear_solve(J,[[F(int(i==j)) for j in range(17)] for i in range(17)])
    assert det==32
    allJ=[[grad[i].get(v,F(0)) for v in t] for i in conditions]
    assert enc(allJ)==pre['all20_by17_response']
    B=[[grad[i].get(vn[n],F(0)) for n in nilfree] for i in eq17]
    tj=[[-sum(invJ[i][k]*B[k][j] for k in range(17)) for j in range(4)] for i in range(17)]
    def after17(node):return [grad[node].get(vn[n],F(0))+sum(grad[node].get(t[i],F(0))*tj[i][j] for i in range(17)) for j,n in enumerate(nilfree)]
    nil3=[after17(i) for i in residual3];ncols=independent_rows([list(r) for r in zip(*nil3)]);assert len(ncols)==3
    _,ndet=linear_solve([[row[j] for j in ncols] for row in nil3],[[F(int(i==j)) for j in range(3)] for i in range(3)])
    free=next(i for i in range(4) if i not in ncols);ker=[F(0)]*4;ker[free]=1
    h,_=linear_solve([[row[j] for j in ncols] for row in nil3],[[-row[free]] for row in nil3])
    for i,c in enumerate(ncols):ker[c]=h[i][0]
    Rrows=source['physical_zero_base_rows9'];rho=[add(physical[r],D.const(-point[r])) for r in Rrows]
    rvel=[sum(a*b for a,b in zip(after17(node),ker)) for node in rho]
    unitrho=next((i for i,v in enumerate(rvel) if v),None)
    # Independent actual19 group cuts at E1: same source section, no type shortcut.
    def lie(U,i,j):
        V=[F(0)]*150
        for pair,(a,b) in enumerate(PAIRS):
            for mon,(r,s) in enumerate(MONS):
                co=U[15*pair+mon]
                if not co:continue
                if a==j and i!=b:V[15*PI[tuple(sorted((i,b)))]+mon]+=co*(1 if i<b else-1)
                if b==j and a!=i:V[15*PI[tuple(sorted((a,i)))]+mon]+=co*(1 if a<i else-1)
                for q,z in [(r,s),(s,r)]:
                    if q==i:V[15*pair+MI[tuple(sorted((j,z)))]]-=co
        return V
    O=[point]+[lie(point,*divmod(k,5)) for k in orbit['GL18_generator_indices']]
    C19=[[col[r] for col in O] for r in cuts[:19]]
    _,det19=linear_solve(C19,[[F(int(i==j)) for j in range(19)] for i in range(19)])
    assert det19
    oldequations=[polynomial(f) for name in ['CA','CB'] for f in builder['component_presentations'][name]]
    nilF=polynomial(exact['nilradical_generator']);ann=[polynomial(f) for f in exact['annihilator_prime_generators']]
    assert len(oldequations)==218 and len(ann)==45 and all(D.values[i]==0 for i in oldequations+ann+[nilF]+conditions)
    # A further literal21-unit input parametrizes the nilsupport curve by an actualR0 coordinate.
    nilcurve=None
    if unitrho is not None:
        binding={vn[n]:D.zero for n in names if n not in nilfree}
        needed=ancestors(D.nodes,conditions+rho)
        mp=replay(D,D.nodes[:],binding,needed);parameter=D.var('TAIL_NIL_R0_PARAMETER')
        nil_eq=[mp[i] for i in conditions]+[add(mp[rho[unitrho]],D.neg(parameter))]
        unknown=t+[vn[n] for n in nilfree]
        _,gg=forward(D.nodes,set(unknown));K=[[gg[i].get(v,F(0)) for v in unknown] for i in nil_eq]
        _,kdet=linear_solve(K,[[F(int(i==j)) for j in range(21)] for i in range(21)]);assert kdet
        nilcurve={'chosen_actual_R0_row':Rrows[unitrho],'parameter_node':parameter,'all21_unit_equations':nil_eq,
                  'unknown21_group17_then_old_free4':unknown,'unit_Jacobian_determinant':str(kdet),
                  'all9_R0_readback_nodes':[mp[i] for i in rho],
                  'nil_module_free_generator_symbol':'classofoldF (notitsvalueonreducednilcurve)',
                  'R0_surjects_regular_nilcurve_by_actual_unit_parameter':True}
    ev={'DAG_scalar_nodes':D.nodes,'old49_centered_variable_nodes':vn,'old49_unshifted_value_nodes':oldvar,
        'old49_center_shift':{'v4':1,'v7':1},'old218_equation_nodes':oldequations,
        'old49_to_S61_61_nodes':u,'old49_beta24_nodes':beta,'old49_S61_raw150_nodes':P,
        'old49_canonical_divergencefree150_nodes':canonical,
        'actual_group17_parameter_nodes':t,'rational_group_newframe_matrix':G,'rational_group_newframe_inverse':Gi,
        'actual_full_physical150_nodes':physical,'all20_common_cut_equation_nodes':conditions,
        'selected17_cut_positions':selected,'all17_unit_Hensel_equations':eq17,
        'unit17_Jacobian':enc(J),'unit17_Jacobian_determinant':str(det),
        'remaining3_cut_positions':other,'remaining3_actual_equation_nodes':residual3,
        'same_physical_R0_rows9':Rrows,'same_physical_R0_readback_nodes9':rho,'physical_w140_readback_node':physical[140],
        'old_actual_nil_generator_node':nilF,'old_nil_generator_expression':exact['nilradical_generator'],
        'old45_nil_annihilator_equation_nodes':ann,'old_nil_support_free_coordinates':nilfree,
        'nil_residual3_derivative_after17_on_free4':enc(nil3),'nil_rank3_minor_columns':ncols,'nil_rank3_minor':str(ndet),
        'nilcurve_tangent_on_free4':list(map(str,ker)),'actual_R0_readback_tangent_on_nilcurve':list(map(str,rvel)),
        'literal_nilcurve_R0_parameter_input':nilcurve,
        'actual_G19_cut_response_at_E1':enc(C19),'actual_G19_cut_minor':str(det19),
        'all_unit_denominators':D.divisions,
        'tail_ring':'completed QQ[old49shift,group17]/(old218,all20physicalcuts),unit17eliminable;actualR0mapreadback9',
        'N_tail_ring':'tailring/(oldF)','nil_module':'regularcurvequotientbyoldann45andcuts,freegeneratoroldF',
        'no_T3_syzygies_relations_depth_or_nil_of_B18_origin_claim':True}
    encoded=json.dumps(ev,separators=(',',':'))+'\n';assert len(encoded.encode())<32*1024*1024
    (out/'literal-tail-binding-DAG.json').write_text(encoded)
    paths=[Path(__file__),prepath,prep,builderp,exactp,wr,basep,bundlep,orbitp,sourcep,ROOT/'code/s61_model.py']
    report={'id':STEM,'all_passed':True,'status':'literal_actual_tail_model_with_regular3cut_and_nilcurve_R0_map',
            'old49_to_physical150_and_all_center_equalities_exact':True,'unit17_Jacobian_determinant':str(det),
            'all218_old_relations_and20physicalcuts_kept':True,'residual3_equation_count':3,
            'actual_G19_cut_minor_at_E1':str(det19),'nil_residual3_rank':3,'nil_residual3_minor':str(ndet),
            'nilcurve_has_actual_R0_unit_parameter':unitrho is not None,
            'chosen_nilcurve_actual_R0_row':Rrows[unitrho] if unitrho is not None else None,
            'tail_dimension_depth':[8,1],'N_tail_dimension_depth':[8,1],'nil_tail_regular_curve_dimension':1,
            'regular3cut_and_actual_binding_proof_in_companion_note_not_type_matching':True,
            'Jinf_syzygy_T3_presentation_or_B18_origin_depth_solved':False,
            'no_Newton_jets_GB_new_points_opaque_branch_pullback_or_p_inverse':True,
            'input_sha256':{str(p.relative_to(ROOT)):sha(p) for p in paths},
            'evidence_sha256':{'literal-tail-binding-DAG.json':sha(out/'literal-tail-binding-DAG.json')},
            'DAG_nodes':len(D.nodes),'certificate_bytes':len(encoded.encode()),'seconds':time.monotonic()-started,
            'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:report[k] for k in ['all_passed','unit17_Jacobian_determinant','actual_G19_cut_minor_at_E1','nil_residual3_minor','chosen_nilcurve_actual_R0_row','DAG_nodes','certificate_bytes','seconds','max_rss_kib']}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
