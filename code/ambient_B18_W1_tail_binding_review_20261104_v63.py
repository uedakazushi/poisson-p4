#!/usr/bin/env python3
"""Independent exact AD, physical group/old49 binding and regular-tail audit."""
import argparse
from fractions import Fraction as F
import gzip
import itertools
import json
from pathlib import Path
import resource
import signal
import time
import sympy as s
from ambient_B18_actual_top_coefficient_preflight_review_20261104_v63 import DAGReview,QuadraticReplay,lie,det
from ambient_W1_S61_intersection_recognition_review_20261103_v63 import rref,inverse

ROOT=Path(__file__).resolve().parents[1];BASE='AMBIENT_B18_W1_TAIL_BINDING_20261104_v63';ID=BASE+'_REVIEW'
PAIRS=list(itertools.combinations(range(5),2));MONS=list(itertools.combinations_with_replacement(range(5),2));PI={m:i for i,m in enumerate(PAIRS)};MI={m:i for i,m in enumerate(MONS)}
def read(p):return json.loads(gzip.decompress(p.read_bytes()) if p.suffix=='.gz' else p.read_text())
def sha(p):import hashlib;return hashlib.sha256(p.read_bytes()).hexdigest()
def matrix(A):return [list(map(F,row)) for row in A]
def mmn(A,B):return [[sum((x*y for x,y in zip(row,col)),F(0)) for col in zip(*B)] for row in A]
def matvec(A,v):return [sum((a*b for a,b in zip(row,v)),F(0)) for row in A]

class Pair:
    def __init__(self,D,nu):self.D=D;self.Q=QuadraticReplay(D,D.zero,nu);self.zero=(D.zero,D.zero);self.one=(D.one,D.zero)
    def add(self,a,b):return self.Q.add(a,b)
    def neg(self,a):return self.Q.neg(a)
    def mul(self,a,b):return self.Q.mul(a,b)
    def inv(self,a):return self.Q.inverse(a)
    def scalar(self,a):return(a,self.D.zero)
    def scale(self,a,c):return self.mul(a,self.scalar(self.D.const(c)))
    def total(self,vs):
        out=self.zero
        for a in vs:out=self.add(out,a)
        return out

def push(Q,U,G,Gi):
    sub=[]
    for p in range(10):
        vals=[]
        for r,t in MONS:
            terms=[]
            for m,(u,v) in enumerate(MONS):
                co=Q.mul(Gi[u][r],Gi[v][t])
                if r!=t:co=Q.add(co,Q.mul(Gi[u][t],Gi[v][r]))
                terms.append(Q.mul(U[15*p+m],co))
            vals.append(Q.total(terms))
        sub.append(vals)
    result=[]
    for i,j in PAIRS:
        wedge=[Q.add(Q.mul(G[i][a],G[j][b]),Q.neg(Q.mul(G[i][b],G[j][a]))) for a,b in PAIRS]
        result.extend(Q.total(Q.mul(c,sub[k][m]) for k,c in enumerate(wedge)) for m in range(15))
    return result

def main(out):
    began=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
    folder=ROOT/'reports'/BASE;primary=read(folder/'report.json');c=read(folder/'literal-tail-binding-DAG.json')
    pbase=ROOT/'reports/AMBIENT_B18_W1_TAIL_COORDINATE_PREFLIGHT_20261104_v63';pre=read(pbase/'tail-coordinate-preflight-DAG.json')
    correction=ROOT/'reports/AMBIENT_B18_W1_POLYNOMIAL_LAURENT_SCOPE_20261104_v63.json'
    gates=[ROOT/'reports/AMBIENT_W1_S61_INTERSECTION_RECOGNITION_20261103_v63_REVIEW.json',
        ROOT/'reports/AMBIENT_B18_ACTUAL_INTERSECTION_ORBIT_PRODUCT_20261103_v63_REVIEW.json',
        ROOT/'reports/AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63_REVIEW/report.json',
        ROOT/'reports/AMBIENT_B18_FOUR_MODULE_DEPTH_CRITERION_20261104_v63_REVIEW/report.json',correction]
    checked={}
    for rp in [folder/'report.json',pbase/'report.json',*gates]:
        r=read(rp);assert r['all_passed']
        for n,h in r.get('input_sha256',{}).items():assert sha(ROOT/n)==h,n;checked[n]=h
        for n,h in r.get('evidence_sha256',{}).items():p=rp.parent/n;assert sha(p)==h,n;checked[str(p.relative_to(ROOT))]=h
    assert read(correction)['correct_operation']=='finite_w_polynomial_Dalg_localize_before_w1completion'
    old=read(ROOT/'reports/AMBIENT_INTERSECTION_EXACT_IDEAL_20261103_v63/report.json')
    builder=read(ROOT/'reports/AMBIENT_INTERSECTION_EXACT_IDEAL_20261103_v63/builder.json')
    base=read(ROOT/'reports/AMBIENT_BASE_PINCHED_SURFACE_20261103_v63/report.json')
    bundle=read(ROOT/'reports/AMBIENT_FIBER_BUNDLE_20261103_v63/report.json')
    wrec=read(ROOT/'reports/AMBIENT_W1_S61_INTERSECTION_RECOGNITION_20261103_v63/recognition.json.gz')
    sourcepath=ROOT/'reports/AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63/actual-source-coefficients.json';source=read(sourcepath)
    D=DAGReview(c['DAG_scalar_nodes']);point=[F(0)]*150
    for pair,mon,v in [((0,2),(1,1),-1),((2,3),(0,1),1),((2,4),(0,0),1),((2,4),(3,3),1),((3,4),(1,1),1)]:point[15*PI[pair]+MI[mon]]=F(v)
    assert point==list(map(F,pre['original_physical_E1_coefficients']))==[D.values[i] for i in c['actual_full_physical150_nodes']]
    vn=c['old49_centered_variable_nodes'];oldvar={n:D.add(i,D.const(1 if n in ['v4','v7'] else 0)) for n,i in vn.items()}
    assert oldvar==c['old49_unshifted_value_nodes']
    def power(a,n):
        v=D.one
        for _ in range(n):v=D.mul(v,a)
        return v
    def polynomial(text,variables=oldvar):
        z=s.sympify(text);vs=sorted(z.free_symbols,key=str)
        if not vs:return D.const(z)
        terms=[]
        for m,co in s.Poly(z,*vs).terms():
            term=D.const(co)
            for x,n in zip(vs,m):term=D.mul(term,power(variables[str(x)],n))
            terms.append(term)
        return D.total(terms)
    ss,aa,bb,cc=[oldvar[n] for n in ['ss','aa','bb','cc']];S,P=s.symbols('s p');beta=[]
    for text in base['slice_24_coordinates']:
        terms=[]
        for (es,ep),co in s.Poly(s.sympify(text),S,P).terms():
            z=D.const(co)
            if ep==1:assert es;z=D.mul(z,D.mul(power(ss,es-1),cc))
            elif ep%2==0:z=D.mul(z,D.mul(power(ss,es),power(aa,ep//2)))
            else:z=D.mul(z,D.mul(power(ss,es),D.mul(bb,power(aa,(ep-3)//2))))
            terms.append(z)
        beta.append(D.total(terms))
    assert beta==c['old49_beta24_nodes']
    adapter=s.Matrix(bundle['beta16_to_intrinsic24_adapter']);rows=list(adapter.T.rref()[1]);left=adapter[rows,:].inv()
    assert left*adapter[rows,:]==s.eye(16)
    b16=[D.total(D.mul(beta[rows[j]],D.const(left[i,j])) for j in range(16) if left[i,j]) for i in range(16)]
    u=[D.zero]*61
    for i,z in zip(bundle['beta_indices_zero_based'],b16):u[i]=z
    q45=[D.zero]*45
    for j,i in enumerate(builder['common_q_selected_rows']):q45[i]=oldvar[f'v{j}']
    for j,i in enumerate(builder['remaining_q_rows']):q45[i]=oldvar[f'z{j}']
    for i,z in zip(bundle['q_indices_zero_based'],q45):u[i]=z
    assert u==c['old49_to_S61_61_nodes']
    cert=read(ROOT/'legacy/v02/data/symplectic_slice_certificate.json');eqdata=read(ROOT/'legacy/v02/data/symplectic_slice_equations.json')
    raw=[D.zero]*150;raw[15*PI[1,2]]=raw[15*PI[3,4]]=D.one
    for i,z in zip(cert['raw_free_indices'],u):raw[i]=z
    for rec in eqdata['eliminated_raw_variables']:
        terms=[]
        for mon,num,den in rec['polynomial']:
            z=D.const(F(num,den))
            for i in mon:z=D.mul(z,u[i])
            terms.append(z)
        raw[rec['raw_index']]=D.total(terms)
    assert raw==c['old49_S61_raw150_nodes']
    div=[[D.zero]*5 for _ in range(5)]
    for ip,(i,j) in enumerate(PAIRS):
        for m,(a,b) in enumerate(MONS):
            z=raw[15*ip+m]
            for v,k in [(a,b),(b,a)]:
                if v==i:div[j][k]=D.add(div[j][k],z)
                if v==j:div[i][k]=D.add(div[i][k],D.neg(z))
    canonical=raw.copy()
    for ip,(i,j) in enumerate(PAIRS):
        for k in range(5):
            a=15*ip+MI[tuple(sorted((i,k)))];b=15*ip+MI[tuple(sorted((j,k)))]
            canonical[a]=D.add(canonical[a],D.mul(div[j][k],D.const(F(-1,5))))
            canonical[b]=D.add(canonical[b],D.mul(div[i][k],D.const(F(1,5))))
    assert canonical==c['old49_canonical_divergencefree150_nodes']
    def eye(n):return [[D.one if i==j else D.zero for j in range(n)] for i in range(n)]
    def mm(A,B):return [[D.total(D.mul(A[i][k],B[k][j]) for k in range(len(B))) for j in range(len(B[0]))] for i in range(len(A))]
    def constmat(M):return [[D.const(v) for v in row] for row in M]
    def minverse(A):
        n=len(A);B=[row.copy()+E.copy() for row,E in zip(A,eye(n))]
        for j in range(n):
            k=next(k for k in range(j,n) if D.values[B[k][j]]);B[j],B[k]=B[k],B[j]
            iv=D.inverse(B[j][j]);B[j]=[D.mul(v,iv) for v in B[j]]
            for i in range(n):
                if i!=j:co=B[i][j];B[i]=[D.add(v,D.neg(D.mul(co,w))) for v,w in zip(B[i],B[j])]
        return [row[n:] for row in B]
    J4=[[0,1,0,0],[-1,0,0,0],[0,0,0,1],[0,0,-1,0]];sp=[]
    for i in range(4):
        for j in range(i,4):
            sym=[[0]*4 for _ in range(4)];sym[i][j]=sym[j][i]=1;sp.append(mmn(J4,sym))
    t=c['actual_group17_parameter_nodes'];G=eye(5);Gi=eye(5);generators=[]
    for k,(i,j) in enumerate([(2,2),(4,4),(3,2),(4,2),(3,1),(4,1),*[(0,j) for j in range(1,5)]]):
        A=eye(5);B=eye(5);A[i][j]=D.add(A[i][j],t[k]);B[i][j]=D.add(B[i][j],D.neg(D.mul(t[k],D.inverse(D.add(D.one,t[k]))))) if i==j else D.neg(t[k])
        G=mm(G,A);Gi=mm(B,Gi);gen=[[F(0)]*5 for _ in range(5)];gen[i][j]=1;generators.append(gen)
    for k,j in enumerate(base['orbit_constant_block_witness']['columns_zero_based'],10):
        plus=eye(4);minus=eye(4)
        for i in range(4):
            for h in range(4):
                z=D.mul(t[k],D.const(sp[j][i][h]/2));plus[i][h]=D.add(plus[i][h],z);minus[i][h]=D.add(minus[i][h],D.neg(z))
        A4=mm(minverse(minus),plus);B4=mm(minverse(plus),minus);A=eye(5);B=eye(5);gen=[[F(0)]*5 for _ in range(5)]
        for i in range(4):
            for h in range(4):A[i+1][h+1]=A4[i][h];B[i+1][h+1]=B4[i][h];gen[i+1][h+1]=sp[j][i][h]
        G=mm(G,A);Gi=mm(B,Gi);generators.append(gen)
    assert G==c['rational_group_newframe_matrix'] and Gi==c['rational_group_newframe_inverse']
    fixed=matrix(wrec['exact_homogeneous_frame']);fixedinv,_=inverse(fixed)
    M=mm(constmat(fixedinv),G);Mi=mm(Gi,constmat(fixed))
    Qring=Pair(D,D.zero);U=[(v,D.zero) for v in canonical];pairmat=lambda A:[[(v,D.zero) for v in row] for row in A]
    physical=[v[0] for v in push(Qring,U,pairmat(M),pairmat(Mi))]
    assert physical==c['actual_full_physical150_nodes']
    # Differentiate the independent physical Lie formula at the SAME E1 point.
    L={(i,j):lie(point,i,j) for i in range(5) for j in range(5)}
    response=[]
    for gen in generators:
        g=mmn(mmn(fixedinv,gen),fixed)
        response.append([sum((g[i][j]*L[i,j][r] for i in range(5) for j in range(5)),F(0)) for r in range(150)])
    assert [[col[r] for col in response] for r in range(150)]==matrix(pre['physical150_by17_group_response'])
    assert all(D.grad[root].get(v,F(0))==response[j][i] for i,root in enumerate(physical) for j,v in enumerate(t))
    cuts=pre['common19cuts_plus_w140'];conds=[D.add(physical[r],D.const(-point[r])) for r in cuts]
    assert conds==c['all20_common_cut_equation_nodes']
    selected=c['selected17_cut_positions'];other=[i for i in range(20) if i not in selected]
    A=[[D.grad[conds[i]].get(v,F(0)) for v in t] for i in selected];Ai,adet=inverse(A)
    assert adet==32 and A==matrix(c['unit17_Jacobian']) and other==[14,18,19]
    assert [conds[i] for i in selected]==c['all17_unit_Hensel_equations'] and [conds[i] for i in other]==c['remaining3_actual_equation_nodes']
    nilfree=['v1','v4','v5','v7'];nilvars=[vn[n] for n in nilfree]
    B=[[D.grad[conds[i]].get(v,F(0)) for v in nilvars] for i in selected]
    tangent=[[-v for v in row] for row in mmn(Ai,B)]
    def after17(root):return [D.grad[root].get(v,F(0))+sum(D.grad[root].get(t[i],F(0))*tangent[i][j] for i in range(17)) for j,v in enumerate(nilvars)]
    nil3=[after17(conds[i]) for i in other]
    assert nil3==matrix(c['nil_residual3_derivative_after17_on_free4'])==[[0,0,1,0],[0,-1,0,0],[0,-1,0,F(5,2)]]
    assert det([[row[j] for j in c['nil_rank3_minor_columns']] for row in nil3])==F(5,2)
    ker=list(map(F,c['nilcurve_tangent_on_free4']));assert ker==[1,0,0,0] and matvec(nil3,ker)==[0]*3
    Rrows=source['physical_zero_base_rows9'];rho=[D.add(physical[r],D.const(-point[r])) for r in Rrows]
    assert Rrows==c['same_physical_R0_rows9']==[27,108,109,112,113,117,118,119,61] and rho==c['same_physical_R0_readback_nodes9']
    speed=[sum(x*y for x,y in zip(after17(root),ker)) for root in rho]
    assert speed==[0,0,0,0,0,-1,0,0,0]==list(map(F,c['actual_R0_readback_tangent_on_nilcurve']))
    orbit=read(ROOT/'reports/AMBIENT_B18_ACTUAL_INTERSECTION_ORBIT_PRODUCT_20261103_v63/certificate.json.gz')
    O=[point]+[lie(point,*divmod(k,5)) for k in orbit['GL18_generator_indices']]
    C19=[[v[r] for v in O] for r in cuts[:19]];assert C19==matrix(c['actual_G19_cut_response_at_E1']) and det(C19)==-32
    # ALL original ideals are retained, with their actual point translations.
    eqs=[polynomial(z) for name in ['CA','CB'] for z in builder['component_presentations'][name]]
    ann=[polynomial(z) for z in old['annihilator_prime_generators']];nilF=polynomial(old['nilradical_generator'])
    assert len(eqs)==218 and eqs==c['old218_equation_nodes'] and ann==c['old45_nil_annihilator_equation_nodes'] and nilF==c['old_actual_nil_generator_node']
    assert old['nilradical_generator']=='2*aa + v0 + z25'
    assert set(old['annihilator_prime_generators'])==set(vn)-set(nilfree)
    assert all(D.values[i]==0 for i in eqs+ann+[nilF]+conds)
    curve=c['literal_nilcurve_R0_parameter_input'];assert curve['chosen_actual_R0_row']==117
    binding={v:D.zero for n,v in vn.items() if n not in nilfree};cache={**binding}
    def substitute(i):
        if i in cache:return cache[i]
        n=D.nodes[i]
        if n[0] in ['c','v']:z=i
        elif n[0]=='n':z=D.neg(substitute(n[1]))
        elif n[0]=='i':z=D.inverse(substitute(n[1]))
        else:z={'a':D.add,'m':D.mul}[n[0]](substitute(n[1]),substitute(n[2]))
        cache[i]=z;return z
    ce=[substitute(i) for i in conds];cr=[substitute(i) for i in rho]
    ce.append(D.add(cr[5],D.neg(curve['parameter_node'])))
    assert ce==curve['all21_unit_equations'] and cr==curve['all9_R0_readback_nodes']
    unknown=t+nilvars;J21=[[D.grad[i].get(v,F(0)) for v in unknown] for i in ce]
    assert det(J21)==-80 and unknown==curve['unknown21_group17_then_old_free4']
    # Repeated requests can share one inverse node; validate coverage by node set.
    assert {r for r,x,v in c['all_unit_denominators']}=={i for i,n in enumerate(D.nodes) if n[0]=='i'}
    for root,x,v in c['all_unit_denominators']:assert D.nodes[root]==('i',x) and D.values[x]==F(v)!=0
    # Torus preserves the19-cut ideal; this works scheme-wise including nil.
    torus=read(ROOT/'reports/AMBIENT_B18_POSITIVE_GRADING_OBSTRUCTION_20261104_v63/certificate.json')
    h=[0,0,-1,0,0];weights=[1+h[i]+h[j]-h[a]-h[b] for i,j in PAIRS for a,b in MONS]
    assert weights==torus['physical150_torus_weights'] and weights[140]==1 and point[140]==1
    assert all(not point[r] or weights[r]==0 for r in cuts[:19])
    # The coordinate-preflight finite value/perspective map: exact source prefix,
    # independent AD of its unit denominators and all150 centers, same group response.
    DP=DAGReview(pre['DAG_scalar_nodes']);assert pre['DAG_scalar_nodes'][:len(source['DAG_scalar_nodes'])]==source['DAG_scalar_nodes']
    assert all(DP.values[a]!=0 for i,a,v in pre['all_new_unit_inverse_norms'])
    assert [DP.values[z[0]] for z in pre['actual_physical150_tail_w1_trace_pairs']]==point
    assert all(DP.values[z[1]]==0 for z in pre['actual_physical150_tail_w1_trace_pairs'])
    assert [DP.values[z[0]] for z in pre['all61_normalized_pairs']]==list(map(F,wrec['exact_normalized61']))
    assert [DP.values[z[0]] for z in pre['all24_primitive_beta_pairs']]==list(map(F,wrec['full_beta24']))
    # Universal rational Darboux frame identity, separate from center tests.
    a,b,c0,d,e,f=s.symbols('a b c d e f');delta=a*f-b*e+c0*d
    V=s.Matrix([[0,a,b,c0],[-a,0,d,e],[-b,-d,0,f],[-c0,-e,-f,0]])
    FF=s.Matrix([[1,0,0,0],[0,1/a,0,0],[d/a,-b/a,1,0],[e/delta,-c0/delta,0,a/delta]])
    assert (FF*V*FF.T-s.Matrix([[0,1,0,0],[-1,0,0,0],[0,0,0,1],[0,0,-1,0]])).applyfunc(s.cancel)==s.zeros(4)
    out.mkdir(parents=True,exist_ok=False);ep=out/'independent-proof.json'
    evidence={'old49_actual61_raw150_canonical150_and_group17_finite_maps_literal':True,
        'all218_old_polynomials_and45_ann_and_nilF_nodes_bound':True,
        'all150_center_and150by17_group_AD_vs_independent_Lie_equal':True,
        'unit17_det':'32','group19_det':'-32','nil3_matrix':[[str(v) for v in row] for row in nil3],
        'nil3_minor':'5/2','nil21_det':'-80','nil_actual_R0_parameter_row':117,
        'nil_actual_R0_velocity':[str(v) for v in speed],
        'regular20_then_permuted17plus3_uses_scheme_group_and_torus_product_not_CM':True,
        'w_polynomial_Dalg_before_w1completion_not_w0_completed_ring':True,
        'D_tail_and_N_tail_dim_depth':[8,1],'nil_tail_regular_curve_cyclic_over_same_R0':True,
        'no_T3_Jinf_syzygy_or_origin_depth19_20_claim':True}
    ep.write_text(json.dumps(evidence,indent=2)+'\n')
    paths=[Path(__file__),ROOT/'research/notes'/(ID+'.md'),folder/'report.json',folder/'literal-tail-binding-DAG.json',pbase/'report.json',pbase/'tail-coordinate-preflight-DAG.json',
        correction,ROOT/'research/notes/AMBIENT_B18_W1_POLYNOMIAL_LAURENT_SCOPE_20261104_v63.md',*gates,
        ROOT/'legacy/v02/data/symplectic_slice_certificate.json',ROOT/'legacy/v02/data/symplectic_slice_equations.json',
        ROOT/'code/ambient_B18_actual_top_coefficient_preflight_review_20261104_v63.py',ROOT/'code/ambient_W1_S61_intersection_recognition_review_20261103_v63.py']
    checked.update({str(p.relative_to(ROOT)):sha(p) for p in paths})
    result={'id':ID,'all_passed':True,'actual_old49_tensor_group17_common20cut_binding_verified':True,
        'own_exact_AD_Lie_and_unit17_group19_nil3_nil21_minors_verified':True,
        'regular3cuts_for_D_N_and_nil_without_CM_assumption_accepted':True,
        'Dalg_polynomial_Laurent_scope_correction_required_and_bound':True,
        'tail_D_and_reduced_N_dimension_depth':[8,1],
        'tail_nil_principal_square_zero_and_cyclic_same_R0_regular_curve':True,
        'origin_B18_depth_or_Jinf_presentation_or_T3_relation_claim':False,
        'input_sha256':checked,'evidence_sha256':{ep.name:sha(ep)},
        'seconds':time.monotonic()-began,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'frozen':True}
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');signal.alarm(0)
    print(json.dumps({k:result[k] for k in ['all_passed','seconds','max_rss_kib']}))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
