#!/usr/bin/env python3
"""ONE bounded actual rank2 coefficient DAG; do not invent a relative graph lift."""
import argparse
from fractions import Fraction as F
import gzip
import hashlib
import itertools as it
import json
from pathlib import Path
import resource
import signal
import time
import sympy as s

ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63'
PAIRS=list(it.combinations(range(5),2));MONS=list(it.combinations_with_replacement(range(5),2))
PI={p:i for i,p in enumerate(PAIRS)};MI={p:i for i,p in enumerate(MONS)}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(gzip.decompress(p.read_bytes())if p.suffix=='.gz'else p.read_text())

class DAG:
    def __init__(self):
        self.nodes=[];self.cache={};self.values=[];self.divisions=[]
        self.zero=self.const(0);self.one=self.const(1)
    def put(self,key,value):
        if key in self.cache:return self.cache[key]
        i=len(self.nodes);self.nodes.append(key);self.values.append(value);self.cache[key]=i;return i
    def const(self,q):
        q=F(str(q));return self.put(('c',q.numerator,q.denominator),q)
    def var(self,name):return self.put(('v',name),F(0))
    def neg(self,a):
        if a==self.zero:return a
        n=self.nodes[a]
        if n[0]=='c':return self.const(-self.values[a])
        if n[0]=='n':return n[1]
        return self.put(('n',a),-self.values[a])
    def add(self,a,b):
        if a==self.zero:return b
        if b==self.zero:return a
        if self.nodes[a]==('n',b) or self.nodes[b]==('n',a):return self.zero
        if self.nodes[a][0]==self.nodes[b][0]=='c':return self.const(self.values[a]+self.values[b])
        a,b=sorted((a,b));return self.put(('a',a,b),self.values[a]+self.values[b])
    def mul(self,a,b):
        if a==self.zero or b==self.zero:return self.zero
        if a==self.one:return b
        if b==self.one:return a
        if self.nodes[a][0]==self.nodes[b][0]=='c':return self.const(self.values[a]*self.values[b])
        a,b=sorted((a,b));return self.put(('m',a,b),self.values[a]*self.values[b])
    def inverse(self,a):
        assert self.values[a],('nonunit denominator at origin',a)
        if self.nodes[a][0]=='c':return self.const(1/self.values[a])
        q=self.put(('i',a),1/self.values[a]);self.divisions.append([q,a,str(self.values[a])]);return q
    def total(self,values):
        a=self.zero
        for b in values:a=self.add(a,b)
        return a

class Quadratic:
    def __init__(self,D,name):
        self.D=D;self.S=D.var(name+'_trace_p');self.T=D.var(name+'_minus_norm_p')
        self.zero=(D.zero,D.zero);self.one=(D.one,D.zero);self.p=(D.zero,D.one)
    def scalar(self,a):return(a,self.D.zero)
    def add(self,a,b):return tuple(self.D.add(a[i],b[i])for i in range(2))
    def neg(self,a):return tuple(self.D.neg(x)for x in a)
    def mul(self,a,b):
        D=self.D;bd=D.mul(a[1],b[1])
        return(D.add(D.mul(a[0],b[0]),D.mul(self.T,bd)),
               D.total([D.mul(a[0],b[1]),D.mul(a[1],b[0]),D.mul(self.S,bd)]))
    def scale(self,a,q):return self.mul(a,self.scalar(self.D.const(q)))
    def power(self,a,n):
        b=self.one
        for _ in range(n):b=self.mul(b,a)
        return b
    def inv(self,a):
        D=self.D
        norm=D.add(D.add(D.mul(a[0],a[0]),D.mul(self.S,D.mul(a[0],a[1]))),D.neg(D.mul(self.T,D.mul(a[1],a[1]))))
        ni=D.inverse(norm)
        return(D.mul(D.add(a[0],D.mul(self.S,a[1])),ni),D.neg(D.mul(a[1],ni)))
    def total(self,values):
        a=self.zero
        for b in values:a=self.add(a,b)
        return a
    def center(self,a):return tuple(self.D.values[i]for i in a)
    def tracepair(self,a):
        D=self.D
        return(D.add(a[0],D.mul(D.const(F(1,2)),D.mul(a[1],self.S))),a[1])

def determinant(Q,M):
    n=len(M);answer=Q.zero
    for perm in it.permutations(range(n)):
        term=Q.scale(Q.one,(-1)**sum(perm[i]>perm[j]for i in range(n)for j in range(i+1,n)))
        for i in range(n):term=Q.mul(term,M[i][perm[i]])
        answer=Q.add(answer,term)
    return answer

def inverse_block(Q,G):
    others=[0,1,3,4];G4=[[G[i][j]for j in others]for i in others]
    det=determinant(Q,G4);invdet=Q.inv(det)
    out=[[Q.zero for _ in range(5)]for _ in range(5)]
    for i in range(4):
        for j in range(4):
            minor=[[G4[r][c]for c in range(4)if c!=i]for r in range(4)if r!=j]
            out[others[i]][others[j]]=Q.scale(Q.mul(determinant(Q,minor),invdet),(-1)**(i+j))
    out[2][2]=Q.inv(G[2][2]);return out

def action(Q,U,G,Gi,scalar):
    transforms={}
    for uv in MONS:
        u,v=uv
        for rs in MONS:
            r,s=rs
            value=Q.mul(Gi[u][r],Gi[v][s])
            if r!=s:value=Q.add(value,Q.mul(Gi[u][s],Gi[v][r]))
            transforms[uv,rs]=value
    substituted=[]
    for pair in PAIRS:
        row=[]
        for rs in MONS:
            row.append(Q.total(Q.mul(U[15*PI[pair]+j],transforms[uv,rs])for j,uv in enumerate(MONS)))
        substituted.append(row)
    result=[]
    for i,j in PAIRS:
        wedge=[Q.add(Q.mul(G[i][a],G[j][b]),Q.neg(Q.mul(G[i][b],G[j][a])))for a,b in PAIRS]
        for k in range(15):
            result.append(Q.mul(scalar,Q.total(Q.mul(wedge[l],substituted[l][k])for l in range(10))))
    return result

def lie(Q,U,i,j):
    result=[Q.zero]*150
    def add(index,value):result[index]=Q.add(result[index],value)
    for pair,(a,b)in enumerate(PAIRS):
        for mon,(u,v)in enumerate(MONS):
            co=U[15*pair+mon]
            if co==Q.zero:continue
            if a==j and i!=b:add(15*PI[tuple(sorted((i,b)))]+mon,Q.scale(co,1 if i<b else-1))
            if b==j and a!=i:add(15*PI[tuple(sorted((a,i)))]+mon,Q.scale(co,1 if a<i else-1))
            for q,other in [(u,v),(v,u)]:
                if q==i:add(15*pair+MI[tuple(sorted((j,other)))],Q.neg(co))
    return result

def solve(Q,M,rhs):
    n=len(M);A=[list(row)+[rhs[i]]for i,row in enumerate(M)];trace=[]
    for j in range(n):
        p=next(i for i in range(j,n)if Q.center(A[i][j])[0])
        A[j],A[p]=A[p],A[j];pivot=A[j][j];trace.append({'column':j,'selected_row':p,'pivot':pivot,'origin':list(map(str,Q.center(pivot)))})
        inverse=Q.inv(pivot)
        A[j]=[Q.mul(x,inverse)for x in A[j]]
        for i in range(n):
            if i==j:continue
            co=A[i][j]
            A[i]=[Q.add(x,Q.neg(Q.mul(co,y)))for x,y in zip(A[i],A[j])]
    return[A[i][-1]for i in range(n)],trace

def main(out):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));out.mkdir(parents=True,exist_ok=False)
    paths=[ROOT/'reports/AMBIENT_B18_QUADRATIC_IMAGE_REDUCTION_20261103_v63/report.json',
           ROOT/'reports/AMBIENT_B18_CORE_GAUGE_PREFLIGHT_20261104_v63/certificate.json',
           ROOT/'reports/AMBIENT_B18_CORE_GAUGE_PREFLIGHT_20261104_v63/implicit-system.json',
           ROOT/'reports/AMBIENT_B18_RELATIVE_WEIGHTED_MODEL_20261104_v63/certificate.json',
           ROOT/'reports/AMBIENT_B18_ACTUAL_INTERSECTION_TANGENT_20261103_v63/certificate.json.gz',
           ROOT/'reports/AMBIENT_B18_ACTUAL_INTERSECTION_ORBIT_PRODUCT_20261103_v63/certificate.json.gz']
    raw,gauge,implicit,weight,tangent,orbit=map(read,paths)
    R9=weight['remaining_actual_even_base10_raw_coordinates'];R9=[i for i in R9 if i!=140]
    gzero=implicit['zero_GL_group_parameter_indices'];gpositive=[i for i in range(19)if i not in gzero]
    cuts=weight['common_physical_group_cuts19'];weights=read(ROOT/'reports/AMBIENT_B18_POSITIVE_GRADING_OBSTRUCTION_20261104_v63/certificate.json')['physical150_torus_weights']
    zerocuts=[i for i in cuts if weights[i]==0];poscuts=[i for i in cuts if weights[i]==1]
    assert len(gzero)==15 and len(gpositive)==4 and len(zerocuts)==15 and len(poscuts)==4
    assert all(orbit['GL18_generator_indices'][i-1]%5==2 for i in gpositive)
    odd=tangent['common10_odd_target_rows'];posodd=[i for i in odd if weights[i]==1];assert len(posodd)==2
    cubics=[r['raw']for r in weight['possible_cubic_discrepancy_raw_coordinates']]
    D=DAG();xids=[D.var(f'x{i}')for i in R9];records={}
    params=s.symbols('p k0:9');sourcezero=[k for k in range(9)if k!=2]
    for family in ['CA','CB']:
        Q=Quadratic(D,family)
        names=[f'gamma{i}'for i in gzero]+[f'k{i}'for i in sourcezero]
        unknown={name:(D.var(f'{family}_{name}_a'),D.var(f'{family}_{name}_b'))for name in names}
        G=[[Q.one if i==j else Q.zero for j in range(5)]for i in range(5)]
        for index in gzero:
            if index==0:continue
            k=orbit['GL18_generator_indices'][index-1];i,j=divmod(k,5)
            G[i][j]=Q.add(G[i][j],unknown[f'gamma{index}'])
        assert all(G[2][i]==G[i][2]==Q.zero for i in range(5)if i!=2)
        scalar=Q.add(Q.one,unknown['gamma0']);Gi=inverse_block(Q,G)
        powers={}
        def exp(var,n):
            key=(var,n)
            if key not in powers:powers[key]=Q.power(Q.p if var=='p'else unknown[var],n)
            return powers[key]
        U0=[];V0=[]
        for row,text in enumerate(raw['families'][family]['core_full150_coordinates']):
            poly=s.Poly(s.sympify(text),*params);parts=[Q.zero,Q.zero]
            for ex,co in poly.terms():
                assert ex[3]<=1
                value=Q.scalar(D.const(co))
                value=Q.mul(value,exp('p',ex[0]))
                for index in sourcezero:value=Q.mul(value,exp(f'k{index}',ex[index+1]))
                parts[ex[3]]=Q.add(parts[ex[3]],value)
                pair=PAIRS[row//15];mon=MONS[row%15]
                if co:
                    assert 2 not in mon
                    assert (2 in pair) == (ex[3]==0)
            U0.append(parts[0]);V0.append(parts[1])
        q0=action(Q,U0,G,Gi,scalar);v0=action(Q,V0,G,Gi,scalar)
        zeroequations=[]
        for row in zerocuts:
            f=Q.add(q0[row],Q.scalar(D.const(-F(orbit['physical_point150'][row]))));zeroequations.extend(f)
        for row,xid in zip(R9,xids):
            f=Q.add(q0[row],Q.scalar(D.neg(xid)));zeroequations.extend(f)
        assert len(zeroequations)==48 and all(D.values[i]==0 for i in zeroequations)
        # The rank2 coefficient Hensel Jacobian has the actual M24 twice on its diagonal.
        O=s.Matrix(orbit['physical150_by19_scalar_then_GL']);L=s.Matrix(raw['families'][family]['linear_J'])[:,19:]
        L=L[:,[i for i in range(9)if i!=2]];W=s.Matrix(raw['families'][family]['raw_second_coefficient'])
        rows=zerocuts+R9
        M24=O[rows,gzero].row_join(L[rows,:]).row_join(W[rows,:])
        det24=s.factor(M24.det(method='domain-ge'))
        expected=s.Rational(gauge['records'][family]['full23_implicit_Jacobian_determinant'])*s.Rational(gauge['records'][family]['fold_quadratic_coefficient_after23_unit_eliminations'])
        assert det24==expected and det24
        columns=[]
        for index in gpositive:
            k=orbit['GL18_generator_indices'][index-1];i,j=divmod(k,5);assert j==2
            columns.append([Q.mul(z,Q.inv(G[2][2]))for z in lie(Q,q0,i,2)])
        columns.append(v0)
        rows5=poscuts+[140];matrix5=[[columns[j][i]for j in range(5)]for i in rows5]
        solution,solve_trace=solve(Q,matrix5,[Q.zero]*4+[Q.one])
        shear=[Q.zero]*5
        for index,value in zip(gpositive,solution[:4]):
            k=orbit['GL18_generator_indices'][index-1];i,j=divmod(k,5)
            shear[i]=Q.mul(value,Q.inv(G[2][2]))
        k2=solution[4]
        # Exact rank-one shear. The source has ey wedge V plus k2*V2,
        # all coefficients independent of y, so degree3 is exact, not a truncation.
        vector=[[Q.zero]*15 for _ in range(5)]
        for i in range(5):
            if i==2:continue
            pair=tuple(sorted((2,i)));sgn=1 if 2<i else-1
            vector[i]=[Q.scale(q0[15*PI[pair]+j],sgn)for j in range(15)]
        def shift(coefficients):
            answer=[[Q.zero]*15 for _ in range(3)]
            for j,(u,v)in enumerate(MONS):
                co=coefficients[j]
                if co==Q.zero:continue
                assert u!=2 and v!=2
                answer[0][j]=Q.add(answer[0][j],co)
                for target,velocity in [(v,shear[u]),(u,shear[v])]:
                    k=MI[tuple(sorted((2,target)))];answer[1][k]=Q.add(answer[1][k],Q.neg(Q.mul(co,velocity)))
                answer[2][MI[2,2]]=Q.add(answer[2][MI[2,2]],Q.mul(co,Q.mul(shear[u],shear[v])))
            return answer
        shifted_vector=[shift(v)for v in vector]
        full=[[Q.zero]*150 for _ in range(4)]
        for i,j in PAIRS:
            for mon in range(15):
                row=15*PI[i,j]+mon
                for n in range(3):
                    if i==2:full[n][row]=Q.add(full[n][row],shifted_vector[j][n][mon])
                    if j==2:full[n][row]=Q.add(full[n][row],Q.neg(shifted_vector[i][n][mon]))
                    value=Q.add(Q.mul(shear[i],shifted_vector[j][n][mon]),Q.neg(Q.mul(shear[j],shifted_vector[i][n][mon])))
                    full[n+1][row]=Q.add(full[n+1][row],value)
                vv=shift(v0[15*PI[i,j]:15*PI[i,j]+15])
                for n in range(3):full[n+1][row]=Q.add(full[n+1][row],Q.mul(k2,vv[n][mon]))
        # Exact source-side coefficient data; no CB-to-CA cover identification.
        records[family]={
            'zero_hensel_unknowns23_order':names,'coefficient_unknown_pairs':unknown,
            'p_quadratic_relation':{'s':Q.S,'t':Q.T,'rule':'p^2=s*p+t'},
            'all48_zero_IFT_quadratic_coefficient_equations':zeroequations,
            'actual_M24_matrix':[[str(z)for z in row]for row in M24.tolist()],
            'M24_determinant':str(det24),
            'unit48_block_determinant_even24_then_odd24_vars_a23_t_b23_s':str(det24**2),
            'positive_group_indices':gpositive,'positive5_equation_physical_rows':rows5,
            'positive5_pair_matrix':matrix5,'positive5_solution':solution,'positive5_unit_Gaussian_trace':solve_trace,
            'actual_tracefree_nu':D.add(Q.T,D.mul(D.const(F(1,4)),D.mul(Q.S,Q.S))),
            'actual_positive_odd_J0_generators':[[row,full[1][row][1]]for row in posodd],
            'actual_zero_odd_I0_generators':[[row,q0[row][1]]for row in odd if row not in posodd],
            'actual_individual_physical_cubic_trace_pairs':[[row,list(Q.tracepair(full[3][row]))]for row in cubics],
            'same_family_trace_basis':'tau=p-s/2; all pairs use(1,tau)',
            'complete_normalized_physical_w_coefficients':[[list(Q.tracepair(z))for z in block]for block in full],
            'positive_J0_origin_values':[str(D.values[full[1][row][1]])for row in posodd],
        }
    certificate={'DAG_scalar_nodes':D.nodes,'DAG_operation_conventions':{'c':'rational constant','v':'unknown/base variable with origin0','a':'sum','m':'product','n':'negation','i':'inverse of recorded origin-unit'},
                 'all_unit_inverse_nodes_and_origin_values':D.divisions,'physical_zero_base_rows9':R9,
                 'positive_base_raw_row':140,'common_positive_odd_physical_rows':posodd,
                 'cubic_graph_raw_rows6':cubics,'families':records,
                 'relative_CB_graph_lift_missing':'Actual individual physical cubic values are NOT CB graph extensions evaluated on CA; b_CB,alpha must first be explicitly expressed in I0_CB+J0_CB with compatible common-target even shifts.'}
    encoded=json.dumps(certificate,separators=(',',':'))+'\n';assert len(encoded.encode())<=32*1024*1024
    (out/'actual-source-coefficients.json').write_text(encoded)
    paths+=[Path(__file__)]
    result={'id':STEM,'all_actual_source_DAG_and_unit_selector_checks_passed':True,
            'status':'actual_J0_and_both_individual_cubic_pairs_constructed_relative_graph_lift_missing',
            'actual_CA_J0_two_generators_as_algebraic_Hensel_DAG':True,
            'actual_nu_and_zero_odd_I0_and_CB_individual_data_also_constructed':True,
            'both_families_48_equations_have_actual_unit_Jacobian':True,
            'all5_positive_variables_solved_by_origin_unit_pivots_without_p_inverse':True,
            'all6_individual_raw_cubic_pairs_in_each_own_trace_basis_constructed':True,
            'requested12_same_CA_basis_CB_minus_CA_graph_pairs_constructed':False,
            'obstruction':'No actual CB cubic graph lift/odd-module membership coefficients serialized; own tau_B pairs cannot be subtracted from tau_A pairs over physical R0.',
            'depth_nil_T3_membership_or_full_D_answer_claimed':False,
            'no_ungauged_R9_rank_rerun_large_GB_formal_series_or_new_point':True,
            'DAG_nodes':len(D.nodes),'certificate_bytes':len(encoded.encode()),
            'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},
            'evidence_sha256':{'actual-source-coefficients.json':sha(out/'actual-source-coefficients.json')},
            'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k]for k in ['status','DAG_nodes','certificate_bytes','seconds','max_rss_kib']}))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True)
    main(p.parse_args().output)
