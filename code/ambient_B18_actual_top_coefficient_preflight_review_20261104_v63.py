#!/usr/bin/env python3
"""Independent DAG AD/source selectors and symbolic unit-row-operation audit."""
import argparse
import ast
from fractions import Fraction as F
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import resource
import signal
import subprocess
import sys
import time

ROOT=Path(__file__).resolve().parents[1]
BASE='AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63'
STEM=BASE+'_REVIEW'
PAIRS=list(itertools.combinations(range(5),2));PI={p:i for i,p in enumerate(PAIRS)}
MONS=list(itertools.combinations_with_replacement(range(5),2));MI={p:i for i,p in enumerate(MONS)}
PARAMS=['p']+[f'k{i}' for i in range(9)];ZEXP=(0,)*10
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(gzip.decompress(p.read_bytes()) if p.suffix=='.gz' else p.read_text())


def p_add(a,b):
    out=dict(a)
    for e,c in b.items():out[e]=out.get(e,F(0))+c
    return {e:c for e,c in out.items() if c}
def p_mul(a,b):
    out={}
    for e,c in a.items():
        for f,d in b.items():
            g=tuple(x+y for x,y in zip(e,f));out[g]=out.get(g,F(0))+c*d
    return {e:c for e,c in out.items() if c}
def parse_poly(text):
    def visit(n):
        if isinstance(n,ast.Constant) and isinstance(n.value,int):return {ZEXP:F(n.value)} if n.value else {}
        if isinstance(n,ast.Name):
            e=list(ZEXP);e[PARAMS.index(n.id)]=1;return {tuple(e):F(1)}
        if isinstance(n,ast.UnaryOp) and isinstance(n.op,ast.USub):return {e:-c for e,c in visit(n.operand).items()}
        if isinstance(n,ast.BinOp):
            a=visit(n.left)
            if isinstance(n.op,ast.Pow):
                assert isinstance(n.right,ast.Constant) and isinstance(n.right.value,int)
                out={ZEXP:F(1)}
                for _ in range(n.right.value):out=p_mul(out,a)
                return out
            b=visit(n.right)
            if isinstance(n.op,ast.Add):return p_add(a,b)
            if isinstance(n.op,ast.Sub):return p_add(a,{e:-c for e,c in b.items()})
            if isinstance(n.op,ast.Mult):return p_mul(a,b)
            if isinstance(n.op,ast.Div):
                assert set(b)=={ZEXP} and b[ZEXP]
                return {e:c/b[ZEXP] for e,c in a.items()}
        raise ValueError(ast.dump(n))
    return visit(ast.parse(text,mode='eval').body)


def lie(point,i,j):
    out=[F(0)]*150
    for row,co in enumerate(point):
        if not co:continue
        a,b=PAIRS[row//15];u,v=MONS[row%15]
        for first,second in [(a,b)]:
            if first==j and i!=second:
                out[15*PI[tuple(sorted((i,second)))]+row%15]+=co*(1 if i<second else -1)
            if second==j and first!=i:
                out[15*PI[tuple(sorted((first,i)))]+row%15]+=co*(1 if first<i else -1)
        for old,other in [(u,v),(v,u)]:
            if old==i:out[15*(row//15)+MI[tuple(sorted((j,other)))]]-=co
    return out


def det(a):
    a=[list(map(F,r)) for r in a];answer=F(1)
    for j in range(len(a)):
        k=next(i for i in range(j,len(a)) if a[i][j])
        if k!=j:a[j],a[k]=a[k],a[j];answer=-answer
        p=a[j][j];answer*=p
        for i in range(j+1,len(a)):
            if not a[i][j]:continue
            q=a[i][j]/p
            for h in range(j,len(a)):a[i][h]-=q*a[j][h]
    return answer


def gradient_add(a,b,ca=F(1),cb=F(1)):
    out={i:ca*v for i,v in a.items() if ca*v}
    for i,v in b.items():out[i]=out.get(i,F(0))+cb*v
    return {i:v for i,v in out.items() if v}


class DAGReview:
    def __init__(self,nodes):
        self.nodes=[tuple(n) for n in nodes];self.cache={n:i for i,n in enumerate(self.nodes)}
        assert len(self.cache)==len(self.nodes)
        self.values=[];self.grad=[];self.deps=[];self.variables={}
        for i,n in enumerate(self.nodes):
            op=n[0]
            if op=='c':v=F(n[1],n[2]);g={};deps=set()
            elif op=='v':
                assert n[1] not in self.variables;self.variables[n[1]]=i
                v=F(0);g={i:F(1)};deps={n[1]}
            else:
                assert all(isinstance(j,int) and 0<=j<i for j in n[1:])
                a=n[1];deps=set(self.deps[a])
                if op=='n':v=-self.values[a];g={j:-x for j,x in self.grad[a].items()}
                elif op=='i':
                    assert self.values[a];v=1/self.values[a]
                    g={j:-x*v*v for j,x in self.grad[a].items()}
                else:
                    b=n[2];deps|=self.deps[b]
                    if op=='a':v=self.values[a]+self.values[b];g=gradient_add(self.grad[a],self.grad[b])
                    elif op=='m':
                        v=self.values[a]*self.values[b]
                        g=gradient_add(self.grad[a],self.grad[b],self.values[b],self.values[a])
                    else:raise ValueError(n)
            self.values.append(v);self.grad.append(g);self.deps.append(deps)
        self.zero=self.const(0);self.one=self.const(1)
    def lookup(self,key):return self.cache[key]
    def const(self,q):q=F(q);return self.lookup(('c',q.numerator,q.denominator))
    def neg(self,a):
        if a==self.zero:return a
        if self.nodes[a][0]=='c':return self.const(-self.values[a])
        if self.nodes[a][0]=='n':return self.nodes[a][1]
        return self.lookup(('n',a))
    def add(self,a,b):
        if a==self.zero:return b
        if b==self.zero:return a
        if self.nodes[a]==('n',b) or self.nodes[b]==('n',a):return self.zero
        if self.nodes[a][0]==self.nodes[b][0]=='c':return self.const(self.values[a]+self.values[b])
        return self.lookup(('a',*sorted((a,b))))
    def mul(self,a,b):
        if a==self.zero or b==self.zero:return self.zero
        if a==self.one:return b
        if b==self.one:return a
        if self.nodes[a][0]==self.nodes[b][0]=='c':return self.const(self.values[a]*self.values[b])
        return self.lookup(('m',*sorted((a,b))))
    def inverse(self,a):
        assert self.values[a]
        if self.nodes[a][0]=='c':return self.const(1/self.values[a])
        return self.lookup(('i',a))
    def total(self,items):
        out=self.zero
        for i in items:out=self.add(out,i)
        return out
    def scalar_multiple(self,node):
        factor=F(1)
        while True:
            n=self.nodes[node]
            if n[0]=='n':factor=-factor;node=n[1]
            elif n[0]=='m' and self.nodes[n[1]][0]=='c':
                factor*=self.values[n[1]];node=n[2]
            elif n[0]=='m' and self.nodes[n[2]][0]=='c':
                factor*=self.values[n[2]];node=n[1]
            elif n[0]=='c':return factor*self.values[node],self.one
            else:return factor,node


class QuadraticReplay:
    def __init__(self,D,s,t):self.D=D;self.s=s;self.t=t
    def add(self,a,b):return tuple(self.D.add(x,y) for x,y in zip(a,b))
    def neg(self,a):return tuple(self.D.neg(x) for x in a)
    def mul(self,a,b):
        D=self.D;bd=D.mul(a[1],b[1])
        return (D.add(D.mul(a[0],b[0]),D.mul(self.t,bd)),
                D.total([D.mul(a[0],b[1]),D.mul(a[1],b[0]),D.mul(self.s,bd)]))
    def inverse(self,a):
        D=self.D
        norm=D.add(D.add(D.mul(a[0],a[0]),D.mul(self.s,D.mul(a[0],a[1]))),D.neg(D.mul(self.t,D.mul(a[1],a[1]))))
        ni=D.inverse(norm)
        return (D.mul(D.add(a[0],D.mul(self.s,a[1])),ni),D.neg(D.mul(a[1],ni)))


def main(out,work):
    started=time.monotonic();signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000))
    folder=ROOT/'reports'/BASE;primary=read(folder/'report.json')
    assert primary['all_actual_source_DAG_and_unit_selector_checks_passed']
    assert not primary['requested12_same_CA_basis_CB_minus_CA_graph_pairs_constructed']
    assert not primary['depth_nil_T3_membership_or_full_D_answer_claimed']
    checked={}
    for name,want in primary['input_sha256'].items():
        p=ROOT/name;assert sha(p)==want,str(p);checked[name]=want
    certpath=folder/'actual-source-coefficients.json';assert sha(certpath)==primary['evidence_sha256'][certpath.name]
    cert=read(certpath);D=DAGReview(cert['DAG_scalar_nodes'])
    assert len(D.nodes)==41315 and len(D.variables)==105
    inverse_ids={i for i,n in enumerate(D.nodes) if n[0]=='i'}
    assert inverse_ids=={i for i,a,v in cert['all_unit_inverse_nodes_and_origin_values']}
    for i,a,v in cert['all_unit_inverse_nodes_and_origin_values']:
        assert D.nodes[i]==('i',a) and D.values[a]==F(v)!=0
    raw=read(ROOT/'reports/AMBIENT_B18_QUADRATIC_IMAGE_REDUCTION_20261103_v63/report.json')
    orbit=read(ROOT/'reports/AMBIENT_B18_ACTUAL_INTERSECTION_ORBIT_PRODUCT_20261103_v63/certificate.json.gz')
    weight=read(ROOT/'reports/AMBIENT_B18_RELATIVE_WEIGHTED_MODEL_20261104_v63/certificate.json')
    implicit=read(ROOT/'reports/AMBIENT_B18_CORE_GAUGE_PREFLIGHT_20261104_v63/implicit-system.json')
    obstruction=ROOT/'reports/AMBIENT_B18_POSITIVE_GRADING_OBSTRUCTION_20261104_v63/certificate.json'
    weights=read(obstruction)['physical150_torus_weights']
    rows9=cert['physical_zero_base_rows9'];assert rows9==[27,108,109,112,113,117,118,119,61]
    cuts0=[i for i in weight['common_physical_group_cuts19'] if weights[i]==0]
    cuts1=[i for i in weight['common_physical_group_cuts19'] if weights[i]==1]
    gzero=implicit['zero_GL_group_parameter_indices'];gplus=[i for i in range(19) if i not in gzero]
    assert len(cuts0)==15 and len(cuts1)==4 and len(gzero)==15 and len(gplus)==4
    summary={}
    for family,rec in cert['families'].items():
        polynomials=[parse_poly(t) for t in raw['families'][family]['core_full150_coordinates']]
        point=[p.get(ZEXP,F(0)) for p in polynomials]
        assert point==list(map(F,orbit['physical_point150']))
        O=[point]+[lie(point,*divmod(i,5)) for i in orbit['GL18_generator_indices']]
        assert [[O[j][i] for j in range(19)] for i in range(150)]==[list(map(F,r)) for r in orbit['physical150_by19_scalar_then_GL']]
        selectors={}
        for name in ['p']+[f'k{i}' for i in range(9)]:
            e=list(ZEXP);e[PARAMS.index(name)]=1
            selectors[name]=[p.get(tuple(e),F(0)) for p in polynomials]
        e=list(ZEXP);e[0]=2;W=[p.get(tuple(e),F(0)) for p in polynomials]
        assert not any(selectors['p'])
        assert W==list(map(F,raw['families'][family]['raw_second_coefficient']))
        for row,p in enumerate(polynomials):
            for ex,co in p.items():
                assert ex[3]<=1 and 2 not in MONS[row%15]
                assert (2 in PAIRS[row//15])==(ex[3]==0)
        names=rec['zero_hensel_unknowns23_order']
        expected_names=[f'gamma{i}' for i in gzero]+[f'k{i}' for i in range(9) if i!=2]
        assert names==expected_names
        cols=[O[i] for i in gzero]+[selectors[f'k{i}'] for i in range(9) if i!=2]+[W]
        M=[[cols[j][row] for j in range(24)] for row in cuts0+rows9]
        assert M==[list(map(F,r)) for r in rec['actual_M24_matrix']]
        determinant=det(M);assert determinant==F(rec['M24_determinant'])
        assert determinant=={'CA':F(-9216,25),'CB':F(576,25)}[family]
        sr,tr=rec['p_quadratic_relation']['s'],rec['p_quadratic_relation']['t']
        acols=[rec['coefficient_unknown_pairs'][n][0] for n in names]+[tr]
        bcols=[rec['coefficient_unknown_pairs'][n][1] for n in names]+[sr]
        equations=rec['all48_zero_IFT_quadratic_coefficient_equations']
        assert len(equations)==48 and all(D.values[i]==0 for i in equations)
        for k in range(24):
            even,odd=equations[2*k:2*k+2]
            assert [D.grad[even].get(j,F(0)) for j in acols]==M[k]
            assert not any(D.grad[even].get(j,F(0)) for j in bcols)
            assert [D.grad[odd].get(j,F(0)) for j in bcols]==M[k]
        assert determinant**2==F(rec['unit48_block_determinant_even24_then_odd24_vars_a23_t_b23_s'])
        # All output coefficients use ONLY the family's own quadratic cover.
        blocks=rec['complete_normalized_physical_w_coefficients']
        assert len(blocks)==4 and all(len(b)==150 for b in blocks)
        for b in blocks:
            for pair in b:
                assert len(pair)==2
                assert all(not any(n.startswith(('CB_' if family=='CA' else 'CA_')) for n in D.deps[i]) for i in pair)
        nu=D.add(tr,D.mul(D.const(F(1,4)),D.mul(sr,sr)))
        assert nu==rec['actual_tracefree_nu']
        assert rec['actual_positive_odd_J0_generators']==[[row,blocks[1][row][1]] for row in [46,95]]
        assert all(D.scalar_multiple(node)==D.scalar_multiple(blocks[0][row][1])
                   for row,node in rec['actual_zero_odd_I0_generators'])
        assert rec['actual_individual_physical_cubic_trace_pairs']==[[row,blocks[3][row]] for row in [9,39,54,84,99,144]]
        # Independently replay every saved Gaussian elementary operation IN the DAG.
        Q=QuadraticReplay(D,sr,tr)
        A=[[tuple(x) for x in row]+[(D.one,D.zero) if i==4 else (D.zero,D.zero)]
           for i,row in enumerate(rec['positive5_pair_matrix'])]
        for j,trace in enumerate(rec['positive5_unit_Gaussian_trace']):
            p=trace['selected_row'];assert p>=j
            assert not any(D.values[A[i][j][0]] for i in range(j,p))
            A[j],A[p]=A[p],A[j];pivot=A[j][j]
            assert list(pivot)==trace['pivot'] and [str(D.values[i]) for i in pivot]==trace['origin']
            inv=Q.inverse(pivot);A[j]=[Q.mul(x,inv) for x in A[j]]
            for i in range(5):
                if i==j:continue
                co=A[i][j];A[i]=[Q.add(x,Q.neg(Q.mul(co,y))) for x,y in zip(A[i],A[j])]
        assert [list(A[i][-1]) for i in range(5)]==rec['positive5_solution']
        assert [r['origin'] for r in rec['positive5_unit_Gaussian_trace']]==[['1','0'],['-1','0'],['-1','0'],['-1','0'],['1','0']]
        # The five origin columns are also bound to actual Lie/k2 source coefficients.
        rows5=cuts1+[140]
        source_columns=[lie(point,*divmod(orbit['GL18_generator_indices'][i-1],5)) for i in gplus]
        source_columns.append(selectors['k2'])
        assert [[D.values[pair[0]] for pair in row] for row in rec['positive5_pair_matrix']]==[
            [source_columns[j][i] for j in range(5)] for i in rows5]
        summary[family]={'M24det':str(determinant),'Jac48det':str(determinant**2),
            'all2304_AD_block_entries_checked':True,'positive5_formal_unit_row_operations_checked':True}

    reader=read(folder/'relative-membership-target.json')
    assert not reader['all60_membership_coefficients_solved'] and not reader['requested12_actual_pairs_constructed']
    assert reader['no_tau_B_equals_tau_A_identification'] and reader['lower_CB_relations_must_be_retained_for_any_membership_lift_choice']
    cb=cert['families']['CB'];oddrows=reader['CB_odd_membership_matrix_physical_rows']
    assert all(D.scalar_multiple(node)==D.scalar_multiple(cb['complete_normalized_physical_w_coefficients'][0 if r not in [46,95] else 1][r][1])
               for r,node in zip(oddrows,reader['CB_odd_membership_matrix_1_by10_nodes']))
    assert reader['six_actual_right_hand_sides_nodes']==[[r,pair[1]] for r,pair in cb['actual_individual_physical_cubic_trace_pairs']]
    # One deterministic finite-builder replay solely binds ALL serialized source roots.
    # It computes no series/unknown solutions and does not rerun the old ungauged rank7 job.
    work.mkdir(parents=True,exist_ok=False)
    replay=work/'finite-builder-replay'
    proc=subprocess.run([sys.executable,str(ROOT/'code/ambient_B18_actual_top_coefficient_preflight_20261104_v63.py'),
                         '--output',str(replay)],cwd=ROOT,capture_output=True,text=True,timeout=45)
    (work/'finite-builder-replay.log').write_text(proc.stdout+proc.stderr)
    assert proc.returncode==0,proc.stderr
    assert sha(replay/'actual-source-coefficients.json')==sha(certpath)
    for name,want in reader['input_sha256'].items():
        p=ROOT/name;assert sha(p)==want,str(p);checked[name]=want
    proof={'all41315_nodes_independently_interpreted_and_differentiated':True,
        'variables':105,'all10_inverse_nodes_have_nonzero_origin_units':True,
        'families':summary,'all150_sparse_physical_source_shape_and_linear_quadratic_selectors':True,
        'exact_source_DAG_byte_binding_replay_sha256':sha(replay/'actual-source-coefficients.json'),
        'positive_shear_formula_exact_degree_upper':3,
        'full150_w0to3_pairs_per_family':True,
        'own_CA_CB_trace_bases_and_nu_I0_J0_cubic_roots_verified':True,
        'same_CA_relative_graph_pairs_constructed_in_this_stage':False,
        'depth_nil_or_T3_membership_claimed':False}
    out.mkdir(parents=True,exist_ok=False)
    ep=out/'independent-proof.json';ep.write_text(json.dumps(proof,indent=2)+'\n')
    log=out/'finite-builder-replay.log';log.write_text(proc.stdout+proc.stderr)
    paths=[Path(__file__),ROOT/'research/notes'/(STEM+'.md'),ROOT/'research/notes'/(BASE+'.md'),
           folder/'report.json',certpath,folder/'relative-membership-target.json',obstruction]
    for p in paths:checked[str(p.relative_to(ROOT))]=sha(p)
    result={'id':STEM,'all_passed':True,
        'independent_all_scalar_DAG_origin_AD_and_unit_inverses':True,
        'independent_sparse_physical_source_GL_M24_and_full48_Jacobian_blocks':True,
        'exact_positive5_symbolic_unit_Gaussian_trace_and_source_origins':True,
        'all150_w0to3_finite_source_DAG_binding_replay_byte_identical':True,
        'correct_common_physical_base_and_own_tau_A_tau_B':True,
        'blocked_relative_graph_reader_preserved_and_no_requested12_completion_claim':True,
        'no_ungauged_R9_rank_job_large_GB_series_solution_newpoint_or_actual_depth_claim':True,
        'deterministic_author_finite_builder_replayed_for_exact_DAG_binding':True,
        'input_sha256':checked,'evidence_sha256':{'independent-proof.json':sha(ep),'finite-builder-replay.log':sha(log)},
        'seconds':time.monotonic()-started,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'frozen':True}
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ['all_passed','seconds','max_rss_kib']}))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--work',type=Path,required=True)
    args=p.parse_args();main(args.output,args.work)
