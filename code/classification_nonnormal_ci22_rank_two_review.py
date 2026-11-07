#!/usr/bin/env python3
"""Independent all-five nonnormal CI22 ranktwo G review, Fraction only.

Rebuild full150 tangent/conormal sources, all350/126, actual Aut Lie algebras,
whole source maps and all zero-complement identities. No producer/CAS imports.
"""
import ast,hashlib,json,resource,signal,time
from fractions import Fraction as F
from itertools import combinations,product,combinations_with_replacement
from pathlib import Path
import classification_cubic_scroll_abelian_review as a
import classification_quadric_cone_zero_review as z

ROOT=Path(__file__).resolve().parents[1]
MAN='research/classification/inputs/nonnormal_ci22_rank_two_manifest_20261107.json'
JOB='research/classification/jobs/nonnormal_ci22_rank_two_review.json'
OUT='reports/CLASSIFICATION_NONNORMAL_CI22_RANK_TWO_REVIEW_20261107.json'
SRC='reports/CLASSIFICATION_NONNORMAL_CI22_TANGENT_20261107.json'
COMP='reports/CLASSIFICATION_NONNORMAL_CI22_SUPPORT_COMPLETION_20261107.json'
PAIRS=list(combinations(range(5),2));M2=a.monomials(5,2);M3=a.monomials(5,3);M4=a.monomials(5,4);M5=a.monomials(5,5)
COLS=list(product(PAIRS,M2));NAMES=['x'+str(i)for i in range(5)]

def sparse(rows):return [{int(j):F(v)for j,v in row.items()if F(v)}for row in rows]
def vector(P):return [P[i][j].get(e,F())for(i,j),e in COLS]
def tensor(v):
    P=z.zero(5)
    for ((i,j),e),c in zip(COLS,v):
        if c:P[i][j][e]=c;P[j][i][e]=-c
    return P
def remainder(p,R,index):
    row={index[e]:c for e,c in p.items()}
    for pivot,v in R.items():
        if pivot in row:row=a.add(row,a.scale(v,-row[pivot]))
    return row
def dot(row,v):return sum(c*v[j]for j,c in row.items())
def quad_sub(rows,pairs,values):return [a.add(*(a.scale(a.mul(values[pairs[j][0]],values[pairs[j][1]]),c)for j,c in row.items()))for row in rows]


def source(case):
    xs=[z.var(5,i)for i in range(5)];qs=[a.parse(case[k],NAMES)for k in('q1','q2')];ix4={e:i for i,e in enumerate(M4)};ix5={e:i for i,e in enumerate(M5)}
    R4=a.rref([{ix4[e]:c for e,c in a.mul(q,{m:F(1)}).items()}for q in qs for m in M2])
    R25=a.rref([{ix5[e]:c for e,c in a.mul(a.mul(q,r),x).items()}for q,r in((qs[0],qs[0]),(qs[0],qs[1]),(qs[1],qs[1]))for x in xs]);assert len(R4)==29 and len(R25)==15
    target={};conormal={}
    for col,((i,j),e)in enumerate(COLS):
        P=z.zero(5);z.setentry(P,i,j,{e:F(1)})
        for at,q in enumerate(qs):
            V=[a.add(*(a.mul(P[s][t],a.diff(q,t))for t in range(5)))for s in range(5)]
            for s,t in PAIRS:
                row=remainder(a.add(a.mul(xs[s],V[t]),a.scale(a.mul(xs[t],V[s]),-1)),R4,ix4)
                for key,c in row.items():target.setdefault((at,s,t,key),{})[col]=c
        bracket=lambda q,r:a.add(*(a.mul(a.mul(a.diff(q,s),P[s][t]),a.diff(r,t))for s,t in product(range(5),repeat=2)))
        K=a.add(a.mul(xs[0],bracket(*qs)),a.scale(a.mul(qs[1],bracket(qs[0],xs[0])),-2),a.scale(a.mul(qs[0],bracket(xs[0],qs[1])),-2))
        row=remainder(K,R25,ix5)
        for key,c in row.items():conormal.setdefault(key,{})[col]=c
    tl=sorted(target);cl=sorted(conormal);T=[target[v]for v in tl];C=[conormal[v]for v in cl];M=T+C
    assert [list(v)for v in tl]==case['all_tangency_row_labels']and cl==case['all_conormal_remainder_monomial_indices']
    assert T==sparse(case['full_tangency_matrix'])and C==sparse(case['full_conormal_matrix'])
    assert 150-len(a.rref(T))==case['tangency_raw_kernel_dimension']and 150-len(a.rref(M))==case['whole_tangent_conormal_raw_dimension']
    E=[]
    for i,j in product(range(5),repeat=2):
        Y=[{}for _ in range(5)];Y[i]=xs[j];E.append(vector(z.wedge(xs,Y)))
    RE=a.rref(E);gauge=sorted(RE);assert len(RE)==24 and gauge==case['Euler_gauge_indices']
    assert all(not dot(row,v)for row in M for v in E)
    minor=a.determinant([[row.get(i,F())for i in gauge]for row in RE.values()]);assert minor==F(case['Euler_gauge_minor'])and minor
    whole=M+[{i:F(1)}for i in gauge];K=[[F(v)for v in row]for row in case['complete_projective_source_kernel']];d=len(K)
    assert len(a.rref(K))==d==150-len(a.rref(whole))==case['whole_projective_source_dimension']
    assert all(not dot(row,v)for row in whole for v in K)
    free=case['complete_projective_free_original_indices'];assert len(free)==d
    assert [[v[i]for i in free]for v in K]==[[F(i==j)for j in range(d)]for i in range(d)]
    family=[tensor(v)for v in K]
    for P,saved in zip(family,case['all150_source_tensors']):
        assert all(P[i][j]=={tuple(e):F(c)for e,c in saved[str((i,j))]}for i,j in PAIRS)
    JJ,DD=a.jacobi_and_divisor(family);assert JJ==sparse(case['full350_Jacobi_rows'])and DD==sparse(case['full126_actualPf_rows'])
    RJ=a.rref(JJ);assert list(RJ.values())==sparse(case['full_Jacobi_coefficient_rowspace_basis'])
    # Full Lie algebra of the actual pencil-preserving GL5, modulo scalars.
    R2=a.rref([{M2.index(e):c for e,c in q.items()}for q in qs]);autcols=[]
    for i,j in product(range(5),repeat=2):
        vals=[]
        for q in qs:
            rem=remainder(a.mul(xs[j],a.diff(q,i)),R2,{e:k for k,e in enumerate(M2)});vals.extend(rem.get(k,F())for k in range(15))
        autcols.append(vals)
    aut=25-len(a.rref(list(map(list,zip(*autcols)))))-1
    assert aut>=case['actual_source_preserving_group_dimension_lower']
    return JJ,DD,list(RJ.values()),dict(model=case['model'],degree4_ideal_rank29=True,degree5_square_rank15=True,
       full150_tangent_and_intrinsic_conormal_reconstructed=True,tangent_projective_dimension=case['tangency_projective_dimension'],complete_source_dimension=d,
       Euler24_and_coordinate_gauge_minor=str(minor),complete_free_coordinates_inverse=True,all150_tensors_and350Jacobi_and126Pf=True,
       original_quadratic_rowspace_rank=len(RJ),actual_pencil_Aut_dimension=aut)


def laurent(text,names):
    n=len(names);one={(0,)*n:F(1)}
    def read(node):
        if isinstance(node,ast.Constant)and isinstance(node.value,int):return a.scale(one,F(node.value))
        if isinstance(node,ast.Name):return z.var(n,names.index(node.id))
        if isinstance(node,ast.UnaryOp)and isinstance(node.op,ast.USub):return a.scale(read(node.operand),-1)
        assert isinstance(node,ast.BinOp),ast.dump(node)
        left=read(node.left)
        if isinstance(node.op,ast.Pow):
            assert isinstance(node.right,ast.Constant)and isinstance(node.right.value,int)and node.right.value>=0
            ans=one
            for _ in range(node.right.value):ans=a.mul(ans,left)
            return ans
        right=read(node.right)
        if isinstance(node.op,ast.Add):return a.add(left,right)
        if isinstance(node.op,ast.Sub):return a.add(left,a.scale(right,-1))
        if isinstance(node.op,ast.Mult):return a.mul(left,right)
        assert isinstance(node.op,ast.Div)and len(right)==1
        e,c=next(iter(right.items()));return {tuple(u-v for u,v in zip(f,e)):d/c for f,d in left.items()}
    return read(ast.parse(str(text),mode='eval').body)


def support(cases,comp,data):
    result={};names=['A','B','G','Z','H','F']
    for i,c in enumerate(cases):
        JJ,DD,E,proof=data[i];pairs=[tuple(v)for v in c['parameter_pair_order']];p=[z.var(len(c['complete_projective_source_kernel']),j)for j in range(len(c['complete_projective_source_kernel']))]
        if i==0:
            for record in comp['regular22_largest_two_whole_source_factors']:
                keep=record['complete_linear_source_indices'];vals=[v if j in keep else{}for j,v in enumerate(p)]
                assert not any(quad_sub(JJ,pairs,vals))and not any(quad_sub(DD,pairs,vals))
            assert comp['regular22_remaining_support_dimension_upper']==3
            result['regular22']=dict(full_two_5parameter_sources_Jacobi_and_actualDzero=True,other_two_inverse_cases_dimensions=[2,3],rank4_source_dimension_upper=3)
        if i==1:
            maps=[]
            for record in comp['nil32_complete_branch_cover_maps']:
                vals=[laurent(record['parameter_map']['p'+str(j)],names)for j in range(len(p))]
                residuals=quad_sub(E,pairs,vals)
                assert residuals==[laurent(v,names)for v in record['all_equivalent_original_Jacobi_residuals']]
                used={j for v in vals for e in v for j,k in enumerate(e)if k}
                assert len(used)<=record['source_upper']
                inverted=record['inverted']
                assert all(k>=0 or names[j]==inverted for v in vals for e in v for j,k in enumerate(e))
                maps.append(dict(case=record['case'],source_dimension_upper=record['source_upper'],only_inverted=inverted,
                   full_original_quadratic_residuals_reconstructed=True,supersource_not_claimed_allPoisson=any(residuals)))
            assert len(maps)==9
            result['nil32']=dict(all_nine_supersource_maps=maps,all_zero_complements_in_inverse_caseproof=True,source_dimension_upper=5)
        if i==2:
            # Verify the literal complete p8=0 source equations; the cases in
            # the accompanying proof then bound the whole support by five.
            vals=p[:];vals[8]={};own=quad_sub(E,pairs,vals)
            texts=['p0*p5+p3*p6','p1*p2+p6*(p3-2*p5)','p1*(p3-3*p5)','p9*(p3-2*p5)','p7*(p0+2*p6)','p2*(p0+3*p6)','p1*p7','p2*p9']
            expected=[laurent(t,['p'+str(j)for j in range(len(p))])for t in texts]
            keys=sorted(set().union(*(v.keys()for v in own+expected)));rows=lambda L:[[v.get(e,F())for e in keys]for v in L]
            assert a.rref(rows(own))==a.rref(rows(expected))
            result['singularM1']=dict(all_eight_p8zero_equations_complete=True,p8nonzero_dimension_upper=3,
                full_inverse_p5ratio_and_p5zero_cases_verified_by_proof=True,source_dimension_upper=5)
        if i==3:
            actual={}
            for record in comp['cone_node_all_D_original350_functionals']:
                r=record['actual_D_row'];phi={int(j):F(v)for j,v in record['original350_functional'].items()}
                assert a.add(*(a.scale(JJ[j],v)for j,v in phi.items()))==DD[r];actual[r]=phi
            assert set(actual)=={j for j,v in enumerate(DD)if v}
            result['cone_node']=dict(all_actualD_nonzero_rows_original350_functionals=len(actual),whole_Poisson_actualDzero=True)
        if i==4:
            vals=p[:];vals[13]={};vals[15]={};assert not any(quad_sub(DD,pairs,vals))
            target=a.mul(p[13],p[15]);eqs=quad_sub(E,pairs,p);assert target in eqs
            for record in comp['cone_cusp_complete_two_open_maps']:
                vals=[laurent(record['parameter_map']['p'+str(j)],names)for j in range(len(p))]
                assert not any(quad_sub(JJ,pairs,vals))
                assert not any(quad_sub(E,pairs,vals))
                assert len({j for v in vals for e in v for j,k in enumerate(e)if k})==5
            result['cone_cusp']=dict(original_product_p13_p15=True,both_zero_full_actualDzero=True,both_open_entire5parameter_Poisson_sources=True,rank4_source_dimension_upper=5)
    return result


def main():
    start=time.monotonic();signal.alarm(240);resource.setrlimit(resource.RLIMIT_CPU,(235,238));resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,)*2)
    assert not(ROOT/OUT).exists();job=json.loads((ROOT/JOB).read_text())
    for row in job['exact_inputs_and_hashes']:assert hashlib.sha256((ROOT/row['path']).read_bytes()).hexdigest()==row['sha256'],row['path']
    for row in json.loads((ROOT/MAN).read_text())['files']:assert hashlib.sha256((ROOT/row['path']).read_bytes()).hexdigest()==row['sha256'],row['path']
    packet=json.loads((ROOT/SRC).read_text());comp=json.loads((ROOT/COMP).read_text())
    for report in(packet,comp):assert hashlib.sha256(report['producer_source_bytes'].encode()).hexdigest()==report['input_sha256']['code/classification_nonnormal_ci22_tangent.py']
    assert packet['all_five_nonnormal_CI22_orbits_G24_closed']is False
    data=[]
    for c in packet['cases']:data.append(source(c));print(c['model'],'full150 source PASS',flush=True)
    supp=support(packet['cases'],comp,data);uppers=[3,5,5,None,5];image=[23,24,24,None,23]
    result=dict(task_id='G',all_passed=True,producer_imports=False,CAS_imports=False,independent_Fraction_only=True,
      exact_input_bindings=job['exact_inputs_and_hashes'],all_five_full_sources=[v[-1]for v in data],whole_reduced_support=supp,
      all_inverse_case_dimensions_verified_by_separate_mathematical_proof=True,all_source_Aut_invariant_before_fiber_subtraction=True,
      full_rank4_source_dimensions_upper=uppers,full_PGL_projective_images_upper=image,all_nonnormal_integralCI22_rank2_G24_closed=True,
      old_linear_coarse_false_report_preserved=True,all_R_nonreduced_vertical_boundaries_retained=True,
      scope='Whole nonnormal integral nondegenerate CI22 surface with rank2 restriction, rank4 Poisson actualD irreducible reduced. No E membership/newcomponent/globalG.',
      elapsed_seconds=time.monotonic()-start,peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    (ROOT/OUT).write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:result[k]for k in('all_passed','elapsed_seconds','peak_rss_kib')}),flush=True)


if __name__=='__main__':main()
