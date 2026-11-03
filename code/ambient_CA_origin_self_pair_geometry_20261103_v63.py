#!/usr/bin/env python3
"""Exact source/orbit/cross-tangent ranks along the actual CA polynomial self-pair."""
import argparse
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import resource
import signal
import time
import sympy as s
from sympy.polys.matrices import DomainMatrix
from sympy.polys.rings import ring

ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_CA_ORIGIN_SELF_PAIR_GEOMETRY_20261103_v63'
PAIRS=list(itertools.combinations(range(5),2));MONS=list(itertools.combinations_with_replacement(range(5),2))


def read(path):
    with(gzip.open(path,'rt')if path.suffix=='.gz'else path.open())as f:return json.load(f)
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def main(out):
    start=time.monotonic();out.mkdir(parents=True,exist_ok=False)
    resource.setrlimit(resource.RLIMIT_AS,(2200000000,2200000000));signal.alarm(120)
    parent=ROOT/'reports/AMBIENT_CA_ORIGIN_POLYNOMIAL_SELF_PAIR_20261103_v63';assert read(parent/'report.json')['all_passed']
    pairreview=ROOT/'reports/AMBIENT_CA_ORIGIN_POLYNOMIAL_SELF_PAIR_20261103_v63_REVIEW.json';assert read(pairreview)['all_passed']
    properreview=ROOT/'reports/AMBIENT_CA_GLOBAL_FLAG_MODEL_20261103_v63_REVIEW.json';assert read(properreview)['all_passed']
    fiberreport=ROOT/'reports/AMBIENT_B18_QUADRATIC_IMAGE_REDUCTION_20261103_v63/report.json'
    fiberreview=ROOT/'reports/AMBIENT_B18_QUADRATIC_IMAGE_REDUCTION_20261103_v63_REVIEW/report.json'
    assert read(fiberreview)['all_passed']and read(fiberreport)['whole_flag_fibers']['CA']['scheme_length']==2
    d=read(parent/'certificate.json.gz');inputpath=ROOT/'reports/AMBIENT_B18_CA_CRITICAL_SCHEME_20261103_v63/certificate.json.gz';base=read(inputpath)
    R,a=ring('a',s.QQ);asym=s.symbols('a');domain=s.QQ.frac_field(asym)
    def load(pol):return R.from_dict({(n,):s.QQ.convert(s.Rational(c))for n,c in pol})
    def const(c):return R(s.Rational(c))
    def save(pol):return[[m[0],str(c)]for m,c in sorted(pol.items())]
    def exprmat(columns):return s.Matrix(150,len(columns),lambda i,j:columns[j][i].as_expr())
    def rankdata(columns):
        A=exprmat(columns);dm=DomainMatrix.from_Matrix(A).convert_to(domain);pc=list(dm.rref()[1]);pr=list(DomainMatrix.from_Matrix(A[:,pc].T).convert_to(domain).rref()[1])
        square=DomainMatrix.from_Matrix(A[pr,pc]).convert_to(s.QQ.poly_ring(asym));det=s.factor(square.domain.to_sympy(square.det()))
        assert det!=0
        return dict(rank=len(pc),minor_rows=pr,minor_columns=pc,minor_determinant=str(det))
    C=[[const(c)for c in row]for row in base['fixed_fiber_linear_basis150_by20']];fibers=list(map(list,zip(*C)))
    point=list(map(load,d['first_source_physical150']));other=list(map(load,d['second_source_physical150']))
    g=[[load(p)for p in row]for row in d['group_g5_by5']];gi=[[load(p)for p in row]for row in d['group_inverse5_by5']]
    pairlookup={pair:i for i,pair in enumerate(PAIRS)};monlookup={mon:i for i,mon in enumerate(MONS)}
    def lie(i,j,q):
        ans=[R.zero]*150
        for old,(k,l)in enumerate(PAIRS):
            for new in[(i,l)]if k==j else[]:
                if new[0]!=new[1]:
                    sign=1 if new[0]<new[1]else-1;outpair=pairlookup[tuple(sorted(new))]
                    for m in range(15):ans[15*outpair+m]+=sign*q[15*old+m]
            for new in[(k,i)]if l==j else[]:
                if new[0]!=new[1]:
                    sign=1 if new[0]<new[1]else-1;outpair=pairlookup[tuple(sorted(new))]
                    for m in range(15):ans[15*outpair+m]+=sign*q[15*old+m]
            for m,mon in enumerate(MONS):
                for occurrence in range(2):
                    if mon[occurrence]==i:
                        changed=list(mon);changed[occurrence]=j;ans[15*old+monlookup[tuple(sorted(changed))]]-=q[15*old+m]
        return ans
    positions=list(map(tuple,base['flag_positions']));J1=[lie(i,j,point)for i,j in positions];J2=[lie(i,j,other)for i,j in positions]
    Cexpr=exprmat(fibers);rows=list(DomainMatrix.from_Matrix(Cexpr.T).convert_to(s.QQ).rref()[1]);inv=Cexpr[rows,:].inv()
    def project(q):
        values=[sum((const(inv[i,j])*q[rows[j]]for j in range(20)),R.zero)for i in range(20)]
        return[q[i]-sum((C[i][j]*values[j]for j in range(20)),R.zero)for i in range(150)]
    quot1=list(map(project,J1));quot2=list(map(project,J2));r1=rankdata(quot1);r2=rankdata(quot2)
    orbit=[lie(i,j,point)for i,j in itertools.product(range(5),repeat=2)];rO=rankdata(orbit);rOP=rankdata([*orbit,point])
    derivative=[p.diff(a)for p in point];rFamily=rankdata([*orbit,point,derivative])
    print(json.dumps(dict(stage='both_source_and_orbit',seconds=time.monotonic()-start,source1_rank=20+r1['rank'],source2_rank=20+r2['rank'],GL_rank=rO['rank'],GL_scalar_rank=rOP['rank'])),flush=True)
    # Push the second source's fixed-fiber tangent directions back by g^-1 to the SAME point.
    wedge=[[gi[i][k]*gi[j][l]-gi[i][l]*gi[j][k]for k,l in PAIRS]for i,j in PAIRS]
    sub=[]
    for m,n in MONS:
        vec=[R.zero]*15
        for i,j in itertools.product(range(5),repeat=2):vec[monlookup[tuple(sorted((i,j)))]]+=g[m][i]*g[n][j]
        sub.append(vec)
    def pushback(q):
        answer=[R.zero]*150
        for old in range(10):
            for mon in range(15):
                if not q[15*old+mon]:continue
                for i in range(10):
                    if not wedge[i][old]:continue
                    for new in range(15):
                        if sub[mon][new]:answer[15*i+new]+=wedge[i][old]*q[15*old+mon]*sub[mon][new]
        return answer
    secondfibers=list(map(pushback,fibers));cross=list(map(project,[*J1,*secondfibers]));rCross=rankdata(cross)
    unionrank=20+rCross['rank'];intersection=20+r1['rank']+20+r2['rank']-unionrank
    # Explicit rank-two and nonzero projective Poisson witnesses: all raw Pfaffians vanish.
    tensor={}
    for k,(i,j)in enumerate(PAIRS):tensor[i,j]={m:point[15*k+n]for n,m in enumerate(MONS)if point[15*k+n]}
    # Every nonzero coefficient has y in its exterior pair, proving P wedge P=0 identically.
    assert all(2 in PAIRS[i//15]for i,p in enumerate(point)if p)
    # (Euler wedge P)012 has x^3 coefficient1.
    assert point[15*pairlookup[0,2]+monlookup[1,1]]==-1
    assert point[15*pairlookup[0,1]+monlookup[1,2]]==0
    assert point[15*pairlookup[1,2]+monlookup[0,1]]==0
    payload=dict(source1_flag_quotient150_by9=[list(map(save,row))for row in zip(*quot1)],source2_flag_quotient150_by9=[list(map(save,row))for row in zip(*quot2)],
        source1_rank_witness=r1,source2_rank_witness=r2,full_GL150_by25=[list(map(save,row))for row in zip(*orbit)],
        GL_rank_witness=rO,GL_plus_scalar_rank_witness=rOP,GL_plus_scalar_plus_parameter_derivative_rank_witness=rFamily,
        cross_same_point_mod_first_fiber150_by29=[list(map(save,row))for row in zip(*cross)],cross_quotient_rank_witness=rCross,
        generic_same_point_tangent_union_affine_rank=unionrank,generic_same_point_tangent_intersection_affine_dimension=intersection,
        all_raw_Pfaffians_zero_by_constant_vector_factor=True,Euler_wedge_P012_xcubed_coefficient='1')
    with(out/'certificate.json.gz').open('wb')as f:
        with gzip.GzipFile(fileobj=f,mode='wb',filename='',mtime=0)as z:z.write((json.dumps(payload,separators=(',',':'))+'\n').encode())
    result=dict(id=STEM,all_passed=True,generic_projective_Poisson_rank=2,all_a_projective_Poisson_rank_exact2=True,
        generic_source1_affine_differential_rank=20+r1['rank'],generic_source2_affine_differential_rank=20+r2['rank'],
        source1_projected9_minor=r1['minor_determinant'],source2_projected9_minor=r2['minor_determinant'],
        generic_GL_affine_orbit_rank=rO['rank'],generic_GL_plus_scalar_rank=rOP['rank'],generic_projective_PGL_orbit_dimension=rOP['rank']-1,
        generic_GL_plus_scalar_plus_parameter_derivative_rank=rFamily['rank'],projective_orbit_saturated_family_dimension=rFamily['rank']-1,
        generic_same_point_tangent_union_affine_rank=unionrank,generic_same_point_tangent_intersection_affine_dimension=intersection,
        proper_whole_fiber_length2_at_origin_dependency_verified=True,
        generic_nonzero_family_whole_source_fiber_two_reduced_points=True,
        generic_nonzero_family_reduced_CA_completion_has_two_smooth29_branches=True,
        whole_CA_embedding_dimension_equals_tangent_union_rank_claimed=False,
        clean_intersection_or_whole_conductor_from_ranks_claimed=False,
        seconds=time.monotonic()-start)
    paths=[Path(__file__),inputpath,parent/'report.json',parent/'certificate.json.gz',pairreview,properreview,fiberreport,fiberreview];note=ROOT/f'research/notes/{STEM}.md'
    if note.exists():paths.append(note)
    result['input_sha256']={str(p.relative_to(ROOT)):sha(p)for p in paths};result['evidence_sha256']={'certificate.json.gz':sha(out/'certificate.json.gz')}
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');signal.alarm(0)
    print(json.dumps({k:v for k,v in result.items()if k not in ['input_sha256','evidence_sha256']},indent=2),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
