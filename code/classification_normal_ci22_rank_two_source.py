#!/usr/bin/env python3
"""Whole normal CI22 tangent sources: all three residual Jordan surface orbits."""
import argparse,hashlib,itertools as it,json,resource,signal,time
from pathlib import Path
import sympy as sp
from classification_cubic_scroll_abelian_review import Q,PAIRS,add,scale,mul,diff,parse,rref,monomials,jacobi_and_divisor,determinant
from classification_quadric_cone_zero_source import empty,kernel
ROOT=Path(__file__).resolve().parents[1]

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',required=True,type=Path);args=ap.parse_args();begun=time.monotonic()
    signal.signal(signal.SIGALRM,lambda s,f:(_ for _ in ()).throw(TimeoutError('three whole normalCI22 sources60s')));signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_CPU,(55,58));resource.setrlimit(resource.RLIMIT_AS,(768*1024**2,)*2)
    names=[f'x{i}'for i in range(5)];xs=[parse(v,names)for v in names];mons2=monomials(5,2);pars=sp.symbols('u0:5 a b c d');models=[('sem221','x0**2+x1**2+x2**2+x3**2+x4**2','x2**2+x3**2+2*x4**2'),('mixed21_11','2*x0*x1+x2**2+x3**2+x4**2','x1**2+x3**2+x4**2'),('nil41','2*x0*x3+2*x1*x2+x4**2','2*x1*x3+x2**2')];out=[]
    def encode(P):return [[list(e),str(v)]for e,v in sorted(P.items(),reverse=True)]
    def tensor(P):return {str((i,j)):encode(P[i][j])for i,j in PAIRS}
    def vec(P):return [P[i][j].get(e,Q(0))for i,j in PAIRS for e in mons2]
    for model,q1text,q2text in models:
        q1,q2=parse(q1text,names),parse(q2text,names);J=[]
        for missing in range(5):
            U=[i for i in range(5)if i!=missing];P=empty()
            for i,j in it.combinations(U,2):
                k,l=[v for v in U if v not in[i,j]];perm=[U.index(v)for v in[i,j,k,l]];sign=(-1)**sum(perm[a]>perm[b]for a in range(4)for b in range(a+1,4))
                P[i][j]=scale(add(mul(diff(q1,k),diff(q2,l)),scale(mul(diff(q1,l),diff(q2,k)),-1)),sign);P[j][i]=scale(P[i][j],-1)
            assert all(not add(*(mul(P[i][j],diff(q,j))for j in range(5)))for q in[q1,q2]for i in range(5));J.append(P)
        raw=[]
        for q in[q1,q2]:
            for i,j in PAIRS:
                P=empty();P[i][j]=q;P[j][i]=scale(q,-1);raw.append(P)
        assert len(rref([vec(P)for P in J+raw]))==25
        vertical=[]
        for a,b in it.product(range(5),repeat=2):
            P=empty()
            for i,j in PAIRS:P[i][j]=add(mul(xs[i],xs[b])if j==a else{},scale(mul(xs[j],xs[b]),-1)if i==a else{})
            vertical.append(vec(P))
        assert len(rref(vertical))==24 and len(rref(vertical+[vec(P)for P in J+raw]))==49
        h=[add(mul(diff(q1,i),diff(q2,j)),scale(mul(diff(q1,j),diff(q2,i)),-1))for i,j in PAIRS]
        rows=[[p.get(e,Q(0))for p in h]+[-q1.get(e,Q(0)),-q2.get(e,Q(0))]for e in mons2];R,free,K=kernel(rows,12);Cbasis=[v[:10]for v in K];assert len(Cbasis)==2 and len(rref(Cbasis))==2
        family=list(J)
        for C in Cbasis:
            for q in[q1,q2]:
                P=empty()
                for value,(i,j)in zip(C,PAIRS):P[i][j]=scale(q,value);P[j][i]=scale(q,-value)
                family.append(P)
        JJ,DD=jacobi_and_divisor(family);pairs=list(it.combinations_with_replacement(range(9),2));assert all(not any(i<5 and j<5 for col in row for i,j in[pairs[col]])for row in JJ+DD)
        RJ=rref(JJ);basis=[[v.get(j,Q(0))for j in range(45)]for v in RJ.values()];polys=[sp.factor(sum(sp.Rational(v.numerator,v.denominator)*pars[i]*pars[j]for v,(i,j)in zip(row,pairs)))for row in basis]
        out.append(dict(model=model,q1=q1text,q2=q2text,whole_tangent_raw25_rank=25,Euler24_combined_rank49=True,full_conormal_matrix=[list(map(str,v))for v in rows],conormal_augmented_nullspace=[list(map(str,v))for v in K],full_two_constant_skew_generators=[list(map(str,v))for v in Cbasis],full9_parameter_order=list(map(str,pars)),all150_source_tensors=[tensor(P)for P in family],all5_Jacobian_Casimirs_verified=True,all_Jacobian_linear_combination_projectivePoisson_and_rank2_verified=True,full350_Jacobi_rows=[{str(j):str(v)for j,v in row.items()}for row in JJ],full126_actualPf_rows=[{str(j):str(v)for j,v in row.items()}for row in DD],full_parameter_pair_order=pairs,Jacobi_coefficient_rowspace_rank=len(RJ),complete_equivalent_quadratic_rowspace_basis=[list(map(str,v))for v in basis],factored_equivalent_Jacobi_basis=list(map(str,polys))))
    paths=[Path(__file__),ROOT/'code/classification_cubic_scroll_abelian_review.py',ROOT/'code/classification_quadric_cone_zero_source.py',ROOT/'research/classification/NORMAL_CI22_RANK_TWO_TANGENT_REDUCTION.md',ROOT/'research/classification/REGULAR_PENCIL_CI22_ABELIAN_CONORMAL_NONDOMINANCE.md',ROOT/'research/classification/ALL_DEGENERATE_CI22_ABELIAN_CONORMAL_EXCLUSION.md']
    result=dict(task_id='G',all_arithmetic_passed=True,cases=out,scope='Full normalCI22 tangent/conormal source reduction and complete original quadratic equations. No wholeG dimension or irreducibleD exclusion inferred from diagnostic rowspace alone.',all_non_normal_CI22_excluded_from_scope=True,input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in paths},producer_source_bytes=Path(__file__).read_text(),seconds=round(time.monotonic()-begun,6),maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    with args.output.open('x')as stream:json.dump(result,stream,indent=2);stream.write('\n')
    print(json.dumps(dict(all_arithmetic_passed=True,cases=[dict(model=v['model'],rank=v['Jacobi_coefficient_rowspace_rank'],equations=v['factored_equivalent_Jacobi_basis'])for v in out],seconds=result['seconds'])))
def factor_proof_main():
    import ast
    ap=argparse.ArgumentParser();ap.add_argument('--factor-proof',action='store_true');ap.add_argument('--input',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args();begun=time.monotonic();src=json.loads(args.input.read_text());names=['u0','u1','u2','u3','u4','a','b','c','d'];parsyms=sp.symbols(' '.join(names));mons4=monomials(5,4);mons5=monomials(5,5)
    frozen={v.name:ast.dump(v,include_attributes=False)for v in ast.parse(src['producer_source_bytes']).body if isinstance(v,ast.FunctionDef)};now={v.name:ast.dump(v,include_attributes=False)for v in ast.parse(Path(__file__).read_text()).body if isinstance(v,ast.FunctionDef)};assert all(now[k]==v for k,v in frozen.items())
    needed=[['a*c','b*c','a*u0','a*u1','d*u2','d*u3','d*(a+b)'],['a*c','b*c','a*u0','a*u2','d*u3','d*u4','d*(a+b)-2*b*u2'],['c**2','2*a*d-5*b*c']];out=[]
    for case,wanted in zip(src['cases'],needed):
        pairs=case['full_parameter_pair_order'];JJ=[{int(j):Q(v)for j,v in row.items()}for row in case['full350_Jacobi_rows']];nonzero=[i for i,v in enumerate(JJ)if v];rows=[[JJ[i].get(j,Q(0))for j in range(45)]+[Q(i==j)for j in nonzero]for i in nonzero];R=rref(rows);certs=[]
        for expr in wanted:
            p=parse(expr,names);v=[p.get(tuple(int(k==i)+int(k==j)for k in range(9)),Q(0))for i,j in pairs];remaining=list(v);coeff=[Q(0)]*len(nonzero)
            for pivot,row in R.items():
                if pivot>=45:continue
                value=remaining[pivot]
                if not value:continue
                for j in range(45):remaining[j]-=value*row.get(j,Q(0))
                for j in range(len(nonzero)):coeff[j]+=value*row.get(45+j,Q(0))
            assert not any(remaining)and all(sum(c*JJ[i].get(j,Q(0))for c,i in zip(coeff,nonzero))==v[j]for j in range(45))
            certs.append(dict(required_equation=expr,full_original350_functional={str(i):str(c)for i,c in zip(nonzero,coeff)if c},literal_full_coefficient_identity=True))
        item=dict(model=case['model'],required_equation_functionals=certs)
        if case['model']!='nil41':
            keep=4 if case['model']=='sem221'else 1;u,a,d=sp.symbols('U A D');sub={v:sp.Integer(0)for v in parsyms};sub.update({parsyms[keep]:u,parsyms[5]:a,parsyms[6]:-a,parsyms[8]:d})
            weights=[sub[v]for v in parsyms]
            def ev(row):return sp.expand(sum(sp.Rational(Q(v).numerator,Q(v).denominator)*weights[pairs[int(j)][0]]*weights[pairs[int(j)][1]]for j,v in row.items()))
            assert all(ev(row)==0 for row in case['full350_Jacobi_rows'])
            pf=[ev(row)for row in case['full126_actualPf_rows']];assert all(not value for e,value in zip(mons5,pf)if not e[keep])
            item.update(rank2_coefficient_open_complete_map={str(k):str(v)for k,v in sub.items()},all_original350_zero_verified=True,actual_D_fixed_linear_factor='x'+str(keep),full126_actualPf_after_complete_map=list(map(str,pf)))
        else:item.update(reduced_support_logic='c²=0 implies c=0; 2ad-5bc=0 then det(ad-bc)=0. No determinant localization retained on Poisson support.')
        out.append(item)
    paths=[Path(__file__),args.input.resolve(),ROOT/'research/classification/NORMAL_CI22_RANK_TWO_EXCLUSION.md',ROOT/'research/classification/NORMAL_CI22_RANK_TWO_TANGENT_REDUCTION.md']
    result=dict(task_id='G',all_arithmetic_passed=True,cases=out,whole_normal_CI22_rank2_irreducible_quintic_empty=True,geometry_scope='Normal CI22 only; k<=1 quadraticfactor, allregular residual3surfaceorbits, smoothellipticquarticcone allmoduli. AllR retained. Non-normalCI22/projectedscroll explicitcomplements.',all_old_function_AST_unchanged=True,input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in paths},producer_source_bytes=Path(__file__).read_text(),seconds=round(time.monotonic()-begun,6))
    with args.output.open('x')as stream:json.dump(result,stream,indent=2);stream.write('\n')
    print(json.dumps(dict(all_arithmetic_passed=True,all_normal_CI22_rank2_empty=True,cases=[v['model']for v in out],seconds=result['seconds'])))

if __name__=='__main__':
    import sys
    factor_proof_main()if '--factor-proof'in sys.argv else main()
