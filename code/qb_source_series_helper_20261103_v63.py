#!/usr/bin/env python3
"""AUTHOR helper: exact point-QQ CA/CB source29 jets to physical150 jets (degree<=3).

Input/output arrays are ORDINARY power-series coefficients (no derivative factorials).
Input parameters are OFFSETS from the frozen Q_b source point; constant terms
must be zero. CB order: binary f3, Hamiltonian16, ratio mu, flag9.
CA order: fixed fiber20, flag9. Group action is exp(sum flag*s_i(t)).
This helper is an author contribution to the cubic probe, not its review.
"""
import argparse
from fractions import Fraction as Q
from functools import lru_cache
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import time
import sympy as s
from ambient_CA_V1_independent_review_20261103_v63 import tensor_factor_operator

ROOT=Path(__file__).resolve().parents[1]
STEM='QB_SOURCE_SERIES_HELPER_20261103_v63'
PAIRS=list(itertools.combinations(range(5),2));MONS=list(itertools.combinations_with_replacement(range(5),2))


def read(path):
    with(gzip.open(path,'rt')if path.suffix=='.gz'else path.open())as f:return json.load(f)
def frac(v):
    if isinstance(v,Q):return v
    v=s.Rational(v);return Q(int(v.p),int(v.q))
def add(*vs):return[sum((v[i]for v in vs),Q(0))for i in range(len(vs[0]))]
def scale(v,c):return[c*x for x in v]
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


@lru_cache(maxsize=1)
def inputs():
    basic=read(ROOT/'reports/AMBIENT_CB_ORIGIN_POLYNOMIAL_SELF_PAIR_20261103_v63/certificate.json.gz')
    cb=read(ROOT/'reports/AMBIENT_CB_GLOBAL_FLAG_MODEL_20261103_v63/report.json')
    ca=read(ROOT/'reports/AMBIENT_B18_CA_CRITICAL_SCHEME_20261103_v63/certificate.json.gz')
    return basic,cb,ca


def constant_map(g,gi):
    wedge=[[g[i][k]*g[j][l]-g[i][l]*g[j][k]for k,l in PAIRS]for i,j in PAIRS];lookup={m:i for i,m in enumerate(MONS)};sub=[]
    for m,n in MONS:
        out=[Q(0)]*15
        for i,j in itertools.product(range(5),repeat=2):out[lookup[tuple(sorted((i,j)))]]+=gi[m][i]*gi[n][j]
        sub.append(out)
    rows=[[]for _ in range(150)]
    for op,ip,old,new in itertools.product(range(10),range(10),range(15),range(15)):
        c=wedge[op][ip]*sub[old][new]
        if c:rows[15*op+new].append((15*ip+old,c))
    return lambda v:[sum((c*v[j]for j,c in row),Q(0))for row in rows]


def physical_push(g,gi,source):
    """Coefficient arrays; g[i][j][degree], source[raw150][degree]."""
    size=len(source[0]);zero=[Q(0)]*size
    def mul(v,w):
        out=zero[:]
        for i,a in enumerate(v):
            if a:
                for j,b in enumerate(w[:size-i]):
                    if b:out[i+j]+=a*b
        return out
    wedge=[[add(mul(g[i][k],g[j][l]),scale(mul(g[i][l],g[j][k]),-1))for k,l in PAIRS]for i,j in PAIRS]
    lookup={m:i for i,m in enumerate(MONS)};sub=[]
    for m,n in MONS:
        row=[zero[:]for _ in range(15)]
        for i,j in itertools.product(range(5),repeat=2):
            k=lookup[tuple(sorted((i,j)))];row[k]=add(row[k],mul(gi[m][i],gi[n][j]))
        sub.append(row)
    result=[zero[:]for _ in range(150)]
    for ip,old in itertools.product(range(10),range(15)):
        if not any(source[15*ip+old]):continue
        for op in range(10):
            if not any(wedge[op][ip]):continue
            v=mul(wedge[op][ip],source[15*ip+old])
            for new in range(15):
                if any(sub[old][new]):result[15*op+new]=add(result[15*op+new],mul(v,sub[old][new]))
    return result


class Source29:
    def __init__(self,kind,parameter,order):
        assert kind in['CB1','CB2','CA']and 0<=order<=3
        self.kind=kind;self.parameter=frac(parameter);self.order=order;self.size=order+1
        if not self.parameter:raise ValueError('CB source charts require a nonzero Q_b parameter')
        basic,cb,ca=inputs();self._basic=basic
        X=s.symbols('h x y z w');h,x,y,z,w=X;b=s.symbols('b');point=s.Rational(self.parameter.numerator,self.parameter.denominator)
        parse=lambda v:s.sympify(v,locals={'b':b,'h':h,'x':x,'y':y,'z':z,'w':w}).subs(b,point)
        if kind=='CA':
            self.basis=list(map(list,zip(*[[frac(c)for c in row]for row in ca['fixed_fiber_linear_basis150_by20']])));self.base=[frac(c)for c in ca['center_fixed20']];self.base[1]-=2*self.parameter;self.base[9]-=self.parameter
            self.flag_positions=list(map(tuple,ca['flag_positions']));self.parameter_order=[f'fiber{i}'for i in range(20)]
            G=[[Q(i==j)for j in range(5)]for i in range(5)];Gi=G
        else:
            enlarged=[[frac(c)for c in v]for v in cb['enlarged_canonical_basis']];self.basis=[*enlarged[:3],*enlarged[6:]];self.g_basis=enlarged[3:6]
            item=basic['actual_two_CB_sources'][int(kind[-1])-1];F=parse(item['f']);H=parse(item['H']);cubics=[s.prod((h,x,z,w)[i]for i in mon)for mon in itertools.combinations_with_replacement(range(4),3)if any(i>=2 for i in mon)]
            self.base=[frac(s.Poly(F,h,x).coeff_monomial(m))for m in[h*h,h*x,x*x]]+[frac(s.Poly(H,h,x,z,w).coeff_monomial(m))for m in cubics]+[Q(0)]
            self.flag_positions=[ij for ij in itertools.product(range(5),repeat=2)if list(ij)not in cb['parabolic_GL_positions']];self.parameter_order=['f_h2','f_hx','f_x2',*[f'H_{m}'for m in cubics],'mu']
            G=[[frac(parse(v))for v in row]for row in item['G']];Gi=[[frac(parse(v))for v in row]for row in item['G_inverse']]
        self.parameter_order.extend([f'flag_{i}_{j}'for i,j in self.flag_positions]);assert len(self.parameter_order)==29
        self.back=constant_map(Gi,G);self.ops=[tensor_factor_operator(*ij)for ij in self.flag_positions]
        # Every length4 path in this opposite-parabolic matrix pattern is absent: A^4=0.
        edges=set(self.flag_positions);paths={(i,j)for i,j in edges}
        for _ in range(3):paths={(i,l)for i,k in paths for j,l in edges if k==j}
        assert not paths
        self.raw_point=[sum((v[i]*c for v,c in zip(self.basis,self.base[:len(self.basis)])),Q(0))for i in range(150)];self.point150=self.back(self.raw_point)
        expected=[frac(parse(v))for v in basic['full_polynomial_Q150']];assert self.point150==expected

    def zero_jets(self):return[[Q(0)]*self.size for _ in range(29)]
    def _prepare(self,jets):
        if len(jets)!=29:raise ValueError('Expected29 source offset jets')
        jets=[[frac(c)for c in row]+[Q(0)]*(self.size-len(row))for row in jets]
        if any(len(row)!=self.size or row[0]for row in jets):raise ValueError('Offsets must have zero constant term and degree <= requested order')
        return jets
    def _mul(self,v,w):return[sum((v[j]*w[i-j]for j in range(i+1)),Q(0))for i in range(self.size)]
    def _source(self,jets):
        params=[row[:]for row in jets[:20]]
        for i,c in enumerate(self.base):params[i][0]+=c
        if self.kind=='CA':basis=self.basis;coeff=params
        else:basis=[*self.basis[:3],*self.g_basis,*self.basis[3:]];coeff=[*params[:3],*[self._mul(params[19],p)for p in params[:3]],*params[3:19]]
        return[[sum((v[i]*c[k]for v,c in zip(basis,coeff)),Q(0))for k in range(self.size)]for i in range(150)]
    def _flag_matrices(self,jets):
        n=self.size;zero=[Q(0)]*n;A=[[zero[:]for _ in range(5)]for _ in range(5)]
        for value,(i,j)in zip(jets[20:],self.flag_positions):A[i][j]=value[:]
        def mm(B,C):return[[add(*[self._mul(B[i][k],C[k][j])for k in range(5)])for j in range(5)]for i in range(5)]
        ident=[[[Q(i==j),*([Q(0)]*(n-1))]for j in range(5)]for i in range(5)];power=ident;E=[[v[:]for v in row]for row in ident];Ei=[[v[:]for v in row]for row in ident]
        factorial=1
        for k in range(1,4):
            power=mm(power,A);factorial*=k
            E=[[add(E[i][j],scale(power[i][j],Q(1,factorial)))for j in range(5)]for i in range(5)]
            Ei=[[add(Ei[i][j],scale(power[i][j],Q((-1)**k,factorial)))for j in range(5)]for i in range(5)]
        assert mm(E,Ei)==ident and mm(Ei,E)==ident
        return E,Ei
    def evaluate(self,jets):
        jets=self._prepare(jets);E,Ei=self._flag_matrices(jets);local=physical_push(E,Ei,self._source(jets));return list(map(list,zip(*[self.back(v)for v in zip(*local)])))
    def _lie(self,op,v):
        W,V=op;out=[Q(0)]*150
        for a,b,c in W:
            for mon in range(15):out[15*a+mon]+=c*v[15*b+mon]
        for old,new,c in V:
            for pair in range(10):out[15*pair+new]+=c*v[15*pair+old]
        return out
    def evaluate_lie(self,jets):
        """Alternative exact exp(L_sum) evaluator; useful for checks or fast Hensel steps."""
        jets=self._prepare(jets);term=list(map(list,zip(*self._source(jets))));out=[v[:]for v in term]
        for power in range(1,self.size):
            following=[[Q(0)]*150 for _ in range(self.size)]
            for op,cofs in zip(self.ops,jets[20:]):
                for source_degree,v in enumerate(term):
                    if any(v):
                        deriv=self._lie(op,v)
                        for delta in range(1,self.size-source_degree):
                            if cofs[delta]:following[source_degree+delta]=add(following[source_degree+delta],scale(deriv,cofs[delta]/power))
            term=following;out=[add(v,w)for v,w in zip(out,term)]
        return list(map(list,zip(*[self.back(v)for v in out])))
    def jacobian(self):
        columns=self.basis[:]if self.kind=='CA'else[*self.basis,add(*[scale(v,c)for v,c in zip(self.g_basis,self.base[:3])])]
        columns=columns+[(self._lie(op,self.raw_point))for op in self.ops]
        return list(map(list,zip(*[self.back(v)for v in columns])))


def make_qb_source(kind='CB1',parameter=Q(1),order=3):return Source29(kind,parameter,order)


def checks(out):
    start=time.monotonic();out.mkdir(parents=True,exist_ok=False);results={}
    for kind in['CB1','CB2','CA']:
        point=make_qb_source(kind,order=3);zero=point.evaluate(point.zero_jets());assert [v[0]for v in zero]==point.point150 and not any(c for v in zero for c in v[1:])
        linear=make_qb_source(kind,order=1);J=linear.jacobian()
        for j in range(29):
            jets=linear.zero_jets();jets[j][1]=1;value=linear.evaluate(jets);assert [v[1]for v in value]==[row[j]for row in J]
        dense=point.zero_jets()
        for i in range(29):
            for d in range(1,4):dense[i][d]=Q((i+2*d)%7-3,(i+d)%4+1)
        before=time.monotonic();physical=point.evaluate(dense);physical_seconds=time.monotonic()-before
        before=time.monotonic();lie=point.evaluate_lie(dense);lie_seconds=time.monotonic()-before;assert physical==lie
        second=make_qb_source(kind,order=2);value2=second.evaluate([row[:3]for row in dense]);assert value2==[row[:3]for row in physical]
        if kind!='CA':
            mixed=point.zero_jets();mixed[0][1]=1;mixed[19][1]=1;q=point.evaluate(mixed);assert [v[2]for v in q]==point.back(point.g_basis[0])and any(v[2]for v in q)
        results[kind]=dict(all29_first_derivatives_match=True,dense_all150_through3_physical_equals_Lie=True,order2_exactly_matches_order3_truncation=True,physical_seconds=physical_seconds,alternative_Lie_seconds=lie_seconds)
    report=dict(id=STEM,all_passed=True,role='AUTHOR helper; not independent review of any cubic probe',input_convention='29 OFFSET jets with zero constant terms',output_convention='150 coefficient arrays degree0..order; Fraction exact',source_kinds=results,seconds=time.monotonic()-start)
    paths=[Path(__file__),ROOT/'reports/AMBIENT_CB_ORIGIN_POLYNOMIAL_SELF_PAIR_20261103_v63/certificate.json.gz',ROOT/'reports/AMBIENT_CB_GLOBAL_FLAG_MODEL_20261103_v63/report.json',ROOT/'reports/AMBIENT_B18_CA_CRITICAL_SCHEME_20261103_v63/certificate.json.gz',ROOT/'code/ambient_CA_V1_independent_review_20261103_v63.py']
    note=ROOT/f'research/notes/{STEM}.md'
    if note.exists():paths.append(note)
    report['input_sha256']={str(p.relative_to(ROOT)):sha(p)for p in paths};(out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);checks(parser.parse_args().output)
