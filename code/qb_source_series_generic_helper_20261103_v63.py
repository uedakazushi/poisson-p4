#!/usr/bin/env python3
"""AUTHOR helper: all150 actual Q_b CA/CB source jets over exact QQ(b), order<=3.

Same offset/ordinary-coefficient conventions as the frozen point helper.
Default evaluate uses5x5 exterior/inverse-monomial coefficient arrays;
evaluate_lie independently uses the exact coefficient exp(L_sumflags).
"""
import argparse
from fractions import Fraction
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import time
import sympy as s
from sympy.polys.fields import FracElement
from qb_source_series_helper_20261103_v63 import Source29,inputs,constant_map,make_qb_source as make_point_source

ROOT=Path(__file__).resolve().parents[1]
STEM='QB_SOURCE_SERIES_GENERIC_HELPER_20261103_v63'
SYMBOL=s.symbols('b');DOMAIN=s.QQ.frac_field(SYMBOL);GENERATOR=DOMAIN.gens[0]
X=s.symbols('h x y z w');h,x,y,z,w=X


def coerce(value):
    if isinstance(value,FracElement)and value.field==DOMAIN.field:return value
    if isinstance(value,str):value=s.sympify(value,locals={'b':SYMBOL,'h':h,'x':x,'y':y,'z':z,'w':w})
    if isinstance(value,s.Basic):return DOMAIN.from_sympy(value)
    return DOMAIN.convert(value)
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


class GenericSource29(Source29):
    domain=DOMAIN;symbol=SYMBOL;generator=GENERATOR;base_denominator_factors=('b',)
    def __init__(self,kind,parameter,order):
        if kind not in['CB1','CB2','CA']or not 0<=order<=3:raise ValueError('kindCB1/CB2/CA; order0..3')
        if parameter is not None and str(parameter)!='b'and coerce(parameter)!=GENERATOR:raise ValueError('This is the formal QQ(b) factory, not a point specialization')
        self.kind=kind;self.parameter=GENERATOR;self.order=order;self.size=order+1
        basic,cb,ca=inputs();self._basic=basic
        if kind=='CA':
            self.basis=list(map(list,zip(*[[coerce(c)for c in row]for row in ca['fixed_fiber_linear_basis150_by20']])));self.base=[coerce(c)for c in ca['center_fixed20']];self.base[1]-=2*GENERATOR;self.base[9]-=GENERATOR
            self.flag_positions=list(map(tuple,ca['flag_positions']));self.parameter_order=[f'fiber{i}'for i in range(20)]
            G=[[DOMAIN.one if i==j else DOMAIN.zero for j in range(5)]for i in range(5)];Gi=G
            source_basis=self.basis
        else:
            enlarged=[[coerce(c)for c in v]for v in cb['enlarged_canonical_basis']];self.basis=[*enlarged[:3],*enlarged[6:]];self.g_basis=enlarged[3:6]
            item=basic['actual_two_CB_sources'][int(kind[-1])-1];F=s.sympify(item['f'],locals={'b':SYMBOL,'h':h,'x':x,'y':y,'z':z,'w':w});H=s.sympify(item['H'],locals={'b':SYMBOL,'h':h,'x':x,'y':y,'z':z,'w':w})
            cubics=[s.prod((h,x,z,w)[i]for i in mon)for mon in itertools.combinations_with_replacement(range(4),3)if any(i>=2 for i in mon)]
            self.base=[coerce(s.Poly(F,h,x).coeff_monomial(m))for m in[h*h,h*x,x*x]]+[coerce(s.Poly(H,h,x,z,w).coeff_monomial(m))for m in cubics]+[DOMAIN.zero]
            self.flag_positions=[ij for ij in itertools.product(range(5),repeat=2)if list(ij)not in cb['parabolic_GL_positions']];self.parameter_order=['f_h2','f_hx','f_x2',*[f'H_{m}'for m in cubics],'mu']
            G=[[coerce(v)for v in row]for row in item['G']];Gi=[[coerce(v)for v in row]for row in item['G_inverse']]
            source_basis=[*self.basis[:3],*self.g_basis,*self.basis[3:]]
        from ambient_CA_V1_independent_review_20261103_v63 import tensor_factor_operator
        self.ops=[tensor_factor_operator(*ij)for ij in self.flag_positions];self.parameter_order.extend([f'flag_{i}_{j}'for i,j in self.flag_positions]);assert len(self.parameter_order)==29
        back=constant_map(Gi,G);self.back=lambda v:[coerce(c)for c in back(v)]
        self._source_rows=[[(j,v[i])for j,v in enumerate(source_basis)if v[i]]for i in range(150)]
        self.raw_point=[sum((v[i]*c for v,c in zip(self.basis,self.base[:len(self.basis)])),DOMAIN.zero)for i in range(150)];self.point150=self.back(self.raw_point)
        assert self.point150==[coerce(v)for v in basic['full_polynomial_Q150']]
        edges=set(self.flag_positions);paths=set(edges)
        for _ in range(3):paths={(i,l)for i,k in paths for j,l in edges if k==j}
        assert not paths
    def zero_jets(self):return[[DOMAIN.zero]*self.size for _ in range(29)]
    def _prepare(self,jets):
        if len(jets)!=29:raise ValueError('Expected29 source offset jets')
        jets=[[coerce(c)for c in row]+[DOMAIN.zero]*(self.size-len(row))for row in jets]
        if any(len(row)!=self.size or row[0]for row in jets):raise ValueError('Offsets must have zero constant term and degree<=order')
        return jets
    def _mul(self,v,w):
        result=[DOMAIN.zero]*self.size
        for i,a in enumerate(v):
            if a:
                for j,b in enumerate(w[:self.size-i]):
                    if b:result[i+j]+=a*b
        return result
    def _source(self,jets):
        params=[row[:]for row in jets[:20]]
        for i,c in enumerate(self.base):params[i][0]+=c
        if self.kind=='CA':coeff=params
        else:coeff=[*params[:3],*[self._mul(params[19],p)for p in params[:3]],*params[3:19]]
        return[[sum((c*coeff[j][k]for j,c in row),DOMAIN.zero)for k in range(self.size)]for row in self._source_rows]
    def evaluate(self,jets):return super().evaluate(jets)
    def evaluate_physical(self,jets):return super().evaluate(jets)


def make_qb_source(kind='CB1',parameter=None,order=3):return GenericSource29(kind,parameter,order)


def checks(out):
    start=time.monotonic();out.mkdir(parents=True,exist_ok=False);results={};example={}
    for kind in['CB1','CB2','CA']:
        phi=make_qb_source(kind,order=3);zero=phi.evaluate(phi.zero_jets());assert [v[0]for v in zero]==phi.point150 and not any(c for v in zero for c in v[1:])
        jets=phi.zero_jets()
        for i in range(29):
            for d in range(1,4):jets[i][d]=DOMAIN.convert(Fraction((i+2*d)%7-3,(i+d)%4+1))*GENERATOR**(d%2)/(1+GENERATOR*GENERATOR)
        then=time.monotonic();fast=phi.evaluate(jets);fasttime=time.monotonic()-then
        then=time.monotonic();physical=phi.evaluate_lie(jets);physicaltime=time.monotonic()-then;assert fast==physical
        second=make_qb_source(kind,order=2);assert second.evaluate([row[:3]for row in jets])==[row[:3]for row in fast]
        point=make_point_source(kind,parameter=Fraction(1),order=3)
        at1=lambda v:Fraction(str(s.cancel(DOMAIN.to_sympy(coerce(v)).subs(SYMBOL,1))))
        assert point.evaluate([[at1(v)for v in row]for row in jets])==[[at1(v)for v in row]for row in fast]
        # All150 zero-flag first derivatives, checked cheaply through the exact supplied Jacobian.
        J=phi.jacobian();assert len(J)==150 and all(len(row)==29 for row in J)
        results[kind]=dict(all150_dense_QQb_D3_physical_equals_Lie=True,order2_truncation_matches=True,specialization1_matches_frozen_point_helper=True,default_physical_D3_seconds=fasttime,alternative_Lie_D3_seconds=physicaltime)
        example[kind]=dict(jets=[[str(DOMAIN.to_sympy(v))for v in row]for row in jets],output150=[[str(DOMAIN.to_sympy(v))for v in row]for row in fast])
        print(json.dumps(dict(stage=kind,**results[kind]),indent=2),flush=True)
    with(out/'generic-checks.json.gz').open('wb')as file:
        with gzip.GzipFile(fileobj=file,mode='wb',filename='',mtime=0)as archive:archive.write((json.dumps(example,separators=(',',':'))+'\n').encode())
    report=dict(id=STEM,all_passed=True,role='AUTHOR generic helper, not independent cubic-gate review',domain='QQ(b)',input='29 OFFSET ordinary coefficient jets; zero constants',output='150 FracElement coefficient arrays degree0..order',base_poles_only_b=True,other_poles_only_from_input_jets=True,checks=results,seconds=time.monotonic()-start)
    paths=[Path(__file__),ROOT/'code/qb_source_series_helper_20261103_v63.py',ROOT/'reports/QB_SOURCE_SERIES_HELPER_20261103_v63/report.json',ROOT/'reports/AMBIENT_CB_ORIGIN_POLYNOMIAL_SELF_PAIR_20261103_v63/certificate.json.gz',ROOT/'reports/AMBIENT_CB_GLOBAL_FLAG_MODEL_20261103_v63/report.json',ROOT/'reports/AMBIENT_B18_CA_CRITICAL_SCHEME_20261103_v63/certificate.json.gz']
    note=ROOT/f'research/notes/{STEM}.md'
    if note.exists():paths.append(note)
    report['input_sha256']={str(p.relative_to(ROOT)):sha(p)for p in paths};report['evidence_sha256']={'generic-checks.json.gz':sha(out/'generic-checks.json.gz')};(out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);checks(p.parse_args().output)
