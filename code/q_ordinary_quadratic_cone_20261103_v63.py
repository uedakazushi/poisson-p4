#!/usr/bin/env python3
"""Exact ordinary two-jet of all350 residuals, then bounded GF(p) cone GB.

The zero-weight parameter is a=2/9+u.  All 35 local variables have ordinary
degree one here.  Rational inverses are expanded through degree two at their
nonzero constant term; this does not expand the 18-million-term weighted model.
"""
import argparse
from fractions import Fraction
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import resource
import subprocess
import time

from q_weighted_finite_model_validation_20261103_v63 import evaluate,physical_equations

ROOT=Path(__file__).resolve().parents[1]
PRIME=32003


class TruncatedQQ:
    zero={}
    one={():Fraction(1)}
    inverse_constants=[]

    @staticmethod
    def q(value):
        value=Fraction(value)
        return {():value} if value else {}

    @staticmethod
    def add(left,right):
        out=left.copy()
        for mon,c in right.items():
            value=out.get(mon,Fraction(0))+c
            if value:out[mon]=value
            else:out.pop(mon,None)
        return out

    @staticmethod
    def neg(value):return {mon:-c for mon,c in value.items()}

    @staticmethod
    def mul(left,right):
        if not left or not right:return {}
        a=[[],[],[]];b=[[],[],[]]
        for mon,c in left.items():a[len(mon)].append((mon,c))
        for mon,c in right.items():b[len(mon)].append((mon,c))
        out={}
        for i in range(3):
            for j in range(3-i):
                for m,c in a[i]:
                    for n,d in b[j]:
                        mon=tuple(sorted(m+n));value=out.get(mon,Fraction(0))+c*d
                        if value:out[mon]=value
                        else:out.pop(mon,None)
        return out

    @classmethod
    def inv(cls,value):
        constant=value.get((),Fraction(0));assert constant
        assert constant.numerator%PRIME and constant.denominator%PRIME
        cls.inverse_constants.append(str(constant))
        linear={mon:c for mon,c in value.items() if len(mon)==1}
        quadratic={mon:c for mon,c in value.items() if len(mon)==2}
        out=cls.q(1/constant)
        out=cls.add(out,{mon:-c/constant**2 for mon,c in linear.items()})
        out=cls.add(out,{mon:-c/constant**2 for mon,c in quadratic.items()})
        out=cls.add(out,{mon:c/constant**3 for mon,c in cls.mul(linear,linear).items()})
        assert cls.mul(value,out)==cls.one
        return out


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def fp(value):
    assert value.denominator%PRIME
    return value.numerator*pow(value.denominator,-1,PRIME)%PRIME


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--model',type=Path,default=ROOT/'reports/Q_WEIGHTED_FINITE_MODEL_20261103_v63/model.json.gz')
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--timeout',type=int,default=300)
    parser.add_argument('--prepare-only',action='store_true')
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=False);start=time.monotonic()
    model=json.load(gzip.open(args.model,'rt'))
    inputs=[{(i,):Fraction(1)} for i in range(35)]
    inputs[0][()]=Fraction(2,9)
    values,inverses=evaluate(model['nodes'],inputs,TruncatedQQ)
    reconstructed=[values[node] for node in model['raw150_reconstruction_nodes']]
    residuals=[values[node] for node in model['all350_residual_nodes']]
    assert all(len(mon)==2 for row in residuals for mon in row), 'All constant/linear terms must vanish'
    _,physical=physical_equations(reconstructed,TruncatedQQ);assert physical==residuals
    monomials=list(itertools.combinations_with_replacement(range(35),2))
    monindex={mon:i for i,mon in enumerate(monomials)}
    # Sparse exact QQ row elimination.  Every pivot used is a p-unit; hence the
    # same selected original rows generate all350 quadrics over Z_(p), as well.
    basis={};selected=[];pivots=[];scales=[]
    for rowindex in sorted(range(350),key=lambda i:(len(residuals[i]),i)):
        row={monindex[mon]:c for mon,c in residuals[rowindex].items()}
        while row:
            pivot=min(row);coefficient=row[pivot]
            if pivot not in basis:
                assert coefficient.numerator%PRIME and coefficient.denominator%PRIME
                basis[pivot]={j:c/coefficient for j,c in row.items()}
                selected.append(rowindex);pivots.append(pivot);scales.append(str(coefficient));break
            for j,c in basis[pivot].items():
                value=row.get(j,Fraction(0))-coefficient*c
                if value:row[j]=value
                else:row.pop(j,None)
    for original in residuals:
        row={monindex[mon]:c for mon,c in original.items()}
        while row:
            pivot=min(row);assert pivot in basis;coefficient=row[pivot]
            for j,c in basis[pivot].items():
                value=row.get(j,Fraction(0))-coefficient*c
                if value:row[j]=value
                else:row.pop(j,None)
    assert all(c.denominator%PRIME for row in residuals for c in row.values())
    variables=['u']+[f'z{i}' for i in range(34)]
    artifact=dict(variables=variables,coefficient_ring='QQ',truncation='ordinary total degree <=2 at a=2/9,z=0',
                  all350_quadratic_rows=[[[*mon,str(c)] for mon,c in sorted(row.items())] for row in residuals],
                  selected_original_rows=selected,pivot_monomial_indices=pivots,
                  exact_QQ_elimination_pivot_values=scales,
                  normalized_row_basis=[dict(pivot=pivot,terms=[[j,str(c)] for j,c in sorted(row.items())])
                                        for pivot,row in sorted(basis.items())])
    encoded=(json.dumps(artifact,separators=(',',':'))+'\n').encode()
    archive=args.output/'quadrics.json.gz'
    with archive.open('wb') as stream:
        with gzip.GzipFile(fileobj=stream,mode='wb',filename='',mtime=0) as z:z.write(encoded)
    summary=(args.output/'summary.txt').resolve();script=args.output/'compute.m2'
    text=f'R=GF({PRIME})['+','.join(variables)+',MonomialOrder=>GRevLex];\n'
    for i in selected:
        terms=[]
        for mon,c in sorted(residuals[i].items()):
            terms.append(str(fp(c))+'*'+'*'.join(variables[j] for j in mon))
        text+=f'q{i}='+ '+'.join(terms)+';\n'
    text+='I=ideal{'+','.join(f'q{i}' for i in selected)+'};\n'
    text+='gbTrace=1; G=gb I; B=gens G;\n'
    text+='f=openOut "'+str(summary)+'";\n'
    text+='f << "DIMENSION " << dim(R/I) << endl;\n'
    text+='f << "GB_GENERATORS " << numgens source B << endl;\n'
    text+='f << "HILBERT " << toString hilbertSeries(R/I,Reduce=>true) << endl;\n'
    text+='f << "LEADING_MONOMIALS " << toString flatten entries leadTerm B << endl;\n'
    text+='close f;\nexit 0;\n';script.write_text(text)
    extraction=time.monotonic()-start
    print('Extracted all350 exact QQ quadrics: rank',len(selected),'nonzero',sum(bool(row) for row in residuals),
          'seconds',extraction,flush=True)
    report=dict(id='Q_ORDINARY_QUADRATIC_CONE_20261103_V63',extraction_completed=True,
                ordinary_local_variables=35,variables=variables,coefficient_basepoint='a=2/9,z=0',
                full350_constant_and_linear_terms_zero=True,
                all350_two_jets_equal_independent_physical_Jacobi=True,
                physical_divergence25_zero=True,
                quadratic_span_rank_QQ=len(selected),nonzero_quadratic_rows=sum(bool(row) for row in residuals),
                full_model_weighted_polynomials_expanded=False,
                truncation_degree=2,all_inverse_Taylor_products_verified_to_degree2=True,
                all_inverse_constants_p_units=True,inverse_nodes_evaluated=inverses,
                all350_quadratic_coefficients_p_integral=True,
                selected_rows_generate_all350_over_Z_localized_p=True,
                prime=PRIME,extraction_seconds=extraction,
                quadratic_archive_sha256=sha(archive),script_sha256=sha(script),
                input_sha256={str(p.resolve().relative_to(ROOT)):sha(p) for p in
                              [Path(__file__),args.model,ROOT/'code/q_weighted_finite_model_validation_20261103_v63.py']},
                full_model_exact_dimension_claimed=False,
                higher_order_equations_or_branch_decomposition_completed=False)
    (args.output/'extraction.json').write_text(json.dumps(report,indent=2)+'\n')
    def limits():
        resource.setrlimit(resource.RLIMIT_AS,(6000000000,6000000000))
        resource.setrlimit(resource.RLIMIT_CPU,(args.timeout,args.timeout+1))
    status='prepared' if args.prepare_only else 'completed'
    with (args.output/'m2.log').open('w') as stream:
        if not args.prepare_only:
            try:
                result=subprocess.run(['M2','--script',str(script)],stdout=stream,stderr=subprocess.STDOUT,
                                      timeout=args.timeout,preexec_fn=limits)
                if result.returncode:status=f'exit_{result.returncode}'
            except subprocess.TimeoutExpired:status='timeout'
    resulttext=summary.read_text() if summary.exists() else ''
    completed=status=='completed' and 'DIMENSION ' in resulttext
    report.update(cone_GB_completed=completed,status=status,
                  memory_limit_bytes=6000000000,wall_timeout_seconds=args.timeout)
    if completed:
        report['finite_field_quadratic_cone_dimension']=int(next(line.split()[1] for line in resulttext.splitlines() if line.startswith('DIMENSION ')))
        report['cone_GB_generators']=int(next(line.split()[1] for line in resulttext.splitlines() if line.startswith('GB_GENERATORS ')))
        report['cone_hilbert_series']=next(line[len('HILBERT '):] for line in resulttext.splitlines() if line.startswith('HILBERT '))
        report['summary_sha256']=sha(summary)
        report['dimension_transfer_requires_proper_projectivization_argument']=True
    report['seconds']=time.monotonic()-start
    (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ('input_sha256','variables')},indent=2))


if __name__=='__main__':main()
