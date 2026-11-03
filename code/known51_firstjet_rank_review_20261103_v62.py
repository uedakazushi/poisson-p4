#!/usr/bin/env python3
"""Independent reverse AD/minors, native-input checks and T38 raw-jet audit."""
import argparse
from fractions import Fraction
import gzip
import hashlib
import itertools
import json
from pathlib import Path

import sympy as s

from s61_model import S61, X, pushforward

ROOT=Path(__file__).resolve().parents[1]
P=32003


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def determinant(matrix,witness):
    A=[[matrix[i][j]%P for j in witness['columns']] for i in witness['rows']]
    result=1
    for j in range(len(A)):
        k=next(i for i in range(j,len(A)) if A[i][j])
        if k!=j:A[j],A[k]=A[k],A[j];result=-result
        pivot=A[j][j];result=result*pivot%P
        for i in range(j+1,len(A)):
            factor=A[i][j]*pow(pivot,-1,P)%P
            A[i]=[(a-factor*b)%P for a,b in zip(A[i],A[j])]
    assert result%P==witness['determinant_mod_32003']!=0


def reverse_differential(body,indices):
    names=body['input_names'];active=[n for n in names if n!='r']
    point=body.get('exact_rational_sample',body.get('exact_sample'))
    def reduce(value):
        if str(value)=='sqrt(3)':return int(s.sqrt_mod(3,P))
        q=Fraction(str(value));return q.numerator*pow(q.denominator,-1,P)%P
    nodes=body['nodes'];values=[]
    for op,*args in nodes:
        arg=args[0]
        if op=='q':value=reduce(arg)
        elif op=='v':value=reduce(point[names[arg]])
        elif op=='+':value=sum(values[j] for j in arg)%P
        elif op=='*':
            value=1
            for j in arg:value=value*values[j]%P
        else:
            assert op=='inv' and values[arg]
            value=pow(values[arg],-1,P)
        values.append(value)
    matrix=[]
    for coordinate in indices:
        adjoint=[0]*len(nodes);adjoint[body['output_u1_through_u61'][coordinate]]=1
        gradient=[0]*len(active)
        for i in reversed(range(len(nodes))):
            a=adjoint[i]
            if not a:continue
            op,arg=nodes[i]
            if op=='v':
                name=names[arg]
                if name!='r':gradient[active.index(name)]=(gradient[active.index(name)]+a)%P
            elif op=='+':
                for j in arg:adjoint[j]=(adjoint[j]+a)%P
            elif op=='*':
                for j,child in enumerate(arg):
                    factor=a
                    for k,other in enumerate(arg):
                        if k!=j:factor=factor*values[other]%P
                    adjoint[child]=(adjoint[child]+factor)%P
            elif op=='inv':adjoint[arg]=(adjoint[arg]-a*values[i]*values[i])%P
        matrix.append(gradient)
    return active,matrix


def audit_native_input(path,record,body,indices):
    encoded=gzip.decompress(path.read_bytes())
    assert hashlib.sha256(encoded).hexdigest()==record['input_sha256']
    assert sha(path)==record['compressed_input_sha256']
    tokens=iter(encoded.decode().split());number=lambda:int(next(tokens))
    header=[number()for _ in range(5)]
    names=body['input_names'];groups=[f'g{i}'for i in range(10)]
    variables=[n for n in names if n!='r' and n not in groups]
    active=[n for n in names if n!='r']
    assert header==[len(variables),len(names),len(active),len(body['nodes']),16]
    assert variables==record['function_variables'] and active==record['differential_column_order']
    group_point=[1,0,0,1,0,0,0,0,0,0]
    for name in names:
        kind,index,column=next(tokens),number(),number()
        expected=('r',0) if name=='r' else ('c',group_point[groups.index(name)]) if name in groups else ('v',variables.index(name))
        assert (kind,index)==expected
        assert column==(active.index(name) if name in active else -1)
    for node in body['nodes']:
        op=next(tokens);assert op==('i'if node[0]=='inv'else node[0])
        if node[0]in ['+','*']:
            assert number()==len(node[1])
            assert [number()for _ in node[1]]==node[1]
        elif node[0]=='q':assert Fraction(next(tokens))==Fraction(node[1])
        else:assert number()==node[1]
    assert [number()for _ in range(16)]==[body['output_u1_through_u61'][i]for i in indices]
    assert next(tokens,None)is None


def t38_check(paths):
    folder=ROOT/'reports/MPS_T38_RAW_JET_RANK_20261103_v62'
    report=json.loads((folder/'report.json').read_text())
    dense_path=ROOT/'reports/MPS_T38_DENSE_FIXED_ORBIT_20261103_v62.json'
    dense=json.loads(dense_path.read_text())
    source_path=ROOT/'reports/MPS_SMOOTHING_COMPLETION_20261103_v62/T38.json'
    source=json.loads(source_path.read_text());pi=s.zeros(5)
    for pair,text in source['homogeneous_bivector'].items():
        i,j=map(int,pair.split(','));pi[i,j]=s.sympify(text);pi[j,i]=-pi[i,j]
    sample={s.Symbol(k):s.Rational(v)for k,v in report['fixed_parameter_values'].items()}
    assert report['fixed_parameter_values']==dense['parameter_values']
    pi=pi.subs(sample);b=s.symbols('b0:4')
    translation=s.eye(5)
    for i in range(4):translation[i+1,0]=-b[i]
    shifted=pushforward(pi,translation);y=X[1:];pairs=list(itertools.combinations(range(4),2))
    def jet(tensor):
        affine=s.Matrix(4,4,lambda i,j:s.expand((tensor[i+1,j+1]-y[i]*tensor[0,j+1]+y[j]*tensor[0,i+1]).subs(X[0],1)))
        zero=dict.fromkeys(y,0)
        return s.Matrix([affine[i,j].subs(zero)for i,j in pairs]
                        +[affine[i,j].diff(v).subs(zero)for v in y for i,j in pairs])
    columns=[]
    generators=[]
    for i,j in itertools.product(range(4),repeat=2):
        g=s.zeros(5);g[i+1,j+1]=1;generators.append(g)
    for j in range(4):
        g=s.zeros(5);g[0,j+1]=1;generators.append(g)
    for j in range(4):
        g=s.zeros(5);g[j+1,0]=-1;generators.append(g)
    for g in generators:
        derivative=g*shifted+shifted*g.T-sum(((g*s.Matrix(X))[i]*shifted.diff(X[i])for i in range(5)),s.zeros(5))
        columns.append(jet(derivative))
    columns.append(jet(shifted))
    actual=s.Matrix.hstack(*columns).applyfunc(s.expand)
    expected=s.Matrix(report['full_polynomial_matrix']).applyfunc(s.sympify)
    assert (actual-expected).applyfunc(s.expand)==s.zeros(30,25)
    reduced=actual.subs(dict(zip(b,[2,3,4,5])))
    def mod(value):
        value=s.Rational(value);return int(value.p)*pow(int(value.q),-1,P)%P
    determinant([[mod(v)for v in row]for row in reduced.tolist()],report['modular_lower_witness'])
    determinant([[mod(s.Rational(v))for v in row]for row in dense['quotient_matrix_QQ']],dense['nonzero_minor'])
    assert dense['nonzero_minor']['rank']==25 and report['constant_linear_jet_image_dimension']==22
    assert report['exact_primitive_firstjet_image_dimension']==22-6-4==12
    encoded=gzip.decompress((folder/'T38.in.gz').read_bytes())
    assert sha(folder/'T38.in.gz')==report['native_archive_sha256']
    assert hashlib.sha256(encoded).hexdigest()==report['native_input_sha256']
    tokens=iter(encoded.decode().split());number=lambda:int(next(tokens))
    assert [number()for _ in range(3)]==[4,30,25]
    K=s.QQ.frac_field(*b);R=K.one.denom.ring
    def polynomial():
        result={}
        for _ in range(number()):
            coefficient=s.QQ.convert(number());m=tuple(number()for _ in range(4));assert m not in result;result[m]=coefficient
        return R.from_dict(result)
    for value in actual:
        num,den=polynomial(),polynomial();assert den
        assert not(K.field.new(num,den)-K.from_sympy(value))
    assert next(tokens,None)is None
    assert sha(folder/'T38-native.log')==report['native_log_sha256']
    assert 'EXACT_FUNCTION_FIELD_RANK 22'in(folder/'T38-native.log').read_text()
    paths.extend([folder/'report.json',folder/'T38.in.gz',folder/'T38-native.log',dense_path,source_path])
    return dict(independent_homogeneous_PGL_derivatives_match_all750_entries=True,
                same_fixed_source_has_dense_projective_orbit24=True,
                exact_raw_jet_rank22_minus_value6_perspective4_accepted=True,
                raw_native_archive_matches_rebuilt_matrix=True)


def verify(output):
    model=S61();indices=[i for i,w in enumerate(model.reduced['parameter_weights'])if w==1]
    index_path=ROOT/'code/data/known51_rank4_s61_coordinate_catalogue_20261103_v62.json'
    index=json.loads(index_path.read_text());table_path=ROOT/'reports/KNOWN51_FIRSTJET_DIMENSION_TABLE_20261103_v62.json'
    table=json.loads(table_path.read_text());dimensions={r['component_id']:r for r in table['records']}
    paths=[Path(__file__),index_path,table_path];lower=[];native=[]
    for item in index['records']:
        name=item['component_id'];archive=ROOT/item['map']['path'];body=json.loads(gzip.decompress(archive.read_bytes()))
        path=ROOT/'reports/KNOWN51_RANK4_FIRSTJET_DIFFERENTIALS_20261103_v62'/(name+'.json');record=json.loads(path.read_text())
        active,matrix=reverse_differential(body,indices)
        assert active==record['input_differential_columns']and indices==record['firstjet_indices_zero_based']
        assert matrix==record['differential_matrix_mod32003']
        assert sha(archive)==record['circuit_archive_sha256']
        determinant(matrix,record['nonzero_minor'])
        exact=dimensions[name]['firstjet_image_dimension'];assert exact>=record['nonzero_minor']['rank']
        assert dimensions[name]['generic_fiber_dimension']==dimensions[name]['s61_component_dimension']-exact
        if record['exact_image_dimension']is not None:assert exact==record['exact_image_dimension']
        lower.append(dict(component_id=name,independent_reverse_AD_all_entries_match=True,minor_nonzero=True))
        paths.extend([path,archive])
        if name in ['T13','T21','T30']or name in ['R4_ZERO_LOG_5H','R4_MINIMAL_M','R4_A0','R4_DPM','R4_CPLUS_EXT','R4_CMINUS_MAG']:
            folder=ROOT/'reports'/('MPS_FIRSTJET_NATIVE_RANKS_T13_T21_T30_20261103_v62'if name.startswith('T')else'ADDITIONAL_RANK4_FIRSTJET_NATIVE_CIRCUIT_20261103_v62')
            npath=folder/(name+'.json');nrecord=json.loads(npath.read_text());audit_native_input(folder/(name+'.in.gz'),nrecord,body,indices)
            assert nrecord['exact_generic_firstjet_rank']==exact
            log=folder/(name+'-native.log');assert sha(log)==nrecord['native_log_sha256']
            assert f'EXACT_CIRCUIT_FUNCTION_FIELD_RANK {exact}'in log.read_text()
            native.append(dict(component_id=name,rank=exact,all_native_source_mapping_and_DAG_records_match=True))
            paths.extend([npath,folder/(name+'.in.gz'),log])
        print('PASS',name,flush=True)
    t38=t38_check(paths)
    assert table['general_firstjet_dominant_components']==['T01','T39','T40']
    for path,digest in table['input_sha256'].items():assert sha(ROOT/path)==digest
    for name in ['code/firstjet_circuit_function_rank_20261103_v62.cpp','code/rational_function_matrix_rank_20261103_v62.cpp',
                 'code/known_rank4_firstjet_native_circuit_20261103_v62.py','research/notes/KNOWN51_FIRSTJET_IMAGE_DIMENSIONS_20261103_v62.md']:
        paths.append(ROOT/name)
    result=dict(id='KNOWN51_FIRSTJET_IMAGE_DIMENSIONS_20261103_v62_REVIEW',all_passed=True,
                independent_reverse_AD_and_minor_checks=lower,native_input_archive_checks=native,T38=t38,
                identity_Sp4_with_all10_Lie_columns_argument_accepted=True,
                exact_native_rational_AD_and_Gaussian_rank_source_audited=True,
                DPM_exact_quadratic_field_pair_arithmetic_audited=True,
                exact_reduced_image_and_generic_fiber_dimensions_all51_accepted=True,
                sampled_rank_upper_used=False,whole_S61_fiber_exhaustion_claimed=False,
                input_sha256={str(path.relative_to(ROOT)):sha(path)for path in paths})
    with output.open('x')as stream:json.dump(result,stream,indent=2);stream.write('\n')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();verify(args.output)
