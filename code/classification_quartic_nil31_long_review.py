#!/usr/bin/env python3
"""Independent nil31-long full support and orbit22 bound, without CAS/producer."""
import argparse
from fractions import Fraction as F
import hashlib
import itertools as it
import json
from pathlib import Path
import re
import resource
import time

from classification_quadratic_cofactor_review import verify_source
from classification_minimal_quartic_review import LaurentRing,determinant
from classification_unimodular_logarithmic_review import RF,rank

ROOT=Path(__file__).resolve().parents[1]

def power(R,p,n):
    if n<0:p=R.inverse_monomial(p);n=-n
    out=R.constant(1)
    for _ in range(n):out=R.mul(out,p)
    return out

def substitute(R,p,values):
    out={}
    for ex,co in p.items():
        term=R.constant(co)
        for i,n in enumerate(ex):
            if n:term=R.mul(term,power(R,values.get(i,R.variable(R.names[i])),n))
        out=R.add(out,term)
    return out

def linear_matrix(R,equations,variables):
    matrix=[]
    for p in equations:
        row=[{}for _ in variables]
        for ex,co in p.items():
            total=sum(ex[i]for i in variables)
            assert total<=1
            if not total:continue
            col=next(j for j,i in enumerate(variables)if ex[i])
            e=list(ex);e[variables[col]]=0;e=tuple(e)
            row[col][e]=row[col].get(e,F())+co
        matrix.append([{e:c for e,c in p.items()if c}for p in row])
    return matrix

def minor(R,matrix,rows,columns):
    return determinant(R,[[matrix[i][j]for j in columns]for i in rows])

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args();start=time.monotonic()
    resource.setrlimit(resource.RLIMIT_CPU,(110,115));resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,)*2)
    srcpath=ROOT/'reports/CLASSIFICATION_HYPERPLANE_QUARTIC_NILPOTENT_SOURCES_20261107.json'
    path=ROOT/'reports/CLASSIFICATION_HYPERPLANE_QUARTIC_NIL31_LONG_SUPPORT_20261107.json'
    report=json.loads(path.read_text());data=report['nil31_long_support']
    assert hashlib.sha256(report['source_bytes'].encode()).hexdigest()==report['source_sha256']
    original=json.loads(srcpath.read_text())['full_nilpotent_sources'][4]
    source_review=verify_source(original,5);assert len(source_review['free_indices'])==15
    R=LaurentRing([f'n{i}'for i in range(15)],['x','y','z','w'],invert=['n14'])
    parse=lambda s:R.parse(re.sub(r'\ba\b','n14',str(s)))
    J=list(map(parse,original['distinct_nonzero_jacobi_coefficients']))
    pf=parse(original['actual_affine_pfaffian']);top=parse(original['highest_quartic_pfaffian'])
    enc=lambda m:[[R.encode(v)for v in row]for row in m]
    # Every necessary boundary forcing is an exact full-coefficient identity.
    assert substitute(R,J[0],{14:{}})==parse('2*n11*(n13+1)/5')
    bad={14:{},13:parse('-1')}
    assert substitute(R,J[2],bad)==parse('-2*n9')
    assert substitute(R,J[13],bad)==parse('2*n7')
    assert substitute(R,J[3],bad|{9:{},7:{}})==parse('-2*n11**2/3')
    assert substitute(R,J[2],{14:{},11:{}})==parse('2*n9*(8*n13+3)/5')
    assert substitute(R,J[13],{14:{},11:{},13:parse('-3/8')})==parse('n7')
    assert substitute(R,J[14],{14:{},11:{},13:parse('-3/8')})==parse('5*n10/8')
    assert substitute(R,J[13],{14:{},11:{},9:{}})==parse('2*n7*(1-4*n13)/5')
    assert not substitute(R,top,{14:{},11:{},9:{},7:{}})
    bounds=data['n14_zero_complement'];boundary_records=[]
    for label,base,variables,rows,expected_det,subkey,pfkey,expected_pf in [
        ('n9_nonzero',{14:'0',13:'-3/8',11:'0',7:'0',10:'0'},[2,6,1,5,0],[3,7,8,9,18],'-45/8192',
         'n9_nonzero_forced_substitution','n9_nonzero_exact_pfaffian',
         '-n9*(3*w-8*n12)**2*(64*n12*n9+8*n8-16*n9*w+2*y-z**2)/18'),
        ('n9_zero_n7_nonzero',{14:'0',13:'1/4',11:'0',9:'0'},[10,2,6,1,0,5],[3,7,8,9,14,18],'-315/1024',
         'n9_zero_n7_nonzero_forced_substitution','exact_pfaffian','-n7*(2*y-z**2+4*n8)**2/4')]:
        base={i:parse(v)for i,v in base.items()}
        eq=[substitute(R,p,base)for p in J]
        matrix=linear_matrix(R,eq,variables)
        assert minor(R,matrix,rows,range(len(variables)))==parse(expected_det)
        values={int(k[1:]):parse(v)for k,v in bounds[subkey].items()}
        assert all(not substitute(R,p,values)for p in J)
        actual=substitute(R,pf,values)
        assert actual==parse(bounds[pfkey])==parse(expected_pf)
        boundary_records.append(dict(name=label,variables=variables,uniform_minor=expected_det,
            all21_Jacobi_coefficients_zero=True,full_Pfaffian_factorization_verified=True,
            full_linear_matrix=enc(matrix)))
    assert J[10]==parse('3*n14*(6*n13+1)/5')
    base={14:parse('n14'),13:parse('-1/6')}
    first={14:parse('n14'),13:parse('-1/6'),11:parse('-9*n14*n4'),9:parse('6*n14*n2'),
           7:parse('-6*n14*(n6+2*n8)'),10:parse('6*n14*(3*n4*n12-n3)')}
    initial=linear_matrix(R,[substitute(R,J[i],base)for i in [0,2,13,14]],[11,9,7,10])
    assert minor(R,initial,range(4),range(4))==parse('2/27')
    assert all(not substitute(R,J[i],first)for i in [0,2,13,14])
    assert substitute(R,J[7],first)==parse('n2*(1-216*n14**2*n2)/2')
    records=[]
    assert len(data['all_two_complete_families'])==2
    for branch,at2,stage1,stage2,free in [
        (0,'0',{6:'0',1:'-27*n14*n4**2'},
         {0:'18*n14*n4*(3*n4*n12-n3)',5:'36*n14*n8**2'},[14,3,4,8,12]),
        (1,'1/(216*n14**2)',{8:'-n12/(6*n14)',1:'-27*n14*n4**2-n6/(12*n14)-n12/(36*n14**2)'},
         {0:'54*n14*n4**2*n12-18*n14*n3*n4+n6*n12/(2*n14)',5:'n12**2/n14-6*n6*n12'},[14,3,4,6,12])]:
        at2=parse(at2);stage1={k:parse(v)for k,v in stage1.items()};stage2={k:parse(v)for k,v in stage2.items()}
        residual=[substitute(R,substitute(R,p,first),{2:at2})for p in J]
        selected=[p for p in residual if all(not ex[0]and not ex[5]for ex in p)]
        matrix1=linear_matrix(R,selected,[1,6,8])
        columns=[0,1]if branch==0 else[0,2]
        det1='1/3'if branch==0 else'-1/3'
        assert minor(R,matrix1,[3,8],columns)==parse(det1)
        assert all(not substitute(R,p,stage1)for p in selected)
        # This affine solution has the remaining one of(n6,n8) free, so rank<=2;
        # the nonzero constant minor proves exact rank2 at EVERY parameter.
        second=[substitute(R,p,stage1)for p in residual]
        matrix2=linear_matrix(R,second,[0,5])
        det2='-5/9'if branch==0 else'-1/2'
        assert minor(R,matrix2,[9,18],[0,1])==parse(det2)
        assert all(not substitute(R,p,stage2)for p in second)
        replacements={2:at2}|stage1|stage2
        complete={j:substitute(R,first.get(j,R.variable(f'n{j}')),replacements)for j in range(15)}
        complete.update(replacements)
        for j in range(15):complete[j]=substitute(R,complete[j],{k:v for k,v in complete.items()if k not in free})
        assert set(i for p in complete.values()for ex in p for i,k in enumerate(ex[:15])if k)==set(free)
        assert all(not substitute(R,p,complete)for p in J)
        claim=data['all_two_complete_families'][branch]
        assert complete=={int(k[1:]):parse(v)for k,v in claim['complete_source_substitutions'].items()}
        assert substitute(R,pf,complete)==parse(claim['actual_pfaffian'])
        assert [substitute(R,parse(p),complete)for p in original['brackets']]==list(map(parse,claim['universal_brackets']))
        records.append(dict(name=claim['name'],free_parameters=free,first_full_matrix=enc(matrix1),
            first_everywhere_rank=2,first_uniform_minor=det1,second_full_matrix=enc(matrix2),
            second_everywhere_rank=2,second_uniform_minor=det2,
            all15_substitutions_and_full21Jacobi_and_sixbrackets_and_Pfaffian_verified=True,
            source_parameter_count=5))
    A=[[F(v)for v in row]for row in original['field_linear_matrix']]
    c=list(map(F,original['field_constant']));columns=[]
    for i,j in it.product(range(4),repeat=2):
        B=[[F(a==i and b==j)for b in range(4)]for a in range(4)]
        columns.append([sum(B[a][k]*A[k][b]-A[a][k]*B[k][b]for k in range(4))for a,b in it.product(range(4),repeat=2)]+[sum(B[a][k]*c[k]for k in range(4))for a in range(4)])
    for j in range(4):columns.append([F()]*16+[-A[a][j]for a in range(4)])
    affine=list(map(list,zip(*columns)))
    assert affine==[[F(v)for v in row]for row in data['dimension_bound']['exact_affine_orbit_tangent_matrix']]
    assert rank([[RF(v)for v in row]for row in affine])==14
    assert rank([[RF(v)for v in row]for row in affine[:16]])==10
    result=dict(task_id='E_reduced_nonnormal',all_passed=True,no_CAS_import=True,no_producer_import=True,
        complete_original_source_review=source_review,boundary_exclusions=boundary_records,
        entire_n14_zero_complement_excluded=True,initial_four_by_four_minor='2/27',
        full_two_families=records,linear_orbit_rank=10,affine_orbit_rank=14,
        projective_saturated_dimension_upper_bound=22,E_membership_proved=False,G_generic_type_excluded=True,
        input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in [path,srcpath]},
        frozen_producer_source_sha256=report['source_sha256'],
        reviewer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),reviewer_bytes=Path(__file__).read_text(),
        seconds=time.monotonic()-start,maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k]for k in ['all_passed','projective_saturated_dimension_upper_bound','seconds','maximum_RSS_KiB']}))

if __name__=='__main__':main()
