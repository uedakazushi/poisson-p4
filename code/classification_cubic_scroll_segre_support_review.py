#!/usr/bin/env python3
"""Independent whole L6+Segre5 scroll support and actual-factor replay."""
import argparse
from fractions import Fraction as Q
import hashlib
import itertools as it
import json
from pathlib import Path
import resource
import time

from classification_cubic_scroll_abelian_review import add, scale, mul, diff, parse, monomials
from classification_cubic_scroll_four_sources_review import substitute

ROOT=Path(__file__).resolve().parents[1]
PAIRS=list(it.combinations(range(5),2))
COUNT=0


def check(a,b,label):
    global COUNT
    assert a==b,label
    COUNT+=1


def audit_support(gens):
    names=[f'r{i}'for i in range(12)]
    R=lambda t:parse(t,names)
    images=[R(f'r{i}')for i in range(12)]
    images[1],images[4],images[10],images[11]=R('2*r5'),{},R('2*r5+r10'),R('2*r2')
    gg=[substitute(g,images,12)for g in gens]
    aa=[R(t)for t in [
        'r2*r9-r5*r10','r2*r7-r8*r10','r8*r9-r5*r7',
        'r5*(r0-2*r8)+r2*(r3+r7)',
        '(r0-2*r8)*r9+r10*(r3+r7)',
        'r8*(r0-2*r8)+r2*r6','r7*(r3+r7)-r6*r9',
        '-((r0-2*r8)*r7+r10*r6)+2*(r8*(r3+r7)-r5*r6)']]
    forward=['0','r0','r0','0','0','-6*r0','0','0','0',
             '2*r3+r1','-2*r3-r1','-2*r4-2*r2','0','r4+r2','r4+3*r2',
             'r5','r7','r6','0','0','6*r0','0','-4*r0','r1']
    for i,t in enumerate(forward):
        check(gg[i],substitute(R(t),aa+[{}]*4,12),f'correct24->A8 g{i+1}')
    back=[gg[1],gg[23],scale(add(gg[14],scale(gg[13],-1)),Q(1,2)),
          scale(add(gg[9],scale(gg[23],-1)),Q(1,2)),
          scale(add(scale(gg[13],3),scale(gg[14],-1)),Q(1,2)),gg[15],gg[17],gg[16]]
    for i in range(8):check(back[i],aa[i],f'A8->24 inverse{i+1}')
    M=[[R('r2'),R('r5'),R('r8')],
       [R('r10'),R('r9'),R('r7')],
       [R('r0-2*r8'),R('-r3-r7'),R('-r6')]]
    minors={(rows,cols):add(mul(M[rows[0]][cols[0]],M[rows[1]][cols[1]]),
                          scale(mul(M[rows[0]][cols[1]],M[rows[1]][cols[0]]),-1))
            for rows in it.combinations(range(3),2)for cols in it.combinations(range(3),2)}
    for rows,cols,index,sign in [((0,1),(0,1),0,1),((0,1),(0,2),1,1),
        ((0,1),(1,2),2,-1),((0,2),(0,1),3,-1),((0,2),(0,2),5,-1),
        ((1,2),(0,1),4,-1),((1,2),(1,2),6,1)]:
        check(minors[rows,cols],scale(aa[index],sign),'allseven known minors')
    x=minors[(0,2),(1,2)];y=minors[(1,2),(0,2)]
    check(add(y,scale(x,2)),aa[7],'remaining y+2x')
    # All earlier unit reductions are checked with literal polynomial identities.
    for i,t in [(19,'(r11-r2)*(r11-2*r2)'),(1,'(r11-2*r2)*r5'),
                (9,'(r11-2*r2)*r8'),(20,'(r11-2*r2)*(r10-r1)'),
                (4,'r4*(-3*(r10-r1)+4*r5)'),(8,'r4*(r10-r1-r5)'),
                (5,'r4*(r4+r9)'),(13,'r4*r7')]:
        check(gens[i-1],R(t),'early whole complement forcing')
    stage=[R(f'r{i}')for i in range(12)];stage[11]=R('2*r2');stage[4]={}
    gstage=[substitute(g,stage,12)for g in gens]
    check(add(gstage[1],scale(gstage[2],-1)),R('(r10-r1+2*r5)*(r1-2*r5)'), 'a2f split1')
    check(add(gstage[5],scale(gstage[2],6)),R('-2*(r1-2*r5)*(r10-r1+r5)'), 'a2f split2')
    check(gstage[21],R('2*r2*(r1-2*r5)'),'a2f split3')
    check(gstage[6],R('-r9*(r1-2*r5)'),'a2f split4')
    return dict(all24_correct_eight_equations_both_directions=True,
        allnine_minors_and_remaining_y_plus_2x=True,
        reduced_support_rankone_hand_proof='row1!=0 => row2=k row1 => y=k minor13cols13=0 =>x0; row1=0 =>x0=>y0',
        all_early_unit_zero_complements_retained=True,
        scheme_radical_or_missing_minor_literal_membership_not_claimed=True)


def actual_source(name,pars,texts,saved,source,B,gens):
    names=pars+[f'x{i}'for i in range(5)];np=len(pars);n=len(names)
    r=[parse(t,names)for t in texts]
    check(pars,saved['all_parameters'],name+' parameter order')
    check(r,[parse(t,names)for t in saved['full_original_r0_to_r11_map']],name+' independent whole map')
    for g in gens:check(substitute(g,r,n),{},name+' all24 generators')
    quadrics=[parse(t,names)for t in ['x0*x3-x1*x2','x0*x4-x1*x3','x2*x4-x3**2']]
    P=[[{}for j in range(5)]for i in range(5)]
    for p,(i,j)in enumerate(PAIRS):
        P[i][j]=add(*(mul(quadrics[g],add(*(scale(r[k],B[10*g+p][k])for k in range(12))))for g in range(3)))
        P[j][i]=scale(P[i][j],-1)
        check(P[i][j],parse(saved['complete_actual_brackets'][str((i,j))],names),name+' full150 tensor')
    deriv=[[[diff(P[i][j],np+l)for l in range(5)]for j in range(5)]for i in range(5)]
    J={tri:add(*(mul(P[i][l],deriv[j][k][l])for i,j,k in
        [tri,(tri[1],tri[2],tri[0]),(tri[2],tri[0],tri[1])]for l in range(5)))
       for tri in it.combinations(range(5),3)}
    for quad in it.combinations(range(5),4):
        check(add(*(scale(mul(parse(f'x{quad[i]}',names),J[quad[:i]+quad[i+1:]]),(-1)**i)
                     for i in range(4))),{},name+' full350 direct')
    rp=[mul(r[i],r[j])for i,j in source['all_reduced_source_quadratic_monomials']]
    ev=lambda row:add(*(scale(rp[int(i)],Q(c))for i,c in row.items()))
    for row in source['all350_reduced_Jacobi_rows']:check(ev(row),{},name+' full350 savedsource')
    D={}
    for m in range(5):
        i,j,k,l=[i for i in range(5)if i!=m]
        pf=add(mul(P[i][j],P[k][l]),scale(mul(P[i][k],P[j][l]),-1),mul(P[i][l],P[j][k]))
        D=add(D,scale(mul(parse(f'x{m}',names),pf),(-1)**m))
    check(D,parse(saved['actual_signed_Euler_Pfaffian'],names),name+' signedPfD')
    for ex,row in zip(monomials(5,5),source['all126_reduced_actual_D_rows']):
        check({e[:np]+(0,)*5:c for e,c in D.items()if e[np:]==ex},ev(row),name+' full126 savedD')
    prod=parse(saved['factor_coefficient'],names);degrees=[];parameter_factors=[]
    for item in saved['full_polynomial_factors']:
        factor=parse(item['factor'],names);assert factor
        degree=item['physical_coordinate_degree']
        check({sum(e[np:])for e in factor},{degree},name+' positive homogeneous degree')
        for _ in range(item['multiplicity']):
            prod=mul(prod,factor)
            if degree>0:degrees.append(degree)
            else:parameter_factors.append(item['factor'])
    check(prod,D,name+' whole actual factor product')
    assert D and len(degrees)>=2 and sum(degrees)==5
    return dict(source=name,all_parameters=pars,all12_original_maps=True,
        all24_original_generators_zero=True,full150_original_tensor_reconstructed=True,
        all350_direct_and_saved_rows_zero=True,all126_signed_Pfaffian_coefficients_match=True,
        full_factor_product_identity=True,physical_factor_degrees=degrees,
        parameter_only_factors=parameter_factors,
        every_nonzero_specialization_reducible=True,all_zero_parameter_complements_retained=True)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--producer',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--job',type=Path,required=True)
    args=ap.parse_args();begun=time.monotonic()
    resource.setrlimit(resource.RLIMIT_CPU,(110,115));resource.setrlimit(resource.RLIMIT_AS,(2048*1024**2,)*2)
    contract=json.loads(args.job.read_text())
    paths=[p['path']for p in contract['exact_inputs_and_hashes']]
    assert all(hashlib.sha256((ROOT/p['path']).read_bytes()).hexdigest()==p['sha256']
               for p in contract['exact_inputs_and_hashes'])
    source=json.load(open(ROOT/'reports/CLASSIFICATION_CUBIC_SCROLL_ABELIAN_CONORMAL_20261107.json'))['cases'][0]
    dim=json.load(open(ROOT/'reports/CLASSIFICATION_CUBIC_SCROLL_ABELIAN_DIMENSION_20261107.json'))['cases'][0]
    oldreview=json.load(open(ROOT/'reports/CLASSIFICATION_CUBIC_SCROLL_ABELIAN_REVIEW_20261107.json'))
    assert oldreview['all_passed'] and oldreview['all24_generators_original_Jacobi_constant_multiples']
    gens=[parse(t,[f'r{i}'for i in range(12)])for t in dim['generators']]
    support=audit_support(gens)
    producer=json.loads(args.producer.read_text());records=producer['cases'];assert len(records)==2
    B=[[Q(v)for v in row]for row in source['complete_original30_by_reduced_source_matrix']]
    out=[]
    for saved in records:
        name=saved['source'];pars=saved['all_parameters'];assert len(pars)==6
        if name=='L6':
            A,P,S,T,V,R=pars
            texts=[A,P,S,T,V,'0',R,'0','0',f'-{V}',P,S]
        else:
            A,Bb,C,D,E,F=pars
            texts=[f'{C}*{D}+2*{A}*{F}',f'2*{A}*{E}',f'{A}*{D}',f'-{C}*{E}-{Bb}*{F}',
                   '0',f'{A}*{E}',f'-{C}*{F}',f'{Bb}*{F}',f'{A}*{F}',f'{Bb}*{E}',
                   f'2*{A}*{E}+{Bb}*{D}',f'2*{A}*{D}']
        out.append(actual_source(name,pars,texts,saved,source,B,gens))
    result=dict(task_id='E_reduced_nonnormal',all_passed=True,
        scope='Whole smooth cubic scroll rankzero with reduced irreducible quintic is empty; all reducible parent membership separate',
        corrected_whole_L6_union_Segre5_support=True,support_review=support,cases=out,
        no_producer_or_CAS_import=True,no_new_rank_minor_operator_point=True,
        positive_degree_factor_specialization_logic='Every specialized nonzero polynomial factor remains homogeneous of the same positive physical degree; otherwise Dzero. At least two physical factors are present at every nonzero D.',
        all_inverted_parameter_complements_retained=True,new_component_count=0,
        independent_identity_count=COUNT,
        input_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},
        reviewer_source_bytes=Path(__file__).read_text(),seconds=time.monotonic()-begun,
        maxrss_kb=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    with args.output.open('x')as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps({k:result[k]for k in ['all_passed','independent_identity_count','seconds','maxrss_kb']}))


if __name__=='__main__':main()
