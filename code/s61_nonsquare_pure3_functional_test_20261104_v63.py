#!/usr/bin/env python3
"""One fixed constant/affine QQ functional test; no higher-degree search."""
import argparse
import copy
import hashlib
import json
import resource
import signal
import time
from fractions import Fraction as Q
from pathlib import Path

from s61_nonsquare_cubic_last_Q9_affine_bezout_20261104_v63 import Span, add, shift, encode

ROOT=Path(__file__).resolve().parents[1]
ID='S61_NONSQUARE_CUBIC_PURE3_FUNCTIONAL_TEST_20261104_v63'


def sha(p):
    h=hashlib.sha256()
    with p.open('rb')as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()


def mul(p,q):
    z={}
    for m,c in p.items():
        for n,d in q.items():
            k=tuple(x+y for x,y in zip(m,n));z[k]=z.get(k,Q())+c*d
    return{k:c for k,c in z.items()if c}


def main(out):
    began=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(57)
    out.mkdir(parents=True,exist_ok=False);read=lambda p:json.loads(p.read_text())
    b=ROOT/'reports/S61_NONSQUARE_CUBIC_PURE3_FUNCTIONAL_INPUT_20261104_v63'
    ip=ROOT/'reports/S61_NONSQUARE_CUBIC_ZW25_AUGMENTED_INPUT_20261104_v63/input.json'
    cp=ROOT/'reports/S61_NONSQUARE_CUBIC_ZW25_UNITS_20261104_v63/checkpoint.json'
    paths=[Path(__file__),ROOT/'code/s61_nonsquare_cubic_last_Q9_affine_bezout_20261104_v63.py',
           b/'report.json',b/'certificate.json',ip,cp]
    ir,inp=read(paths[2]),read(paths[3]);initial=read(ip);checkpoint=read(cp)
    assert ir['frozen']and sha(paths[3])==ir['certificate_sha256']
    assert ir['input_sha256'][str(ip.relative_to(ROOT))]==sha(ip)
    assert ir['actual_initial_and_post_G11_heavy_zero_gate_passed']and ir['G11_cleanup_exactly_zero_not_omitted_by_assumption']
    dec=lambda p:{tuple(m):Q(c)for m,c in p if Q(c)}
    quadrics=[[dec(p)for p in row]for row in inp['Q10x26_polynomials']]
    sources=[[dec(p)for p in row]for row in inp['cleared3x26_polynomials']]
    C={(1,0,0):Q(100),(0,2,0):Q(18),(0,1,0):Q(-9),(0,0,0):Q(-2)}
    targets=[(name,dec(p))for name,p in inp['known_targets']]
    clearedtargets=[(name,mul(C,p))for name,p in targets]
    rows=inp['old3_source_labels'];assert rows==[10,13,16]
    assert all(max(exps)==1 for exps in inp['normalized_denominator_C_exponents'])
    for source,row in zip(rows,sources):
        assert all(not initial['matrix25x32'][source][j]for j in range(26,32))
        assert all(dec(initial['matrix25x32'][source][j])=={m:-64000*c for m,c in row[j].items()}for j in range(26))
    shifts=[(0,0,0),(1,0,0),(0,1,0),(0,0,1)]
    labels=[(j,mon)for mon in shifts for j in range(26)]
    source_results=[];success=[]
    p11=dec(initial['matrix25x32'][11][26])
    expected_p=mul(mul({(0,0,0):Q(2,15625)},dict(dec(inp['known_targets'][3][1]))),mul(dec(inp['known_targets'][4][1]),dec(inp['known_targets'][4][1])))
    # p11 also has C^2; verify the known source pivot directly.
    expected_p=mul(expected_p,mul(C,C));assert p11==expected_p
    pure11=[dec(p)for p in initial['matrix25x32'][11][:26]]
    def verify_w(w,sourceindex,target):
        zero={}
        for row in quadrics:
            z={}
            for q,p in zip(row,w):z=add(z,mul(q,p))
            assert z==zero
        z={}
        for q,p in zip(sources[sourceindex],w):z=add(z,mul(q,p))
        assert z==mul(C,target)
        raw={}
        for q,p in zip(initial['matrix25x32'][rows[sourceindex]][:26],w):raw=add(raw,mul(dec(q),p))
        lam_u=mul({(0,0,0):Q(-64000)},mul(C,target));assert raw==lam_u
        numerator={}
        for q,p in zip(pure11,w):numerator=add(numerator,mul(q,p))
        # Lift only the deleted H*L component; no other heavy coordinate.
        lift={'pure26': [encode(p)for p in w],
              'deleted_HL26_component_num_den':[encode({m:-c for m,c in numerator.items()}),encode(p11)],
              'remaining5_heavy_components':['0']*5,
              'G11_annihilation_cleared_identity_verified':True,
              'Q10_annihilation_all_coefficients_verified':True,
              'actual_Gi_value':encode(lam_u),
              'unit_normalization':'Divide by lambda_i*u to obtain value1; both are fixed known FULLU units.'}
        # The cleared G11 pairing is p11*numerator + p11*(-numerator).
        assert not add(mul(p11,numerator),mul(p11,{m:-c for m,c in numerator.items()}))
        unitcomponents=[]
        for j,p in enumerate(w):
            if len(p)==1:
                (mon,c),=p.items()
                if mon[1]==0:
                    unitcomponents.append({'target_pure_column':j,'constant':str(c),'a_power':mon[0],'d_power':mon[2]})
        lift['known_unit_functional_components']=unitcomponents
        return lift
    for sourceindex,row in enumerate(sources):
        span=None;previous_mons=None;stages=[]
        for n,name in [(26,'constant'),(104,'affine_1_a_e_d')]:
            ss=shifts[:1]if n==26 else shifts
            equations=quadrics+[row]
            rowmons=[]
            for eq,entryrow in enumerate(equations):
                mons=set().union(*(set(shift(p,h))for h in ss for p in entryrow))
                if eq==10:mons|=set().union(*(set(p)for _,p in clearedtargets))
                rowmons.extend((eq,m)for m in sorted(mons,key=lambda m:(-sum(m),m)))
            ri={item:i for i,item in enumerate(rowmons)}
            if span is None:
                span=Span();start=0
            else:
                # Embed the constant span without re-running its pivots.
                mapping={i:ri[item]for i,item in enumerate(previous_mons)}
                span.basis=[(mapping[pivot],{mapping[k]:v for k,v in col.items()},rep.copy())for pivot,col,rep in span.basis]
                start=26
            columns=[]
            for j,h in labels[:n]:
                column={}
                for eq,entryrow in enumerate(equations):
                    for m,c in shift(entryrow[j],h).items():column[ri[eq,m]]=c
                columns.append(column)
            for index in range(start,n):span.extend(columns[index],index)
            for pivot,v,rep in span.basis:
                z={}
                for index,c in rep.items():z=add(z,columns[index],c)
                assert z==v and v[pivot]==1
            assert all(not span.query(col)[0]for col in columns)
            wanted=inp['equation_profiles'][sourceindex]['profiles'][0 if n==26 else 1]
            assert len(rowmons)==wanted['total_QQ_row_equations']
            results=[]
            for (targetname,target),(_,cleared)in zip(targets,clearedtargets):
                vector={ri[10,m]:c for m,c in cleared.items()}
                remainder,rep=span.query(vector);candidate=not remainder
                w=[{}for _ in range(26)];lift=None
                if candidate:
                    for index,c in rep.items():
                        j,h=labels[index];w[j][h]=w[j].get(h,Q())+c
                    lift=verify_w(w,sourceindex,target)
                    success.append([rows[sourceindex],name,targetname])
                else:
                    z=vector.copy()
                    for index,c in rep.items():z=add(z,columns[index],-c)
                    assert z==remainder
                results.append({'target_knownunit':targetname,'candidate':candidate,
                                'pure_functional26_polynomials':[encode(p)for p in w]if candidate else None,
                                'candidate_full_lift_and_identities':lift,
                                'exact_reduced_residual':encode(remainder),
                                'QQ_reduction_representation':encode(rep)})
            stages.append({'cofactor_class':name,'unknowns':n,'QQ_row_count':len(rowmons),
                           'column_rank':len(span.basis),'coefficient_row_labels':[[eq,list(m)]for eq,m in rowmons],
                           'basis':[{'pivot':p,'vector':encode(v),'source_witness':encode(rep)}for p,v,rep in span.basis],
                           'all_source_columns_and_basis_witnesses_checked':True,'target_results':results})
            previous_mons=rowmons
        source_results.append({'source_label':rows[sourceindex],'raw135':inp['raw_genuine_source135_indices'][sourceindex],
                               'constant_then_affine_stages':stages})
    cert={'cofactor_basis':['1','a','e','d'],'unknown_labels':[[j,list(m)]for j,m in labels],
          'equations':'Q10*w=0 and Gi_clear*w=C*u for each correctly source-boundGi',
          'fixed_targets':[name for name,_ in targets],
          'all3_source_results':source_results,'successful_source_class_unit_triples':success,
          'candidate_lift_recipe':'Only deletedH*L component=-G11pure*w/G11HL. Its denominator is a certified FULLU unit; sourceGi initialH6=0 so its value remainslambda_i*u. No full source/target transport.',
          'failure_scope_only_constant_and_affine_functional_classes':True}
    (out/'certificate.json').write_text(json.dumps(cert,separators=(',',':'))+'\n')
    result={'id':ID,'all_fixed_QQ_linear_algebra_and_witness_checks_passed':True,
            'source_labels':rows,'unknown_counts':[26,104],
            'QQ_ranks':[[x['column_rank']for x in r['constant_then_affine_stages']]for r in source_results],
            'successful_source_class_unit_triples':success,
            'all_candidate_Q10_and_originalGi_identities_verified':True,
            'all_candidate_G11_HL_lifts_verified':True,
            'known_lambda_C_clearing_preserved':True,
            'no_fulltransport_newpivot_higherdegree_GB_point_or_unknowninverse':True,
            'failure_not_global_rank_or_ideal_proper_statement':True,
            'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},
            'certificate_sha256':sha(out/'certificate.json'),
            'seconds':time.monotonic()-began,'frozen':True,'independent_review_pending':True}
    assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');signal.alarm(0)
    print(json.dumps(result,indent=2),flush=True)


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True)
    main(ap.parse_args().output.resolve())
