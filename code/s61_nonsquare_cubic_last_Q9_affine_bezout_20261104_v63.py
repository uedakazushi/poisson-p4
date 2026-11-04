#!/usr/bin/env python3
"""Fixed constant/affine QQ coefficient-span tests on nine frozen polynomials."""
import argparse
import hashlib
import json
import resource
import signal
import time
from fractions import Fraction as Q
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ID = 'S61_NONSQUARE_CUBIC_LAST_Q9_AFFINE_BEZOUT_20261104_v63'
ZERO = (0,0,0)


def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()


def add(v,w,scale=Q(1)):
    z=v.copy()
    for m,c in w.items():
        z[m]=z.get(m,Q())+scale*c
        if not z[m]:del z[m]
    return z


def shift(p,m):
    return {tuple(x+y for x,y in zip(k,m)):c for k,c in p.items()}


def encode(p):
    return [[list(m) if isinstance(m,tuple) else m,str(c)] for m,c in sorted(p.items())]


class Span:
    def __init__(self):self.basis=[]
    def extend(self,column,index):
        v,rep=column.copy(),{index:Q(1)}
        for pivot,b,br in self.basis:
            co=v.get(pivot,Q())
            if co:v=add(v,b,-co);rep=add(rep,br,-co)
        if not v:return False
        pivot=min(v);lead=v[pivot]
        self.basis.append((pivot,{k:c/lead for k,c in v.items()},{k:c/lead for k,c in rep.items()}))
        return True
    def query(self,column):
        v,rep=column.copy(),{}
        for pivot,b,br in self.basis:
            co=v.get(pivot,Q())
            if co:v=add(v,b,-co);rep=add(rep,br,co)
        return v,rep


def main(out):
    began=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(57)
    out.mkdir(parents=True,exist_ok=False)
    read=lambda p:json.loads(p.read_text())
    b=ROOT/'reports/S61_NONSQUARE_CUBIC_LAST_Q9_PREFLIGHT_20261104_v63'
    paths=[Path(__file__),b/'report.json',b/'certificate.json']
    report,cert=read(paths[1]),read(paths[2]);assert report['frozen']and sha(paths[2])==report['certificate_sha256']
    frozen=cert['actual_nine_coefficients'];raw=[];normalized=[]
    for i,v in enumerate(frozen):
        nn,dd=v['actual_coefficient'];assert len(dd)==1 and dd[0][0]==[0,0,0]
        denom=Q(dd[0][1]);p={tuple(m):Q(c)/denom for m,c in nn};raw.append(p)
        content=1 if i<8 else 0
        assert v['named_factor_powers_a_d_F_B_C_J']==[content,0,0,0,0,0]
        assert all(m[0]>=content for m in p)
        g={(m[0]-content,m[1],m[2]):c for m,c in p.items()};normalized.append(g)
        assert shift(g,(content,0,0))==p
    targets=[('1',{ZERO:Q(1)}),('a',{(1,0,0):Q(1)}),('d',{(0,0,1):Q(1)}),
             ('F',{(1,0,0):Q(500),(0,2,0):Q(90),(0,1,0):Q(3),ZERO:Q(-24)}),
             ('B',{(1,0,1):Q(50),(0,2,1):Q(18),(0,1,1):Q(-9),(0,0,1):Q(-2),(0,1,0):Q(60),ZERO:Q(-40)}),
             ('C',{(1,0,0):Q(100),(0,2,0):Q(18),(0,1,0):Q(-9),ZERO:Q(-2)}),
             ('J',{(1,0,1):Q(5),(0,1,0):Q(-6),ZERO:Q(-4)}),
             ('a*d',{(1,0,1):Q(1)})]
    shifts=[ZERO,(1,0,0),(0,1,0),(0,0,1)]
    labels=[(i,m)for m in shifts for i in range(9)]
    polys=[shift(normalized[i],m)for i,m in labels]
    mon=sorted(set().union(*(set(p)for p in polys),*(set(p)for _,p in targets)),key=lambda m:(-sum(m),m))
    assert len(mon)==117 and sum(map(len,polys))+sum(len(p)for _,p in targets)==1513
    mi={m:i for i,m in enumerate(mon)}
    cols=[{mi[m]:c for m,c in p.items()}for p in polys]
    tc=[{mi[m]:c for m,c in p.items()}for _,p in targets]
    span=Span();stages=[]
    for stop,name in [(9,'constant'),(36,'affine_1_a_e_d')]:
        begin=0 if stop==9 else 9
        for i in range(begin,stop):span.extend(cols[i],i)
        # Independent witness reconstruction: every saved basis equals A*rep.
        for pivot,v,rep in span.basis:
            check={}
            for i,c in rep.items():check=add(check,cols[i],c)
            assert check==v and v[pivot]==1
        for column in cols[:stop]:assert not span.query(column)[0]
        results=[]
        for (target,p),column in zip(targets,tc):
            remainder,rep=span.query(column)
            candidate=not remainder
            cofs=[{}for _ in range(9)]
            if candidate:
                for index,c in rep.items():
                    source,m=labels[index];cofs[source][m]=cofs[source].get(m,Q())+c
                identity={}
                for i,cof in enumerate(cofs):
                    for m,c in cof.items():identity=add(identity,shift(normalized[i],m),c)
                assert identity==p
                # Verify directly in ALL original QQ coefficients, clearing
                # only the certified a content rather than evaluating points.
                original_identity={}
                for i,cof in enumerate(cofs):
                    for m,c in cof.items():
                        extra=(0,0,0)if i<8 else(1,0,0)
                        original_identity=add(original_identity,shift(raw[i],tuple(x+y for x,y in zip(m,extra))),c)
                assert original_identity==shift(p,(1,0,0))
            else:
                # Save the exact reduced obstruction in this finite class.
                check=column.copy()
                for i,c in rep.items():check=add(check,cols[i],-c)
                assert check==remainder
            results.append({'target_known_unit':target,'in_fixed_span':candidate,
                            'normalized_source_cofactors9': [encode(v)for v in cofs]if candidate else None,
                            'actual_source_cofactor_denominators':['a']*8+['1']if candidate else None,
                            'all_original_QQ_coefficients_identity_verified':candidate,
                            'reduced_remainder':encode(remainder),
                            'reduction_representation_on_fixed_columns':encode(rep)})
        stages.append({'cofactor_class':name,'unknowns':stop,'QQ_coefficient_rank':len(span.basis),
                       'basis':[{'pivot_monomial_row':pv,'column':encode(v),'source_witness':encode(r)}for pv,v,r in span.basis],
                       'all_fixed_source_columns_in_saved_basis_span':True,
                       'target_results':results})
    good=[(st['cofactor_class'],r['target_known_unit'])for st in stages for r in st['target_results']if r['in_fixed_span']]
    certout={'normalization':'g_i=f_i/a for first8; g_8=f_8; every literal a factor verified',
             'source9_original_target_labels':[v['original_target_column']for v in frozen],
             'normalized9_polynomials':[encode(p)for p in normalized],
             'cofactor_monomials':[list(m)for m in shifts],
             'QQ_row_monomials117':[list(m)for m in mon],
             'fixed36_column_labels':[[i,list(m)]for i,m in labels],
             'known_unit_targets':[[name,encode(p)]for name,p in targets],
             'all_constant_and_affine_stages':stages,
             'success_class_target_pairs':good,
             'no_higherdegree_or_idealproper_inference':True,
             'all_pure_and_heavy_tails_remain_in_full_actual_row':True}
    (out/'certificate.json').write_text(json.dumps(certout,separators=(',',':'))+'\n')
    result={'id':ID,'fixed_QQ_coefficient_linear_algebra_completed':True,
            'input_expanded_terms':1513,'coefficient_row_equations':117,
            'cofactor_unknowns':[9,36],'QQ_ranks':[st['QQ_coefficient_rank']for st in stages],
            'successful_target_pairs':good,
            'all_candidates_verified_in_ALL_original_QQ_coefficients':True,
            'only_exact_recorded_a_contents_stripped':True,
            'no_higherdegree_GB_point_saturation_unknown_factor_or_automatic_transport':True,
            'failure_scope_only_constant_and_affine_cofactor_classes':True,
            'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},
            'certificate_sha256':sha(out/'certificate.json'),
            'seconds':time.monotonic()-began,'frozen':True,'independent_review_pending':True}
    assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n')
    signal.alarm(0);print(json.dumps(result,indent=2),flush=True)


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True)
    main(ap.parse_args().output.resolve())
