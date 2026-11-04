#!/usr/bin/env python3
"""One fixed Q10-only degree2 QQ right-kernel computation; no targets."""
import argparse
import hashlib
import json
import resource
import signal
import time
from fractions import Fraction as Q
from pathlib import Path

from s61_nonsquare_cubic_last_Q9_affine_bezout_20261104_v63 import add, shift, encode

ROOT=Path(__file__).resolve().parents[1]
ID='S61_NONSQUARE_CUBIC_Q10_DEGREE2_RIGHT_KERNEL_20261104_v63'


def sha(p):
    h=hashlib.sha256()
    with p.open('rb')as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()


def main(out):
    began=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(58)
    out.mkdir(parents=True,exist_ok=False);read=lambda p:json.loads(p.read_text())
    ib=ROOT/'reports/S61_NONSQUARE_CUBIC_PURE3_FUNCTIONAL_INPUT_20261104_v63'
    ab=ROOT/'reports/S61_NONSQUARE_CUBIC_PURE3_FUNCTIONAL_TEST_20261104_v63'
    paths=[Path(__file__),ROOT/'code/s61_nonsquare_cubic_last_Q9_affine_bezout_20261104_v63.py',
           ib/'report.json',ib/'certificate.json',ab/'report.json',ab/'certificate.json']
    ir,inp,ar,affine=[read(p)for p in paths[2:]]
    assert sha(paths[3])==ir['certificate_sha256']and sha(paths[5])==ar['certificate_sha256']
    assert ir['frozen']and ar['frozen']
    dec=lambda p:{tuple(m):Q(c)for m,c in p if Q(c)}
    matrix=[[dec(q)for q in row]for row in inp['Q10x26_polynomials']]
    assert sum(len(p)for row in matrix for p in row)==317
    mons=[(0,0,0),(1,0,0),(0,1,0),(0,0,1),(2,0,0),(1,1,0),(1,0,1),(0,2,0),(0,1,1),(0,0,2)]
    labels=[(j,m)for m in mons for j in range(26)]
    rowlabels=[];counts=[]
    for eq,row in enumerate(matrix):
        ss=set().union(*(set(shift(p,m))for m in mons for p in row))
        counts.append(len(ss));rowlabels.extend((eq,m)for m in sorted(ss,key=lambda m:(-sum(m),m)))
    rowid={v:i for i,v in enumerate(rowlabels)}
    assert len(rowlabels)==347 and len(labels)==260 and 317*len(mons)==3170
    columns=[]
    for j,m in labels:
        col={}
        for eq,row in enumerate(matrix):
            for ex,c in shift(row[j],m).items():col[rowid[eq,ex]]=c
        columns.append(col)
    previous=affine['all3_source_results'][0]['constant_then_affine_stages'][1]
    assert previous['unknowns']==104 and previous['column_rank']==104
    oldlabels=[(eq,tuple(m))for eq,m in previous['coefficient_row_labels']]
    basis=[]
    for b in previous['basis']:
        assert oldlabels[b['pivot']][0]<10
        pv=rowid[oldlabels[b['pivot']]]
        vector={rowid[oldlabels[int(i)]]:Q(c)for i,c in b['vector']if oldlabels[int(i)][0]<10}
        rep={int(i):Q(c)for i,c in b['source_witness']}
        assert all(i<104 for i in rep)
        reconstruct={}
        for i,c in rep.items():reconstruct=add(reconstruct,columns[i],c)
        assert reconstruct==vector and vector[pv]==1
        assert all(not vector.get(oldpivot,Q())for oldpivot,_,_ in basis)
        basis.append((pv,vector,rep))
    assert len(basis)==104
    kernels=[];free=[];complete=False;stop=None;processed=104
    def guard():
        if time.monotonic()-began>50:raise TimeoutError('50s_fixed_linear_algebra_boundary')
    def reconstruct(rep):
        value={}
        for i,c in rep.items():value=add(value,columns[i],c)
        return value
    def polyvector(rep):
        ww=[{}for _ in range(26)]
        for index,c in rep.items():
            j,mon=labels[index];ww[j][mon]=ww[j].get(mon,Q())+c
        return ww
    try:
        for index in range(104,260):
            guard();v=columns[index].copy();rep={index:Q(1)}
            for pv,b,br in basis:
                c=v.get(pv,Q())
                if c:v=add(v,b,-c);rep=add(rep,br,-c)
            if v:
                pv=min(v);lead=v[pv]
                vv={k:c/lead for k,c in v.items()};rr={k:c/lead for k,c in rep.items()}
                assert reconstruct(rr)==vv
                basis.append((pv,vv,rr))
            else:
                assert not reconstruct(rep)
                assert rep[index]==1 and all(i<index or i==index for i in rep)
                kernels.append(rep);free.append(index)
            processed=index+1
            if processed%16==0 or processed==260:
                (out/'progress.json').write_text(json.dumps({'processed_columns':processed,'QQ_rank':len(basis),'right_kernel_vectors':len(kernels),'seconds':time.monotonic()-began},indent=2)+'\n')
        complete=True
    except(TimeoutError,MemoryError)as exc:
        stop=str(exc)or type(exc).__name__
    # Every returned vector is checked against ALL original coefficient rows
    # and polynomial equations; no parameter specialization is used.
    polys=[]
    for index,rep in zip(free,kernels):
        assert not reconstruct(rep)
        ww=polyvector(rep);assert sum(labels[index][1])==2
        for row in matrix:
            value={}
            for p,w in zip(row,ww):
                for mon,c in w.items():value=add(value,shift(p,mon),c)
            assert not value
        polys.append(ww)
    # Distinct free-coordinate pivots make the returned vectors independent.
    assert all(rep.get(free[j],Q())==int(i==j)for i,rep in enumerate(kernels)for j in range(len(free)))
    assert len(basis)+len(kernels)==processed
    certificate={'coefficient_row_labels347':[[eq,list(m)]for eq,m in rowlabels],
                 'coefficient_rows_by_Qproduct':counts,
                 'right_unknown_labels260':[[j,list(m)]for j,m in labels],
                 'degree2_monomial_basis':mons,
                 'Q10x26_polynomial_input':[[encode(p)for p in row]for row in matrix],
                 'reused104_affine_Q_pivots_verified_by_exact_projection':True,
                 'processed_columns':processed,'complete_all260':complete,
                 'QQ_column_rank':len(basis),'QQ_right_kernel_dimension':len(kernels),
                 'column_basis':[{'pivot':p,'vector':encode(v),'source_witness':encode(r)}for p,v,r in basis],
                 'free_unknown_indices':free,
                 'right_kernel_vectors260':[encode(v)for v in kernels],
                 'actual26_polynomial_functionals_degree2':[[encode(p)for p in row]for row in polys],
                 'all_kernel_identities_verified_in_all347_QQ_rows_and10_polynomial_equations':True,
                 'all_kernel_free_coordinates_identity':True,
                 'all_kernel_vectors_have_degree2_component':True,
                 'exactness_proof':'Each successive column is either a new independent normalized pivot with an exact input witness, or has a saved dependence with distinctfree coefficient1. Rank+nullity equals all260 iffcomplete; no targets or genuine rows enter.',
                 'no_module_generator_or_unit_functional_claim':True,
                 'no_degree3_target_Gi_or_field_rank_probe':True}
    text=json.dumps(certificate,separators=(',',':'))+'\n';assert len(text)<31*1024*1024
    (out/'certificate.json').write_text(text)
    result={'id':ID,'complete_Q10_degree2_QQ_right_kernel':complete,
            'coefficient_equations':347,'unknowns':260,'input_coefficient_terms':3170,
            'QQ_rank':len(basis),'QQ_right_kernel_dimension':len(kernels),
            'reused_lower104_Q_pivots_no_lowerdegree_rerun':True,
            'all_returned_vectors_exactly_verified':True,
            'no_Gi_targets_degree3_GB_point_or_207expansion':True,
            'no_global_unit_or_freecomplement_or_full57_claim':True,
            'stop_reason':stop,
            'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},
            'certificate_sha256':sha(out/'certificate.json'),
            'seconds':time.monotonic()-began,'frozen':True,'independent_review_pending':True}
    assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');signal.alarm(0)
    print(json.dumps(result,indent=2),flush=True)


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True)
    main(ap.parse_args().output.resolve())
