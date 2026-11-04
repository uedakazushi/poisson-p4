#!/usr/bin/env python3
"""Seventy fixed scaling identities; no minors, Gaussian or new B-block work."""
import argparse
import hashlib
import json
import resource
import signal
import time
from fractions import Fraction as Q
from pathlib import Path

import sympy as s
from sympy.polys.rings import ring

ROOT=Path(__file__).resolve().parents[1]
ID='S61_NONSQUARE_CUBIC_GENUINE14_H5_PREFLIGHT_20261104_v63'


def sha(p):
    h=hashlib.sha256()
    with p.open('rb')as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()


def main(out):
    began=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(57)
    out.mkdir(parents=True,exist_ok=False);read=lambda p:json.loads(p.read_text())
    final=ROOT/'reports/S61_NONSQUARE_CUBIC_LAST_Q9_UNIT_COMPLETION_20261104_v63'
    parent=ROOT/'reports/S61_NONSQUARE_CUBIC_Q2_L_UNIT_COMPLETION_20261104_v63'
    old=ROOT/'reports/S61_NONSQUARE_CUBIC_ZW25_UNITS_20261104_v63'
    raw=ROOT/'reports/S61_NONSQUARE_RAW_H_LINEAR_UNITS_20261104_v63'
    generic=ROOT/'reports/S61_NONSQUARE_RAW_H_GENERIC_RANK_20261104_v63'
    inputs=ROOT/'reports/S61_NONSQUARE_CUBIC_ZW25_AUGMENTED_INPUT_20261104_v63/input.json'
    peers=[ROOT/'reports'/name/'report.json'for name in[
        'S61_NONSQUARE_RAW_H_LINEAR_UNITS_20261104_v63_REVIEW',
        'S61_NONSQUARE_RAW_H_GENERIC_RANK_20261104_v63_REVIEW',
        'S61_NONSQUARE_RAW_H_MINOR_OPEN_20261104_v63_REVIEW']]
    paths=[Path(__file__),final/'report.json',final/'certificate.json',parent/'report.json',parent/'certificate.json',
           old/'report.json',old/'checkpoint.json',raw/'report.json',raw/'preflight.json',raw/'checkpoint.json',
           generic/'report.json',generic/'certificate.json',inputs,*peers]
    fr,fc,pr,pc,orr,oc,rr,rp,rc,gr,gc,data=[read(p)for p in paths[1:13]]
    assert sha(paths[2])==fr['certificate_sha256']and sha(paths[4])==pr['certificate_sha256']
    assert sha(paths[6])==orr['evidence_sha256']['checkpoint.json']
    assert sha(paths[9])==rr['evidence_sha256']['checkpoint.json']
    assert sha(paths[11])==gr['certificate_sha256']
    assert all(read(p)['all_passed']for p in peers)
    rows=fc['new14_source_labels'];columns=[27,28,29,30,31]
    assert rows==[10,12,13,14,15,16,17,18,19,20,21,22,23,24]
    traces=[]
    for i,original in enumerate(rows):
        for j in columns:
            node=fc['actual14x21_scalar_roots'][i][fc['new21_target_labels'].index(j)]
            if not node:
                assert not oc['matrix25x32'][original][j][0];trace=[0]
            else:
                n=fc['scalar_DAG']['nodes'][node];assert n[0]=='matrix_parent'
                leaf=pc['scalar_DAG']['nodes'][n[1]];assert leaf==['checkpoint',original,j]
                trace=[node,n[1],leaf]
            traces.append([original,j,trace])
    R,a,e,d=ring('a,e,d',s.QQ);K=R.to_field()
    decp=lambda ts:R.from_dict({tuple(m):s.QQ(c)for m,c in ts})
    dec=lambda v:K(decp(v[0]))/K(decp(v[1]))
    enc=lambda p:[[list(m),str(c)]for m,c in p.items()]
    known=[a,d,500*a+90*e**2+3*e-24,(50*a+18*e**2-9*e-2)*d+60*e-40,
           100*a+18*e**2-9*e-2,5*a*d-6*e-4]
    F,B,C=known[2:5]
    def u(q):
        z=K.from_expr(s.sympify(q['constant']))
        for f,n in zip(known,q['a_d_F_B_C_J_powers']):z*=K(f)**n
        return z
    def mulunit(v,w,divide=False):
        return {'constant':str(Q(v['constant'])/Q(w['constant'])if divide else Q(v['constant'])*Q(w['constant'])),
                'a_d_F_B_C_J_powers':[x-y if divide else x+y for x,y in zip(v['a_d_F_B_C_J_powers'],w['a_d_F_B_C_J_powers'])]}
    def factorunit(poly):
        powers=[]
        for f in known:
            n=0
            while not poly.rem(f):poly=poly.exquo(f);n+=1
            powers.append(n)
        assert poly.is_ground and poly
        return {'constant':str(poly),'a_d_F_B_C_J_powers':powers}
    rawM=[[decp(q)for q in row]for row in rp['matrix15x6']]
    N=[[decp(q)for q in row]for row in rc['matrix15x6']]
    keep=rp['remaining15_genuine_indices'];assert keep==[0,1,2,3,4,5,7,8,9,10,11,12,13,14,15]
    sunit={'constant':'1/12500000','a_d_F_B_C_J_powers':[0,0,1,2,1,0]}
    gammaunit={'constant':'-4/625','a_d_F_B_C_J_powers':[0,0,1,1,2,0]}
    assert rawM[1][1]==u(sunit)
    initial=rc['initial_row_divisors'];st=rc['steps'];assert len(st)==1 and st[0]['pivot_row']==1
    elim={x['row']:x for x in st[0]['eliminations']}
    colmap=[0,2,3,4,5]
    rowunits=[];mapping=[];literal=[];compact=[]
    for original in rows:
        label=data['source25_labels'][original];assert label[0]=='ZW_genuine'
        genuine=label[1];rawindex=keep.index(genuine)
        assert rawindex in rc['remaining_rows']
        lambdap=decp(data['genuine16_source_clearing_units'][genuine]);lambdaunit=factorunit(lambdap)
        if rawindex in elim:
            tau=mulunit(mulunit(initial[1],initial[rawindex]),elim[rawindex]['row_divisor_after'])
        else:
            assert not rawM[rawindex][1]
            tau=mulunit(sunit,initial[rawindex])
        rho=mulunit(mulunit(lambdaunit,tau),sunit,True)
        rowunits.append(rho);mapping.append([original,genuine,rawindex])
        lr=[];nr=[]
        for jj,j in enumerate(columns):
            value=dec(oc['matrix25x32'][original][j]);expected=u(rho)*K(N[rawindex][colmap[jj]])
            if jj==4:expected*=u(gammaunit)
            assert value==expected
            lr.append(oc['matrix25x32'][original][j]);nr.append(N[rawindex][colmap[jj]])
        literal.append(lr);compact.append(nr)
    assert sum(bool(q[0])for row in literal for q in row)==43
    assert sum(len(q[0])for row in literal for q in row)==13854
    cert={'final14_source_labels':rows,'final5_heavy_target_labels':columns,
          'final_target_positions':[fc['new21_target_labels'].index(j)for j in columns],
          'all70_literal_external_root_traces':traces,
          'old_checkpoint70_rationals':literal,
          'actual_genuine_raw_checkpoint_index_map':mapping,
          'known_row_units_rho14':rowunits,'known_column_units':['1','1','1','1',gammaunit],
          'scaling_identity':'FinalH14x5=diag(rho)*RawCompact14x5*diag(1,1,1,1,gamma). All70 identities exact.',
          'raw_compact14x5':[[enc(q)for q in row]for row in compact],
          'scale_derivation':'Let s=X1_x1,lambda_i sourceclearing; rawunitcleanup gives tau_i. Then rho_i=lambda_i*tau_i/s. The L coordinate correction cancels in s*Xi_h-Xi_x1*X1_h. Psi also has columnunitgamma.',
          'generic_rank5_only_on_known_extra_open':'Accepted rawHrank7 implies finalheavy5 rank5 overQQ(a,e,d); FULLU D(pQ) explicit5minor is nonzero. No allpoint rank5 inferred.',
          'previously_accepted_generic_sparse_pair':{'raw_rows':gc['sparse_pair_rows'],'raw_columns':gc['sparse_pair_columns'],'determinant_formula':gc['pairdet_formula']},
          'parameter_scope':'FULLU; generic heavy rank5 and explicit FULLU D(pQ) rank5 are accepted dependencies only. The unknown p,Q remain nonunits on FULLU.',
          'B7x3_frontier_owned_by_other_agent_not_recomputed':True,
          'concrete_next_strategy':'This exact knownunit preconditioner identifies the actual finalheavyblock with the already-reviewed compactraw14x5. Its sparsepair and7x3block may be studied under separately authorized bounds. No repeat named-entry Gaussian or assumption of full H7 globally.'}
    (out/'certificate.json').write_text(json.dumps(cert,separators=(',',':'))+'\n')
    result={'id':ID,'all70_literal_heavyroots_bind_unchanged_old_checkpoint':True,
            'all70_exact_knownunit_scaling_identities_passed':True,
            'expanded_checkpoint_terms':13854,'compact70_terms':sum(len(q)for row in compact for q in row),
            'compact70_nonzero':sum(bool(q)for row in compact for q in row),
            'compact70_max_entry_terms':max(len(q)for row in compact for q in row),
            'generic_final_H5_rank':5,'global_H5_rank5_not_claimed':True,
            'no_B7x3_content_minor_or_cost_work_duplicated':True,
            'no_POWER_expansion_newminor_Gaussian_GB_point_or_rank_probe':True,
            'no_global_H7_full57_full239_or_automatic_next_operation':True,
            'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},
            'certificate_sha256':sha(out/'certificate.json'),
            'seconds':time.monotonic()-began,'frozen':True}
    assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');signal.alarm(0)
    print(json.dumps(result,indent=2),flush=True)


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True)
    main(ap.parse_args().output.resolve())
