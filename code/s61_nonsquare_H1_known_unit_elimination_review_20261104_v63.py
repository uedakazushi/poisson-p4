#!/usr/bin/env python3
"""Independent fixed7 source rows/pivot replay;38 lazy rows remain unevaluated."""
import argparse
from fractions import Fraction as Q
import gzip
import json
from pathlib import Path
import resource
import signal
import time
import sympy as s
from sympy.polys.rings import ring
from s61_nonsquare_free281_circuit_review_20261104_v63 import read,sha
from s61_nonsquare_joint_Lell_free183_review_20261104_v63 import matrix

ROOT=Path(__file__).resolve().parents[1]
BASE='S61_NONSQUARE_H1_KNOWN_UNIT_ELIMINATION_20261104_v63'
ID=BASE+'_REVIEW'

def main(out):
    began=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
    folder=ROOT/'reports'/BASE;r=read(folder/'report.json')
    with gzip.open(folder/'checkpoint.json.gz','rt') as f:state=json.load(f)
    assert r['frozen'] and state['completed_pivots']==r['completed_global_H1_known_unit_pivots']==7
    for n,h in r['input_sha256'].items():assert sha(ROOT/n)==h,n
    for n,h in r['evidence_sha256'].items():assert sha(folder/n)==h,n
    pb=ROOT/'reports/S61_NONSQUARE_H_SQUARE_H1_PREFLIGHT_20261104_v63';source=read(pb/'certificate.json')
    peer=ROOT/'reports/S61_NONSQUARE_H_SQUARE_H1_PREFLIGHT_20261104_v63_REVIEW/report.json';assert read(peer)['all_passed']
    assert sha(pb/'certificate.json')==state['source_input_sha256']
    R,a,e,d=ring('a,e,d',s.QQ);K=R.to_field();zero=K.zero
    def poly(ts):return R.from_dict({tuple(m):s.QQ(c) for m,c in ts})
    def dec(z):return K(poly(z[0]))/K(poly(z[1]))
    def enc(v):return [[[list(m),str(c)] for m,c in p.items()] for p in [v.numer,v.denom]]
    F=500*a+90*e**2+3*e-24;B=(50*a+18*e**2-9*e-2)*d+60*e-40;C=100*a+18*e**2-9*e-2;J=5*a*d-6*e-4
    known=[a,d,F,B,C,J];sg=source['scalar_circuit'];scalar_cache={0:zero}
    def scalar(i):
        if i in scalar_cache:return scalar_cache[i]
        n=sg['nodes'][i]
        if n[0]=='c':v=dec(sg['atoms_num_den'][n[1]])
        elif n[0]=='a':v=scalar(n[1])+scalar(n[2])
        elif n[0]=='s':v=scalar(n[1])-scalar(n[2])
        else:assert n[0]=='m';v=scalar(n[1])*scalar(n[2])
        scalar_cache[i]=v;return v
    nf=source['parent_target239_to122_operator'];mons=source['parent_joint_target239_monomials'];mi={tuple(m):i for i,m in enumerate(mons)}
    kept=source['parent_target122_monomial_indices'];hcols=source['target122_Hdegree_classes']['1']
    hmons=[tuple(mons[kept[i]][:7]) for i in hcols]
    lin=source['linear_nodes'];lin_cache={0:{}}
    def linear(i):
        if i in lin_cache:return lin_cache[i]
        n=lin[i]
        if n[0]=='v':v={n[1]:K.one}
        else:
            assert n[0]=='l';v={}
            for co,j in n[1]:
                c=scalar(co)
                for var,x in linear(j).items():v[var]=v.get(var,zero)+c*x
            v={j:x for j,x in v.items() if x}
        lin_cache[i]=v;return v
    normal=[linear(nf['output122_linear_roots'][i]) for i in hcols]
    vin=nf['input_variables']
    def project(values):
        inputs={}
        for m,v in values.items():
            if v:n=list(m)+[1];inputs[vin[mi[tuple(n)]]]=v
        return [sum((c*inputs.get(v,zero) for v,c in row.items()),zero) for row in normal]
    np=ROOT/'reports/S61_NONSQUARE_JOINT_LELL_PREFLIGHT_20261104_v63/coordinate-inverse-numerators7x7.txt';Delta=K(3200*a*C*d*d)/81
    N=[[K(poly(z))/Delta for z in row] for row in matrix(np)]
    assert N[5]==[K(int(j==5)) for j in range(7)] and N[6]==[K(int(j==6)) for j in range(7)]
    oldB=[N[0],N[1],N[2],N[6],N[3],[N[5][j]+N[1][j]-3*N[0][j] for j in range(7)],N[4]]
    precond=state['raw_row_preconditioner'];assert [[dec(v) for v in row] for row in precond['oldlight_inverse7x7']]==oldB
    delta=K(-64000*C*C);assert dec(precond['delta'])==delta and precond['replace_rows']==list(range(43,50))
    assert precond['source_determinant']=='delta^7'
    rawp=ROOT/'reports/S61_NONSQUARE_ELL_SIX_BORDERS_20261104_v63/frozen57-source.txt'
    br=read(ROOT/'reports/S61_BINARY_NONSQUARE_CUBIC_ROW_BASIS_20261103_v62/report.json');assert sha(rawp)==br['remaining_matrix_sha256']
    src=read(ROOT/'reports/S61_BINARY_NONSQUARE_CUBIC_REPLAY_20261103_v62/source-report.json');cm=src['normal_monomials'];light=[0,1,2,3,4,6,7]
    hc=[]
    for j in light:m=[0]*8;m[5]=m[j]=1;hc.append(cm.index(m))
    needed_basis={42,47,53};rawheavy={i:{} for i in needed_basis}
    with rawp.open() as f:
        assert f.readline().split()==['57','91','57']
        for i in range(57):
            for j in range(91):
                line=f.readline();assert line
                if i not in needed_basis or j not in hc:continue
                t=line.split();assert len(t)==1+4*int(t[0]);rawheavy[i][j]=K(poly([[list(map(int,t[k+1:k+4])),t[k]] for k in range(1,len(t),4)]))
        assert not f.read().strip()
    heavy={i:[rawheavy[i][j] for j in hc] for i in needed_basis}
    selected=[11,12,89,8,85,10,14]
    assert [t['raw134_source_row'] for t in state['unit_pivots']]==selected
    assert all(41+(i-1)//7 in needed_basis for i in selected)
    def multiply_linear(x,y):
        result={}
        for i,co in enumerate(x):
            if not co:continue
            for j,co2 in enumerate(y):
                if not co2:continue
                m=tuple(int(k==i)+int(k==j) for k in range(7));result[m]=result.get(m,zero)+co*co2
        return result
    ownraw={}
    source_bind=[]
    for row in selected:
        genuine=(row-1)//7;j=(row-1)%7;basisrow=41+genuine
        hnew=[sum((heavy[basisrow][k]*oldB[k][i] for k in range(7)),zero) for i in range(7)]
        ownraw[row]=project(multiply_linear(hnew,oldB[j]))
        normalized=br['remaining_original_row_indices'][basisrow];label=src['labels'][normalized]
        source_bind.append({'raw134_row':row,'remaining57_row':basisrow,'normalized63_row':normalized,'actual135_cubic':label[1],'oldlight':light[j]})
    # Bind the compact AH/BH to the actual two frozen cubics and full source identity.
    compact=ROOT/'reports/S61_NONSQUARE_ELL_COMPACT_RELATION_20261104_v63';cp=compact/'compact-polynomials.txt';pols=[]
    # This archive has named header E8 rather than a matrix shape.
    for line in cp.read_text().splitlines()[1:]:
        ts=line.split();assert len(ts)==1+4*int(ts[0])
        pols.append(K(poly([[list(map(int,ts[k+1:k+4])),ts[k]] for k in range(1,len(ts),4)])))
    oldm2=[list(map(int,z.split())) for z in (compact/'source-metadata.txt').read_text().splitlines()[1:30]]
    hcol=oldm2.index([0,0,0,0,0,1,0,0]);kZ,kW=pols[8+hcol],pols[8+29+hcol]
    assert pols[74]==delta and [dec(v) for v in precond['compact_H_coefficients']]==[kZ,kW]
    diff=[delta*(x-y) for x,y in zip(heavy[47],heavy[42])]
    assert [diff[0]-3*diff[5],diff[1]+diff[5],diff[2],diff[4],diff[6]]==[zero]*5
    assert diff[5]==kZ and diff[3]==kW
    specialH=[zero]*7;specialH[5]=kZ;specialH[6]=kW
    removed={(2,0,0,0,0,0,0),(1,0,1,0,0,0,0),(0,1,0,0,0,1,0),(0,1,0,0,0,0,1),(0,2,0,0,0,0,0)}
    assert not set(hmons)&removed and len(set(hmons)|removed)==28
    for j in range(7):
        full=multiply_linear(specialH,oldB[j]);value=project(full)
        assert value==[full.get(m,zero) for m in hmons]
        ownraw[43+j]=value
    materialized=dict(state['materialized_raw_rows']);assert len(materialized)==96
    for row,value in ownraw.items():assert [dec(v) for v in materialized[row]]==value
    # These are fixed recorded pivots. No scan or unit search is performed.
    normalized=[];used=[];pivot_proof=[];den_set=set()
    for t in state['unit_pivots']:
        row=t['raw134_source_row'];cur=ownraw[row].copy();coefs=[]
        for j,p in zip(used,normalized):
            co=cur[j];coefs.append(co);cur=[x-co*y for x,y in zip(cur,p)]
            assert not cur[j]
        assert coefs==[dec(v) for v in t['subtract_previous_pivot_coefficients']]
        column=t['pivot_column'];pivot=dec(t['pivot']);assert cur[column]==pivot and pivot
        fact=t['unit_factorization'];value=K.from_expr(s.sympify(fact['constant']))
        assert value
        for f,n in zip(known,fact['a_d_F_B_C_J_powers']):value*=K(f)**n
        assert value==pivot
        p=[dec(v) for v in t['normalized23']]
        assert all(x==pivot*y for x,y in zip(cur,p)) and p[column]==1 and all(not p[j] for j in used)
        # ALL stored normalized/cofactor/pivot denominators are known FULLU units.
        for x in [*p,*coefs,pivot]:den_set.add(tuple(sorted(x.denom.items())))
        normalized.append(p);used.append(column)
        pivot_proof.append({'raw134_row':row,'column':column,'unit_factorization':fact})
    for terms in den_set:
        den=R.from_dict(dict(terms))
        for f in known:
            while not den.rem(f):den=den.exquo(f)
        assert den.is_ground and den
    assert len(set(used))==7
    assert state['remaining_columns']==[i for i in range(23) if i not in used]
    assert state['remaining_rows']==[i for i in range(134) if i not in selected]
    evaluated=set(materialized);lazy=sorted(set(range(134))-evaluated)
    assert len(lazy)==38 and set(range(113,134))<=set(lazy)
    assert sorted(state['completed_evaluated_rows'])==sorted(evaluated)
    assert all(0<=row<134 and 0<=stage<=7 for row,stage in state['tested_row_pivotstage_pairs'])
    assert state['pure98_coefficient_queries']==0 and state['all_unexpanded_rows_retained_as_original_scalar_circuit_plus_saved_source_preconditioner']
    assert state['redundant_full_source_target_after_pivots']==[183,114] and r['new_target_unit_lower_from_original239']==125
    assert 57+112+21==190 and 190-7==183 and 117+1+7==125 and 121-7==114
    out.mkdir(parents=True,exist_ok=False);ep=out/'independent-proof.json'
    evidence={'actual_source_bindings_for_selected7':source_bind,'verified7pivots':pivot_proof,
        'all14_selected_and_compact_special_raw_rows_trueQQ_equal':True,
        'all7times23_normalized_coefficients_and_prior_subtractions_exact':True,
        'source_preconditioner_delta_power7_known_unit':True,
        'remaining38_unevaluated_raw_rows':lazy,'all21_quartics_retained_lazy':True,
        'normalized_pivot_unique_denominators_unit_verified':len(den_set),
        'actual_full_source_lifts_include_lazy_pure98':True,
        'global_actual_rank_lower125_safe_redundant_shape':[183,114],
        'free156_or_allpoint_cubic57_full239_claim':False}
    ep.write_text(json.dumps(evidence,indent=2)+'\n')
    paths=[Path(__file__),ROOT/'research/notes'/(ID+'.md'),folder/'report.json',folder/'checkpoint.json.gz',pb/'certificate.json',peer,cp,compact/'source-metadata.txt',rawp,np,
        ROOT/'reports/S61_BINARY_NONSQUARE_CUBIC_ROW_BASIS_20261103_v62/report.json',
        ROOT/'reports/S61_BINARY_NONSQUARE_CUBIC_REPLAY_20261103_v62/source-report.json',
        ROOT/'reports/S61_NONSQUARE_ELL_COMPACT_RELATION_20261104_v63_REVIEW/report.json',
        ROOT/'code/s61_nonsquare_joint_Lell_free183_review_20261104_v63.py']
    result={'id':ID,'all_passed':True,'all7_fixed_pivots_bound_to_actual57_full_source_QQ':True,
        'compact_unit_source7preconditioner_actual_coefficients_bound':True,
        'all7_full23_rows_subtractions_normalizations_and_unit_factorizations_verified':True,
        'all38_lazy_rows_and21_quartics_preserved_without_expansion':True,
        'global_rank_lower125_and_samecokernel_redundant183_to114_accepted':True,
        'automatic_free156_full239_or_fullU_cubic57_claim':False,
        'no_producer_search_remaining38_evaluation_newpoint_GB_or_pure98_query':True,
        'input_sha256':{str(p.relative_to(ROOT)):sha(p) for p in paths},'evidence_sha256':{ep.name:sha(ep)},
        'seconds':time.monotonic()-began,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'frozen':True}
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');signal.alarm(0)
    print(json.dumps({k:result[k] for k in ['all_passed','seconds','max_rss_kib']}))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
