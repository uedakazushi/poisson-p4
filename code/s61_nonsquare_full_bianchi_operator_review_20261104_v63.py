#!/usr/bin/env python3
"""Independent core-polynomial source/Schouten/weight audit; no C122 expansion."""
import argparse
import collections
from fractions import Fraction as F
import itertools
import json
from pathlib import Path
import resource
import signal
import time

import ambient_B18_zero_linear_contact4_review_20261103_v63 as own

ROOT=Path(__file__).resolve().parents[1]
STEM='S61_NONSQUARE_FULL_BIANCHI_OPERATOR_20261104_v63_REVIEW'


def add(D,key,value):
    value=own.Poly(value)
    if not value.terms:return
    value=D.get(key,own.Poly(0))+value
    if value.terms:D[key]=value
    else:D.pop(key,None)


def core(encoded):return own.Poly({tuple(e):F(c) for e,c in encoded})


def parity(indices):
    return (-1)**sum(a>b for i,a in enumerate(indices) for b in indices[i+1:])


def main(out):
    start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000))
    signal.alarm(60);assert not out.exists();out.mkdir(parents=True)
    ap=ROOT/'reports/S61_NONSQUARE_FULL_BIANCHI_OPERATOR_20261104_v63/report.json'
    author=own.read(ap);assert author['all_full_raw_map_and_literal_identity_checks_passed']
    for n,h in author['input_sha256'].items():assert own.sha(ROOT/n)==h,n
    opath=ap.parent/'operator.json';checkpoint=ap.parent/'source-verified-checkpoint.json'
    assert own.sha(opath)==author['operator_sha256'] and own.sha(checkpoint)==author['source_checkpoint_sha256']
    op,cp=own.read(opath),own.read(checkpoint)
    for k in ['raw_Jacobi_maps','homogeneous_bivector_q_t_coefficients','full_homogeneous135_equations','full175_Schouten_rows','extra_tG6_row']:assert op[k]==cp[k]
    universal_review=ROOT/'reports/S61_BINARY_UNIFORM_CUBIC_REVIEW_20261101.json'
    ur=own.read(universal_review);assert ur['accepted'] and ur['exact_rerun_passed']
    for n,h in ur['source_sha256'].items():assert own.sha(ROOT/n)==h,n
    chartpath=ROOT/'reports/S61_BINARY_HYPERPLANE_WEIGHTED_CHART_20261102.json'
    oldpath=ROOT/'reports/S61_BINARY_DEEPER_TRIANGULAR_CHART_20261101.json'
    extrapath=ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json'
    certpath=ROOT/'legacy/v02/data/symplectic_slice_certificate.json'
    datapath=ROOT/'legacy/v02/data/symplectic_slice_equations.json'
    reducedpath=ROOT/'legacy/v02/data/symplectic_slice_reduced_equations.json'
    idpath=ROOT/'legacy/v02/data/symplectic_bianchi_identities.json'
    chart,old,extra,cert,data,reduced,ids=[own.read(p) for p in [chartpath,oldpath,extrapath,certpath,datapath,reducedpath,idpath]]
    old=old['linear'];variables=own.symbols(['a','e','d']);one=own.Poly(1)
    weights=reduced['parameter_weights'];two=[i for i,w in enumerate(weights) if w==2]
    beta={int(k[1:])-1:own.Poly(own.evaluate(v,variables)) for k,v in chart['beta'].items()}
    def encoded(poly):
        result={}
        for mon,n,d in poly:
            value=own.Poly(F(n,d));qm=[]
            for i in mon:
                if weights[i]==1:value*=beta[i]
                else:qm.append(two.index(i))
                if not value.terms:break
            add(result,tuple(sorted(qm)),value)
        return result
    # Reconstruct the affine bivector directly from immutable encoded150 data.
    phys=[{} for _ in range(150)]
    for pair in [(1,2),(3,4)]:phys[15*own.PI[pair]+own.MI[0,0]]={():one}
    for raw,i in zip(cert['raw_free_indices'],range(61)):
        phys[raw]={(two.index(i),):one} if weights[i]==2 else {():beta[i]}
    for row in data['eliminated_raw_variables']:phys[row['raw_index']]=encoded(row['polynomial'])
    assert all(len(mon)<=1 for p in phys for mon in p)
    P={pair:{} for pair in itertools.combinations(range(4),2)}
    def put_affine(pair,rawpair,extra_variable=None,sign=1):
        for m,mon in enumerate(own.MONS):
            spatial=tuple(mon.count(i) for i in range(1,5))
            if extra_variable is not None:spatial=tuple(a+int(i==extra_variable) for i,a in enumerate(spatial))
            for qm,c in phys[15*own.PI[rawpair]+m].items():add(P[pair],(qm[0] if qm else 45,spatial),sign*c)
    for i,j in P:
        put_affine((i,j),(i+1,j+1))
        put_affine((i,j),(0,j+1),i,-1)
        put_affine((i,j),(0,i+1),j,1)
    savedP={(i,j):{(v,tuple(m)):core(c) for v,m,c in row} for i,j,row in op['homogeneous_bivector_q_t_coefficients']}
    assert P==savedP
    allF=reduced['equations']+[extra['universal_new_equation']];Fh=[]
    for equation in allF:
        row={}
        for qm,c in encoded(equation).items():
            assert len(qm)<=2
            add(row,tuple(sorted((*qm,*([45]*(2-len(qm)))))),c)
        Fh.append(row)
    assert len(Fh)==135
    assert Fh==[{tuple(m):core(c) for m,c in row} for row in op['full_homogeneous135_equations']]
    (out/'progress.json').write_text(json.dumps({'stage':'immutable P and135 source bound','seconds':time.monotonic()-start})+'\n')
    pairs={}
    for i in range(4):
        for j in range(4):pairs[i,j]=P[i,j] if i<j else ({k:-v for k,v in P[j,i].items()} if i>j else {})
    triples=list(itertools.combinations(range(4),3));jac={tr:{} for tr in triples}
    for tr in triples:
        a,b,c=tr
        for i,j,k in [(a,b,c),(b,c,a),(c,a,b)]:
            for h in range(4):
                for (v,m),pc in pairs[i,h].items():
                    for (w,n),qc in pairs[j,k].items():
                        if n[h]:
                            spatial=tuple(x+y-int(l==h) for l,(x,y) in enumerate(zip(m,n)))
                            add(jac[tr],(spatial,tuple(sorted((v,w)))),n[h]*pc*qc)
    oldmap={g:i for i,g in enumerate(reduced['old_generator_indices'])}
    redundancy={r['eliminated_old_generator']:r for r in ids['generator_redundancy_identities']}
    def old_generator(j):
        if j in oldmap:return {oldmap[j]:one}
        r=redundancy[j];out={}
        for k,n,d in r['constant_coefficients_on_retained_old_generators']:add(out,oldmap[k],F(n,d))
        for k,u,n,d in r['linear_coefficients_on_weight4_old_generators']:
            assert weights[u]==1;add(out,oldmap[k],F(n,d)*beta[u])
        return out
    T={tr:{} for tr in triples};rawcounts=[]
    for entries,oldids,D in zip(op['raw_Jacobi_maps'],[op['old154_weight4_generator_indices'],op['old154_weight5_generator_indices']],[3,4]):
        assert len(entries)==(80 if D==3 else 140)
        seen=set()
        for rec in entries:
            tr,m=tuple(rec['triple']),tuple(rec['spatial_monomial']);assert sum(m)==D and (tr,m) not in seen;seen.add((tr,m))
            co={}
            for k,c in rec['old_basis_coordinates']:
                for g,v in old_generator(oldids[k]).items():add(co,g,F(c)*v)
            assert co=={g:core(v) for g,v in rec['reduced134_coordinates']}
            predicted={}
            for g,v in co.items():
                assert g<134
                add(T[tr],(g,m),v)
                for qm,c in Fh[g].items():add(predicted,qm,v*c)
            actual={qm:v for (n,qm),v in jac[tr].items() if n==m}
            assert predicted==actual
        assert len(seen)==len(triples)*(20 if D==3 else 35);rawcounts.append(len(seen))
    assert all(sum(m) in [3,4] for row in jac.values() for m,qm in row)
    (out/'progress.json').write_text(json.dumps({'stage':'direct full core Jacobi matches220 raw maps','seconds':time.monotonic()-start})+'\n')
    # Alternating/unshuffle component calculation of [hP,T], with free source symbols.
    operator={}
    def term(m,g,v,c):
        if m not in operator:operator[m]={}
        add(operator[m],(g,v),c)
    for i in range(4):
        others=tuple(j for j in range(4) if j!=i)
        for h in range(4):
            for (v,pm),pc in pairs[i,h].items():
                for (g,tm),tc in T[others].items():
                    if tm[h]:term(tuple(x+y-int(j==h) for j,(x,y) in enumerate(zip(pm,tm))),g,v,(-1)**i*tm[h]*pc*tc)
    for i,k in itertools.combinations(range(4),2):
        others=tuple(j for j in range(4) if j not in [i,k])
        for h in range(4):
            tr=(h,*others)
            if len(set(tr))<3:continue
            sign=(-1)**(i+k)*parity(tr)
            for (v,pm),pc in P[i,k].items():
                if not pm[h]:continue
                for (g,tm),tc in T[tuple(sorted(tr))].items():term(tuple(x+y-int(j==h) for j,(x,y) in enumerate(zip(pm,tm))),g,v,sign*pm[h]*pc*tc)
    operator={m:row for m,row in operator.items() if row}
    records=op['full175_Schouten_rows'];assert len(records)==175
    nonzero=collections.Counter()
    for r in records:
        D,m=r['spatial_degree'],tuple(r['spatial_monomial']);assert sum(m)==D
        actual=operator.get(m,{})
        assert actual=={(g,v):core(c) for g,v,c in r['coefficients']}
        nonzero[D]+=bool(actual)
    assert nonzero=={4:35,5:56,6:0}
    assert all(sum(m) in [4,5,6] for m in operator)
    # All175 literal source identities now follow from this independent Jac(P)
    # correspondence and graded Jacobi [P,[P,P]]=0, without cubic re-expansion.
    for monindices,stored in zip(extra['affine_quartic_monomial_order'],extra['canonical_35_Bianchi_vectors']):
        m=tuple(monindices.count(i) for i in range(4));wanted={}
        for index,c in stored:v,g=divmod(index,70);add(wanted,(g,v),F(c))
        assert wanted=={key:c for key,c in operator.get(m,{}).items() if key[1]<45}
    erow={(134,45):one}
    for index,c in extra['selected_non_Bianchi_syzygy']:
        v,g=divmod(index,70);add(erow,(g,v),-F(c))
    assert erow=={(g,v):core(c) for g,v,c in op['extra_tG6_row']}
    literal={}
    for (g,v),c in erow.items():
        for qm,f in Fh[g].items():add(literal,tuple(sorted((v,*qm))),c*f)
    assert not literal
    fixed,cols=old['fixed'],old['cols'];cw=chart['local_coordinate_weights'];rw=chart['equation_weights']
    assert rw==op['equation_weights'] and cw==op['local_coordinate_weights'] and len(rw)==135
    rows=[operator.get(tuple(r['spatial_monomial']),{}) for r in records]+[erow]
    pieces=[];counts=collections.Counter()
    for idx,row in enumerate(rows):
        groups={}
        for (g,v),c in row.items():
            local=-1 if v==fixed else (44 if v==45 else cols.index(v))
            wt=rw[g]+(0 if local==-1 else cw[local]);groups.setdefault(wt,{})
            add(groups[wt],(g,local),c)
        for w,group in sorted(groups.items()):
            D=records[idx]['spatial_degree'] if idx<175 else 'extra_tG6'
            pieces.append(dict(operator_row=idx,spatial_degree=D,weight=w,coefficients=group));counts[str(D),w]+=1
    assert len(pieces)==len(op['local_weight_pieces'])
    for a,b in zip(pieces,op['local_weight_pieces']):
        assert a['operator_row']==b['operator_row'] and a['spatial_degree']==b['spatial_degree'] and a['weight']==b['weight']
        assert a['coefficients']=={(g,v):core(c) for g,v,c in b['coefficients']}
    assert [[D,w,n] for (D,w),n in sorted(counts.items())]==op['local_weight_piece_counts']==author['actual_local_weight_piece_counts']
    zero=[i for i,w in enumerate(cw) if w==0];assert zero==[19,23]
    assert all(chart['T34'][i][j]=='0' for i in zero for j in chart['normal34'])
    assert all(v=='0' for i in zero for v in chart['U34'][i])
    quartic=[i for i,w in enumerate(rw) if w==4 and i not in chart['selected34_rows']]
    assert quartic==op['genuine_quartic_original_rows'] and len(quartic)==33
    block=[];labels=[];discard=0
    for i,p in enumerate(pieces):
        if p['weight']!=4:continue
        row={}
        for (g,v),c in p['coefficients'].items():
            if g in quartic:
                if v==-1:row[g]=c
                else:assert v in zero;discard+=1
        block.append([row.get(g,own.Poly(0)) for g in quartic]);labels.append(i)
    saved=[[core(c) for c in row] for row in op['weight4_quartic_constant_block']]
    assert block==saved and labels==op['weight4_quartic_constant_block_row_piece_indices']
    assert discard==132==op['zero_weight_nuisance_terms_dropped_after_exact_graph_zero_proof']
    assert len(block)==15 and sum(any(c.terms for c in row) for row in block)==11
    assert sum(len(c.terms) for row in block for c in row)==88
    paths=[Path(__file__),Path(own.__file__),ROOT/f'research/notes/{STEM}.md',ap,opath,checkpoint,universal_review,chartpath,oldpath,extrapath,certpath,datapath,reducedpath,idpath]
    result=dict(id=STEM,all_passed=True,
        independent_immutable150_reconstruction_and_full_affine_P_q_t_core_binding=True,
        independent_all135_q_t_homogeneous_equations_including_G6=True,
        old_universal_full61_Jacobi_map_gate_accepted_not_rerun=True,
        independent_direct_full_core_Jacobi_reconstructs_all80_and140_correspondences=True,
        independent_old154_to134_redundancy_maps_exact=True,
        independent_Schouten_operator_matches_all175_serialized_rows=True,
        all175_literal_syzygies_proved_by_verified_Jacobi_correspondence_and_graded_Jacobi=True,
        cubic_polynomial_substitution_all175_not_duplicated=True,
        old35_canonical_q_parts_and_literal_tG6_identity_verified=True,
        spatial_degree_slots=[35,56,84],nonzero_spatial_degree_rows=[35,56,0],
        actual_weight_piece_counts=op['local_weight_piece_counts'],
        exact_zero_normal_graph_coordinates=[19,23],discarded_only_block_terms=132,
        weight4_constant_quartic_block_shape=[15,33],nonzero_block_rows=11,nonzero_core_terms=88,
        block_has_not_been_pivoted_or_transported_to_normal_DAG_or319=True,
        no80_directions_rankclaim_fullquartic_or_C122_job=True,
        verification_kind='Independent sparse QQ[a,e,d] source/Jacobi/Schouten/weight reconstruction; no universal61 replay, C122, GB, point, field or rank job.',
        input_sha256={str(p.relative_to(ROOT)):own.sha(p) for p in paths},seconds=time.monotonic()-start)
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');signal.alarm(0)
    print(json.dumps({k:result[k] for k in ['all_passed','nonzero_spatial_degree_rows','weight4_constant_quartic_block_shape','nonzero_block_rows','seconds']},indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);args=p.parse_args()
    try:main(args.output)
    except BaseException as error:
        if args.output.exists() and not (args.output/'report.json').exists():
            (args.output/'incomplete.json').write_text(json.dumps({'all_passed':False,'error':type(error).__name__,'detail':str(error),'no_automatic_larger_retry':True},indent=2)+'\n')
        raise
