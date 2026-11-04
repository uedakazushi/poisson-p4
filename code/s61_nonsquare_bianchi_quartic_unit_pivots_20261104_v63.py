#!/usr/bin/env python3
"""Bounded known-unit pivots and exact source transport for the full Bianchi block.

The source135 identities are preserved as finite substitution circuits.  No
239-column quartic coefficient matrix, new point, field or GB is computed.
"""
import argparse
import hashlib
import json
from pathlib import Path
import resource
import signal
import tarfile
import time
import sympy as s
from s61_nonsquare_oracle_frame_adapter_20261103_v63 import FrameAdapter

ROOT=Path(__file__).resolve().parents[1]
ID='S61_NONSQUARE_BIANCHI_QUARTIC_UNIT_PIVOTS_20261104_v63'
UP=ROOT/'reports/S61_NONSQUARE_FULL_BIANCHI_OPERATOR_20261104_v63'

def read(p):return json.loads(p.read_text())
def sha(p):
    h=hashlib.sha256()
    with p.open('rb')as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()

def main(out,resume=False):
    start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(53 if resume else 60)
    out.mkdir(parents=True,exist_ok=resume)
    def progress(stage,**kw):
        p=dict(stage=stage,seconds=time.monotonic()-start,**kw)
        (out/'progress.json').write_text(json.dumps(p,indent=2)+'\n');print(json.dumps(p),flush=True)
    operator_path=UP/'operator.json';upreport=read(UP/'report.json')
    assert sha(operator_path)==upreport['operator_sha256']
    assert upreport['all175_literal135_syzygies_verified']and upreport['homogenized_t_G6_identity_verified']
    O=read(operator_path);f=FrameAdapter.load();o=f.oracle;K,R=o.K,o.R;initial_nodes=len(o.nodes)
    enc=lambda p:[[list(ex),str(c)]for ex,c in p.items()]
    frac=lambda ts:K(R.from_dict({tuple(ex):s.QQ(c)for ex,c in ts}))
    expr=lambda x:K.from_expr(s.sympify(x))
    def unit(c,refuse=False):
        c=K(c)
        if not c:return None if refuse else (_ for _ in ()).throw(AssertionError('zero unit'))
        powers=[];constants=[]
        for p in[c.numer,c.denom]:
            row=[]
            for q in o.known:
                n=0
                while not p.rem(q):p=p.exquo(q);n+=1
                row.append(n)
            if not(p.is_ground and p):
                if refuse:return None
                raise AssertionError('unknown factor '+str(c))
            constants.append(p);powers.append(row)
        return dict(constant=str(K(constants[0])/K(constants[1])),a_F_B_C_powers=[a-b for a,b in zip(*powers)])
    cols=O['local_other_q_order'];fixed=O['local_fixed_q_index'];cw=O['local_coordinate_weights'];rw=O['equation_weights']
    pieces=O['weight4_quartic_constant_block_row_piece_indices'];quartics=O['genuine_quartic_original_rows']
    assert len(pieces)==15 and len(quartics)==33
    prior=None
    if resume:
        prior=read(out/'first-attempt-frontier.json')
        assert prior['operator_sha256']==sha(operator_path)
        assert prior['block_sha256']==sha(out/'block-pivots.json')
        assert prior['full135_binding_and15_literal_identities_passed']and prior['rank11_and_remaining_zero']
        literal_checks=15
    else:
        # Bind the published full source135 literally to the oracle's J/B atoms.
        cols=O['local_other_q_order'];fixed=O['local_fixed_q_index'];cw=O['local_coordinate_weights'];rw=O['equation_weights']
        assert cols==f.data['other_original_coordinate_order']and fixed==f.data['fixed_original_coordinate']
        assert cw==f.data['local_coordinate_weights']and rw==f.data['equation_weights']
        def add(t,k,v):
            v=t.get(k,K.zero)+v
            if v:t[k]=v
            else:t.pop(k,None)
        def loc(v):return -1 if v==fixed else(44 if v==45 else cols.index(v))
        Fh=[[(tuple(mon),frac(c))for mon,c in row]for row in O['full_homogeneous135_equations']]
        for i,row in enumerate(Fh):
            actual={}
            for mon,c in row:add(actual,tuple(sorted(loc(v)for v in mon if v!=fixed)),c)
            expected={}
            for j,k in f.data['J_sparse_atom_ids'][i]:add(expected,(j,),o.atoms[k])
            for j,k,c in f.data['B_sparse_atom_ids'][i]:add(expected,tuple(sorted((j,k))),o.atoms[c])
            assert actual==expected,(i,'source135 binding')
            for mon,c in actual.items():
                if any(j in[19,23]for j in mon):continue
                assert sum(cw[j]for j in mon)==rw[i],(i,mon)
        assert all(not f.data[key][j]for key in ['linear_graph_roots','local_Z2_roots','local_Z3_roots','local_Z4_roots']for j in[19,23])
        pieces=O['weight4_quartic_constant_block_row_piece_indices'];quartics=O['genuine_quartic_original_rows']
        assert len(pieces)==15 and len(quartics)==33
        literal_checks=0
        for j in pieces:
            p=O['local_weight_pieces'][j];raw=O['full175_Schouten_rows'][p['operator_row']];check={};local={}
            for g,v,c0 in raw['coefficients']:
                c=frac(c0)
                for mon,b in Fh[g]:add(check,tuple(sorted((*mon,v))),c*b)
                l=loc(v)
                if rw[g]+(0 if l==-1 else cw[l])==4:add(local,(g,l),c)
            assert not check
            assert local=={(g,l):frac(c)for g,l,c in p['coefficients']}
            literal_checks+=1
        progress('full135 binding and15 literal source identities passed')
    # Only known units may be inverted; record every elementary row operation.
    block=[[frac(c)for c in row]for row in O['weight4_quartic_constant_block']]
    if resume:
        pivot_result=read(out/'block-pivots.json');r=pivot_result['rank_by_global_known_unit_pivots']
        U=[[expr(c)for c in row]for row in pivot_result['elementary_transform']]
        B=[[expr(c)for c in row]for row in pivot_result['transformed_block']]
        pivot_cols=pivot_result['pivot_columns'];all_remaining_zero=pivot_result['all_remaining_rows_exactly_zero']
        assert [[str(c)for c in row]for row in block]==pivot_result['original_block']
        assert r==11 and all_remaining_zero
        progress('resume from frozen11-pivot checkpoint; no source or pivot rerun',rank=r)
    else:
        B=[row[:]for row in block];U=[[K.one if i==j else K.zero for j in range(15)]for i in range(15)]
        steps=[];pivot_cols=[];r=0
        while r<15:
            candidates=[]
            for i in range(r,15):
                for j,c in enumerate(B[i]):
                    if j in pivot_cols or not c:continue
                    proof=unit(c,True)
                    if proof is not None:candidates.append((len(c.numer)+len(c.denom),i,j,proof))
            if not candidates:break
            _,i,j,proof=min(candidates,key=lambda t:t[:3]);p=B[i][j]
            if i!=r:B[r],B[i]=B[i],B[r];U[r],U[i]=U[i],U[r]
            B[r]=[c/p for c in B[r]];U[r]=[c/p for c in U[r]]
            elim=[]
            for k in range(15):
                if k==r or not B[k][j]:continue
                c=B[k][j];elim.append([k,str(c)])
                B[k]=[a-c*b for a,b in zip(B[k],B[r])];U[k]=[a-c*b for a,b in zip(U[k],U[r])]
            steps.append(dict(pivot_row=r,source_row_before_swap=i,quartic_column=j,original_quartic_row=quartics[j],pivot=str(p),unit=proof,eliminations=elim))
            pivot_cols.append(j);r+=1
        assert [[sum((U[i][k]*block[k][j]for k in range(15)),K.zero)for j in range(33)]for i in range(15)]==B
        for i,j in enumerate(pivot_cols):assert[B[k][j]for k in range(r)]==[K.one if k==i else K.zero for k in range(r)]
        all_remaining_zero=all(not c for row in B[r:]for c in row)
        pivot_result=dict(rank_by_global_known_unit_pivots=r,pivot_columns=pivot_cols,steps=steps,
                          all_remaining_rows_exactly_zero=all_remaining_zero,
                          original_block=[[str(c)for c in row]for row in block],
                          elementary_transform=[[str(c)for c in row]for row in U],
                          transformed_block=[[str(c)for c in row]for row in B])
        (out/'block-pivots.json').write_text(json.dumps(pivot_result,indent=2)+'\n')
        progress('known-unit block pivots completed',rank=r,remaining_zero=all_remaining_zero)
    # Extend the already accepted exact row-covariance adapter to the five
    # deleted cubic rows. Short anchors only FIX the scalar; covariance is the
    # accepted full63 physical-source reconstruction, not sample evidence.
    source_path=ROOT/'reports/S61_BINARY_NONSQUARE_CUBIC_REPLAY_20261103_v62/source-report.json';source=read(source_path)
    row_units={i:c for i,c in zip(f.frame['basis57_original_row_indices'],f.basis_units)if i>=42}
    norm_path=ROOT/'reports/S61_BINARY_NONSQUARE_CUBIC_ROW_BASIS_20261103_v62/normalization.txt'
    def native_poly(line):
        t=line.split();n=int(t[0]);assert len(t)==1+4*n
        return K(R.from_dict({tuple(map(int,t[2+4*k:5+4*k])):s.QQ(t[1+4*k])for k in range(n)}))
    nl=norm_path.read_text().splitlines();assert nl[0]=='63';norms=[native_poly(line)for line in nl[1:]]
    missing=[i for i in range(42,63)if i not in row_units];assert missing==[48,55,57,59,61]
    archive=ROOT/'reports/S61_BINARY_NONSQUARE_CUBIC_REPLAY_20261103_v62/matrices.tar.gz';anchors={};h=hashlib.sha256()
    with tarfile.open(archive)as tf:
        stream=tf.extractfile('full');line=stream.readline();h.update(line);assert line.decode().strip()=='63 91 57'
        for i in range(63):
            choices=[]
            for k,mon in enumerate(source['normal_monomials']):
                line=stream.readline();h.update(line)
                if i in missing:
                    n=int(line.split(maxsplit=1)[0])
                    if n:choices.append((0 if mon[5]else 1,n,k,line.decode()))
            if i in missing:anchors[i]=min(choices,key=lambda x:x[:3])
        assert not stream.read()
    assert h.hexdigest()==source['full_matrix_sha256']
    # Immutable writer sets the target metadata57; its Schur parser reads r
    # but checks only dimensions63/91. It is not an observed rank premise.
    header_writer=ROOT/'code/s61_binary_hyperplane_cubic_matrix_20261102.py'
    header_parser=ROOT/'code/s61_binary_hyperplane_cubic_unit_schur_20261102.cpp'
    assert "f.write('63 91 57"+chr(92)+"n')" in header_writer.read_text()
    assert 'in>>n>>m>>r;assert(n==63&&m==91)' in header_parser.read_text()

    supplemental=[]
    for i in missing:
        _,_,k,line=anchors[i];mon=tuple(source['normal_monomials'][k]);g=source['labels'][i][1]
        saved=native_poly(line)/norms[i]
        value=f.core_value(f.equation_coefficient(g,mon))/saved;proof=unit(value);row_units[i]=value
        supplemental.append(dict(original63_row=i,physical135_row=g,anchor_monomial=list(mon),unit=str(value),unit_factorization=proof,
                                 normalized_full63_anchor=str(saved),actual_physical_anchor=str(value*saved)))
    assert len(row_units)==21
    progress('all21 cubic physical row units fixed',supplemental5=supplemental)
    cert_path=ROOT/'reports/S61_NONSQUARE_QUARTIC_319_ROWS_20261103_v63/certificate.json';C=read(cert_path)
    rows=C['source_rows'];assert len(rows)==354
    q_lookup={(row[1],tuple(row[2])):i for i,row in enumerate(rows)if row[0]=='quadric'}
    cubic_lookup={(row[2],row[3]):i for i,row in enumerate(rows)if row[0]=='cubic_normalized'}
    qrows={g:i for i,g in enumerate(f.frame['quadric_original_rows'])}
    crows={source['labels'][i][1]:i for i in range(42,63)}
    q4lookup={row[1]:i for i,row in enumerate(rows)if row[0]=='quartic'}
    assert list(q4lookup)==quartics
    weight2=sorted({mon for i,mon in q_lookup});assert len(weight2)==29
    light=[0,1,2,3,4,6,7];one=[0]*8;zeros=(0,)*8
    # Build only weight1/2 graph coefficients: later terms have too high weight.
    local_graph={}
    for j,w in enumerate(cw):
        if w==1:local_graph[j]=f.data['linear_graph_roots'][j]
        elif w==2:local_graph[j]=o.lin([(1,f.data['linear_graph_roots'][j]),(1,f.data['local_Z2_roots'][j])])
    cache={}
    def graph_coefficient(j,mon):
        key=(j,mon)
        if key not in cache:cache[key]=o.lin([(f.scalar(mon),o.coefficient(local_graph[j],mon))])
        return cache[key]
    selected=set(f.data['selected34_rows']);lifts=[];vectors354=[];vectors319=[];dropsteps=[]
    for bi,piece_index in enumerate(pieces):
        p=O['local_weight_pieces'][piece_index];terms={};skipped=[]
        for g,j,co in p['coefficients']:
            c=frac(co)
            if g in selected:skipped.append([g,j]);continue
            w=rw[g];assert w in[2,3,4]
            if w==4:
                if j in[19,23]:continue
                assert j==-1
                terms.setdefault(q4lookup[g],[]).append((c*f.D[5]**2,o.const(1)))
            elif w==3:
                assert cw[j]==1;orig=crows[g]
                for k in light:
                    mon=tuple(int(i==k)for i in range(8));n=graph_coefficient(j,mon)
                    if n:terms.setdefault(cubic_lookup[g,k],[]).append((c*row_units[orig],n))
            else:
                assert cw[j]==2;qi=qrows[g]
                for mon in weight2:
                    n=graph_coefficient(j,mon)
                    if n:terms.setdefault(q_lookup[qi,mon],[]).append((c*f.q_units[qi],n))
        vector=[o.lin(terms.get(i,[]))for i in range(354)]
        assert all(o.weights[n]==0 for n in vector if n)
        for k,g in enumerate(quartics):assert o.core_value(vector[q4lookup[g]])==block[bi][k]*f.D[5]**2
        vectors354.append(vector[:]);deletions=[]
        for rel in C['lifted_relations']:
            i=rel['deleted_row'];pivot=expr(rel['pivot']);unit(pivot);n=vector[i]
            if n:
                deletions.append(dict(deleted_row=i,coefficient_node=n,pivot=str(pivot)))
                for k,c in rel['coefficients']:
                    if k!=i:vector[k]=o.lin([(1,vector[k]),(-expr(c)/pivot,n)])
                vector[i]=0
        assert all(not vector[j]for j in C['deleted_row_indices'])
        retained=[vector[i]for i in C['retained_row_indices']];vectors319.append(retained);dropsteps.append(deletions)
        lifts.append(dict(block_row=bi,source_piece_index=piece_index,full_operator_row=p['operator_row'],
                          discarded_selected_equation_terms=skipped,source354_coefficients=vector if False else vectors354[-1],
                          after35_deletions_coefficients319=retained))
    progress('all15 exact relations transported through354to319',nodes=len(o.nodes))
    # The same elementary row operations combine actual relation vectors.  A
    # common D_H^2 normalization makes their selected quartic block exactly I.
    directions=[]
    for i in range(r):
        v=[o.lin([(U[i][k]/f.D[5]**2,vectors319[k][j])for k in range(15)if U[i][k]])for j in range(319)]
        for k,col in enumerate(pivot_cols):
            target=C['retained_row_indices'].index(q4lookup[quartics[col]])
            assert o.core_value(v[target])==(K.one if i==k else K.zero)
        directions.append(v)
    assert len(o.nodes)<200000
    result=dict(id=ID,all_passed=True,base_open_for_new_relations='a*F*B*C!=0',
        literal_source135_oracle_binding=True,literal_full135_identity_checks=literal_checks,
        direct_weight4_relation_count=15,nonzero_quartic_block_rows=11,quartic_block_shape=[15,33],
        uniform_new_independent_kernel_directions=r,known_unit_pivot_rank=r,
        remaining_quartic_block_rows_exactly_zero=all_remaining_zero,
        all15_relations_transported_to_aligned354_and_through35_unit_deletions=True,
        quartic_convention='Qhat_g(v)=F_g(D*v)/D_H^2; D=(a,a,a,1,a,B/50,a,a)',
        old22_scope='FULLU: a*d*F*B*C*J!=0',combined_uniform_directions_on_FULLU=22+r,
        missing_directions_for_80_gate=80-22-r,
        no_80_kernel_full239_rank64_or_fullquartic_claim=True,
        no_new_point_field_GB_or_expanded239_column_matrix=True,
        initial_DAG_nodes=initial_nodes,final_DAG_nodes=len(o.nodes),core_atoms=len(o.atoms),
        resume_checkpoint_used=resume,first_attempt_seconds=prior['first_attempt_seconds']if resume else 0,
        seconds=time.monotonic()-start,frozen=True,independent_review_pending=True)
    packet=dict(source_operator_path=str(operator_path.relative_to(ROOT)),source_operator_sha256=sha(operator_path),
        row_unit_covariance_proof='Accepted full63 native source reconstruction is T34*D with the unique selected34 graph; extraction has common ud^3. The five saved anchors fix only its already-proved whole-row scalars.',
        source_full63_header=[63,91,57],header_third_integer='declared target metadata; ignored by saved Schur parser, not a rank premise',
        supplemental5_cubic_units=supplemental,all21_cubic_units={str(i):str(c)for i,c in row_units.items()},
        raw15_source_lifts=lifts,deletion_steps=dropsteps,
        reduced_pivot_relation_coefficients319=directions,pivot_quartic_columns=pivot_cols,
        source319_row_indices=C['retained_row_indices'],actual319_source_definitions=[rows[i]for i in C['retained_row_indices']],
        arithmetic=o.encoded())
    (out/'certificate.json').write_text(json.dumps(packet,separators=(',',':'))+'\n')
    paths=[Path(__file__),operator_path,UP/'report.json',ROOT/'reports/S61_NONSQUARE_LAZY_QUARTIC_ORACLE_20261103_v63/oracle.json',
        ROOT/'reports/S61_NONSQUARE_ORACLE_FRAME_ADAPTER_20261103_v63/adapter.json',
        ROOT/'reports/S61_NONSQUARE_ORACLE_FRAME_ADAPTER_20261103_v63_REVIEW/report.json',cert_path,
        ROOT/'reports/S61_NONSQUARE_QUARTIC_319_ROWS_20261103_v63/report.json',source_path,norm_path,archive,
        ROOT/'reports/S61_BINARY_NONSQUARE_CUBIC_ROW_BASIS_20261103_v62_REVIEW/report.json',
        ROOT/'reports/S61_NONSQUARE_UNIVERSAL_QUAD22_20261103_v63/report.json',header_writer,header_parser]
    result['input_sha256']={str(p.relative_to(ROOT)):sha(p)for p in paths}
    result['evidence_sha256']={name:sha(out/name)for name in ['certificate.json','block-pivots.json']}
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');signal.alarm(0)
    print(json.dumps(result,indent=2),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--resume-checkpoint',action='store_true');args=p.parse_args();main(args.output,args.resume_checkpoint)
