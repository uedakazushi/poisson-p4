#!/usr/bin/env python3
"""Bounded independent binding/recorded-pivot/DAG-transport audit, no239 matrix."""
import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import resource
import signal
import tarfile
import time

import sympy as s
from s61_nonsquare_lazy_quartic_oracle_20261103_v63 import Oracle
from s61_nonsquare_oracle_frame_adapter_20261103_v63 import FrameAdapter

ROOT=Path(__file__).resolve().parents[1]
STEM='S61_NONSQUARE_BIANCHI_QUARTIC_UNIT_PIVOTS_20261104_v63_REVIEW'


def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def accepted(p):
    r=read(p);assert r['all_passed']
    for n,h in r.get('input_sha256',{}).items():assert sha(ROOT/n)==h,n
    for n,h in r.get('evidence_sha256',{}).items():
        q=ROOT/n if (ROOT/n).exists() else p.parent/n
        assert sha(q)==h,n
    return r


def main(out):
    start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60)
    assert not out.exists();out.mkdir(parents=True)
    base=ROOT/'reports/S61_NONSQUARE_BIANCHI_QUARTIC_UNIT_PIVOTS_20261104_v63'
    ap=base/'report.json';author=accepted(ap)
    upp=ROOT/'reports/S61_NONSQUARE_FULL_BIANCHI_OPERATOR_20261104_v63_REVIEW/report.json'
    framereview=ROOT/'reports/S61_NONSQUARE_ORACLE_FRAME_ADAPTER_20261103_v63_REVIEW/report.json'
    basisreview=ROOT/'reports/S61_BINARY_NONSQUARE_CUBIC_ROW_BASIS_20261103_v62_REVIEW/report.json'
    quadreview=ROOT/'reports/S61_NONSQUARE_UNIVERSAL_QUAD22_20261103_v63_REVIEW/report.json'
    up,fr,br,qr=[accepted(p) for p in [upp,framereview,basisreview,quadreview]]
    assert up['all175_literal_syzygies_proved_by_verified_Jacobi_correspondence_and_graded_Jacobi']
    assert fr['all16_genuine_global_row_units_proved_by_native_covariance_then_fixed_by_anchors']
    assert qr['pointwise_22_directions_independent_on_FULLU']
    cp=base/'certificate.json';c=read(cp);bp=base/'block-pivots.json';pivot=read(bp)
    pubpath=base/'publish-receipt.json';pub=read(pubpath)
    assert pub['all_passed'] and pub['report_sha256']==sha(ap) and pub['source_sha256']==sha(ROOT/'code/s61_nonsquare_bianchi_quartic_unit_pivots_20261104_v63.py')
    assert pub['note_sha256']==sha(ROOT/'research/notes/S61_NONSQUARE_BIANCHI_QUARTIC_UNIT_PIVOTS_20261104_v63.md')
    for n,h in pub['artifacts_sha256'].items():assert sha(base/n)==h
    opath=ROOT/c['source_operator_path'];O=read(opath);assert sha(opath)==c['source_operator_sha256']
    oraclepath=ROOT/'reports/S61_NONSQUARE_LAZY_QUARTIC_ORACLE_20261103_v63/oracle.json'
    oldoracle=read(oraclepath)
    framepath=ROOT/'reports/S61_NONSQUARE_ORACLE_FRAME_ADAPTER_20261103_v63/adapter.json';frame=read(framepath)
    o,_=Oracle.load(cp);f=FrameAdapter(o,oldoracle,frame);K,R=o.K,o.R
    assert len(o.nodes)==10710 and len(o.atoms)==4228
    assert o.nodes[:1373]==oldoracle['arithmetic']['nodes']
    assert c['arithmetic']['atoms'][:len(oldoracle['arithmetic']['atoms'])]==oldoracle['arithmetic']['atoms']
    for i,node in enumerate(o.nodes):
        if i==0:assert node==['z'];continue
        tag=node[0]
        if tag=='v':assert o.weights[i]==o.variable_weights[node[1]] and o.degrees[i]==1
        elif tag=='c':assert o.weights[i]==0 and o.degrees[i]==0
        elif tag=='m':
            assert 0<node[1]<i and 0<node[2]<i
            assert o.weights[i]==o.weights[node[1]]+o.weights[node[2]]<=4
        else:
            assert tag=='l' and all(0<n<i and o.atoms[a] for a,n in node[1])
            assert all(o.weights[n]==o.weights[i] for a,n in node[1])
    expr=lambda v:K.from_expr(s.sympify(v))
    poly=lambda terms:R.from_dict({tuple(e):s.QQ(v) for e,v in terms})
    frac=lambda terms:K(poly(terms))
    def unit(value):
        value=K(value);assert value
        constant=K.one;powers=[]
        for sign,pol in [(1,value.numer),(-1,value.denom)]:
            counts=[]
            for p in o.known:
                n=0
                while not pol.rem(p):pol=pol.exquo(p);n+=1
                counts.append(sign*n)
            assert pol.is_ground and pol
            constant=constant*K(pol)**sign;powers.append(counts)
        return constant,[sum(v) for v in zip(*powers)]
    # Source135 compatibility follows the accepted operator proof; rebind only
    # the local J/B atom interface, without re-expanding any175 identities.
    fixed=O['local_fixed_q_index'];cols=O['local_other_q_order'];rw=O['equation_weights'];cw=O['local_coordinate_weights']
    assert fixed==oldoracle['fixed_original_coordinate'] and cols==oldoracle['other_original_coordinate_order']
    assert rw==oldoracle['equation_weights'] and cw==oldoracle['local_coordinate_weights']
    def add(d,k,v):
        v=d.get(k,K.zero)+v
        if v:d[k]=v
        else:d.pop(k,None)
    def local(v):return -1 if v==fixed else (44 if v==45 else cols.index(v))
    for i,row in enumerate(O['full_homogeneous135_equations']):
        actual={}
        for mon,c0 in row:add(actual,tuple(sorted(local(v) for v in mon if v!=fixed)),frac(c0))
        expected={}
        for j,a in oldoracle['J_sparse_atom_ids'][i]:add(expected,(j,),o.atoms[a])
        for j,k,a in oldoracle['B_sparse_atom_ids'][i]:add(expected,tuple(sorted((j,k))),o.atoms[a])
        assert actual==expected
    for j in [19,23]:
        assert all(oldoracle[key][j]==0 for key in ['linear_graph_roots','local_Z2_roots','local_Z3_roots','local_Z4_roots'])
    quartics=O['genuine_quartic_original_rows'];pieces=O['weight4_quartic_constant_block_row_piece_indices']
    block=[[frac(v) for v in row] for row in O['weight4_quartic_constant_block']]
    assert len(block)==15 and len(quartics)==33
    assert [[str(v) for v in row] for row in block]==pivot['original_block']
    B=[row[:] for row in block];U=[[K(int(i==j)) for j in range(15)] for i in range(15)]
    for r,step in enumerate(pivot['steps']):
        assert step['pivot_row']==r
        i,j=step['source_row_before_swap'],step['quartic_column'];p=expr(step['pivot'])
        assert B[i][j]==p and quartics[j]==step['original_quartic_row']
        ground,powers=unit(p);assert str(ground)==step['unit']['constant'] and powers==step['unit']['a_F_B_C_powers']
        B[r],B[i]=B[i],B[r];U[r],U[i]=U[i],U[r]
        B[r]=[a/p for a in B[r]];U[r]=[a/p for a in U[r]]
        elim=[[k,str(B[k][j])] for k in range(15) if k!=r and B[k][j]];assert elim==step['eliminations']
        for k,v in elim:
            c0=expr(v);B[k]=[a-c0*b for a,b in zip(B[k],B[r])];U[k]=[a-c0*b for a,b in zip(U[k],U[r])]
    assert len(pivot['steps'])==11 and all(not v for row in B[11:] for v in row)
    assert [[str(v) for v in row] for row in U]==pivot['elementary_transform']
    assert [[str(v) for v in row] for row in B]==pivot['transformed_block']
    for i,j in enumerate(pivot['pivot_columns']):assert [B[k][j] for k in range(11)]==[K(int(k==i)) for k in range(11)]
    (out/'progress.json').write_text(json.dumps({'stage':'source interface and saved11-unit elementary steps verified','seconds':time.monotonic()-start})+'\n')
    sourcepath=ROOT/'reports/S61_BINARY_NONSQUARE_CUBIC_REPLAY_20261103_v62/source-report.json';source=read(sourcepath)
    normpath=ROOT/'reports/S61_BINARY_NONSQUARE_CUBIC_ROW_BASIS_20261103_v62/normalization.txt'
    def native(line):
        t=line.split();n=int(t[0]);assert len(t)==1+4*n
        return K(R.from_dict({tuple(map(int,t[2+4*k:5+4*k])):s.QQ(Fraction(t[1+4*k])) for k in range(n)}))
    lines=normpath.read_text().splitlines();assert lines[0]=='63';norms=[native(v) for v in lines[1:]]
    for v in norms:unit(v)
    rowunits={i:c0 for i,c0 in zip(frame['basis57_original_row_indices'],f.basis_units) if i>=42}
    supplements=c['supplemental5_cubic_units'];anchors={}
    wanted={(r['original63_row'],source['normal_monomials'].index(r['anchor_monomial'])) for r in supplements}
    archive=ROOT/'reports/S61_BINARY_NONSQUARE_CUBIC_REPLAY_20261103_v62/matrices.tar.gz';digest=hashlib.sha256()
    with tarfile.open(archive) as tar:
        stream=tar.extractfile('full');header=stream.readline();digest.update(header);assert header.decode().strip()=='63 91 57'
        for i in range(63):
            for j in range(91):
                line=stream.readline();assert line;digest.update(line)
                if (i,j) in wanted:anchors[i,j]=native(line.decode())
        assert not stream.read()
    assert digest.hexdigest()==source['full_matrix_sha256']
    writer=ROOT/'code/s61_binary_hyperplane_cubic_matrix_20261102.py';parser=ROOT/'code/s61_binary_hyperplane_cubic_unit_schur_20261102.cpp'
    assert "f.write('63 91 57"+chr(92)+"n')" in writer.read_text()
    assert 'in>>n>>m>>r;assert(n==63&&m==91)' in parser.read_text()
    for rec in supplements:
        i,g=rec['original63_row'],rec['physical135_row'];mon=tuple(rec['anchor_monomial']);j=source['normal_monomials'].index(list(mon))
        assert source['labels'][i][1]==g
        saved=anchors[i,j]/norms[i];assert saved==expr(rec['normalized_full63_anchor'])
        factor=expr(rec['unit']);ground,powers=unit(factor)
        assert str(ground)==rec['unit_factorization']['constant'] and powers==rec['unit_factorization']['a_F_B_C_powers']
        actual=o.core_value(f.equation_coefficient(g,mon))
        assert actual==factor*saved==expr(rec['actual_physical_anchor'])
        rowunits[i]=factor
    assert len(rowunits)==21 and {str(i):str(v) for i,v in rowunits.items()}==c['all21_cubic_units']
    (out/'progress.json').write_text(json.dumps({'stage':'five unit anchors and full63 header/memberSHA verified','seconds':time.monotonic()-start})+'\n')
    rows_path=ROOT/'reports/S61_NONSQUARE_QUARTIC_319_ROWS_20261103_v63/certificate.json';C=read(rows_path);rows=C['source_rows']
    assert c['source319_row_indices']==C['retained_row_indices'] and c['actual319_source_definitions']==[rows[i] for i in C['retained_row_indices']]
    qlookup={(r[1],tuple(r[2])):i for i,r in enumerate(rows) if r[0]=='quadric'}
    clookup={(r[2],r[3]):i for i,r in enumerate(rows) if r[0]=='cubic_normalized'}
    q4={r[1]:i for i,r in enumerate(rows) if r[0]=='quartic'};assert list(q4)==quartics
    qrows={g:i for i,g in enumerate(frame['quadric_original_rows'])};crows={source['labels'][i][1]:i for i in range(42,63)}
    mons2=sorted({m for i,m in qlookup});light=[0,1,2,3,4,6,7]
    E=oldoracle['linear_graph_roots'];Z2=oldoracle['local_Z2_roots']
    graph={j:E[j] if w==1 else o.lin([(1,E[j]),(1,Z2[j])]) for j,w in enumerate(cw) if w in [1,2]}
    def coefficient(j,mon):return o.lin([(f.scalar(mon),o.coefficient(graph[j],mon))])
    selected=set(oldoracle['selected34_rows']);v319=[]
    for bi,pi in enumerate(pieces):
        piece=O['local_weight_pieces'][pi];record=c['raw15_source_lifts'][bi]
        assert record['source_piece_index']==pi and record['full_operator_row']==piece['operator_row']
        terms={};skipped=[]
        for g,j,cs in piece['coefficients']:
            a=frac(cs)
            if g in selected:skipped.append([g,j]);continue
            if rw[g]==4:
                if j in [19,23]:continue
                assert j==-1;terms.setdefault(q4[g],[]).append((a*f.D[5]**2,o.const(1)))
            elif rw[g]==3:
                assert cw[j]==1
                for k in light:
                    mon=tuple(int(i==k) for i in range(8));n=coefficient(j,mon)
                    if n:terms.setdefault(clookup[g,k],[]).append((a*rowunits[crows[g]],n))
            else:
                assert rw[g]==2 and cw[j]==2
                for mon in mons2:
                    n=coefficient(j,mon)
                    if n:terms.setdefault(qlookup[qrows[g],mon],[]).append((a*f.q_units[qrows[g]],n))
        vector=[o.lin(terms.get(i,[])) for i in range(354)]
        assert skipped==record['discarded_selected_equation_terms'] and vector==record['source354_coefficients']
        deletions=[]
        for rel in C['lifted_relations']:
            i=rel['deleted_row'];p=expr(rel['pivot']);unit(p);n=vector[i]
            if n:
                deletions.append(dict(deleted_row=i,coefficient_node=n,pivot=str(p)))
                for k,cs in rel['coefficients']:
                    if k!=i:vector[k]=o.lin([(1,vector[k]),(-expr(cs)/p,n)])
                vector[i]=0
        assert deletions==c['deletion_steps'][bi]
        assert all(vector[i]==0 for i in C['deleted_row_indices'])
        reduced=[vector[i] for i in C['retained_row_indices']]
        assert reduced==record['after35_deletions_coefficients319'];v319.append(reduced)
    for i in range(11):
        direction=[o.lin([(U[i][k]/f.D[5]**2,v319[k][j]) for k in range(15) if U[i][k]]) for j in range(319)]
        assert direction==c['reduced_pivot_relation_coefficients319'][i]
        for k,col in enumerate(pivot['pivot_columns']):
            position=C['retained_row_indices'].index(q4[quartics[col]])
            assert o.core_value(direction[position])==K(int(i==k))
    assert len(o.nodes)==10710 and len(o.atoms)==4228
    frontier=base/'first-attempt-frontier.json';fr0=read(frontier)
    assert fr0['block_sha256']==sha(bp) and fr0['rank11_and_remaining_zero']
    assert pub['no_second_pivot_or_source_verification'] and pub['cumulative_author_job_seconds']<60
    paths=[Path(__file__),ROOT/f'research/notes/{STEM}.md',ap,cp,bp,pubpath,frontier,upp,framereview,basisreview,quadreview,
        opath,oraclepath,framepath,sourcepath,normpath,archive,writer,parser,rows_path,
        ROOT/'code/s61_nonsquare_lazy_quartic_oracle_20261103_v63.py',ROOT/'code/s61_nonsquare_oracle_frame_adapter_20261103_v63.py']
    result=dict(id=STEM,all_passed=True,
        upstream_full175_source_proof_reused_without_reexpansion=True,
        full135_local_J_B_literal_binding_and_exact_zero_nuisance_checked=True,
        DAG_prefix_atoms_denominators_and_topology_checked=True,
        saved11_pivot_steps_replayed_no_search_and_all_pivots_known_units=True,
        full63_archive_member_SHA_and_header63_91_57_independently_checked=True,
        header57_is_unused_writer_target_metadata_not_rank_evidence=True,
        five_supplemental_cubic_anchors_and_all21_whole_row_units_bound=True,
        whole_row_covariance_proof_used_not_anchor_extrapolation=True,
        all15_exact354_vectors_and35unit_deletions_match_all319_coefficients=True,
        all11_final_vectors_have_literal_identity_quartic_block=True,
        exact_DAG_nodes=10710,exact_core_atoms=4228,
        new_pointwise_independent_kernel_directions_on_aFBC=11,
        old22_independence_scope='FULLU=a*d*F*B*C*J!=0',combined_independent_on_FULLU=33,remaining_for80=47,
        no_C122_rank64_full239_or80_claim=True,
        no_expanded239_column_matrix_GB_new_field_point_or_pivot_search=True,
        verification_kind='Independent interface/unit/pivot certificate and exact finite-DAG transport verification; no producer invocation or fulloperator replay.',
        input_sha256={str(p.relative_to(ROOT)):sha(p) for p in paths},seconds=time.monotonic()-start)
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');signal.alarm(0)
    print(json.dumps({k:result[k] for k in ['all_passed','new_pointwise_independent_kernel_directions_on_aFBC','combined_independent_on_FULLU','remaining_for80','seconds']},indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);args=p.parse_args()
    try:main(args.output)
    except BaseException as e:
        if args.output.exists() and not (args.output/'report.json').exists():
            (args.output/'incomplete.json').write_text(json.dumps({'all_passed':False,'error':type(e).__name__,'detail':str(e),'no_automatic_larger_retry':True},indent=2)+'\n')
        raise
