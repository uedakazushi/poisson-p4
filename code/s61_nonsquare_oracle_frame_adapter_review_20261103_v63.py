#!/usr/bin/env python3
"""Independent DAG coefficient arithmetic and native frame-covariance audit.

Does not rebuild all148858 genuine cubic core terms, probe a point, or compute rank.
"""
import argparse
from functools import lru_cache
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import resource
import signal
import tarfile
import time
import sympy as s
from sympy.polys.rings import ring

ROOT=Path(__file__).resolve().parents[1]
ID='S61_NONSQUARE_ORACLE_FRAME_ADAPTER_20261103_v63_REVIEW'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main(output):
    start=time.monotonic();signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
    output.mkdir(parents=True,exist_ok=False)
    parent=ROOT/'reports/S61_NONSQUARE_ORACLE_FRAME_ADAPTER_20261103_v63'
    arpath=parent/'report.json';framepath=parent/'adapter.json';ar,frame=read(arpath),read(framepath)
    assert ar['all_exact_coordinate_and_row_unit_checks_passed']
    oraclepath=ROOT/frame['oracle_path'];data=read(oraclepath);assert sha(oraclepath)==frame['oracle_sha256']
    oracle_review=ROOT/'reports/S61_NONSQUARE_LAZY_QUARTIC_ORACLE_20261103_v63_REVIEW/report.json'
    row_review=ROOT/'reports/S61_BINARY_NONSQUARE_CUBIC_ROW_BASIS_20261103_v62_REVIEW/report.json'
    replay=ROOT/'reports/S61_BINARY_NONSQUARE_CUBIC_REPLAY_20261103_v62/report.json'
    old=ROOT/'reports/S61_BINARY_NONSQUARE_CUBIC_REPLAY_20261103_v62/source-report.json'
    source=read(old);assert source['all_passed']
    checked={}
    for p in [arpath,oracle_review,row_review,replay]:
        obj=read(p);assert obj.get('all_passed',obj.get('all_exact_coordinate_and_row_unit_checks_passed',False))
        for key in ['input_sha256','evidence_sha256','generated_sha256']:
            for name,want in obj.get(key,{}).items():
                q=next((q for q in [ROOT/name,p.parent/name]if q.is_file()and sha(q)==want),None)
                assert q is not None,(p,name);checked[str(q.relative_to(ROOT))]=want
    for name,want in source['input_sha256'].items():assert sha(ROOT/name)==want;checked[name]=want
    chartpath=ROOT/'reports/S61_BINARY_HYPERPLANE_WEIGHTED_CHART_20261102.json';chart=read(chartpath)
    quadpath=ROOT/'reports/S61_BINARY_HYPERPLANE_QUADRATIC_CODIMENSION_20261102.json';quad=read(quadpath)
    baseroot=ROOT/'reports/S61_BINARY_NONSQUARE_CUBIC_ROW_BASIS_20261103_v62';basis=read(baseroot/'report.json')
    hpath=ROOT/'reports/S61_BINARY_HYPERPLANE_HEAVY_SQUARE_PIVOT_20261102.json';heavy=read(hpath)
    R,a,e,d=ring('a,e,d',s.QQ);K=R.to_field();frac=lambda x:K.from_expr(s.sympify(x))
    known=[a,500*a+90*e*e+3*e-24,(50*a+18*e*e-9*e-2)*d+60*e-40,100*a+18*e*e-9*e-2]
    def unit(value):
        assert value
        powers=[];constants=[]
        for polynomial in [value.numer,value.denom]:
            exponents=[]
            for f in known:
                count=0
                while polynomial.rem(f)==0:polynomial=polynomial.exquo(f);count+=1
                exponents.append(count)
            assert polynomial and polynomial.is_ground;powers.append(exponents);constants.append(polynomial)
        return {'constant':str(K(constants[0])/K(constants[1])),
                'a_F_B_C_powers':[i-j for i,j in zip(*powers)]}
    D=list(map(frac,frame['Tden']));assert D==[K(a),K(a),K(a),K.one,K(a),K(known[2])/50,K(a),K(a)]
    assert D==list(map(frac,source['Tden']))
    ud=frac(source['Uden']);assert ud==K(a*a*known[1]*known[2])/25000;unit(ud)
    record=data['arithmetic'];atoms=[]
    for numerator,denominator in record['atoms']:
        polynomial=lambda ts:R.from_dict({tuple(ex):s.QQ(c)for ex,c in ts})
        atoms.append(K(polynomial(numerator))/K(polynomial(denominator)))
    T=[[atoms[i]for i in row]for row in data['T_atom_ids']]
    U=[[atoms[i]for i in row]for row in data['U_atom_ids']]
    assert T==[[frac(row[j])for j in chart['normal34']]for row in chart['T34']]
    assert U==[[frac(c)for c in row]for row in chart['U34']]
    assert data['selected34_rows']==chart['selected34_rows']
    assert data['source_original_equations']==135
    oldchart=read(ROOT/'reports/S61_BINARY_DEEPER_TRIANGULAR_CHART_20261101.json')['linear']
    assert data['fixed_original_coordinate']==oldchart['fixed'] and data['other_original_coordinate_order']==oldchart['cols']
    archivepath=replay.parent/'matrices.tar.gz'
    with tarfile.open(archivepath)as archive:
        archivedweights=list(map(int,archive.extractfile('weights').read().split()))
    assert archivedweights==chart['local_coordinate_weights']+chart['equation_weights']+chart['normal_weights34']
    #Independent direct coefficient recursion in the frozen DAG, not Oracle.coefficient/core_value.
    nodes=record['nodes'];weights=record['node_weights'];degrees=record['node_ordinary_degrees'];vw=record['variable_weights']
    @lru_cache(None)
    def coefficient(node,mon):
        if not node or sum(i*j for i,j in zip(vw,mon))!=weights[node]:return K.zero
        if degrees[node]is not None and degrees[node]>=0 and sum(mon)!=degrees[node]:return K.zero
        nd=nodes[node]
        if nd[0]=='c':return atoms[nd[1]]if not any(mon)else K.zero
        if nd[0]=='v':return K.one if mon==tuple(int(j==nd[1])for j in range(8))else K.zero
        if nd[0]=='l':return sum((atoms[c]*coefficient(child,mon)for c,child in nd[1]),K.zero)
        assert nd[0]=='m';value=K.zero
        for left in itertools.product(*(range(k+1)for k in mon)):
            if sum(i*j for i,j in zip(vw,left))!=weights[nd[1]]:continue
            right=tuple(i-j for i,j in zip(mon,left))
            value+=coefficient(nd[1],left)*coefficient(nd[2],right)
        return value
    def scalar(mon):
        answer=K.one
        for c,power in zip(D,mon):answer*=c**power
        return answer
    def native(line):
        ts=line.split();n=int(ts[0]);assert len(ts)==1+4*n
        return K(R.from_dict({tuple(map(int,ts[2+4*j:5+4*j])):s.QQ(ts[1+4*j])for j in range(n)}))
    normalization=baseroot/'normalization.txt';norms=list(map(native,normalization.read_text().splitlines()[1:]));assert len(norms)==63
    for c in norms:unit(c)
    nativepath=baseroot/'remaining.txt.gz';raw=gzip.decompress(nativepath.read_bytes());assert hashlib.sha256(raw).hexdigest()==basis['remaining_matrix_sha256']
    lines=raw.decode().splitlines();assert lines[0]=='57 91 57'and len(lines)==5188
    savedq=[{tuple(mon):K(R.from_dict({tuple(ex):s.QQ(c)for ex,c in terms}))for mon,terms in row}for row in quad['quadratic_rows_sparse']]
    q_units=list(map(frac,frame['quadratic_encoded_source_units']));assert len(q_units)==6
    light=[0,1,2,3,4,6,7]
    mon2=[tuple(pair.count(i)for i in range(8))for pair in itertools.combinations_with_replacement(light,2)]
    mon2.append((0,0,0,0,0,1,0,0));qcount=0
    for qi,original in enumerate(quad['quadratic_rows']):
        assert frame['quadric_original_rows'][qi]==original and unit(q_units[qi])==frame['quadratic_unit_factorizations'][qi]
        for mon in mon2:
            assert coefficient(data['residual_weight_le4_roots'][original],mon)*scalar(mon)==q_units[qi]*savedq[qi].get(mon,K.zero)
            qcount+=1
    definitions=[];product_count=0;anchors=[];bunits=list(map(frac,frame['basis57_encoded_source_units']))
    assert frame['basis57_original_row_indices']==basis['remaining_original_row_indices']
    for position,original in enumerate(basis['remaining_original_row_indices']):
        label=source['labels'][original]
        if label[0]=='product':
            _,qi,k=label;definitions.append({'kind':'product','original_equation':quad['quadratic_rows'][qi],'light_index':k})
            assert bunits[position]==q_units[qi]*norms[original]
            for j,mon in enumerate(source['normal_monomials']):
                lower=list(mon);lower[k]-=1
                expected=savedq[qi].get(tuple(lower),K.zero)if min(lower)>=0 else K.zero
                assert native(lines[1+position*91+j])*norms[original]==expected;product_count+=1
        else:
            assert label[0]=='cubic';equation=label[1]
            definitions.append({'kind':'genuine_cubic','original_equation':equation})
            assert equation not in data['selected34_rows'] and data['equation_weights'][equation]==3
            seed=next(item for item in frame['sixteen_genuine_anchor_checks']if item['basis_position']==position)
            mon=tuple(seed['monomial']);j=source['normal_monomials'].index(list(mon))
            factor=coefficient(data['residual_weight_le4_roots'][equation],mon)*scalar(mon)/native(lines[1+position*91+j])
            assert factor==bunits[position] and unit(factor)==seed['known_unit_factor']
            anchors.append({'basis_position':position,'original_equation':equation,'factor':str(factor),'unit':unit(factor)})
        unit(bunits[position])
    assert definitions==frame['basis57_source_definitions'] and product_count==3731 and len(anchors)==16 and qcount==174
    #Old heavy formulas are T34, after dividing native A2 by ud*D_H^2.
    hm=(0,0,0,0,0,2,0,0)
    for row,oldvalue in zip([77,131],heavy['coefficient_formulas']):
        value=coefficient(data['residual_weight_le4_roots'][row],hm)
        assert value==frac(oldvalue) and value*scalar(hm)/D[5]**2==frac(oldvalue)
    assert frac(frame['common_quartic_unit_rescale'])==1/D[5]**2
    assert frame['actual_quartic_original_rows']==data['genuine_quartic_original_rows']and len(frame['actual_quartic_original_rows'])==33
    #API checks only after the independent comparisons.
    from s61_nonsquare_oracle_frame_adapter_20261103_v63 import FrameAdapter
    api=FrameAdapter.load();before=len(api.oracle.nodes)
    qfixture=(2,0,0,0,0,0,0,0)
    assert api.core_value(api.quad_coefficient(0,qfixture))==savedq[0].get(qfixture,K.zero)
    fixtures=[]
    for position in [0,41]:
        j=next(j for j in range(91)if int(lines[1+position*91+j].split()[0]))
        mon=tuple(source['normal_monomials'][j]);value=api.core_value(api.basis57_coefficient(position,mon))
        assert value==native(lines[1+position*91+j]);fixtures.append([position,list(mon)])
    assert api.core_value(api.quartic_coefficient(77,hm))==frac(heavy['coefficient_formulas'][0])
    evidence={'all174_quadratic_coefficients':qcount,'all3731_product_basis_entries':product_count,
        'sixteen_genuine_exact_anchors':anchors,'all57_source_definitions_checked':True,
        'native_covariance':'W_d=ud^(d-1)*Z_d; A_d=ud^(d-1)*F_d(T34*D*v); extraction A2*ud^2+A3*ud=ud^3*(F2+F3)',
        'no_row_combination_in_extract_or_original63_to_remaining57':True,
        'four_API_fixtures':[list(qfixture),*fixtures,list(hm)],'API_query_nodes_added':len(api.oracle.nodes)-before,
        'all33_quartic_common_unit':'D_H^-2; preserves literal H² coefficients and old u,v',
        'direct_coefficient_cache_entries':coefficient.cache_info().currsize,
        'no_all148858_genuine_cubic_term_expansion_or_actual_rank_claim':True}
    ev=output/'independent-certificate.json';ev.write_text(json.dumps(evidence,indent=2)+'\n')
    note=ROOT/'research/notes'/(ID+'.md')
    files=[Path(__file__),note,arpath,framepath,oraclepath,oracle_review,row_review,replay,old,archivepath,chartpath,quadpath,
           baseroot/'report.json',normalization,nativepath,hpath,ROOT/'code/s61_nonsquare_oracle_frame_adapter_20261103_v63.py',
           ROOT/'code/parametric_implicit_jets_20261101.cpp',ROOT/'code/s61_binary_hyperplane_cubic_extract_20261102.cpp',
           ROOT/'code/s61_binary_hyperplane_cubic_matrix_20261102.py',ROOT/'code/s61_binary_hyperplane_heavy_square_pivot_20261102.py',
           ROOT/'research/notes/S61_NONSQUARE_ORACLE_FRAME_ADAPTER_20261103_v63.md']
    result={'id':ID,'all_passed':True,'same_original135_selected34_T_U_source_and_diagonal_known_units_accepted':True,
        'all174_quadratic_and3731_product_entries_independently_checked':True,
        'all16_genuine_global_row_units_proved_by_native_covariance_then_fixed_by_anchors':True,
        'no_unproved_anchor_to_whole_polynomial_inference_or_product41_row_combination':True,
        'all33_quartic_Qhat_convention_and_literal_old_Hsquare_data_accepted':True,
        'four_public_API_fixtures_independently_match':True,
        'full_genuine_cubic_core_expansion_new_rank_point_field_or_GB':False,
        'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in files},'checked_dependency_sha256':checked,
        'evidence_sha256':{ev.name:sha(ev)},'peak_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        'seconds':time.monotonic()-start}
    (output/'report.json').write_text(json.dumps(result,indent=2)+'\n')
    signal.alarm(0);print(json.dumps({'all_passed':True,'quadratic_entries':qcount,'product_entries':product_count,'seconds':result['seconds']}))

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,default=ROOT/'reports'/ID)
    main(parser.parse_args().output.resolve())
