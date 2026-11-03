#!/usr/bin/env python3
"""Lift five accepted exact cubic relations to35 unit quartic row deletions."""
import argparse
from fractions import Fraction
import hashlib
import itertools
import json
from pathlib import Path
import time
import sympy as s
from sympy.polys.rings import ring

ROOT=Path(__file__).resolve().parents[1]
ID='S61_NONSQUARE_QUARTIC_319_ROWS_20261103_v63'


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main(output):
    start=time.monotonic()
    base=ROOT/'reports/S61_BINARY_NONSQUARE_CUBIC_ROW_BASIS_20261103_v62'
    review_path=ROOT/'reports/S61_BINARY_NONSQUARE_CUBIC_ROW_BASIS_20261103_v62_REVIEW/report.json'
    basis_path=base/'report.json';review=json.loads(review_path.read_text());basis=json.loads(basis_path.read_text())
    assert review['all_passed']and review['exact_all6_QQ_relations_all91_columns']
    assert review['exact_all63_known_unit_reconstructions']and review['normalization_uses_only_a_F_B_C']
    for path in [review_path,basis_path]:
        for name,want in json.loads(path.read_text())['input_sha256'].items():assert sha(ROOT/name)==want,name
    cert_path=ROOT/'code/data/s61_binary_nonsquare_cubic_row_basis_20261103_v62.json'
    cert=json.loads(cert_path.read_text());source_path=ROOT/'reports/S61_BINARY_NONSQUARE_CUBIC_REPLAY_20261103_v62/source-report.json'
    source=json.loads(source_path.read_text());labels=source['labels'];assert len(labels)==63
    normal_path=base/'normalization.txt';assert sha(normal_path)==basis['generated_sha256']['normalization.txt']
    R,a,e,d=ring('a,e,d',s.QQ);K=R.to_field()
    known=[a,500*a+90*e**2+3*e-24,(50*a+18*e**2-9*e-2)*d+60*e-40,100*a+18*e**2-9*e-2]
    def unit(p):
        assert p
        for f in known:
            while p.rem(f)==0:p=p.exquo(f)
        assert p.is_ground
    def parse(line):
        ts=line.split();n=int(ts[0]);assert len(ts)==1+4*n
        return R.from_dict({tuple(map(int,ts[2+4*k:5+4*k])):s.QQ(Fraction(ts[1+4*k]))for k in range(n)})
    lines=normal_path.read_text().splitlines();assert lines[0]=='63'and len(lines)==64
    factors=[parse(line)for line in lines[1:]]
    for p in factors:unit(p)
    light=[0,1,2,3,4,6,7]
    exponent=lambda *indices:tuple(indices.count(i)for i in range(8))
    weight2=sorted([exponent(i,j)for i,j in itertools.combinations_with_replacement(light,2)]+[exponent(5)])
    weight4=sorted([exponent(*p)for p in itertools.combinations_with_replacement(light,4)]+
                   [exponent(5,i,j)for i,j in itertools.combinations_with_replacement(light,2)]+[exponent(5,5)])
    assert len(weight2)==29 and len(weight4)==239
    chart_path=ROOT/'reports/S61_BINARY_HYPERPLANE_WEIGHTED_CHART_20261102.json'
    chart=json.loads(chart_path.read_text());q4rows=[i for i,w in enumerate(chart['equation_weights'])if w==4 and i not in chart['selected34_rows']]
    assert len(q4rows)==33
    rows=[['quadric',i,list(mon)]for i in range(6)for mon in weight2]
    rows += [['cubic_normalized',i,labels[i][1],k]for i in range(42,63)for k in light]
    rows += [['quartic',i]for i in q4rows]
    assert len(rows)==354
    lookup={(i,k):174+(i-42)*7+light.index(k)for i in range(42,63)for k in light}
    relations=[];deleted=[];cubic_deleted=[]
    expected=[R.one,R.one,R.one,a,-a]
    for relation,pivot in zip(cert['relations'][:5],expected):
        original=relation['deleted_original_row'];assert original in [48,55,59,57,61]
        coefficients={item['original_row']:R.from_dict({tuple(map(int,term[1:])):s.QQ(Fraction(term[0]))for term in item['terms']})for item in relation['coefficients']}
        assert coefficients[original]==pivot
        assert all(not coefficients.get(i,R.zero)for i in cubic_deleted)
        for k in light:
            terms={}
            for i,coefficient in coefficients.items():
                if i<42:
                    _,qi,j=labels[i]
                    mon=exponent(j,k);row=qi*29+weight2.index(mon)
                    c=K(coefficient)/K(factors[i])
                else:
                    row=lookup[i,k];c=K(coefficient)
                terms[row]=terms.get(row,K.zero)+c
            terms={i:c for i,c in terms.items()if c}
            drop=lookup[original,k]
            assert terms[drop]==K(pivot)and all(not terms.get(i,K.zero)for i in deleted)
            for c in terms.values():unit(c.denom)
            relations.append({'deleted_row':drop,'deleted_source':rows[drop],
                              'pivot':str(pivot),
                              'coefficients':[[i,str(c)]for i,c in sorted(terms.items())]})
            deleted.append(drop)
        cubic_deleted.append(original)
    assert len(deleted)==len(set(deleted))==35
    retained=[i for i in range(354)if i not in deleted]
    assert len(retained)==319
    triangle=R.one
    for p in expected:triangle*=p**7
    assert triangle==-a**14
    output.mkdir(parents=True,exist_ok=False)
    evidence={'original_shape':[354,239],'retained_shape':[319,239],
              'source_rows':rows,'weight4_columns':weight4,
              'lifted_relations':relations,'deleted_row_indices':deleted,'retained_row_indices':retained,
              'triangular35minor':'-a^14','genuine_quartic_original_rows_retained':q4rows,
              'source_cubic_identity_scope':'All91 columns exact in independent accepted review',
              'new_full135_quartic_expansion_performed':False}
    (output/'certificate.json').write_text(json.dumps(evidence,indent=2)+'\n')
    note=ROOT/'research/notes'/(ID+'.md')
    result={'id':ID,'all_passed':True,'open':'a*F*B*C!=0; no d or J inversion',
            'actual_full_weight4_row_span_preserved':True,'old_shape':[354,239],'new_shape':[319,239],
            'all35_lifted_cubic_relations_and_unit_pivots_checked':True,
            'row_span_equality_both_directions':True,'all33_genuine_quartic_rows_retained':True,
            'triangular_deleted_minor':'-a^14','unknown_core_factor_inverted':False,
            'quartic_surjectivity_length65_or_ambient_claim':False,
            'no_new_GB_field_or_parametric_matrix_expansion':True,
            'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in [Path(__file__),note,review_path,basis_path,cert_path,source_path,normal_path,chart_path]},
            'evidence_sha256':{'certificate.json':sha(output/'certificate.json')},
            'seconds':time.monotonic()-start}
    (output/'report.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'all_passed':True,'new_shape':[319,239],'seconds':result['seconds']}))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=ROOT/'reports'/ID)
    main(parser.parse_args().output.resolve())
