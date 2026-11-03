#!/usr/bin/env python3
"""Exact actual 9x43 cubic Schur circuit and a specialization-rank obstruction.

Only immutable existing matrices are read. No parametric matrix is expanded,
and no new finite field, Groebner basis, or generic-function-field solver is run.
"""
import argparse
from fractions import Fraction
import gzip
import hashlib
import json
from pathlib import Path
import resource
import subprocess
import time

ROOT=Path(__file__).resolve().parents[1]
ID='S61_NONSQUARE_CUBIC_NINE_ROW_CRITERION_20261103_v63'


def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()


def main(output):
    start=time.monotonic();output.mkdir(parents=True,exist_ok=False)
    base=ROOT/'reports/S61_BINARY_NONSQUARE_CUBIC_ROW_BASIS_20261103_v62'
    report_path=base/'report.json';report=json.loads(report_path.read_text())
    review_path=ROOT/'reports/S61_BINARY_NONSQUARE_CUBIC_ROW_BASIS_20261103_v62_REVIEW/report.json'
    review=json.loads(review_path.read_text());assert report['all_passed'] and review['all_passed']
    assert review['global_free_row_module_rank']==57
    assert review['quadric_product_row_submodule_free_rank']==41
    assert review['abstract_row_module_quotient_free_rank']==16
    for path in [report_path,review_path]:
        data=json.loads(path.read_text())
        for name,want in data['input_sha256'].items():assert sha(ROOT/name)==want,name
    data=gzip.decompress((base/'remaining.txt.gz').read_bytes())
    assert hashlib.sha256(data).hexdigest()==report['remaining_matrix_sha256']
    rows=report['remaining_original_row_indices']
    assert len(rows)==57 and sum(i<42 for i in rows)==41
    source_path=ROOT/'reports/S61_BINARY_NONSQUARE_CUBIC_REPLAY_20261103_v62/source-report.json'
    source=json.loads(source_path.read_text());mons=source['normal_monomials']
    assert all(source['labels'][i][0]=='product' for i in rows[:41])
    assert all(source['labels'][i][0]=='cubic' for i in rows[41:])
    assert len(mons)==91
    heavy=[j for j,mon in enumerate(mons)if mon[5]]
    pure=[j for j,mon in enumerate(mons)if not mon[5]]
    assert len(heavy)==7 and len(pure)==84
    assert all(mons[j][5]==1 and sum(mons[j])==2 for j in heavy)
    assert all(sum(mons[j])==3 for j in pure)
    lines=data.decode().splitlines();assert lines[0]=='57 91 57' and len(lines)==1+57*91
    point=(1,2,3);matrix=[];max_degree=0;coefficient_terms=0
    aa,ee,dd=point
    assert all([aa,500*aa+90*ee**2+3*ee-24,
                (50*aa+18*ee**2-9*ee-2)*dd+60*ee-40,
                100*aa+18*ee**2-9*ee-2])
    for i in range(57):
        values=[]
        for j in range(91):
            tokens=lines[1+i*91+j].split();count=int(tokens[0]);assert len(tokens)==1+4*count
            assert not (i<41 and j in heavy and count)
            q=Fraction(0)
            for k in range(count):
                c=Fraction(tokens[1+4*k]);ex=tuple(map(int,tokens[2+4*k:5+4*k]));assert min(ex)>=0
                q+=c*point[0]**ex[0]*point[1]**ex[1]*point[2]**ex[2]
                max_degree=max(max_degree,sum(ex))
            coefficient_terms+=count;values.append(q)
        matrix.append(values)
    inp=output/'point-matrix.txt'
    inp.write_text('57 91\n'+' '.join(map(str,heavy))+'\n'+
                   '\n'.join(' '.join(map(str,row))for row in matrix)+'\n')
    cpp=ROOT/'code/s61_nonsquare_cubic_nine_row_20261103_v63.cpp';exe=output/'numeric-qq'
    flags=subprocess.check_output(['pkg-config','--cflags','--libs','flint'],text=True).split()
    subprocess.run(['c++','-O2','-std=c++17',str(cpp),*flags,'-o',str(exe)],check=True,capture_output=True,text=True)
    cp=subprocess.run([str(exe),str(inp),str(output/'point-schur.txt')],check=True,capture_output=True,text=True,timeout=20)
    assert 'ALL387_ENTRIES_PASS' in cp.stdout
    (output/'native.log').write_text(cp.stdout+cp.stderr)
    out=(output/'point-schur.txt').read_text().splitlines()
    assert out[0]=='EXACT_QQ_NINE_ROW_CRITERION 41 7 48 9 57'
    ip,g7,g9,remaining,schur9=[list(map(int,line.split()))for line in out[1:6]]
    assert [len(x)for x in [ip,g7,g9,remaining,schur9]]==[41,7,9,43,9]
    dp,dh,d48,ds=map(Fraction,out[6:10]);assert dp and dh and ds and dp*dh==d48
    numeric=[Fraction(x)for x in out[10:]];assert len(numeric)==387
    circuit={
        'source':str((base/'remaining.txt.gz').relative_to(ROOT)),
        'source_uncompressed_sha256':report['remaining_matrix_sha256'],
        'source_row_original_ids':rows,
        'product_rows_basis_positions':list(range(41)),
        'heavy_columns':heavy,'product_pivot_columns':ip,
        'heavy_pivot_rows_basis_positions':g7,'remaining_rows_basis_positions':g9,
        'remaining_columns':remaining,
        'pivot48_rows':list(range(41))+g7,'pivot48_columns':ip+heavy,
        'entry_definition':'S_ij=det A[pivot48_rows+[remaining_rows[i]], pivot48_columns+[remaining_columns[j]]]',
        'entry_count':387,'entry_determinant_size':49,
        'no_entry_expansion_performed':True,
        'domain':'D(DeltaP*DeltaH), aFBC!=0; chosen minors are additional tracked chart units, not known global units',
        'exact_rank_identity':'rank(actual full cubic)=48+rank(S) on this chart',
        'full_rank57_criterion':'rank(S)=9',
        'finite_cover_rule':'All 41-column choices of the fixed41 product rows, and all7-row choices of the16 genuine heavy rows',
        'point':list(point),'point_DeltaP':str(dp),'point_DeltaH':str(dh),
        'point_Delta48':str(d48),'point_Schur9_determinant':str(ds),
        'point_border9_determinant':str(ds*d48**9),
        'point_nonzero_9minor_local_columns':schur9,
        'bordered_entry_total_degree_upper':49*max_degree,
    }
    (output/'circuit.json').write_text(json.dumps(circuit,indent=2)+'\n')
    # A polynomial-module countermodel to the proposed rank implication:
    # I41 in pure84, I7 in heavy7, I8 in the remaining pure coordinates,
    # last row x*e49+y*e50. No claim that this toy satisfies Jacobi.
    toy={
        'base':'QQ[x,y,z]','ambient_columns':91,'pure_columns':84,'heavy_columns':7,
        'product_rank_everywhere':41,'heavy_projection_rank_everywhere':7,
        'abstract_row_module_free_rank':57,'abstract_quotient_by_products_free_rank':16,
        'last_row':'x*e49+y*e50 in two unused pure columns',
        'full_rank_on_D_x_or_y':57,'full_rank_on_x_y_zero':56,
        'full_rank_failure_locus':'V(x,y), dimension1',
        'Jacobi_solution_or_actual_counterexample':False,
    }
    assert 41+7+8+1==57 and 41+7+8==56 and 41+7==48
    (output/'module-countermodel.json').write_text(json.dumps(toy,indent=2)+'\n')
    exe.unlink()
    note=ROOT/'research/notes'/(ID+'.md')
    result={
        'id':ID,'all_passed':True,
        'actual_original57_by91_matrix_used':True,
        'actual_product41_heavy7_block_has_exact9_by43_remainder':True,
        'all387_entries_specified_as_actual49_bordered_minor_circuits':True,
        'chosen_chart_nonempty_by_exact_QQ_point':True,
        'point_rank41_7_48_9_57_verified':True,
        'actual_coefficient_terms_read':coefficient_terms,'maximum_entry_core_degree':max_degree,
        'global_product41_plus_heavy7_does_not_formally_imply_full57':True,
        'toy_countermodel_not_actual_Jacobi':True,
        'full_cubic_failure_finite_not_proved':True,
        'quartic_coverage_or_ambient_removal_not_proved':True,
        'd_zero_not_covered_by_new_pair_product_theorem':True,
        'no_GB_no_new_field_no_parametric_matrix_expansion':True,
        'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in [Path(__file__),cpp,note,report_path,review_path,source_path,base/'remaining.txt.gz']},
        'evidence_sha256':{p.name:sha(p)for p in output.iterdir()if p.is_file()},
        'peak_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        'seconds':time.monotonic()-start,
    }
    (output/'report.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'all_passed':True,'shape':[9,43],'seconds':result['seconds']}))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=ROOT/'reports'/ID)
    main(parser.parse_args().output.resolve())
