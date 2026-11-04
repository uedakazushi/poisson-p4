#!/usr/bin/env python3
"""READ ONLY Laurent degree-width obstruction to nonmonomial slice content."""
from pathlib import Path
from math import prod
import hashlib,json,resource,signal,time
ROOT=Path(__file__).resolve().parents[1]
ID='S61_SPARSE_GENERIC5_CONSTRAINED83_CONTENT_WIDTH_SCREEN_20261104_v63'
PARENT='S61_SPARSE_GENERIC5_CONSTRAINED83_LOW_GROWTH_CONTINUATION_20261104_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def run():
    started=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
    out=ROOT/'reports'/ID;out.mkdir(parents=True,exist_ok=False)
    base=ROOT/'reports'/PARENT;rp=base/'report.json';r=read(rp)
    assert sha(rp)=='0b4bd1b5530d13e58ea0674753c7c60777ce3278c5022c269a5af0038dcd86c7' and r['all_passed']
    mp=base/'residual-matrix.json';assert sha(mp)==r['evidence_sha256'][mp.name]
    m=read(mp);assert m['total_terms']==143215
    A={(i,j):p for i,j,p in m['sparse_entries']}
    def width(p):return [max(t[v+1]for t in p['terms'])-min(t[v+1]for t in p['terms'])for v in range(5)]
    def record(index,keys):
        if not keys:return {'index':index,'zero':True}
        ws=[width(A[k])for k in keys]
        minimum=[min(w[v]for w in ws)for v in range(5)]
        witnesses=[next(list(k)for k,w in zip(keys,ws)if w[v]==minimum[v])for v in range(5)]
        return {'index':index,'zero':False,'nonzero_entries':len(keys),'coefficient_terms':sum(len(A[k]['terms'])for k in keys),
                'monomial_entry':next((list(k)for k in keys if len(A[k]['terms'])==1),None),
                'common_factor_coordinate_degree_bounds':minimum,
                'coordinate_minimum_witness_source_positions':witnesses,
                'nonmonomial_common_factor_structurally_impossible':not any(minimum),
                'normalized_common_factor_dense_support_bound':prod(w+1 for w in minimum),
                'max_input_polynomial_terms':max(len(A[k]['terms'])for k in keys),
                'max_input_polynomial_total_width':max(sum(w)for w in ws)}
    rows=[record(i,[(i,j)for j in m['column_original_ids']if(i,j)in A])for i in m['row_original_ids']]
    cols=[record(j,[(i,j)for i in m['row_original_ids']if(i,j)in A])for j in m['column_original_ids']]
    possible_rows=[z for z in rows if not z['zero']and not z['nonmonomial_common_factor_structurally_impossible']]
    possible_cols=[z for z in cols if not z['zero']and not z['nonmonomial_common_factor_structurally_impossible']]
    evidence={'parameters':['H','L','s0','q0','t0'],'exact_source_row_width_screen':rows,'exact_source_column_width_screen':cols,
              'remaining_possible_rows':possible_rows,'remaining_possible_columns':possible_cols,
              'proof':'Normalize a Laurent common factor and complementary factor by their coordinate monomial minima. For each variable, polynomial degrees of a product add over the domain of the other variables. Hence deg_v(g)<=min width_v(P) over all nonzero slice entries. All bounds zero forces the normalized g constant. No evaluation, gcd or polynomial division occurs.',
              'positive_bounds_are_only_possible_factor_support_bounds':True,
              'no_polynomial_coefficient_arithmetic_or_gcd':True}
    ep=out/'exact-width-obstructions.json';ep.write_text(json.dumps(evidence,indent=2)+'\n')
    report={'id':ID,'all_passed':True,'source_coefficient_count_unchanged':143215,
            'nonzero_rows':sum(not z['zero']for z in rows),'nonzero_columns':sum(not z['zero']for z in cols),
            'rows_with_possible_nonmonomial_content':len(possible_rows),'columns_with_possible_nonmonomial_content':len(possible_cols),
            'possible_row_ids':[z['index']for z in possible_rows],'possible_column_ids':[z['index']for z in possible_cols],
            'no_gcd_pivot_division_series_or_new_source_job':True,
            'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in [rp,mp,ROOT/f'code/{Path(__file__).name}',ROOT/f'research/notes/{ID}.md']},
            'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-started,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items()if k not in('input_sha256','evidence_sha256')}))
if __name__=='__main__':run()
