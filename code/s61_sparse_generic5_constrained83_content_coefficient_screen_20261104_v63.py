#!/usr/bin/env python3
"""Read-only coefficient-slice refinement of Laurent common-factor widths."""
from pathlib import Path
import hashlib,json,resource,signal,time
ROOT=Path(__file__).resolve().parents[1]
ID='S61_SPARSE_GENERIC5_CONSTRAINED83_CONTENT_COEFFICIENT_SCREEN_20261104_v63'
PARENT='S61_SPARSE_GENERIC5_CONSTRAINED83_LOW_GROWTH_CONTINUATION_20261104_v63'
WIDTH='S61_SPARSE_GENERIC5_CONSTRAINED83_CONTENT_WIDTH_SCREEN_20261104_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def run():
    started=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
    out=ROOT/'reports'/ID;out.mkdir(parents=True,exist_ok=False)
    base=ROOT/'reports'/PARENT;rp=base/'report.json';r=read(rp);assert r['all_passed'] and sha(rp)=='0b4bd1b5530d13e58ea0674753c7c60777ce3278c5022c269a5af0038dcd86c7'
    mp=base/'residual-matrix.json';assert sha(mp)==r['evidence_sha256'][mp.name]
    m=read(mp);A={(i,j):p for i,j,p in m['sparse_entries']}
    wp=ROOT/'reports'/WIDTH/'report.json';wr=read(wp);assert wr['all_passed']
    wep=wp.parent/'exact-width-obstructions.json';assert sha(wep)==wr['evidence_sha256'][wep.name];we=read(wep)
    def refine(initial,axis):
        index=initial['index'];bounds=initial['common_factor_coordinate_degree_bounds'][:]
        keys=([(index,j)for j in m['column_original_ids']if(index,j)in A]if axis=='row'else[(i,index)for i in m['row_original_ids']if(i,index)in A])
        stages=[]
        for step in range(5):
            active=[v for v in range(5)if bounds[v]];outside=[v for v in range(5)if not bounds[v]]
            if not active:break
            improved=bounds[:];witness={}
            for key in keys:
                p=A[key];groups={}
                for term_index,t in enumerate(p['terms']):
                    e=tuple(p['shift'][v]+t[v+1]for v in range(5));groups.setdefault(tuple(e[v]for v in outside),[]).append((term_index,e))
                assert sum(len(g)for g in groups.values())<=190
                for group_key,g in groups.items():
                    for v in active:
                        width=max(e[v]for _,e in g)-min(e[v]for _,e in g)
                        if width<improved[v]:
                            improved[v]=width;witness[str(v)]={'source_position':list(key),'outside_variables':outside,'outside_exponents':list(group_key),'source_term_indices':[k for k,e in g],'coefficient_width':width}
            stages.append({'before':bounds,'after':improved,'strict_decrease_witnesses':witness})
            if improved==bounds:break
            bounds=improved
        return {'axis':axis,'index':index,'initial_bounds':initial['common_factor_coordinate_degree_bounds'],'stages':stages,'remaining_factor_degree_bounds':bounds,'nonmonomial_common_factor_excluded':not any(bounds),'no_coefficient_arithmetic':True}
    results=[refine(z,'row')for z in we['remaining_possible_rows']]+[refine(z,'column')for z in we['remaining_possible_columns']]
    remaining=[z for z in results if not z['nonmonomial_common_factor_excluded']]
    evidence={'proof':'A common factor independent of outside variables divides every coefficient polynomial after grouping by those outside exponents. Each coefficient support therefore gives the same coordinate width obstruction. Iterate when new coordinate bounds become zero. Original coefficients are merely indexed; no additions, products, gcd, or division are computed.',
              'all_refinements':results,'remaining_possible_slices':remaining,
              'positive_bounds_do_not_assert_a_factor_exists':True,'no_pivot_gcd_division_or_series':True}
    ep=out/'coefficient-support-width-witnesses.json';ep.write_text(json.dumps(evidence,indent=2)+'\n')
    report={'id':ID,'all_passed':True,'source143215_terms_unchanged':True,'initial_possible_rows':39,'initial_possible_columns':10,
            'remaining_possible_row_ids':[z['index']for z in remaining if z['axis']=='row'],
            'remaining_possible_column_ids':[z['index']for z in remaining if z['axis']=='column'],
            'no_coefficient_arithmetic_gcd_pivot_or_series':True,
            'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in [rp,mp,wp,wep,ROOT/f'code/{Path(__file__).name}',ROOT/f'research/notes/{ID}.md']},
            'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-started,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items()if k not in('input_sha256','evidence_sha256')}))
if __name__=='__main__':run()
