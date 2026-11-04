#!/usr/bin/env python3
"""ONE Ka13 17x24 coefficient preparation; no matrix-rank or prime job."""
import argparse,hashlib,json,math,resource,signal,sys,time
from pathlib import Path
from s61_nonsquare_raw_H_B_K13_named_saturation_20261104_v63 import Ring,Stop,pa,pm
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_RAW_H_B_K13_CUBIC17_COEFFICIENTS_20261104_v63'
sys.set_int_max_str_digits(0)
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def ip(p):return {int(i):int(v)for i,v in p}
def ser(p):return [[i,str(v)]for i,v in sorted(p.items(),reverse=True)]
def exactdiv(p,g):
    p=dict(p);q={};deg=max(g);lc=g[deg]
    while p and max(p)>=deg:
        k=max(p)-deg;v=p[max(p)];assert v%lc==0;c=v//lc;q[k]=c;p=pa(p,{i+k:c*v for i,v in g.items()},-1)
    assert not p;return q
def main(out):
    start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000))
    signal.signal(signal.SIGALRM,lambda *_:(_ for _ in()).throw(TimeoutError('60s total cap')));signal.alarm(60)
    out=out.resolve();out.mkdir(parents=True,exist_ok=False)
    pb=ROOT/'reports/S61_NONSQUARE_RAW_H_B_K13_FULL_CUBIC_PREFLIGHT_20261104_v63'
    assert sha(pb/'report.json')=='879cf16978e498f79a44a8aba4d0bc1d6ef99069249d43ae0f73497fe01756c3'
    pr=read(pb/'report.json');assert pr['all_passed']
    for name in ['source-and-types.json','monomial-factor-recipes.json','height-and-output-bounds.json']:
        assert sha(pb/name)==pr['evidence_sha256'][name]
    bind=read(pb/'source-and-types.json');recipes=read(pb/'monomial-factor-recipes.json');height=read(pb/'height-and-output-bounds.json')
    kb=ROOT/'reports/S61_NONSQUARE_RAW_H_B_K13_ACTUAL_SOURCE_SCHEME_20261104_v63'
    dg=ROOT/'reports/S61_NONSQUARE_RAW_H_B_K13_D_EUCLID_20261104_v63';g=ip(read(dg/'checkpoint.json')['g13']);assert g[13]>0
    ring=Ring(g);old=[json.loads(s)for s in(kb/'monomial-nodes.jsonl').read_text().splitlines()]
    flat=ROOT/bind['existing_flat_universal_input'];assert sha(flat)==pr['input_sha256'][str(flat.relative_to(ROOT))]
    units=read(kb/'source-open-units.json');a=ring.make({1:1});B=ring.dec(units['FULLU_product_all6_factors'][3]);aB=ring.mul(a,B)
    scaleLCM=int(height['all_constant_denominators_LCM']);baseD=int(height['base_denominator_LCM']);D=baseD**3*g[13]**24
    assert D.bit_length()==height['uniform_monomial_denominator_bits']
    pre=dict(source_shape=[17,24],source_entries=408,source_numerator_terms=84253,source_is_EXTERNAL_no_input_copy=True,
        known_base_nodes=272,monomial_cache_entries=1113,common_integer_denominator_bits=D.bit_length(),
        matrix_RAW_bound=height['cleared_matrix_RAW_JSON_bytes_upper'],hard_JOB_TOTAL_artifact_bytes=32*1024*1024,
        hard_AS_bytes=2000000000,hard_seconds=60,hard_explicit_terms=150000,
        no_serialized_monomial_cache_or_old_witness_duplication=True,
        source_identity_recipe='c*Cij=scaleLCM*(aB/u)*Pij inKa13; sourcePij and monomialfactors are pinned arrays. Exact integer quotient is checked and hashed, not duplicated; all operands/longdivision uniquely specified.',
        no_rank_inverse_search_modularprime_or_next_job=True)
    (out/'preflight.json').write_text(json.dumps(pre,indent=2)+'\n')
    completed=0;mon={};scaled={};peak_terms=0;maxbits=0;identities=[];profiles=[];status='pending'
    mpath=out/'matrix17x24.jsonl'
    def guard(extra=0):
        n=sum(p.stat().st_size for p in out.iterdir()if p.is_file())
        if n+extra+2*1024*1024>32*1024*1024:raise Stop('TOTAL_artifact_size_guard',dict(existing=n,extra=extra))
        if peak_terms>150000:raise Stop('live_explicit_terms_guard',dict(terms=peak_terms))
    def checkpoint(status,extra=None):
        (out/'checkpoint.json').write_text(json.dumps(dict(status=status,completed_entries=completed,completed_rows=completed//24,
            monomial_cache_entries=len(mon),max_canonical_coefficient_bits=maxbits,peak_live_explicit_terms=peak_terms,
            field_integer_operations=ring.operations,seconds=time.monotonic()-start,extra=extra),indent=2)+'\n')
    try:
        bases={j:ring.dec(old[j]['value'])for rec in recipes for j in rec['accepted_node_factors']}
        assert len(bases)==272
        for rec in recipes:
            ids=rec['accepted_node_factors'];value=bases[ids[0]]
            for j in ids[1:]:value=ring.mul(value,bases[j])
            assert D%value[1]==0;m=tuple(rec['exponent']);mon[m]=value;scaled[m]={i:v*(D//value[1])for i,v in value[0].items()}
            peak_terms=max(peak_terms,sum(len(v[0])for v in bases.values())+sum(len(v[0])for v in mon.values())+sum(map(len,scaled.values())))
        checkpoint('all1113_monomials_ready_INTERNAL_cache_only');guard()
        multipliers={'1':aB,'a':B,'B':a};unitkind={0:'1',1:'a',2:'B'}
        with flat.open()as f:
            assert f.readline().split()==['17','24']
            for index in range(408):
                kind,c,count=map(int,f.readline().split());expected=bind['entries'][index];tag=unitkind[kind]
                assert tag==expected['denominator_unit']and c==int(expected['denominator_integer_scale'])and count==expected['numerator_terms']
                raw={};source=[]
                for _ in range(count):
                    v,i,j,k=map(int,f.readline().split());m=(i,j,k);assert m in scaled;raw=pa(raw,scaled[m],v);source.append([[i,j,k],str(v)])
                assert hashlib.sha256(json.dumps(source,separators=(',',':')).encode()).hexdigest()==expected['numerator_coefficients_SHA256']
                p=ring.make(raw,D);mult=multipliers[tag];z=ring.scale(ring.mul(p,mult),scaleLCM//c);assert scaleLCM%c==0
                # Integer proof of c*z = scaleLCM*mult*p. All operands are
                # explicit saved/pinned arrays; the exact quotient is replayable.
                lhs={i:v*c*p[1]*mult[1]for i,v in z[0].items()}
                rhs={i:v*scaleLCM*z[1]for i,v in pm(p[0],mult[0]).items()};q=exactdiv(pa(lhs,rhs,-1),g)
                rec=dict(index=index,row17=index//24,column24=index%24,source_row25=expected['source_row25'],
                    target_column32=expected['target_column32'],value=ring.enc(z))
                text=json.dumps(rec,separators=(',',':'))+'\n';guard(len(text.encode()))
                with mpath.open('a')as o:o.write(text)
                identities.append(dict(index=index,source_numerator_SHA256=expected['numerator_coefficients_SHA256'],denominator_type=tag,
                    source_integer_denominator=str(c),clearing_integer_scale=str(scaleLCM//c),
                    source_identity='c*Cnum*den(Peval)*den(mult)-scaleLCM*den(C)*num(Peval)*num(mult)=g13*Q',
                    Peval_recipe='raw=sum(integer_source_coefficient*scaled_monomial);P=canonical(raw/D),degree<13',
                    multiplier_recipe={'1':'a*B','a':'B','B':'a'}[tag],
                    exact_integer_division_remainder_zero=True,integer_quotient_degree=max(q,default=-1),
                    integer_quotient_SHA256=hashlib.sha256(json.dumps(ser(q),separators=(',',':')).encode()).hexdigest()))
                nb=max((abs(v).bit_length()for v in z[0].values()),default=0);db=z[1].bit_length();maxbits=max(maxbits,nb,db)
                profiles.append(dict(nonzero=bool(z[0]),QQ_terms=len(z[0]),numerator_bits=nb,denominator_bits=db))
                completed+=1
                peak_terms=max(peak_terms,sum(len(v[0])for v in bases.values())+sum(len(v[0])for v in mon.values())+sum(map(len,scaled.values()))+count+len(q)+len(z[0]))
                checkpoint('entry_source_identity_and_value_saved');guard()
            assert not f.read().strip()
        payload=json.dumps(dict(source_flat_input=str(flat.relative_to(ROOT)),source_SHA256=sha(flat),
            shared_monomial_factor_recipes=str((pb/'monomial-factor-recipes.json').relative_to(ROOT)),
            uniform_D=str(D),scale_LCM=str(scaleLCM),g13=ser(g),multipliers={k:ring.enc(v)for k,v in multipliers.items()},
            all408_source_identities=identities,no_unknown_scalar_solve_or_cofactor_choice=True),indent=2)+'\n'
        guard(len(payload.encode()));(out/'source-identities.json').write_text(payload)
        status='all408_actual_coefficients_source_checked_NO_rank_test';checkpoint(status)
    except Stop as e:status=e.status;checkpoint(status,e.data)
    complete=status=='all408_actual_coefficients_source_checked_NO_rank_test'
    proof=dict(status=status,coefficient_algebra='Ka13=QQ[a]/g13',matrix='C=cLCM*a*B*M17x24',
        scalar_known_unit='cLCM is nonzeroQQ; a andB are units by accepted source-open-units FULLU inverse. Thus C and actualM have identical rank/module-surjectivity atALLKa13 components.',
        unit8_minor='-a^4*F*B²*C²/4608, existingverified8prefix remainsinvertible.',
        rank_criterion='Surjectivity Ka13^24->Ka13^17 ofthismatrix wouldimplyfullcubic57 byaccepteduniversalcriterion, but NOTtestedhere.',
        no_field_irreducibility_etale_assumption_in_preparation=True,
        no_sourcepoint_prime_rank_Gaussian_POWER_or_Poisson_failure_claim=True)
    (out/'rank-scope.json').write_text(json.dumps(proof,indent=2)+'\n')
    paths=[Path(__file__),ROOT/'code/s61_nonsquare_raw_H_B_K13_named_saturation_20261104_v63.py',ROOT/f'research/notes/{ID}.md',pb/'report.json',pb/'source-and-types.json',pb/'monomial-factor-recipes.json',pb/'height-and-output-bounds.json',flat,kb/'report.json',kb/'monomial-nodes.jsonl',kb/'source-open-units.json',dg/'report.json',dg/'checkpoint.json',ROOT/'reports/S61_NONSQUARE_CUBIC_ZW25_UNITS_20261104_v63_REVIEW/report.json']
    report=dict(id=ID,all_passed=complete,status=status,entries=completed,rows=17,columns=24,
        nonzero_entries=sum(p['nonzero']for p in profiles),QQ_terms=sum(p['QQ_terms']for p in profiles),
        max_numerator_bits=max((p['numerator_bits']for p in profiles),default=0),max_denominator_bits=max((p['denominator_bits']for p in profiles),default=0),
        field_integer_operations=ring.operations,peak_live_explicit_terms=peak_terms,
        no_input_or_cache_artifact_duplication=True,no_rank_or_prime_or_newpoint_or_next_job=True,
        independent_review_pending=True,input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},
        evidence_sha256={p.name:sha(p)for p in out.iterdir()if p.is_file()},seconds=time.monotonic()-start,
        max_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    assert sum(p.stat().st_size for p in out.iterdir()if p.is_file())<32*1024*1024
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');signal.alarm(0)
    print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    try:main(a.output)
    except Exception as e:
        if a.output.exists():(a.output/'failure.json').write_text(json.dumps(dict(all_passed=False,error=str(e),completed_matrix_rows_preserved=(a.output/'matrix17x24.jsonl').exists(),no_auto_retry=True),indent=2)+'\n')
        raise
