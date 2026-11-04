#!/usr/bin/env python3
"""Independent exact408-source audit on full Ka13, no rank/prime job."""
import hashlib,importlib.util,json,math,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_RAW_H_B_K13_CUBIC17_COEFFICIENTS_20261104_v63_REVIEW'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main(out):
    started=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));out.mkdir(exist_ok=False,parents=True)
    helper=ROOT/'code/s61_nonsquare_raw_H_B_K13_d_Euclid_review_20261104_v63.py';spec=importlib.util.spec_from_file_location('own_Ka13_fixed_table',helper);ar=importlib.util.module_from_spec(spec);spec.loader.exec_module(ar)
    pb=ROOT/'reports/S61_NONSQUARE_RAW_H_B_K13_CUBIC17_COEFFICIENTS_20261104_v63';pr=read(pb/'report.json');assert pr['all_passed']
    paths={Path(__file__),helper,ROOT/'research/notes'/(ID+'.md'),pb/'report.json'}
    for rel,h in pr['input_sha256'].items():assert sha(ROOT/rel)==h,rel;paths.add(ROOT/rel)
    for rel,h in pr['evidence_sha256'].items():assert sha(pb/rel)==h,rel;paths.add(pb/rel)
    c=read(pb/'source-identities.json');g=ar.ip(c['g13']);A=ar.OwnAlgebra(g)
    pref=ROOT/'reports/S61_NONSQUARE_RAW_H_B_K13_FULL_CUBIC_PREFLIGHT_20261104_v63';bind=read(pref/'source-and-types.json');recipes=read(pref/'monomial-factor-recipes.json')
    flat=ROOT/c['source_flat_input'];assert flat==ROOT/bind['existing_flat_universal_input'] and sha(flat)==c['source_SHA256']
    identical=ROOT/'reports/S61_NONSQUARE_THETA_P100_CUBIC_GENERIC_20261104_v63/source17x24-input.txt';assert flat.read_bytes()==identical.read_bytes();paths.add(identical)
    # Reconstruct the universal pre-curve input from the actual25x32 unit checkpoint.
    checkpoint=ROOT/bind['source_checkpoint'];cp=read(checkpoint);paths.add(checkpoint)
    rows,cols=cp['remaining_rows'],cp['remaining_columns'];assert rows==bind['source_rows25'] and cols==bind['target_columns32'] and(len(rows),len(cols))==(17,24)
    unitpolys={'1':{(0,0,0):Q(1)},'a':{(1,0,0):Q(1)},'B':{(1,0,1):Q(50),(0,2,1):Q(18),(0,1,1):Q(-9),(0,0,1):Q(-2),(0,1,0):Q(60),(0,0,0):Q(-40)}}
    lines=['17 24'];numerators=[];denominators=[];termcount=0
    for index,(i,j)in enumerate((i,j)for i in rows for j in cols):
        p,q=cp['matrix25x32'][i][j];b=bind['entries'][index];tag=b['denominator_unit'];scale=int(b['denominator_integer_scale']);assert scale
        assert(i,j)==(b['source_row25'],b['target_column32']) and len(p)==b['numerator_terms']
        assert {tuple(m):Q(v)for m,v in q}=={m:v*scale for m,v in unitpolys[tag].items()}
        lines.append(f"{ {'1':0,'a':1,'B':2}[tag]} {scale} {len(p)}")
        for m,v in p:
            vv=Q(v);assert vv.denominator==1 and len(m)==3 and min(m)>=0
            lines.append(str(vv.numerator)+' '+' '.join(map(str,m)))
        assert hashlib.sha256(json.dumps(p,separators=(',',':')).encode()).hexdigest()==b['numerator_coefficients_SHA256']
        numerators.append(p);denominators.append(scale);termcount+=len(p)
    assert ('\n'.join(lines)+'\n').encode()==flat.read_bytes() and termcount==84253
    oldbase=ROOT/'reports/S61_NONSQUARE_RAW_H_B_K13_ACTUAL_SOURCE_SCHEME_20261104_v63';old=[json.loads(line)for line in(oldbase/'monomial-nodes.jsonl').read_text().splitlines()]
    oldpeer=ROOT/'reports/S61_NONSQUARE_RAW_H_B_K13_ACTUAL_SOURCE_SCHEME_20261104_v63_REVIEW/report.json';assert read(oldpeer)['all_passed'];paths.add(oldpeer)
    exps=[]
    initial={'one':(0,0,0),'a':(1,0,0),'e=t':(0,1,0),'d=s':(0,0,1)}
    for k,n in enumerate(old):
        assert n['id']==k
        if n['op']=='literal_accepted_input':e=initial[n['name']]
        else:e=tuple(a+b for a,b in zip(exps[n['left']],exps[n['right']]))
        exps.append(e)
    mon={};D=int(c['uniform_D']);assert D>0
    for rec in recipes:
        ids=rec['accepted_node_factors'];assert 1<=len(ids)<=3
        exponent=tuple(sum(exps[i][j]for i in ids)for j in range(3));assert exponent==tuple(rec['exponent'])
        value=A.dec(old[ids[0]]['value'])
        for i in ids[1:]:value=A.mul(value,A.dec(old[i]['value']))
        assert D%value[1]==0;mon[exponent]={j:v*(D//value[1])for j,v in enumerate(value[0])if v}
    assert len(mon)==1113
    units=read(oldbase/'source-open-units.json');a=A.make({1:1});B=A.dec(units['FULLU_product_all6_factors'][3]);aB=A.mul(a,B)
    multipliers={'1':aB,'a':B,'B':a}
    for tag,v in multipliers.items():assert v==A.dec(c['multipliers'][tag])
    L=int(c['scale_LCM']);assert L==math.lcm(*map(abs,denominators))
    values=[json.loads(line)for line in(pb/'matrix17x24.jsonl').read_text().splitlines()];assert len(values)==408
    for index,(p,cc,zrec,identity)in enumerate(zip(numerators,denominators,values,c['all408_source_identities'])):
        raw={}
        for m,v in p:raw=ar.plus(raw,mon[tuple(m)],int(v))
        P=A.make(raw,D);tag=bind['entries'][index]['denominator_unit'];mult=multipliers[tag];assert L%cc==0
        z=A.mul(A.mul(P,mult),A.c(Q(L,cc)));assert z==A.dec(zrec['value'])
        assert(zrec['index'],zrec['row17'],zrec['column24'])==(index,index//24,index%24)
        assert(zrec['source_row25'],zrec['target_column32'])==(rows[index//24],cols[index%24])
        lhs={i:v*cc*P[1]*mult[1]for i,v in enumerate(z[0])if v}
        pnum={i:v for i,v in enumerate(P[0])if v};mnum={i:v for i,v in enumerate(mult[0])if v}
        rhs={i:v*L*z[1]for i,v in ar.times(pnum,mnum).items()}
        quotient=ar.exactdiv(ar.plus(lhs,rhs,-1),g)
        ser=[[i,str(v)]for i,v in sorted(quotient.items(),reverse=True)]
        assert hashlib.sha256(json.dumps(ser,separators=(',',':')).encode()).hexdigest()==identity['integer_quotient_SHA256']
        assert max(quotient,default=-1)==identity['integer_quotient_degree']
    prefix=ROOT/'reports/S61_NONSQUARE_CUBIC_ZW25_UNITS_20261104_v63_REVIEW/report.json';assert read(prefix)['all_passed'];paths.add(prefix)
    proof={'all_passed':True,'arithmetic':'Own fixed-table Ka13 arithmetic and ordinary integer long division, no author Ring/producer.',
      'all84253_universal_numerator_coefficients_direct_checkpoint_to_flat_byte_match':True,
      'P80_and_P100_named_flat_files_identical_UNSPECIALIZED_a_e_d_not_zeta':True,
      'all408_denominators_exact_c_times_1_a_B':True,'all1113_monomial_factor_exponents_and_values':True,
      'all408_Ka13_cleared_values_and_integer_quotient_hashes':True,
      'matrix_relation':'C=cLCM*a*B*M17x24. cLCM nonzeroQQ and a,B are whole-Ka13 units by accepted source/open proof. The old8prefix determinant−a⁴FB²C²/4608 is a QQ/FULLU unit. Thus rank17 of this actual C is sufficient for actual fullcubic57 on this entire finite bad algebra.',
      'no_rank_or_modular_inference_in_this_source_gate':True,'new_rank_inverse_GCD_point_prime_GB_jobs':0,'quotient_arithmetic_operations':A.ops}
    (out/'independent-proof.json').write_text(json.dumps(proof,indent=2)+'\n')
    report={'id':ID,'all_passed':True,'all408_actual_coefficients_independently_source_bound':True,'source_universal_not_curve_specialized':True,'no_rank_job':True,
      'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in sorted(paths)},'evidence_sha256':{'independent-proof.json':sha(out/'independent-proof.json')},
      'seconds':time.monotonic()-started,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'frozen':True}
    assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024;signal.alarm(0)
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2),flush=True)
if __name__=='__main__':
    import argparse
    pa=argparse.ArgumentParser();pa.add_argument('--output',type=Path,default=ROOT/'reports'/ID);main(pa.parse_args().output.resolve())
