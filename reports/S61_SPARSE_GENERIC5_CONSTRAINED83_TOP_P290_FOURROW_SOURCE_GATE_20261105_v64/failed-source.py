#!/usr/bin/env python3
"""SOURCE only:16 changed top cells, exact support/height/whole-job bounds."""
from pathlib import Path
from fractions import Fraction as Q
from math import gcd
import hashlib,json,resource,signal,time
from s61_sparse_generic5_primitive_binary_codec_20261105_v64 import decode,encode
ROOT=Path(__file__).resolve().parents[1]
ID='S61_SPARSE_GENERIC5_CONSTRAINED83_TOP_P290_FOURROW_SOURCE_GATE_20261105_v64'
def sha(p):
    h=hashlib.sha256()
    with p.open('rb')as f:
        for b in iter(lambda:f.read(65536),b''):h.update(b)
    return h.hexdigest()
def read(p):return json.loads(p.read_text())
def enc(v):return(json.dumps(v,separators=(',',':'))+'\n').encode()
def numeric(v):
    if isinstance(v,bool):return 0
    if isinstance(v,(int,float)):return 1
    if isinstance(v,list):return sum(map(numeric,v))
    if isinstance(v,dict):return sum(map(numeric,v.values()))
    return 0
def supp(p):return{tuple(p['shift'][j]+t[j+1]for j in range(5))for t in p['terms']}if p else set()
def stats(p):
    if not p or not p['terms']:return 0,Q(0),0
    return len(p['terms']),Q(p['content']),max(abs(int(t[0]))for t in p['terms'])
def main():
    begin=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000))
    out=ROOT/'reports'/ID;out.mkdir(parents=True,exist_ok=False);(out/'supports').mkdir()
    def gate(stem,h):
        p=ROOT/f'reports/{stem}/report.json';assert sha(p)==h;r=read(p);assert r['all_passed'];return p,r
    cp,c=gate('S61_SPARSE_GENERIC5_CONSTRAINED83_CACHED_KERNEL_COMMON_SHAPE_20261105_v64','3bc324e022f1b5f0666f9eb2082c5e265335e06f6b27ea36b90be52d302d55ea')
    sr,s=gate('S61_SPARSE_GENERIC5_CONSTRAINED83_STREAM49_ACTIVE_FRAGMENT_MIXED_BINARY_COMPACTION_20261105_v64','62a95473e52e2a85054586bd63e28b79ae58850bebf49d6e5a0a4daeb69e7970')
    tr,t=gate('S61_SPARSE_GENERIC5_CONSTRAINED83_FACTORED_CONSTRAINT_SCHUR_PREFLIGHT_20261105_v64','8d0c175945ac170d8275da8fb2ecd1b7673d40dfd6c3bdc035f54d58b75912bf')
    gp=cp.parent/'fourteen-source-denominators-and-factored-P290.json';pp=cp.parent/'actual-common-denominator-P290.json'
    for p in(gp,pp):assert sha(p)==c['evidence_sha256'][p.name]
    common=read(gp);P=read(pp);assert len(P['terms'])==290
    ip=sr.parent/s['current_index'];matrix=sr.parent/s['current_snapshot'];assert sha(ip)==s['current_index_sha256']and sha(matrix)==s['current_snapshot_sha256'];ix=read(ip)
    costp=tr.parent/'row-stream-cost-source-only.json';assert sha(costp)==t['evidence_sha256'][costp.name];oldcost=read(costp)
    eligible=sorted([a for a in oldcost['per_row']if a['row']<134],key=lambda a:(a['lazy_numerator_pair_majorant_only'],a['current_input_terms'],a['row']))
    chosen=[a['row']for a in eligible[:4]];assert chosen==[52,11,44,51]
    rp=ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_EXACT49_FIXED_CONSTRAINT_MINOR_PREFLIGHT_20261105_v64/two-actual-fixed-constraint-minors.json';small=read(rp)
    paths={cp,sr,tr,gp,pp,ip,matrix,costp,rp,Path(__file__),ROOT/f'research/notes/{ID}.md',ROOT/'code/s61_sparse_generic5_primitive_binary_codec_20261105_v64.py'}
    # Charge all actual proof inputs read by a future provenance checker too.
    for p,h in c['input_sha256'].items():assert sha(ROOT/p)==h;paths.add(ROOT/p)
    numerators={};recipes={};source_rows={}
    for a in common['all14_fixed_denominator_equalities']:
        rec=a['common_numerator_recipe'];d=rec['source_row_payload'];path=ROOT/d['path'];paths.add(path)
        if d['row_original0']not in source_rows:
            with path.open('rb')as f:f.seek(d['offset']);payload=f.read(d['bytes'])
            assert hashlib.sha256(payload).hexdigest()==d['payload_sha256'];source_rows[d['row_original0']]=json.loads(payload)
        v=source_rows[d['row_original0']][a['component']]['numerator']
        numerators[(a['row_original0'],a['component'])]={'content':rec['QQ_scalar'],'shift':rec['parameter_exponent'],'terms':v}
        recipes[(a['row_original0'],a['component'])]=a
    del source_rows
    cells=ix['all_literal_cells'];lookup={(v['row'],v['column']):(k,v)for k,v in enumerate(cells)}
    source_base=numeric(ix)+numeric(numerators)+numeric(recipes)+numeric(common)+numeric(small)+3*len(lookup)+4313+290+12288
    source_peak=source_base+5*4060;assert source_peak<150000;future_peak=0;plans=[];caps=0;coeff_pair_ops=0;METADATA=25000
    for row in chosen:
        if time.monotonic()-begin>45:raise RuntimeError('SOURCE45s stop, complete prior cell certificates preserved, no coefficients evaluated')
        cur={};refs={}
        needed={v for b in small['blocks']for v in b['pivot_columns']+b['free_columns']}
        with matrix.open('rb')as f:
            for col in sorted(needed):
                if(row,col)not in lookup:continue
                k,v=lookup[(row,col)];f.seek(v['payload_offset']);payload=f.read(v['payload_bytes']);assert hashlib.sha256(payload).hexdigest()==v['payload_sha256'];p=decode(payload);assert encode(p)==payload;cur[col]=p;refs[col]=k
        row_terms=sum(len(p['terms'])for p in cur.values());retained=4313+290+row_terms
        for block in small['blocks']:
            origin=1 if block['name']=='alpha'else 21
            for j,col in enumerate(block['free_columns']):
                operands=[{'kind':'P_times_oldfree','a':cur.get(col),'b':P,'source_column':col}]
                for pivot in block['pivot_columns']:
                    operands.append({'kind':'oldpivot_times_cached_numerator','a':cur.get(pivot),'b':numerators[(pivot-origin,j)],'source_column':pivot,'numerator_key':[pivot-origin,j]})
                union=set();base=[0]*5;inter=[];D=1;M=0;active=[];resident=retained+METADATA
                for op in operands:
                    na,qa,ma=stats(op['a']);nb,qb,mb=stats(op['b'])
                    if na and nb:active.append((op,na,nb,qa*qb,ma,mb))
                for op,na,nb,qq,ma,mb in active:
                    D=D*qq.denominator//gcd(D,qq.denominator)
                for step,(op,na,nb,qq,ma,mb)in enumerate(active):
                    aa,bb=supp(op['a']),supp(op['b']);product={tuple(x[k]+y[k]for k in range(5))for x in aa for y in bb};union.update(product)
                    source_peak=max(source_peak,source_base+numeric(cur)+row_terms+numeric(plans)+numeric(inter)+5*(len(aa)+len(bb)+len(product)+2*len(union))+8192);assert source_peak<150000
                    bnd=resident+2*len(product)+3*len(union)+256;future_peak=max(future_peak,bnd)
                    for k in range(5):base[k]=min(base[k],op['a']['shift'][k]+op['b']['shift'][k])
                    M+=abs(qq.numerator)*(D//qq.denominator)*min(na,nb)*ma*mb;coeff_pair_ops+=na*nb
                    inter.append({'step':step,'kind':op['kind'],'source_column':op['source_column'],'source_cell_index':refs[op['source_column']],'source_cell_reference':lookup[(row,op['source_column'])][1],
                        'cached_numerator_key':op.get('numerator_key'),'n_a':na,'n_b':nb,'product_support_upper':len(product),'union_support_upper':len(union),
                        'whole_owned_scratch_upper':bnd,'QQ_product_content':str(qq),'integer_accumulator_coefficient_abs_upper':str(M)})
                    assert bnd<150000
                    del aa,bb,product
                exps=[max((e[k]-base[k]for e in union),default=0)for k in range(5)];cap=2048+len(str(M))+len(str(D))+len(union)*(len(str(M))+sum(len(str(v))for v in exps)+48)
                source_peak=max(source_peak,source_base+numeric(cur)+row_terms+numeric(plans)+numeric(inter)+10*len(union)+8704);assert source_peak<150000
                sp=out/'supports'/f'row-{row}-col-{col}.json';sp.write_bytes(enc(sorted(union)));caps+=cap
                plans.append({'row':row,'column':col,'block':block['name'],'common_denominator_source':str(pp.relative_to(ROOT)),
                    'exact_integer_scalar_common_denominator':str(D),'integer_numerator_coefficient_abs_upper':str(M),
                    'actual_accumulator_shift':base,'output_support_upper':len(union),'support_file':str(sp.relative_to(ROOT)),'support_sha256':sha(sp),
                    'RAW_literal_polynomial_upper':cap,'per_cell_input_terms':retained,'source_row_input_terms':row_terms,'steps':inter,
                    'formula':'P290*old_free+sum5(old_pivot*rescaled_cached_numerator)','coefficient_products_evaluated':False})
                source_peak=max(source_peak,source_base+numeric(cur)+row_terms+numeric(plans)+10*len(union)+8192);assert source_peak<150000
                del union,inter,active
        del cur
    assert len(plans)==16
    evidence={'all_passed':True,'source_only':True,'selected_rows':chosen,'free_columns':[7,13,27,33],'all16_cell_plans':plans,
        'literal_coefficient_source_file':str(matrix.relative_to(ROOT)),'literal_coefficient_source_sha256':sha(matrix),'all20_other_columns_immutable_source_references':True,'unchanged66_constraint_rows':'Retain the actual formal roots and rank-basis zero proof from the177x24 operator. No zero-fill or row deletion.',
        'column_scaling':'Multiply ONLY the four free columns of the Schur residual by nonzero P290. This scales those unknowns inversely and is a K5 target unit operation. Other20 coefficient columns remain literal.',
        'cokernel':'Full177x24 residual is retained; dropping its66zero rows removes K5^66. Kernel/column rank can use the111top. No complete inverse.',
        'chart':'Prior parameter-monomial units and P290*detA*detB nonzero; complements remain open.',
        'no_all111_output_unsafe_claim':True,'no_scalar_product_worker_or_new_rank_search':True}
    ep=out/'sixteen-source-support-height-cell-plans.json';ep.write_bytes(enc(evidence))
    active=[{'path':str(p.relative_to(ROOT)),'sha256':sha(p),'bytes':p.stat().st_size}for p in sorted(paths)]
    actual_input=sum(p['bytes']for p in active);sourcefiles=sum(p.stat().st_size for p in out.rglob('*')if p.is_file());reserve=3*1024*1024
    total=actual_input+sourcefiles+caps+reserve;safe=future_peak<=150000 and total<=32*1024*1024
    budget={'all_passed':True,'source_only':True,'FULL_safe_sixteen_cell_coefficient_job':safe,'actual_READ_inputs':active,
        'actual_READ_input_RAW_bytes':actual_input,'all_source_plan_and_support_files_RAW_bytes':sourcefiles,'all16_output_RAW_upper':caps,
        'worker_code_executable_receipt_control_reserve':reserve,'fresh_RAW_upper':total,'whole_owned_with_fixed_metadata_upper':future_peak,
        'fixed_future_metadata_upper':METADATA,'future_wrapper_obligation':'Use the thin row/cell references and old-kernel payload offsets in these source plans; hash complete provenance files by streaming. Do not parse the full current index or whole old kernel in the coefficient worker. Release source-wrapper arrays before native.','source_owned_peak':source_peak,'coefficient_term_pair_operations_upper':coeff_pair_ops,
        'hard_seconds':60,'soft_cell_boundary_seconds':45,'hard_address_space_bytes':2000000000,
        'only_one_current_row_one_output_plus4313cachedNumeratorTerms290P':True,
        'no_coefficient_job_before_parent_GO':True,'all_prior_source_files_and_fullmatrix_versions_retained':True}
    bp=out/'whole-budget.json';bp.write_bytes(enc(budget))
    report={'id':ID,'all_passed':True,'source_only':True,'selected_rows':chosen,'cells':16,'FULL_safe_coefficient_job':safe,
        'source_owned_peak':source_peak,'future_owned_upper':future_peak,'future_RAW_upper':total,'coefficient_pair_upper':coeff_pair_ops,
        'no_coefficient_products_evaluated':True,'no_all111_or_fullinverse_claim':True,'input_sha256':{p['path']:p['sha256']for p in active},
        'evidence_sha256':{str(p.relative_to(out)):sha(p)for p in out.rglob('*')if p.is_file()},'seconds':time.monotonic()-begin}
    rp2=out/'report.json';rp2.write_bytes(enc(report));assert actual_input+sum(p.stat().st_size for p in out.rglob('*')if p.is_file())<32*1024*1024
    signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in('input_sha256','evidence_sha256')}));print('report_sha256',sha(rp2))
if __name__=='__main__':main()
