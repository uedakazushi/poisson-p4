#!/usr/bin/env python3
"""SOURCE only:16 changed top cells, exact support/height/whole-job bounds."""
from pathlib import Path
from fractions import Fraction as Q
from math import gcd
import argparse,hashlib,json,resource,signal,time
from s61_sparse_generic5_primitive_binary_codec_20261105_v64 import decode,encode
ROOT=Path(__file__).resolve().parents[1]

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
def main(a):
    ID=a.id
    begin=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000))
    out=ROOT/'reports'/ID;out.mkdir(parents=True,exist_ok=False);(out/'supports').mkdir();(out/'cell-plans').mkdir()
    cursor={'phase':'source-input'}
    def guard(v,kind='source_owned'):
        if v>=150000:
            (out/'guard-stop.json').write_bytes(enc({'all_passed':False,'source_only':True,'guard':kind,'attempted_bound':v,'limit':150000,'cursor':dict(cursor),'complete_prior_cell_forecasts_preserved':True,'no_coefficient_products_evaluated':True}));raise RuntimeError('exact source guard stop')
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
    hr=ROOT/a.prior_report;assert sha(hr)==a.prior_sha;h=read(hr);assert h['all_passed']
    hp=hr.parent/'evaluated16-cells-and-full177-frame.json';assert sha(hp)==h['evidence_sha256'][hp.name];history=read(hp);history_cells=history.get('all_evaluated_cells',history.get('evaluated_cells16',[]));done={z['row']for z in history_cells};assert len(history_cells)==4*len(done)
    chosen=[a['row']for a in eligible if a['row']not in done][:4];assert len(chosen)==4
    rp=ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_EXACT49_FIXED_CONSTRAINT_MINOR_PREFLIGHT_20261105_v64/two-actual-fixed-constraint-minors.json';small=read(rp)
    paths={cp,sr,tr,hr,hp,gp,pp,ip,matrix,costp,rp,Path(__file__),ROOT/f'research/notes/{ID}.md',ROOT/'code/s61_sparse_generic5_primitive_binary_codec_20261105_v64.py'}
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
    parse_phase_owned=numeric(ix)+numeric(source_rows)+numeric(common)+numeric(small)+4313+4060+8192;assert parse_phase_owned<150000
    del source_rows
    cells=ix['all_literal_cells'];lookup={(v['row'],v['column']):(k,v)for k,v in enumerate(cells)}
    needed_columns={v for b in small['blocks']for v in b['pivot_columns']+b['free_columns']}
    thin={key:value for key,value in lookup.items()if key[0]in chosen and key[1]in needed_columns}
    del ix,cells,lookup,recipes,eligible,oldcost
    source_base=numeric(thin)+numeric(numerators)+numeric(common)+numeric(small)+3*len(thin)+4313+290+12288
    source_peak=max(parse_phase_owned,source_base+5*4060);guard(source_peak);future_peak=0;plans=[];caps=0;coeff_pair_ops=0;METADATA=25000
    for row in chosen:
        cursor={'phase':'row','row':row}
        if time.monotonic()-begin>45:raise RuntimeError('SOURCE45s stop, complete prior cell certificates preserved, no coefficients evaluated')
        cur={};refs={}
        needed={v for b in small['blocks']for v in b['pivot_columns']+b['free_columns']}
        with matrix.open('rb')as f:
            for col in sorted(needed):
                if(row,col)not in thin:continue
                k,v=thin[(row,col)];f.seek(v['payload_offset']);payload=f.read(v['payload_bytes']);assert hashlib.sha256(payload).hexdigest()==v['payload_sha256'];p=decode(payload);assert encode(p)==payload;cur[col]=p;refs[col]=k
        row_terms=sum(len(p['terms'])for p in cur.values());retained=4313+290+row_terms
        for block in small['blocks']:
            origin=1 if block['name']=='alpha'else 21
            for j,col in enumerate(block['free_columns']):
                cursor={'phase':'cell','row':row,'column':col}
                operands=[{'kind':'P_times_oldfree','a':cur.get(col),'b':P,'source_column':col}]
                for pivot in block['pivot_columns']:
                    operands.append({'kind':'oldpivot_times_cached_numerator','a':cur.get(pivot),'b':numerators[(pivot-origin,j)],'source_column':pivot,'numerator_key':[pivot-origin,j]})
                union=set();base=[0]*5;inter=[];D=1;M=0;active=[];rowbase=290+row_terms;onecache=0
                for op in operands:
                    na,qa,ma=stats(op['a']);nb,qb,mb=stats(op['b'])
                    if na and nb:active.append((op,na,nb,qa*qb,ma,mb))
                for op,na,nb,qq,ma,mb in active:
                    D=D*qq.denominator//gcd(D,qq.denominator)
                ranges=[]
                for op,*_ in active:
                    lo=[];hi=[]
                    for k in range(5):
                        lo.append(min(op['a']['shift'][k]+t[k+1]for t in op['a']['terms'])+min(op['b']['shift'][k]+t[k+1]for t in op['b']['terms']))
                        hi.append(max(op['a']['shift'][k]+t[k+1]for t in op['a']['terms'])+max(op['b']['shift'][k]+t[k+1]for t in op['b']['terms']))
                    ranges.append((lo,hi))
                lower=[min(z[0][k]for z in ranges)for k in range(5)];upper=[max(z[1][k]for z in ranges)for k in range(5)];width=[upper[k]-lower[k]+1 for k in range(5)];weights=[1]*5
                for k in range(3,-1,-1):weights[k]=weights[k+1]*width[k+1]
                volume=width[0]*weights[0];assert volume<2**63;bias=sum(lower[k]*weights[k]for k in range(5))
                code_abs=sum(max(abs(lower[k]),abs(upper[k]))*weights[k]for k in range(5))
                for op,*_ in active:
                    pair_abs=0
                    for poly in(op['a'],op['b']):pair_abs+=sum(max(abs(poly['shift'][k]+t[k+1])for t in poly['terms'])*weights[k]for k in range(5))
                    code_abs=max(code_abs,pair_abs)
                assert code_abs<2**63
                def packed(p):return{sum((p['shift'][k]+t[k+1])*weights[k]for k in range(5))for t in p['terms']}
                def unpack(code):
                    u=code-bias;assert 0<=u<volume;e=[]
                    for k in range(5):v,u=divmod(u,weights[k]);assert v<width[k];e.append(v+lower[k])
                    assert u==0;return e
                for step,(op,na,nb,qq,ma,mb)in enumerate(active):
                    cursor={'phase':'support-step','row':row,'column':col,'step':step}
                    aa,bb=packed(op['a']),packed(op['b']);assert len(aa)==na and len(bb)==nb;product={x+y for x in aa for y in bb};union.update(product)
                    source_peak=max(source_peak,source_base+numeric(cur)+row_terms+numeric(plans)+numeric(inter)+(len(aa)+len(bb)+len(product)+2*len(union))+8192);guard(source_peak)
                    cache=0 if op['kind']=='P_times_oldfree'else nb;onecache=max(onecache,cache);resident=rowbase+cache+METADATA;update_bound=resident+2*len(product)+3*len(union)+256;parse_bound=rowbase+len(union)+4*cache+256+METADATA;bnd=max(update_bound,parse_bound);future_peak=max(future_peak,bnd)
                    for k in range(5):base[k]=min(base[k],op['a']['shift'][k]+op['b']['shift'][k])
                    M+=abs(qq.numerator)*(D//qq.denominator)*min(na,nb)*ma*mb;coeff_pair_ops+=na*nb
                    inter.append({'step':step,'kind':op['kind'],'source_column':op['source_column'],'source_cell_index':refs[op['source_column']],'source_cell_reference':thin[(row,op['source_column'])][1],
                        'cached_numerator_key':op.get('numerator_key'),'n_a':na,'n_b':nb,'product_support_upper':len(product),'union_support_upper':len(union),
                        'whole_owned_scratch_upper':bnd,'single_cache_input_terms':rowbase+cache,'cached_factor_parser_upper':parse_bound,'QQ_product_content':str(qq),'integer_accumulator_coefficient_abs_upper':str(M)})
                    guard(bnd,'future_owned')
                    del aa,bb,product
                exps=[upper[k]-base[k]for k in range(5)];cap=2048+len(str(M))+len(str(D))+len(union)*(len(str(M))+sum(len(str(v))for v in exps)+48)
                source_peak=max(source_peak,source_base+numeric(cur)+row_terms+numeric(plans)+numeric(inter)+2*len(union)+8704);guard(source_peak)
                sp=out/'supports'/f'row-{row}-col-{col}.json';tmp=sp.with_suffix('.json.partial')
                with tmp.open('wb')as f:
                    f.write(b'[');first=True
                    for code in sorted(union):
                        if not first:f.write(b',')
                        first=False;f.write(json.dumps(unpack(code),separators=(',',':')).encode())
                    f.write(b']\n')
                tmp.rename(sp);caps+=cap
                binary_phase=rowbase+4*len(union)+256+METADATA;future_peak=max(future_peak,binary_phase)
                plans.append({'row':row,'column':col,'block':block['name'],'common_denominator_source':str(pp.relative_to(ROOT)),
                    'exact_integer_scalar_common_denominator':str(D),'integer_numerator_coefficient_abs_upper':str(M),
                    'actual_accumulator_shift':base,'output_support_upper':len(union),'support_file':str(sp.relative_to(ROOT)),'support_sha256':sha(sp),
                    'RAW_literal_polynomial_upper':cap,'per_cell_input_terms':rowbase+onecache,'binary_roundtrip_phase_upper':binary_phase,'packed_support_box':{'lower':lower,'upper':upper,'width':width,'weights':weights,'volume':volume,'bias':bias,'packed_signed_abs_upper':code_abs,'packed_signed_bits_upper':code_abs.bit_length()},'source_row_input_terms':row_terms,'steps':inter,
                    'formula':'P290*old_free+sum5(old_pivot*rescaled_cached_numerator)','coefficient_products_evaluated':False})
                (out/'cell-plans'/f'row-{row}-col-{col}.json').write_bytes(enc(plans[-1]))
                source_peak=max(source_peak,source_base+numeric(cur)+row_terms+numeric(plans)+2*len(union)+8192);guard(source_peak)
                del union,inter,active
        del cur
    assert len(plans)==16
    evidence={'all_passed':True,'source_only':True,'selected_rows':chosen,'free_columns':[7,13,27,33],'all16_cell_plans':plans,
        'literal_coefficient_source_file':str(matrix.relative_to(ROOT)),'literal_coefficient_source_sha256':sha(matrix),'all20_other_columns_immutable_source_references':True,'unchanged66_constraint_rows':'Retain the actual formal roots and rank-basis zero proof from the177x24 operator. No zero-fill or row deletion.',
        'column_scaling':'Multiply ONLY the four free columns of the Schur residual by nonzero P290. This scales those unknowns inversely and is a K5 target unit operation. Other20 coefficient columns remain literal.',
        'cokernel':'Full177x24 residual is retained; dropping its66zero rows removes K5^66. Kernel/column rank can use the111top. No complete inverse.',
        'chart':'Prior parameter-monomial units and P290*detA*detB nonzero; complements remain open.',
        'previous16_coefficient_cells_unchanged':str(hr.relative_to(ROOT)),'no_all111_output_unsafe_claim':True,'no_scalar_product_worker_or_new_rank_search':True}
    ep=out/'sixteen-source-support-height-cell-plans.json';ep.write_bytes(enc(evidence))
    active=[{'path':str(p.relative_to(ROOT)),'sha256':sha(p),'bytes':p.stat().st_size}for p in sorted(paths)]
    actual_input=sum(p['bytes']for p in active);sourcefiles=sum(p.stat().st_size for p in out.rglob('*')if p.is_file());reserve=3*1024*1024
    total=actual_input+sourcefiles+caps+reserve;prep_upper=numeric(plans)+numeric(common)+numeric(P)+6*max(z['output_support_upper']for z in plans)+8192;future_peak=max(future_peak,prep_upper);safe=future_peak<=56473 and total<=32*1024*1024
    budget={'all_passed':True,'source_only':True,'FULL_safe_sixteen_cell_coefficient_job':safe,'actual_READ_inputs':active,
        'actual_READ_input_RAW_bytes':actual_input,'all_source_plan_and_support_files_RAW_bytes':sourcefiles,'all16_output_RAW_upper':caps,
        'worker_code_executable_receipt_control_reserve':reserve,'fresh_RAW_upper':total,'whole_owned_with_fixed_metadata_upper':future_peak,'coefficient_owned_hard_cap':56473,
        'fixed_future_metadata_upper':METADATA,'future_wrapper_obligation':'Use the thin row/cell references and old-kernel payload offsets in these source plans; hash complete provenance files by streaming. Do not parse the full current index or whole old kernel in the coefficient worker. Release source-wrapper arrays before native.','source_owned_peak':source_peak,'coefficient_term_pair_operations_upper':coeff_pair_ops,
        'hard_seconds':60,'soft_cell_boundary_seconds':45,'hard_address_space_bytes':2000000000,
        'cache_policy':'Keep290P terms and one current row. Read ONE cached numerator by authenticated controller offset per product, then release it. All14 cached values/fingerprints are unchanged; pair counts unchanged. No4313-term resident cache.','future_preparation_owned_upper':prep_upper,'packed_support_representation_exact':'Linear mixed-radix map injective on explicit cell bounding box; product exponents map to sums of codes. volume<2^63; semantic support files stream exact unpacked5-vectors.',
        'prior_source_guard_stop_preserved':'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_TOP_P290_FOURROW_SOURCE_GATE_20261105_v64/report.json','thin_cell_descriptor_count':len(thin),'no_coefficient_job_before_parent_GO':True,'all_prior_source_files_and_fullmatrix_versions_retained':True}
    bp=out/'whole-budget.json';bp.write_bytes(enc(budget))
    report={'id':ID,'all_passed':True,'source_only':True,'selected_rows':chosen,'cells':16,'FULL_safe_coefficient_job':safe,
        'source_owned_peak':source_peak,'future_owned_upper':future_peak,'future_RAW_upper':total,'coefficient_pair_upper':coeff_pair_ops,
        'no_coefficient_products_evaluated':True,'no_all111_or_fullinverse_claim':True,'input_sha256':{p['path']:p['sha256']for p in active},
        'evidence_sha256':{str(p.relative_to(out)):sha(p)for p in out.rglob('*')if p.is_file()},'seconds':time.monotonic()-begin}
    rp2=out/'report.json';rp2.write_bytes(enc(report));assert actual_input+sum(p.stat().st_size for p in out.rglob('*')if p.is_file())<32*1024*1024
    signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in('input_sha256','evidence_sha256')}));print('report_sha256',sha(rp2))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--id',required=True);p.add_argument('--prior-report',required=True);p.add_argument('--prior-sha',required=True);main(p.parse_args())
