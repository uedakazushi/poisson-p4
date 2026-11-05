#!/usr/bin/env python3
"""SOURCE only: two fixed invertible constraint minors and adjugate bounds."""
from pathlib import Path
from fractions import Fraction as Q
from math import gcd,factorial
import gc,hashlib,itertools,json,resource,signal,time
ROOT=Path(__file__).resolve().parents[1]
ID='S61_SPARSE_GENERIC5_CONSTRAINED83_EXACT49_FIXED_CONSTRAINT_MINOR_PREFLIGHT_20261105_v64'
PROFILE='S61_SPARSE_GENERIC5_CONSTRAINED83_EXACT49_CONSTRAINT_PROFILE_20261105_v64'
def sha(p):
    h=hashlib.sha256()
    with p.open('rb')as f:
        for b in iter(lambda:f.read(65536),b''):h.update(b)
    return h.hexdigest()
def read(p):return json.loads(p.read_text())
def enc(v):return(json.dumps(v,separators=(',',':'))+'\n').encode()
def num(v):
    if isinstance(v,bool):return 0
    if isinstance(v,(int,float)):return 1
    if isinstance(v,list):return sum(map(num,v))
    if isinstance(v,dict):return sum(map(num,v.values()))
    return 0
def main():
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000))
    out=ROOT/'reports'/ID;out.mkdir(parents=True,exist_ok=False);(out/'support-only').mkdir()
    pd=ROOT/'reports'/PROFILE;pr=pd/'report.json';assert sha(pr)=='1c97a068769399f85198e5d47943c8f925bfaf7757accee5fb874259ab0c9a62'
    p=read(pr);assert p['all_passed']and p['constraint_ranks_from_old_L18']==[5,5]
    ep=pd/'actual-constraint-blocks-and-frontier.json';assert sha(ep)==p['evidence_sha256'][ep.name];profile=read(ep)
    od=ROOT/'reports/S61_SPARSE_GENERIC5_SMALL_TAIL_FIELD_20261104_v63';orr=od/'report.json';orrj=read(orr);assert orrj['all_passed']and orrj['kernel_complete']
    cd=ROOT/'reports/S61_SPARSE_GENERIC5_CONSTANT_UNIT_PREFIX_20261104_v63';cr=cd/'report.json';cj=read(cr);assert cj['all_passed']
    cp=cd/'constant-prefix.json';fp=od/'field-prefix-checkpoint.json';kp=od/'generic-field-kernel.json'
    for path,rr in [(cp,cj),(fp,orrj),(kp,orrj)]:assert sha(path)==rr['evidence_sha256'][path.name]
    C=read(cp);F=read(fp);assert C['constant_QQ_pivots']==11 and F['completed_pivots']==7
    basis=C['row_order'][:11]+[C['row_order'][11+j]for j in F['row_order'][:7]];assert len(set(basis))==18
    source_peak=num(C)+num(F)+num(profile)+8192;del C,F;gc.collect()
    # Whole kernel has159249 numeric slots. Do not parse it: reuse its accepted
    # all51/free-I2 certificate, read only its small literal leading header.
    with kp.open('rb')as f:head=f.read(240).decode()
    compact=''.join(head.split());assert '"all_passed":true' in compact and '"full_L_rank":18' in compact and '"free_original_columns":[6,12]' in compact
    gate=ROOT/'reports/S61_SPARSE_GENERIC5_PLANE_RANK3_CIRCUIT_20261104_v63_REVIEW/report.json';q=read(gate);assert q['all_passed']and q['rank_L18_plane3_normal42'];del q
    paths={pr,ep,cp,fp,kp,orr,cr,gate,Path(__file__),ROOT/f'research/notes/{ID}.md',ROOT/'code/s61_sparse_generic5_small_tail_field_20261104_v63.py',ROOT/'code/s61_sparse_generic5_small_tail_field_20261104_v63.cpp'}
    base_slots=num(profile)+8192;support_pool_peak=0;future_peak=0;retained=0;plans=[];inputs=[];output_caps=0;all_proofs=[]
    INPUT_COEFFICIENTS=1150;FUTURE_METADATA=25000
    for group in profile['actual_constraint_blocks']:
        name=group['name'];off=134 if name=='alpha'else 185;free=[7,13]if off==134 else[27,33]
        cols=[c for c in group['columns']if c not in free];cells={(r,c):v for r,c,k,h,v in group['entries']};survive=[off+r for r in basis if off+r in group['rows']]
        nonzero=[r for r in survive if any((r,c)in cells for c in group['columns'])];zero=[r for r in survive if r not in nonzero]
        rows=sorted(nonzero);assert rows==[off+j for j in[6,14,16,41,42]]and zero==[off+1,off+26]
        assert len(survive)==7 and len(rows)==5
        matrix=[];B=[];scales=[];supports=[];rowL1=[]
        for r in rows:
            selected=[cells.get((r,c))for c in cols];bs=[cells.get((r,c))for c in free];polys=[v for v in selected+bs if v]
            den=1
            for v in polys:den=den*Q(v['content']).denominator//gcd(den,Q(v['content']).denominator)
            shift=[min(v['shift'][j]+t[j+1]for v in polys for t in v['terms'])for j in range(5)]
            scales.append({'row':r,'positive_integer_denominator':str(den),'parameter_monomial_exponent':[-v for v in shift]})
            ss=[];l1=[]
            for v in selected:
                if not v:ss.append(set());l1.append(0);continue
                c=Q(v['content']);assert den%c.denominator==0
                l1.append(abs(c.numerator)*(den//c.denominator)*sum(abs(int(t[0]))for t in v['terms']))
                ss.append({tuple(v['shift'][j]+t[j+1]-shift[j]for j in range(5))for t in v['terms']});assert all(min(e)>=0 for e in ss[-1])
            supports.append(ss);rowL1.append(l1);matrix.append(selected);B.append(bs)
        inputs.append({'name':name,'rows':rows,'pivot_columns':cols,'free_columns':free,'literal_H5x5':matrix,'literal_B5x2':B,'row_clearings':scales,
                       'definition':'Htilde=diag(den_i*parameter^monomial_i)*H, Btilde uses the identical row scale. All inverses are K5 base functions only.'})
        local_support={};local_records=[]
        # Every output is evaluated sequentially. Store only support certificates
        # on disk; never retain all43847*5 exponent metadata integers in RAM.
        tasks=[('det',list(range(5)),list(range(5)))]+[(f'minor-{r}-{c}',[i for i in range(5)if i!=r],[j for j in range(5)if j!=c])for r in range(5)for c in range(5)]
        for label,rr,cc in tasks:
            dp={0:{(0,0,0,0,0)}};nodes=[];peak=0
            for mask in range(1,1<<len(cc)):
                k=mask.bit_count();union=set();terms=[]
                for j,c in enumerate(cc):
                    if not mask>>j&1:continue
                    a=supports[rr[k-1]][c];b=dp[mask^(1<<j)];prod={tuple(x[v]+y[v]for v in range(5))for x in a for y in b}
                    union.update(prod);pool=sum(map(len,dp.values()));scratch=pool+len(a)+len(b)+len(prod)+len(union)+96
                    support_pool_peak=max(support_pool_peak,5*scratch);peak=max(peak,scratch)
                    source_peak=max(source_peak,base_slots+num(plans)+num(local_records)+5*scratch);assert source_peak<150000
                    future_peak=max(future_peak,INPUT_COEFFICIENTS+FUTURE_METADATA+retained+pool+2*len(prod)+3*len(union)+256)
                    terms.append([c,mask^(1<<j),len(a),len(b),len(prod),len(union)])
                dp[mask]=union;nodes.append([mask,len(union),terms])
            value=dp[(1<<len(cc))-1];bound=len(value);degree=max(map(sum,value),default=0)
            # Integer coefficient L1 bound is the sum over the fixed Leibniz
            # permutations of products of source integer row-cleared L1 norms.
            height=sum(__import__('functools').reduce(lambda a,b:a*b,(rowL1[r][c]for r,c in zip(rr,perm)),1)for perm in itertools.permutations(cc))
            cap=1024+bound*(len(str(height))+5*max(1,len(str(degree)))+48)
            output_caps+=cap;op=out/'support-only'/f'{name}-{label}.json';op.write_bytes(enc(sorted(value)));local_support[label]=op
            local_records.append({'output':label,'rows':rr,'columns':cc,'support_upper':bound,'support_path':str(op.relative_to(ROOT)),'support_sha256':sha(op),'degree_upper':degree,
                                  'integer_coefficient_abs_upper':str(height),'literal_output_RAW_upper':cap,'DP_nodes':nodes,'DP_pool_scratch_bound':peak})
            retained+=bound;del dp,value,union,prod
        # Fixed100 adjugate identities, both orders for each matrix. Signs do
        # not alter supports. No coefficients or determinants are evaluated.
        ident=[]
        for order in ['HAdj','AdjH']:
            for i in range(5):
                for j in range(5):
                    u=set();record=[]
                    for k in range(5):
                        if order=='HAdj':a=supports[i][k];path=local_support[f'minor-{j}-{k}']
                        else:a=supports[k][j];path=local_support[f'minor-{k}-{i}']
                        b={tuple(e)for e in read(path)};prod={tuple(x[v]+y[v]for v in range(5))for x in a for y in b};u.update(prod)
                        source_peak=max(source_peak,base_slots+num(plans)+num(local_records)+5*(len(a)+len(b)+len(prod)+len(u)+96));assert source_peak<150000
                        future_peak=max(future_peak,INPUT_COEFFICIENTS+FUTURE_METADATA+retained+2*len(prod)+3*len(u)+256)
                        record.append([k,len(a),len(b),len(prod),len(u)]);del b,prod
                    ident.append([order,i,j,record])
        plans.append({'name':name,'nonzero_det_from_source_basis_quotient':True,'accepted_old_row_basis18':basis,'surviving_old_rows7':survive,'literal_zero_survivors2':zero,
                      'selected_rows5':rows,'selected_columns5':cols,'free_columns2':free,'outputs':local_records,'fixed_adjugate_identity_support_plans':ident})
        all_proofs.append({'name':name,'det_support_upper':local_records[0]['support_upper'],'all26_output_support_sum':sum(z['support_upper']for z in local_records)})
        del local_records,supports,matrix,B,cells
    assert retained==43847 and future_peak<150000
    metadata_slots=num(plans)+num(inputs)+2048;assert metadata_slots<=FUTURE_METADATA
    inp=out/'two-actual-fixed-constraint-minors.json';inp.write_bytes(enc({'variable_order':['H','L','s0','q0','t0'],'original_free_coordinates':[6,12],'blocks':inputs,'no_evaluated_determinant_or_adjugate':True}))
    plan=out/'support-and-worker-plan.json';plan.write_bytes(enc({'source_only':True,'plans':plans,'all_output_support_sum':retained,'future_fixed_metadata_upper':FUTURE_METADATA,'future_coefficient_and_scratch_upper':future_peak,'actual_metadata_slots':metadata_slots,
        'source_owned_peak_upper':source_peak,'support_metadata_peak_upper':support_pool_peak,'all26_outputs_per_block':all_proofs,'no_coefficient_products_computed':True}))
    base=sum(x.stat().st_size for x in paths);newbytes=sum(x.stat().st_size for x in out.rglob('*')if x.is_file());rawupper=base+newbytes+output_caps+3*1024*1024
    assert rawupper<32*1024*1024
    proof={'all_passed':True,'source_only':True,'coefficient_worker_authorized':False,'future_hard_seconds':60,'future_memory_bytes':2000000000,'future_owned_including_metadata_upper':future_peak,
           'future_RAW_including_actual_source_support_receipts_all52_outputs_and3MiB_other_upper':rawupper,'source_RAW_current_input_output_upper':base+newbytes+1024*1024,
           'literal52_coefficient_outputs_RAW_upper':output_caps,'runtime_obligation':'A future separately authorized worker must enforce the same60s and save every complete minor; no automatic enlarged retry.',
           'unit_scope':'Both fixed determinants are nonzero K5 base Laurent functions by the accepted18-row basis modulo13 pure unit rows and surviving free6,12 projection. They are field units, not units at all complex parameter values. No normal unknown occurs.',
           'localization_and_complement':'The future explicit chart requires both cleared determinant numerators nonzero and every previous parameter-monomial denominator nonzero. The complementary base divisors stay outside that chart; no uniform coverage or exhaustion claim.',
           'inverse_instruction':'x_pivot=Adj(Htilde)*(scaled_rhs-Btilde*x_free)/det(Htilde). All38 zero-RHS constraint rows follow from rank5 and the same kernel. Do not expand the111-row top block without a new whole-support gate.'}
    bp=out/'whole-budget-and-scope.json';bp.write_bytes(enc(proof))
    report={'id':ID,'all_passed':True,'source_only':True,'source_generation':49,'fixed_minors_nonzero_by_accepted_source_linear_algebra':True,'new_determinant_coefficients_computed':False,
            'fixed_original_L_rows':[6,14,16,41,42],'source_owned_peak':source_peak,'future_owned_upper':future_peak,'future_RAW_upper':rawupper,'output_support_sums':[17738,26109],
            'input_sha256':{str(x.relative_to(ROOT)):sha(x)for x in paths},'evidence_sha256':{str(x.relative_to(out)):sha(x)for x in out.rglob('*')if x.is_file()},'seconds':time.monotonic()-start}
    rp=out/'report.json';rp.write_bytes(enc(report));assert base+sum(x.stat().st_size for x in out.rglob('*')if x.is_file())<rawupper
    signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in('input_sha256','evidence_sha256')}));print('report_sha256',sha(rp))
if __name__=='__main__':main()
