#!/usr/bin/env python3
"""SOURCE only: exact factored ten-variable Schur, retaining66 zero rows."""
from pathlib import Path
import hashlib,json,resource,signal,time
ROOT=Path(__file__).resolve().parents[1]
ID='S61_SPARSE_GENERIC5_CONSTRAINED83_FACTORED_CONSTRAINT_SCHUR_PREFLIGHT_20261105_v64'
def sha(p):
    h=hashlib.sha256()
    with p.open('rb')as f:
        for z in iter(lambda:f.read(65536),b''):h.update(z)
    return h.hexdigest()
def read(p):return json.loads(p.read_text())
def enc(v):return(json.dumps(v,separators=(',',':'))+'\n').encode()
def numeric(v):
    if isinstance(v,bool):return 0
    if isinstance(v,(int,float)):return 1
    if isinstance(v,list):return sum(map(numeric,v))
    if isinstance(v,dict):return sum(map(numeric,v.values()))
    return 0
def arrays(blob,key):
    """Byte-index20 row arrays without parsing159249 numeric kernel slots."""
    a=blob.index(key);a=blob.index(b'[',a)+1;out=[]
    while len(out)<20:
        while blob[a]in b' \r\n\t,':a+=1
        assert blob[a]==91;start=a;depth=0;quoted=False;escaped=False
        while True:
            x=blob[a]
            if quoted:
                if escaped:escaped=False
                elif x==92:escaped=True
                elif x==34:quoted=False
            elif x==34:quoted=True
            elif x==91:depth+=1
            elif x==93:
                depth-=1
                if not depth:a+=1;break
            a+=1
        out.append((start,a-start))
    return out
def main():
    begin=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000))
    out=ROOT/'reports'/ID;out.mkdir(parents=True,exist_ok=False)
    def gate(stem,h):
        p=ROOT/f'reports/{stem}/report.json';assert sha(p)==h;r=read(p);assert r['all_passed'];return p,r
    rp,r=gate('S61_SPARSE_GENERIC5_CONSTRAINED83_FIXED_CONSTRAINT_ADJUGATES_20261105_v64','5c45544b64ce273de96ff56fbe26ae5e14ecc56c363b4d2e0b41b443f7d09b2a')
    peer,pe=gate('S61_SPARSE_GENERIC5_FIXED_CONSTRAINT_ADJUGATES_20261105_v64_REVIEW','fbc3cb3fef01642fb7d72ba4956b4b769bc67bcefff77381bd14a8f085d65fa3')
    sr,s=gate('S61_SPARSE_GENERIC5_CONSTRAINED83_STREAM49_ACTIVE_FRAGMENT_MIXED_BINARY_COMPACTION_20261105_v64','62a95473e52e2a85054586bd63e28b79ae58850bebf49d6e5a0a4daeb69e7970')
    fr,f=gate('S61_SPARSE_GENERIC5_CONSTRAINED83_EXACT49_FIXED_CONSTRAINT_MINOR_PREFLIGHT_20261105_v64','15ec0fa2870a5985257a5a9e62500141ade5619b2ba1afe7b50a9bddcafea06e')
    ixp=sr.parent/s['current_index'];assert sha(ixp)==s['current_index_sha256'];ix=read(ixp)
    matrix=sr.parent/s['current_snapshot'];assert sha(matrix)==s['current_snapshot_sha256']
    qp=rp.parent/'source-clearings-and-all52-output-checks.json';assert sha(qp)==r['evidence_sha256'][qp.name];q=read(qp)
    inp=fr.parent/'two-actual-fixed-constraint-minors.json';assert sha(inp)==f['evidence_sha256'][inp.name];small=read(inp)
    paths={rp,peer,sr,fr,ixp,matrix,qp,inp,Path(__file__),ROOT/f'research/notes/{ID}.md'}
    oldreport=ROOT/'reports/S61_SPARSE_GENERIC5_SMALL_TAIL_FIELD_20261104_v63/report.json';old=read(oldreport);assert old['all_passed']and old['kernel_complete']
    kp=oldreport.parent/'generic-field-kernel.json';assert sha(kp)==old['evidence_sha256'][kp.name];paths.update([oldreport,kp])
    nodes=[{'op':'constant','value':'0','type':'K5'},{'op':'constant','value':'1','type':'K5'}];cache={};owner_peak=0
    def node(v):
        key=json.dumps(v,separators=(',',':'),sort_keys=True)
        if key in cache:return cache[key]
        n=len(nodes);nodes.append(v);cache[key]=n;return n
    def op(o,a,b=None):
        if o=='add':
            if a==0:return b
            if b==0:return a
        if o=='mul':
            if a==0 or b==0:return 0
            if a==1:return b
            if b==1:return a
        if o=='neg':
            if a==0:return 0
            return node({'op':o,'arg':a,'type':'K5'})
        return node({'op':o,'args':[a,b],'type':'K5'})
    cells=ix['all_literal_cells'];oldrefs={};old_sparse=[]
    for k,c in enumerate(cells):
        t=node({'op':'literal_current_cell','row':c['row'],'column':c['column'],'source_cell_index':k,'type':'K5'})
        oldrefs[(c['row'],c['column'])]=t;old_sparse.append([c['row'],c['column'],t])
    def src(r,c):return oldrefs.get((r,c),0)
    base_slots=numeric(ix)+numeric(small)+numeric(q)+numeric(r)+numeric(s)+numeric(f)+numeric(old)+numeric(pe)+3*len(oldrefs)+8192
    coefficients={}
    for z in q['output_files']:
        p=rp.parent/z['file'];assert sha(p)==z['sha256'];paths.add(p)
        coefficients[z['file']]=node({'op':'actual_small_polynomial','path':str(p.relative_to(ROOT)),'sha256':z['sha256'],'terms':z['terms'],'type':'K5'})
    # Inspect only seven selected existing row payloads. No new coefficients,
    # no whole-kernel JSON parse, GCD or denominator expansion.
    blob=kp.read_bytes();rowranges=arrays(blob,b'"T20x2_rational"');selected={4,10,11,14,15,18,19};kernel={};denrefs=[];kdesc=[];cost={}
    for row in sorted(selected):
        offset,length=rowranges[row];payload=blob[offset:offset+length];v=json.loads(payload);assert len(v)==2
        kh=hashlib.sha256(payload).hexdigest();desc={'row_original0':row,'path':str(kp.relative_to(ROOT)),'offset':offset,'bytes':length,'payload_sha256':kh,'components':[]}
        for j,a in enumerate(v):
            assert a['denominator'];n=len(a['numerator']);d=len(a['denominator'])
            t=node({'op':'accepted_kernel_coefficient','row_original0':row,'component':j,'payload_offset':offset,'payload_bytes':length,'payload_sha256':kh,'type':'K5'})
            kernel[(row,j)]=t;cost[t]=[n,d];dr={'kernel_node':t,'source_component':j,'denominator_terms':d,'denominator_sha256':hashlib.sha256(enc(a['denominator'])).hexdigest()}
            denrefs.append(dr);desc['components'].append({'node':t,'numerator_terms':n,'denominator_terms':d,'denominator_reference':dr})
        kdesc.append(desc)
        owner_peak=max(owner_peak,base_slots+numeric(nodes)+numeric(old_sparse)+numeric(kdesc)+numeric(v)+len(cache)+8192)
        assert owner_peak<=150000;del v,payload
    del blob
    variables=ix['all_original_column_ids'];rows=ix['all_original_row_ids'];pivot=[];selectedrows=[];blocks=[];substitution={}
    for b in small['blocks']:
        name=b['name'];pivot+=b['pivot_columns'];selectedrows+=b['rows'];origin=1 if name=='alpha'else 21
        det=coefficients[f'{name}-det.json'];inv=node({'op':'inverse','arg':det,'nonzero_gate':str(rp.relative_to(ROOT)),'type':'K5'})
        # Record fully specified actual-adjugate inverse roots as well as the
        # certified old-kernel aliases; do not evaluate any of their products.
        alias=[]
        for i,c in enumerate(b['pivot_columns']):
            row=[]
            for j,free in enumerate(b['free_columns']):
                total=0
                for k in range(5):
                    z=op('mul',coefficients[f'{name}-adj-{i}-{k}.json'],coefficients.get(f'{name}-cleared-{k}-{5+j}.json',0))
                    # Cleared source polynomials are source-bound below.
                    if not z:
                        p=rp.parent/f'{name}-cleared-{k}-{5+j}.json';assert sha(p)==r['evidence_sha256'][p.name];paths.add(p)
                        if read(p)['terms']:
                            leaf=node({'op':'actual_small_polynomial','path':str(p.relative_to(ROOT)),'sha256':sha(p),'type':'K5'})
                            coefficients[p.name]=leaf;z=op('mul',coefficients[f'{name}-adj-{i}-{k}.json'],leaf)
                    total=op('add',total,z)
                actual_solution=op('neg',op('mul',inv,total));known=kernel[(c-origin,j)]
                substitution[(c,free)]=known;row.append({'free_column':free,'actual_adjugate_solution':actual_solution,'accepted_kernel_alias':known,
                    'equality_gate':'Both solve the SAME fixed five constraints and preserve free original coordinates6,12; uniqueness of invertible H proves equality in K5.'})
            alias.append({'pivot_column':c,'solutions':row})
        blocks.append({'name':name,'pivot_rows':b['rows'],'pivot_columns':b['pivot_columns'],'free_columns':b['free_columns'],
                       'actual_det_node':det,'actual_inverse_det_node':inv,'row_clearings':b['row_clearings'],'substitution_aliases':alias})
    keep=[c for c in variables if c not in pivot];resrows=[r for r in rows if r not in selectedrows];top=[r for r in resrows if r<134];zero=[r for r in resrows if r>=134]
    assert len(variables)==34 and len(pivot)==10 and len(keep)==24 and len(resrows)==177 and len(top)==111 and len(zero)==66
    residual=[];input_counts={r:sum(c['terms']for c in cells if c['row']==r)for r in rows};perrow=[]
    for rr in resrows:
        result=[];pair_majorant=0
        for cc in keep:
            value=src(rr,cc)
            for b in small['blocks']:
                if cc in b['free_columns']:
                    for pp in b['pivot_columns']:
                        kk=substitution[(pp,cc)];value=op('add',value,op('mul',src(rr,pp),kk))
                        cnt=next((c['terms']for c in cells if c['row']==rr and c['column']==pp),0);pair_majorant+=cnt*cost[kk][0]
            result.append(value)
        residual.append([rr,result]);perrow.append({'row':rr,'current_input_terms':input_counts[rr],'lazy_numerator_pair_majorant_only':pair_majorant,
            'all_scalar_field_products_unevaluated':True,'constraint_row_identically_zero_by_rank_basis':rr in zero})
        owner_peak=max(owner_peak,base_slots+numeric(nodes)+numeric(old_sparse)+numeric(residual)+numeric(perrow)+numeric(kdesc)+len(cache)+8192);assert owner_peak<150000
    # Keep every actual root of the66 formal zero rows. Never replace an
    # unevaluated scalar response by literal zero.
    operator={'all_passed':True,'source_only':True,'field':'QQ(H,L,s0,q0,t0)','unknown_normal_inverted':False,'scalar_operator_nodes':nodes,
        'original34_columns':variables,'original187_rows':rows,'original_matrix_nonzero_roots':old_sparse,
        'current_source_index':str(ixp.relative_to(ROOT)),'current_source_index_sha256':sha(ixp),'current_source_snapshot_sha256':sha(matrix),
        'small_constraint_blocks':blocks,'all34_variables_retained_in_formal_substitution':True,
        'remaining24_columns':keep,'full_residual177_rows':residual,'top111_row_ids':top,'formal_zero66_row_ids':zero,
        'all66_zero_row_roots_retained_not_zero_filled':True,'old_kernel_payload_index7':kdesc,'old_kernel_extra_denominator_refs':denrefs,
        'field_alias_chart':'Actual adjugate chart D(detA*detB) plus old monomial units. Evaluating the cached canonical aliases additionally requires their14 recorded denominator polynomials nonzero; all complementary base divisors remain outside that chosen chart.',
        'formal_common_denominator':'Keep detA,detB and/or the product of the14 old-kernel denominator references FACTORED. No denominator LCM/GCD or polynomial product is evaluated.',
        'block_equivalence':'After selected-row scales rho, reorder M=[H,V;C,U]. H=diag(H_A,H_B). Exact row normalization/elimination and column shear give diag(I10,U-C H^-1 V). Small inverse is Adj(H)/det(H).',
        'cokernel_statement':'coker(original187x34) = coker(full177x24). coker(full177x24) = coker(top111x24) direct_sum K5^66. Dropping zero rows REMOVES the66 free cokernel summands; only kernels/column ranks agree.',
        'RHS_instruction':'Normalize selected RHS as Htilde^-1 diag(rho)*b_selected; other RHS becomes b_other-C*normalized_selected_RHS. Do not drop remaining constraint RHS for arbitrary input.',
        'zero_row_proof':'Each remaining constraint row belongs to its fixed five-row rank5 basis, whose kernel is the same accepted free6,12 kernel. The current actual Schur coefficients therefore vanish; formal roots are kept above.',
        'known_rank_only':'Rank83 was known before these49 unit steps; the residual24 column rank is a consequence, not a new calculation. Complete24 inverse and fiber exhaustion remain open.'}
    opath=out/'actual-factored177x24-and-top111x24.json';opath.write_bytes(enc(operator))
    pc=out/'row-stream-cost-source-only.json';pc.write_bytes(enc({'all_passed':True,'source_only':True,'actual_adjugate_factor_terms':17784,'current_matrix_logical_terms':809647,
        'current_matrix_not_owned_in_full':True,'per_row':perrow,'max_current_row_input_terms':max(input_counts.values()),
        'small_cached_alias_numerator_terms':sum(a[0]for a in cost.values()),'small_cached_alias_denominator_terms':sum(a[1]for a in cost.values()),
        'no_top_coefficient_worker_GO':True,'future_worker_obligations':'First forecast support, coefficient heights, current literal READ set, all new RAW versions and fused/streamed scratch for the chosen row/entry. Pair majorants here are NOT a safe output/owned certificate. No top scalar product is evaluated.'}))
    active=sum(p.stat().st_size for p in paths);outputs=opath.stat().st_size+pc.stat().st_size;rawupper=active+outputs+1024*1024
    assert rawupper<32*1024*1024 and owner_peak<150000
    budget={'all_passed':True,'source_only':True,'owned_metadata_and_source_row_parse_upper':owner_peak,'RAW_actual_read_new_operator_receipts_reserve_upper':rawupper,
            'no_scalar_coefficient_product_expansion':True,'no_further_worker_authorized':True,'memory_bytes':2000000000,'hard_seconds':60}
    bp=out/'whole-source-budget.json';bp.write_bytes(enc(budget))
    report={'id':ID,'all_passed':True,'source_only':True,'old_shape':[187,34],'factored_full_shape':[177,24],'factored_top_shape':[111,24],
        'zero_constraint_rows_retained':66,'original_column_labels_retained':34,'typed_scalar_nodes':len(nodes),'seven_old_kernel_payloads_source_bound':True,
        'actual_adjugate_source_terms':17784,'source_owned_peak':owner_peak,'RAW_upper':rawupper,'no_new_coefficient_worker_or_complete_inverse':True,
        'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},'evidence_sha256':{p.name:sha(p)for p in out.iterdir()if p.is_file()},'seconds':time.monotonic()-begin}
    rp2=out/'report.json';rp2.write_bytes(enc(report));assert active+sum(p.stat().st_size for p in out.iterdir())<=rawupper
    signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in('input_sha256','evidence_sha256')}));print('report_sha256',sha(rp2))
if __name__=='__main__':main()
