#!/usr/bin/env python3
"""Independent literal-source and ONE coefficient checker; no author imports."""
import argparse, collections, gc, hashlib, json, math, pathlib, resource, signal, time
from fractions import Fraction as F

ROOT = pathlib.Path(__file__).resolve().parents[1]
BASE = 'S61_NONSQUARE_PREPOWER_M12_AXIS_'
GATE = ROOT/'reports'/f'{BASE}WORKER_GATE_20261105_v63'
PRIMARY = ROOT/'reports'/f'{BASE}ENTRY_20261105_v63'
SOURCE_ID = f'{BASE}ONE_COEFFICIENT_SOURCE_20261105_v63_REVIEW'
VALUE_ID = f'{BASE}ONE_COEFFICIENT_20261105_v63_REVIEW'
EXPECTED = F(460045760701204897047977860413366725648080149283864576,
             1756614506664276123046875)
PINS = {
 'reports/S61_NONSQUARE_PREPOWER_M12_AXIS_ENTRY_20261105_v63/report.json': '1cd0ddbac43a9657d22b951befef4d5834e823e6cc01880a531482a7ad164602',
 'reports/S61_NONSQUARE_PREPOWER_M12_AXIS_WORKER_GATE_20261105_v63/report.json': '68e5d8a0cde431f8faf258533f1ffc9dcefedde2a19bd26c9e0870ed3ada379c',
 'reports/S61_NONSQUARE_PREPOWER_M12_AXIS_WORKER_GATE_20261105_v63/recipe.json': '192fdcd9bee80fe26bbe368494a67140e314aa5cf23e3dcbbceacd2b3972a398',
 'reports/S61_NONSQUARE_PREPOWER_M12_AXIS_WORKER_GATE_20261105_v63/literal-atoms67.jsonl': 'ce8d552d7e7dbf20b43e46ee784da94786bec81e082a3516f99583f56a7191d8',
 'reports/S61_NONSQUARE_PREPOWER_M12_AXIS_WORKER_GATE_20261105_v63/old-axis579.json': '6f0b77960137c303cbafffc47a4d55f0403c4e7ce19ff4408026c22f421f4157',
 'reports/S61_NONSQUARE_CUBIC_Q2_L_UNIT_COMPLETION_20261104_v63/certificate.json': 'b1521bc07df9b31c8d23b2bdc519108e80aa859b9582a0213d72d14661db2516',
 'reports/S61_NONSQUARE_CUBIC_ZW25_UNITS_20261104_v63/checkpoint.json': '9723b4430d1bf8cf59cb98b5392bba13af83ce1413ee2f2b17205d25af3d58c3',
 'reports/S61_NONSQUARE_CUBIC_GENUINE14_H5_PREFLIGHT_20261104_v63/certificate.json': 'f7b495ef7b47f74c6c5e3067209f8c2fd4d9294a30a4371048948137ea42a2c7',
 'reports/S61_NONSQUARE_PURE_KERNEL12_SOURCE_PREFLIGHT_20261104_v63/prePOWER-h9-source.json': 'aa7969ea2c606d375d4c9a078752998ab47afb1759ee5ddd8586b7a3a473b612',
 'reports/S61_NONSQUARE_PREPOWER_K92_WHOLE_AXIS_SOURCE_V2_20261105_v63/whole-axis-literal-source.json': '4bcad246e83aea056c4e84622b16299e016bf12ec18a4da94d97c6782d6eb40a',
 'reports/S61_NONSQUARE_PREPOWER_DIRECT_M12_FULL15_SOURCE_20261105_v63_REVIEW/report.json': 'f946ca6312c0b5167aa5cc257434d9338727a0ff3a6142496bc2ffa55d3890c8',
}

def sha(p):
    h=hashlib.sha256()
    with pathlib.Path(p).open('rb') as f:
        for b in iter(lambda:f.read(65536), b''): h.update(b)
    return h.hexdigest()
def load(p): return json.loads(pathlib.Path(p).read_text())
def save(p,x): pathlib.Path(p).write_text(json.dumps(x,indent=2)+'\n')
def numeric(x):
    if isinstance(x,dict): return sum(numeric(k)+numeric(v) for k,v in x.items())
    if isinstance(x,(list,tuple)): return sum(map(numeric,x))
    return int(isinstance(x,(int,F)) and not isinstance(x,bool) or isinstance(x,str) and x.lstrip('-').replace('/','').isdigit())
def field(text,key):
    """Decode exactly one named field; do not decode unrelated matrices."""
    pos=text.index(json.dumps(key)); pos=text.index(':',pos)+1
    return json.JSONDecoder().raw_decode(text,pos+len(text[pos:])-len(text[pos:].lstrip()))[0]
def array_spans(text,start):
    """Independent bracket/string scanner, yielding immediate element spans."""
    assert text[start]=='['
    i=start+1
    while True:
        while text[i].isspace() or text[i]==',': i+=1
        if text[i]==']': return
        lo=i; depth=0; quoted=False; escaped=False
        while i<len(text):
            c=text[i]
            if quoted:
                if escaped: escaped=False
                elif c=='\\': escaped=True
                elif c=='"': quoted=False
            elif c=='"': quoted=True
            elif c in '[{': depth+=1
            elif c in ']}':
                if depth==0: break
                depth-=1
                if depth==0: i+=1; break
            elif c==',' and depth==0: break
            i+=1
        yield lo,i
def norm3(p):
    d={}
    for m,c in p:
        m=tuple(m); assert len(m)==3 and m not in d
        if F(c): d[m]=F(c)
    return d
def literal_equal(a,b): return norm3(a[0])==norm3(b[0]) and norm3(a[1])==norm3(b[1])

def bind_source(out):
    pins={p:sha(ROOT/p) for p in PINS}
    assert pins==PINS
    gate=load(GATE/'report.json'); recipe=load(GATE/'recipe.json'); old=load(GATE/'old-axis579.json')
    assert gate['all_passed']
    qtext=(ROOT/'reports/S61_NONSQUARE_CUBIC_Q2_L_UNIT_COMPLETION_20261104_v63/certificate.json').read_text()
    dag=field(qtext,'scalar_DAG'); roots=field(qtext,'remaining15x22_scalar_roots')
    assert roots==recipe['current15_roots']
    checkpoint=(ROOT/dag['checkpoint_reference']).read_text()
    start=checkpoint.index('[',checkpoint.index('"matrix25x32"'))
    rows=list(array_spans(checkpoint,start)); assert len(rows)==25
    spans={}
    for r,(a,b) in enumerate(rows):
        cells=list(array_spans(checkpoint,a)); assert len(cells)==32
        for c,s in enumerate(cells): spans[r,c]=s
    hfile=ROOT/'reports/S61_NONSQUARE_PURE_KERNEL12_SOURCE_PREFLIGHT_20261104_v63/prePOWER-h9-source.json'
    h=load(hfile); h9=h['source_h9'][h['original_target_labels'].index(17)]
    heavytext=(ROOT/'reports/S61_NONSQUARE_CUBIC_GENUINE14_H5_PREFLIGHT_20261104_v63/certificate.json').read_text()
    hs=heavytext.index('[',heavytext.index('"raw_compact14x5"'))
    hr=list(array_spans(heavytext,hs)); hc=list(array_spans(heavytext,hr[13][0]))
    assert len(hr)==14 and len(hc)==5
    source_checks=[]; atom_terms=0; max_literal_leaves=0
    bp={(1,0,1):F(50),(0,2,1):F(18),(0,1,1):F(-9),(0,0,1):F(-2),(0,1,0):F(60),(0,0,0):F(-40)}
    with (GATE/'literal-atoms67.jsonl').open() as f:
        for line in f:
            a=json.loads(line); ref=a['source_reference']; lit=a['literal_num_den']; meta=a['expected']
            max_literal_leaves=max(max_literal_leaves,numeric(a)); atom_terms+=len(lit[0])+len(lit[1])
            if 'constant' in ref: expected=[[[[0,0,0],'1']],[[[0,0,0],'1']]]
            elif 'text_span' in ref:
                s=spans[ref['row'],ref['column']]; assert list(s)==ref['text_span']
                expected=json.loads(checkpoint[s[0]:s[1]])
            elif 'json_pointer' in ref: expected=dag['new_unit_atoms_num_den'][int(ref['json_pointer'].split('/')[-1])]
            elif 'compact' in ref:
                s=hc[ref['heavy_column']]
                expected=[json.loads(heavytext[s[0]:s[1]]),[[[0,0,0],'16119/8000000']]]
                assert ref['compact']==13 and ref['known_exact_rho13_over_L']=='8000000/16119'
            else:
                assert ref['target']==17 and ref['path']==str(hfile.relative_to(ROOT))
                expected=[h9,[[[0,0,0],'1']]]
            assert literal_equal(lit,expected),a['atom']
            num,den=map(norm3,lit); aa=meta['denominator_a_power']; bb=meta['denominator_B_power']; unit=F(meta['denominator_QQ_unit'])
            assert aa>=0 and bb in (0,1) and unit
            assert all(m[1]>=0 and m[2]>=0 for m in num)
            if bb: assert aa==0 and den=={m:unit*c for m,c in bp.items()}
            else: assert den=={(aa,0,0):unit}
            source_checks.append({'atom':a['atom'],'literal_original_source_equal':True,'legal_a_B_localized_face':True,'literal_terms':len(num)+len(den)})
    assert len(source_checks)==67
    seen=set()
    for n in recipe['nodes']:
        i=n['node']; d=n['definition']; assert dag['nodes'][i]==d
        if n['atom'] is None and d[0] in ('a','s','m'): assert d[1] in seen and d[2] in seen
        seen.add(i)
    assert len(seen)==214 and recipe['current22_original_labels']==list(field(qtext,'new22_target_original32_columns'))
    oldsource=load(ROOT/old['old_reference']); assert sha(ROOT/old['old_reference'])==old['old_sha256']
    assert list(map(int,old['old_row_LCD']))==list(map(int,oldsource['axis_row_denominators']))
    assert len(old['axis_coefficients'])==121
    for x,y in zip(old['axis_coefficients'],oldsource['axis_coefficients']):
        assert (x['row'],x['column'],x['polynomial'])==(y['row'],y['column'],y['axis_integer_polynomial'])
    assert sum(len(x['polynomial']) for x in old['axis_coefficients'])==579
    rs=recipe['rowdual']; cs=recipe['coldual']; shift=-sum(rs)-sum(cs)
    assert shift==5; target=15-shift; assert target==10
    rowD=list(map(int,old['old_row_LCD'])); es=recipe['extra_row_QQ_integer_scales']; D12=int(recipe['new_row_LCD'])
    norms=[0]*12
    for x in old['axis_coefficients']: norms[x['row']]+=sum(abs(v) for _,v in x['polynomial'])*es[x['row']]
    for i,b in enumerate(recipe['source_and_extra_axis_bounds']['addedcol']):
        assert b['lo']+rs[i]>=0 and b['hi']+rs[i]<=9
        norms[i]+=(b['L1']*rowD[i]*es[i]+b['D']-1)//b['D']
    for j,b in enumerate(recipe['source_and_extra_axis_bounds']['bottom']):
        assert b['lo']+cs[j]>=0 and b['hi']+cs[j]<=9
        norms[11]+=(b['L1']*D12+b['D']-1)//b['D']
    b=recipe['source_and_extra_axis_bounds']['corner']; assert b['lo']>=0 and b['hi']<=9
    norms[11]+=(b['L1']*D12+b['D']-1)//b['D']
    H=math.prod(norms); assert H<=int(recipe['all_minor_H']) and H.bit_length()<=1096
    masks=max(math.comb(12,r)+math.comb(12,r+1) for r in range(12))
    dpcoeff=masks*(target+1); numeric_upper=60000
    digits=len(str(H))+1
    dpraw=4096*(40+11*(digits+4)); raw_input=sum((ROOT/p).stat().st_size for p in PINS)
    # The large checkpoint is only a raw string plus <=one decoded literal; never the full rational matrix.
    raw_bound=raw_input+dpraw+2000000
    assert raw_bound<32*1024*1024 and dpcoeff==18876
    source_numeric=numeric(recipe)+numeric(old)+numeric(dag)+numeric(h)+max_literal_leaves+3000
    assert source_numeric<60000
    proof={
      'all67_original_literal_bindings':source_checks,'all214_original_DAG_definitions_equal':True,
      'all_dependencies_topological':True,'old579_exact_frozen_input_reuse':True,
      'old579_source_arithmetic_not_repeated':True,'face_subring':'QQ[a,e,d,1/a,1/B] -> QQ[a,1/a], e=d=0, B=-40',
      'axis_is_not_a_FULLU_parameter_point':True,'new_integer_target_degree':target,'D12_to_integer_C_exponent_shift':shift,
      'det_normalization':'detM12_axis = a^5 detC12 / (new_LCD * product(old_LCD_i * extra_i))',
      'row_L1_upper':list(map(str,norms)),'integer_DP_absolute_upper':str(H),
      'integer_DP_absolute_bits':H.bit_length(),'source_Fraction_bit_upper':238,
      'adjacent_mask_layers_upper':masks,'DP_coefficient_slots_upper':dpcoeff,
      'DP_plus_input_metadata_numeric_upper':numeric_upper,'source_numeric_upper':source_numeric,
      'source_phase_parser_holds_one_checkpoint_literal_only':True,'RAW_input_plus_output_upper':raw_bound,
      'DP_output_upper':dpraw,'hard_numeric_cap':150000,'hard_RAW_cap':32*1024*1024,
      'hard_bit_cap':2048,'hard_seconds':60,'hard_memory_bytes':2000000000,
      'single_coefficient_only':True,'no_values_evaluated_in_this_SOURCE_phase':True,
      'DP_identity':'D_(r+1,J union {j}) += (-1)^#{i in J:i>j} D_(r,J) C_(r,j) mod a^11',
      'truncation_justified_by_all_integer_entry_exponents_nonnegative':True,
      'coefficient_bound_proof':'Each partial signed sum is bounded by the product of prefix row L1 norms, itself <= H; all row norms >=1.',
      'source_frame':'M11 row1+2row2 also on target17; G134 H3 = (8000000/16119)*compact13[1,3,2]; no kappa3',
      'scope':'One coefficient proves detD12 nonzero. On oldOmega intersect D(D12), accepted DIRECT identity gives full57; no global FULLU assertion.'}
    save(out/'source-and-budget-proof.json',proof)
    report={'id':SOURCE_ID,'all_passed':True,'frozen':True,'input_sha256':pins,
            'code_sha256':sha(__file__),'evidence_sha256':{'source-and-budget-proof.json':sha(out/'source-and-budget-proof.json')},
            'no_signed_source_or_determinant_evaluation':True}
    save(out/'report.json',report)
    return recipe,old,proof

def run_values(out):
    t=time.monotonic(); resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000))
    resource.setrlimit(resource.RLIMIT_CPU,(58,58)); signal.alarm(60)
    src=ROOT/'reports'/SOURCE_ID
    sr=load(src/'report.json'); assert sr['all_passed'] and sr['code_sha256']==sha(__file__)
    for p,h in sr['input_sha256'].items(): assert sha(ROOT/p)==h
    proof=load(src/'source-and-budget-proof.json'); recipe=load(GATE/'recipe.json'); old=load(GATE/'old-axis579.json')
    atoms={}; values={}; bits=0; peak=0
    def boundbits(q):
        nonlocal bits
        bits=max(bits,abs(q.numerator).bit_length(),q.denominator.bit_length()); assert bits<=238
    def guard(scratch=0):
        nonlocal peak
        # All scalar coefficient, exponent, dictionary-key, recipe/input and audit metadata slots.
        n=2*(sum(map(len,atoms.values()))+sum(map(len,values.values())))+scratch+22000
        peak=max(peak,n); assert n<=60000
        assert time.monotonic()-t<57
    def simplify(d): return {k:v for k,v in d.items() if v}
    def plus(a,b,sign=1):
        c=dict(a)
        for k,v in b.items(): c[k]=c.get(k,F(0))+sign*v; boundbits(c[k])
        return simplify(c)
    def times(a,b):
        c={}
        for i,x in a.items():
            for j,y in b.items():
                z=x*y; boundbits(z); c[i+j]=c.get(i+j,F(0))+z; boundbits(c[i+j])
        return simplify(c)
    def scale(a,q=F(1),shift=0):
        c={}
        for i,x in a.items():
            z=x*q; boundbits(z)
            if z:c[i+shift]=z
        return c
    def check(a,b):
        if not a:return
        assert min(a)>=b['lo'] and max(a)<=b['hi'] and len(a)<=b['terms']
        assert all((x*b['D']).denominator==1 for x in a.values())
        assert sum(abs(x*b['D']) for x in a.values())<=b['L1']
    def serial(a):return [[i,str(x)] for i,x in sorted(a.items())]
    for line in (GATE/'literal-atoms67.jsonl').open():
        a=json.loads(line); meta=a['expected']; aa=meta['denominator_a_power']; bb=meta['denominator_B_power']
        den=F(meta['denominator_QQ_unit'])*(-40)**bb
        p={}
        for m,c in a['literal_num_den'][0]:
            if not m[1] and not m[2]:
                q=F(c)/den; boundbits(q)
                if q:p[m[0]-aa]=q
        check(p,meta['bound']['axis']); atoms[a['atom']]=p; guard(numeric(a)+300)
    del a,p,meta
    with (PRIMARY/'source214-axis-witness.jsonl').open() as f, (out/'all214-fixed-source-checks.jsonl').open('w') as w:
        for n in recipe['nodes']:
            d=n['definition']; i=n['node']
            if n['atom'] is not None:p=atoms[n['atom']]
            elif d[0]=='z':p={}
            elif d[0]=='one':p={0:F(1)}
            elif d[0]=='m':p=times(values[d[1]],values[d[2]])
            else:p=plus(values[d[1]],values[d[2]],-1 if d[0]=='s' else 1)
            check(p,n['axis_bound']); values[i]=p
            expected=json.loads(next(f)); assert expected['node']==i and expected['definition']==d and expected['axis_polynomial']==serial(p)
            w.write(json.dumps({'node':i,'own_axis_polynomial':serial(p),'original_definition_and_primary_fixed_value_match':True})+'\n'); guard(500)
        assert not f.read().strip()
    labels=recipe['current22_original_labels']; rt=recipe['current15_roots']; c={int(k):F(v) for k,v in recipe['constants'].items()}
    def m(i,j):return values[rt[i][labels.index(j)]]
    col=[atoms[recipe['h9_atom']],scale(plus(m(1,17),scale(m(3,17),F(2))),shift=-1),scale(m(3,17),shift=-1),m(6,17),m(2,17),
         scale(m(11,17),1/c[10]),scale(plus(m(12,17),m(10,17),-1),1/c[11]),scale(m(9,17),1/c[8]),scale(m(7,17),1/c[6]),
         scale(m(13,17),1/c[12],1),scale(m(10,17),1/c[9])]
    bottom=[m(14,j) for j in [22,11,13,25,15,19,14,18]]+[atoms[j] for j in recipe['H3_atoms_in_order']]
    corner=m(14,17)
    extra=load(PRIMARY/'added23-axis-values.json')
    assert extra['column17']==list(map(serial,col)) and extra['G134row']==list(map(serial,bottom)) and extra['corner']==serial(corner)
    D=list(map(int,old['old_row_LCD'])); E=recipe['extra_row_QQ_integer_scales']; Dn=int(recipe['new_row_LCD'])
    C=[[{} for j in range(12)] for i in range(12)]
    for x in old['axis_coefficients']: C[x['row']][x['column']]={k:v*E[x['row']] for k,v in x['polynomial']}
    def integers(a,f,shift):
        z={}
        for k,v in a.items():
            q=v*f; assert q.denominator==1 and 0<=k+shift<=9
            z[k+shift]=q.numerator
        return z
    for i,p in enumerate(col):C[i][11]=integers(p,D[i]*E[i],recipe['rowdual'][i])
    for j,p in enumerate(bottom):C[11][j]=integers(p,Dn,recipe['coldual'][j])
    C[11][11]=integers(corner,Dn,0)
    with (PRIMARY/'native-input.txt').open() as f:
        assert next(f).split()[0]=='M12AXIS001'; assert int(next(f))==int(recipe['all_minor_H'])
        for i in range(12):
            for j in range(12):
                v=list(map(int,next(f).split())); assert v[:2]==[i,j] and v[2]==len(C[i][j])
                assert dict(zip(v[3::2],v[4::2]))==C[i][j]
        assert not f.read().strip()
    for i in range(12):assert sum(sum(abs(x) for x in p.values()) for p in C[i])<=int(proof['row_L1_upper'][i])
    save(out/'own-integer-C12-normalization.json',{'matrix':[[serial(p) for p in row] for row in C],
         'all144_primary_integer_inputs_match':True,'all23_actual_source_adapter_values_match':True,
         'target_degree':proof['new_integer_target_degree'],'LCD_product':str(Dn*math.prod(a*b for a,b in zip(D,E)))})
    denom=Dn*math.prod(a*b for a,b in zip(D,E)); K=proof['new_integer_target_degree']; H=int(proof['integer_DP_absolute_upper'])
    del atoms,values,col,bottom,corner,extra,old,recipe,p,n,d,expected; gc.collect()
    dp={0:[1]+[0]*K}; stage=[]; products=0; dp_peak=0
    with (out/'truncated-subset-DP-states.jsonl').open('w') as w:
        for r in range(12):
            nxt={}
            for mask,pol in dp.items():
                for j in range(12):
                    if mask>>j&1:continue
                    nm=mask|(1<<j); ans=nxt.setdefault(nm,[0]*(K+1)); sign=-1 if (mask>>(j+1)).bit_count()%2 else 1
                    for a,x in enumerate(pol):
                        if not x:continue
                        for b,y in C[r][j].items():
                            if a+b<=K:
                                ans[a+b]+=sign*x*y; products+=1
                                assert abs(ans[a+b])<=H
            states=len(dp)+len(nxt); owned=states*(K+2)+10000
            dp_peak=max(dp_peak,owned); assert owned<=60000 and time.monotonic()-t<57
            digest=hashlib.sha256()
            for mask,pol in sorted(nxt.items()):
                rec=json.dumps({'rows':r+1,'column_mask':mask,'coefficients0through10':list(map(str,pol))},separators=(',',':'))+'\n'
                w.write(rec);digest.update(rec.encode())
            stage.append({'rows':r+1,'masks':len(nxt),'state_stream_SHA256':digest.hexdigest()});dp=nxt
    integer=dp[(1<<12)-1][K]; answer=F(integer,denom); assert answer==EXPECTED and answer
    ev={'single_coefficient_a15':str(answer),'integer_C12_coefficient_a10':str(integer),'LCD_product':str(denom),
        'derived_degree10_from_normalization_not_assumed':True,'all67_original_sources_and214_own_values_checked':True,
        'all144_integer_inputs_and23_added_source_adapter_values_checked':True,'DP_stage_receipts':stage,
        'all_integer_accumulators_bounded_by_source_absolute_row_L1_product':True,
        'source_fraction_observed_bits':bits,'source_numeric_upper_observed':peak,'DP_numeric_upper_observed':dp_peak,
        'integer_products':products,'one_coefficient_only':True,'all29_axis_coefficients_NOT_replayed':True,
        'char0_actual_D12_nonzero':True,'no_numeric_parameter_point_or_prime':True,
        'scope':'On oldOmega intersect D(D12), the accepted DIRECT actual-source identity implies full57. Complements remain open; no global FULLU57 or quartic239 claim.'}
    save(out/'one-coefficient-proof.json',ev)
    report={'id':VALUE_ID,'all_passed':True,'frozen':True,'source_report_sha256':sha(src/'report.json'),
      'primary_report_sha256':PINS[str((PRIMARY/'report.json').relative_to(ROOT))],
      'code_sha256':sha(__file__),'elapsed_seconds':time.monotonic()-t,
      'evidence_sha256':{p.name:sha(p) for p in out.iterdir() if p.is_file()},
      'single_a15_coefficient_independently_verified':True,'full_QQ_source_scope':True,
      'producer_or_native_imports':False,'no_old572_replay':True,'no_full29_coefficient_claim':True}
    assert sum(p.stat().st_size for p in out.iterdir() if p.is_file())<proof['DP_output_upper']+2000000
    save(out/'report.json',report); print(json.dumps({'all_passed':True,'report':str(out/'report.json'),'sha256':sha(out/'report.json'),'coefficient':str(answer),'elapsed_seconds':report['elapsed_seconds']}))

def main():
    p=argparse.ArgumentParser();p.add_argument('phase',choices=['source','value']);args=p.parse_args()
    out=ROOT/'reports'/(SOURCE_ID if args.phase=='source' else VALUE_ID)
    assert not out.exists(),f'Frozen output exists: {out}'
    out.mkdir(parents=True)
    try:
        if args.phase=='source':
            bind_source(out);print(json.dumps({'all_passed':True,'report':str(out/'report.json'),'sha256':sha(out/'report.json')}))
        else:run_values(out)
    except Exception as e:
        save(out/'failure.json',{'type':type(e).__name__,'message':str(e),'no_positive_claim':True});raise
if __name__=='__main__':main()
