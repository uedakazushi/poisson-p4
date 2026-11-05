#!/usr/bin/env python3
"""Known-action linear G2 profile; no physical/source coefficient query."""
import argparse,hashlib,itertools,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
from ambient_B18_ell2_raw33_quotient_20261105_v63 import add,scale,mv,srec
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_ELL2_C2_G2_LINEAR_SOURCE_PROFILE_20261105_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def parse(v):return {i:F(c) for i,c in v}
def guard(v):
    assert all(max(abs(x.numerator).bit_length(),x.denominator.bit_length())<=384 for x in v.values());return v
def main(out):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024))
    out=out.resolve();out.mkdir(parents=True,exist_ok=False);hs={}
    def read(p):hs[p]=sha(ROOT/p);return json.loads((ROOT/p).read_text())
    t=read('reports/AMBIENT_B18_ELL2_HONEST_T2_ACTIVE_UNIT_TRANSPORT_20261105_v63/actual-active-unit-bases-and26-frame-transport.json')
    pref=read('reports/AMBIENT_B18_ELL2_C2_SECOND_BETTI_ACTIVE_PREFLIGHT_20261105_v63/active-five-d3-source-split-and-budget.json')
    peer=read('reports/AMBIENT_B18_ELL2_HONEST_T2_ACTIVE_UNIT_TRANSPORT_20261105_v63_REVIEW/report.json');assert peer['all_passed']
    active=pref['active_variables'];assert active==[0,1,4,5,7]
    selected=t['active56_source_columns'];free=t['active_F2_free42'];Pi=list(map(parse,t['active_Pi17_columns']));D2=list(map(parse,pref['active_d2_columns']))
    fixed=list(map(parse,t['active_F1_fixed67_columns']));norm={};reps={}
    for r in t['active_F1_fixed67_receipts']:
        j=r['source_column'];v=dict(fixed[j]);w={j:F(1)}
        for p,c in r['subtract']:v=add(v,norm[p],-F(c));w=add(w,reps[p],-F(c))
        d=F(r['divide_by']);p=r['pivot'];norm[p]=guard(scale(v,1/d));reps[p]=guard(scale(w,1/d));assert norm[p]==parse(r['normalized'])
    alpha=[];I2=[]
    for j in free:
        v=dict(D2[j]);w={}
        while v:
            p=min(v);assert p in norm;c=v[p];w=add(w,reps[p],c);v=add(v,norm[p],-c)
        assert all(i>=11 for i in w)
        a={i-11:c for i,c in w.items()};assert mv([D2[k] for k in selected],a)==D2[j]
        i2=add({j:F(1)},{selected[k]:c for k,c in a.items()},-1)
        assert not mv(D2,i2);alpha.append(a);I2.append(i2)
    blocks1=[(1,())]+[(0,(i,)) for i in active];target={w:k for k,(p,w) in enumerate(blocks1)}
    blocks2=[(2-q,w) for q in range(3) for w in itertools.combinations(active,q)]
    varidx={i:j for j,i in enumerate(active)};linear=[]
    for p,w in blocks2:
        for j in range(14):
            v={}
            if p:v[(14*target[w]+j)*6+5]=F(1)
            for s,i in enumerate(w):
                row=14*target[w[:s]+w[s+1:]]+j;v=add(v,{row*6+varidx[i]:F((-1)**(p+s))})
            linear.append(v)
    assert len(linear)==224
    cols=[]
    for i2 in I2:
        v=mv(linear,i2);q={}
        for k,c in v.items():
            r,var=divmod(k,6)
            q=add(q,{j*6+var:c*a for j,a in Pi[r].items()})
        cols.append(guard(q))
    assert len(cols)==42
    preflight={'source':'actual canonical d2linear, fixedunit constants andPi17 only','matrix_shape':[102,42],
      'upper_dense_slots':102*42,'no_physical_source_jet_or_rank_field_sample':True,
      'total_live_bound':20000,'limits':{'seconds':60,'bytes':2_000_000_000,'artifact_bytes':32*1024*1024,'coefficients':150000,'bits':384}}
    (out/'preflight.json').write_text(json.dumps(preflight,separators=(',',':'))+'\n')
    data={'free42':free,'alpha42':[srec(v) for v in alpha],'I2_origin42':[srec(v) for v in I2],
      'linear_G2_columns':[srec(v) for v in cols],'source_formula':'Pi17*d2linear*I2origin',
      'source_zero_checks':'all42 d2origin*I2origin=0; allinputsourcecoefficientsalreadyfixed',
      'source_inputs':hs}
    ep=out/'known-G2-linear-input-and-provenance.json';ep.write_text(json.dumps(data,separators=(',',':'))+'\n')
    count=sum(map(len,cols))+sum(map(len,alpha))+sum(map(len,I2))+sum(map(len,Pi))+sum(map(len,norm.values()))+sum(map(len,reps.values()))
    assert count+102*42+42*42+42*41//2<20000
    report={'id':'AMBIENT_B18_ELL2_C2_G2_LINEAR_INPUT_GATE_20261105_v63','all_passed':True,
      'rank_started':False,'input_source_coefficients_queried':False,'active_linear_shape':[102,42],
      'whole_named_slot_upper':count+102*42+42*42+42*41//2,'input_sha256':hs,'code_sha256':sha(Path(__file__)),
      'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-start}
    (out/'report.json').write_text(json.dumps(report,separators=(',',':'))+'\n');print(json.dumps(report))

def rank_only(out,gate_path):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024))
    out=out.resolve();out.mkdir(parents=True,exist_ok=False)
    gp=ROOT/gate_path;g=json.loads(gp.read_text());assert g['all_passed'] and not g['rank_started']
    ip=gp.parent/'known-G2-linear-input-and-provenance.json';assert sha(ip)==g['evidence_sha256'][ip.name]
    data=json.loads(ip.read_text());cols=list(map(parse,data['linear_G2_columns']));alpha=list(map(parse,data['alpha42']));I2=list(map(parse,data['I2_origin42']));free=data['free42']
    hs={str(gp.relative_to(ROOT)):sha(gp),str(ip.relative_to(ROOT)):sha(ip)}
    piv={};receipts=[];reps2={};kernel=[]
    for j,orig in enumerate(cols):
        v=dict(orig);w={j:F(1)};steps=[]
        while v:
            p=min(v)
            if p not in piv:
                d=v[p];piv[p]=guard(scale(v,1/d));reps2[p]=guard(scale(w,1/d));receipts.append({'source_column':j,'pivot':p,'divide_by':str(d),'subtract':[[i,str(c)] for i,c in steps],'normalized':srec(piv[p])});break
            d=v[p];steps.append((p,d));v=guard(add(v,piv[p],-d));w=guard(add(w,reps2[p],-d))
        if not v:kernel.append(w)
    assert all(not mv(cols,k) for k in kernel)
    porder=sorted(piv);pidx={p:i for i,p in enumerate(porder)};factors=[]
    for orig in cols:
        v=dict(orig);c={}
        while v:
            p=min(v);a=v[p];c[pidx[p]]=a;v=add(v,piv[p],-a)
        factors.append(c)
    B=[piv[p] for p in porder];assert all(mv(B,v)==c for v,c in zip(factors,cols))
    sel=[x['source_column'] for x in receipts];sequence=[x['pivot'] for x in receipts]
    sign=(-1)**sum(sequence[i]>sequence[j] for i in range(len(sequence)) for j in range(i+1,len(sequence)))
    evidence={'variable_order6':['zA0','zA1','zA4','zA5','rA','ell'],
      'source42_free_indices':free,'constant_middle56_coefficients_alpha42':[srec(v) for v in alpha],
      'actual_I2_source_origin42_columns':[srec(v) for v in I2],
      'actual_linear_G2_active102_by42_columns':[srec(v) for v in cols],
      'rank':len(piv),'fixed_rank_receipts':receipts,'normalized_basis':[srec(v) for v in B],
      'all42_column_factorizations':[srec(v) for v in factors],
      'fixed_minor_rows':porder,'fixed_minor_columns':sel,
      'fixed_minor_sparse_columns':[srec({pidx[i]:c for i,c in cols[j].items() if i in pidx}) for j in sel],
      'fixed_minor_nonzero_exact_product':{'sign':sign,'QQ_factors':[x['divide_by'] for x in receipts]},
      'exact_kernel_on_source42':[srec(v) for v in kernel],
      'full102_linear_rank':len(piv)+60,
      'linear_injectivity_full102':len(piv)==42,
      'formula_proof':'G2lin=Pi17*d2lin*I2origin; constantd2I2origin=0 killsprojectionderivative andPi17*d2origin=0 killsI2derivative',
      'tensor_injectivity_proof':'targetactive17 givesactiveG2lin*a -sum(zerozi*bi); variable separation forcesall51bi0, thenzeroKoszul2forces9ci0',
      'T2_origin_values_not_computed':True,'no_actual_X1_or_new_function_values_queried':True}
    ep=out/'actual-G2-linear-fixed-source-profile.json';ep.write_text(json.dumps(evidence,separators=(',',':'))+'\n')
    slots=sum(map(len,cols))+sum(map(len,alpha))+sum(map(len,I2))+sum(map(len,piv.values()))+sum(map(len,factors))+sum(map(len,kernel));assert slots<20000
    report={'id':STEM,'all_passed':True,'active_G2_linear_rank':len(piv),'active_source_dimension':42,
      'full_G2_linear_rank':len(piv)+60,'full_G2_source_dimension':102,'full_linear_injective':len(piv)==42,
      'physical_source_or_X1_values_queried':False,'actual_T2_or_H0_length_computed':False,
      'input_sha256':hs,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),
      'evidence_sha256':{p.name:sha(p) for p in out.iterdir()},'named_QQ_slots':slots,'seconds':time.monotonic()-start}
    (out/'report.json').write_text(json.dumps(report,separators=(',',':'))+'\n');print(json.dumps({k:report[k] for k in ['all_passed','active_G2_linear_rank','full_G2_linear_rank','full_linear_injective','named_QQ_slots','seconds']}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);p.add_argument('--rank-source-gate');a=p.parse_args()
    if a.rank_source_gate:rank_only(a.output,a.rank_source_gate)
    else:main(a.output)
