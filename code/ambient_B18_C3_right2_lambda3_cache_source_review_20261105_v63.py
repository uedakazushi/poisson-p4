#!/usr/bin/env python3
"""Independent fixed six-root SOURCE bounds; no lambda composition."""
import hashlib,importlib.util,json,math,pathlib,resource,signal,time
from collections import Counter
from fractions import Fraction as Q
ROOT=pathlib.Path(__file__).resolve().parents[1]
BASE='AMBIENT_B18_C3_RIGHT2_LAMBDA3_CACHE_SOURCE_20261105_v63'
ID=BASE+'_REVIEW'
PIN='ae2b9c8162f53a396dcfb52288edbb53035569742b384de9b34286cf8cdc0d53'
def sha(p):return hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(pathlib.Path(p).read_text())
def save(p,d):pathlib.Path(p).write_text(json.dumps(d,indent=2)+'\n')
def main():
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30))
    out=ROOT/'reports'/ID;out.mkdir(exist_ok=False);pins={}
    helper=ROOT/'code/ambient_B18_C3_alpha4_cached_source_review_20261105_v63.py'
    assert sha(helper)=='430cc27fc33b263cdcaec0fd616b87e88824bef702a3ff49fbc7136eaf428313'
    sp=importlib.util.spec_from_file_location('own_envelope',helper);own=importlib.util.module_from_spec(sp);sp.loader.exec_module(own)
    def pin(p,h=None):
        p=pathlib.Path(p);s=sha(p)
        if h is not None:assert s==h,(str(p),s,h)
        pins[str(p.relative_to(ROOT))]=s
    pin(ROOT/'reports'/BASE/'report.json',PIN);pr=read(ROOT/'reports'/BASE/'report.json')
    for p,h in pr['input_sha256'].items():pin(ROOT/p,h)
    for p,h in pr['evidence_sha256'].items():pin(ROOT/'reports'/BASE/p,h)
    c=read(ROOT/'reports'/BASE/'compiled-six-even-roots-and-provided-joint105-plane3-cache.json')
    profile=read(ROOT/'reports'/BASE/'lambda3-cache-and-right2-whole-source-profile.json')
    records=c['operator_records_originalID_op_args'];assert len(records)==7095
    oldpath=ROOT/'reports/AMBIENT_B18_C3_PURE_ELL_READBACK_SOURCE_20261105_v63/compiled-original164-DAG-and-provided-pure-ell-inputs.json'
    original=own.field(oldpath,'operator_records_originalID_op_args');k=0
    for r in original:
        if k<len(records) and r[0]==records[k][0]:assert r==records[k];k+=1
    assert k==7095;del original
    roots=[m[k] for m in c['root_pairs_physical3'] for k in ['aA_node','aB_node']]
    assert roots==[20114,40892,20138,40916,20153,40931]
    expected_meta=[r for r in own.field(oldpath,'physical12_coefficient_roots') if r['physical_row'] in [18,60,67]]
    assert c['root_pairs_physical3']==expected_meta
    formal=read(ROOT/'reports/AMBIENT_B18_C3_RIGHT2_QUARTIC_SOURCE_20261105_v63/two-coherent-HB-quartic-sources-and-unpaired-affine-operator.json')
    assert formal['pivot_original_roots']==expected_meta
    formalpeer=ROOT/'reports/AMBIENT_B18_C3_RIGHT2_QUARTIC_SOURCE_20261105_v63_REVIEW/report.json'
    pin(formalpeer,'05d8c4eb681961627e8d3d6210856052b52a6c6945efebfc055e2621b0b975af');assert read(formalpeer)['all_passed']
    needed=set(roots)
    for r in reversed(records):
        if r[0] in needed and r[1] not in ['input','c']:needed.update(r[2:])
    assert needed=={r[0] for r in records};del needed
    oldalpha=read(ROOT/'reports/AMBIENT_B18_C3_ALPHA4_CACHED_SOURCE_20261105_v63/compiled-eight-even-roots-and-provided-joint105-plane3-cache.json')
    assert c['provided105_originalID_r_s_coefficients10']==oldalpha['provided105_originalID_r_s_coefficients10']
    assert c['provided_monomials10']==oldalpha['provided_monomials10'];del oldalpha
    qpath=ROOT/'reports/AMBIENT_B18_CPRIME_QUADRATIC20_VALUES_20261105_v63/phase1-common105-and-114-source-equations.json'
    qm=own.field(qpath,'base27_monomials');qv=dict(own.field(qpath,'all105_cached_composite_values'))
    qi=[qm.index(m) for m in [[],[7],[8],[7,7],[7,8],[8,8]]]
    for i,values in c['provided105_originalID_r_s_coefficients10']:
        a=list(map(Q,values));lower=[a[0],a[1]+a[2],a[2]/5,a[3]+a[4]+a[5],(a[4]+2*a[5])/5,a[5]/25]
        assert lower==[Q(qv[i][j]) for j in qi]
    del qv,qm
    provided=own.input_envelopes(c['provided105_originalID_r_s_coefficients10'],c['provided_monomials10'])
    env=own.Majorants();uses=Counter();cache={};hasher=hashlib.sha256();peak=0;inverse_pre=0;inverses=[];fixtures=[]
    for r in records:
        if r[1] not in ['input','c']:uses.update(r[2:])
    uses.update(roots)
    for r in records:
        i,op=r[:2]
        if op=='input':v=provided[i]
        elif op=='c':assert len(r)==4 and r[3]>0;v=[Q(r[2],r[3])]+[(0,1)]*5;fixtures.append([i,r[2],r[3],str(v[0])])
        elif op=='a':v=env.vector_add(cache[r[2]],cache[r[3]])
        elif op=='n':v=[-cache[r[2]][0]]+cache[r[2]][1:]
        elif op=='m':v=env.conv(cache[r[2]],cache[r[3]])
        elif op=='i':
            v=env.inverse(cache[r[2]]);inverses.append([i,r[2],'1'])
            test=own.Majorants();test.conv(v,cache[r[2]]);test.conv(cache[r[2]],v);inverse_pre=max(inverse_pre,test.pre)
        else:raise AssertionError(op)
        env.stored(v);hasher.update((own.encode([i,env.serialize(v)])+'\n').encode());cache[i]=v;peak=max(peak,len(cache))
        if op not in ['input','c']:
            for j in r[2:]:
                uses[j]-=1
                if not uses[j]:del cache[j]
        if not uses[i]:del cache[i]
    assert peak==162 and inverses==c['both_original_inverse_centers']
    assert hasher.hexdigest()==profile['every_node_majorant_sha256']
    assert env.node_coarse==39 and [[i,env.serialize(cache[i])] for i in roots]==profile['root10_final_majorants']
    diff=own.Majorants();lambda_bounds=[]
    for meta in expected_meta:
        a,b=cache[meta['aA_node']],cache[meta['aB_node']];d=diff.vector_add(a,b)
        lambda_bounds.append([meta['physical_row'],list(d[5])])
    assert lambda_bounds==[[18,[1680,3456000]],[60,[726,1728000]],[67,[1094,6912000]]]
    assert max(env.pre,inverse_pre,diff.pre,40)<=79
    # The30 lower expected values remain the previously evaluated C2/ell
    # data, rather than values from this bound-only source computation.
    old=dict((i,dict(v)) for i,v in own.field(ROOT/'reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63/actual-source105-ell2-raw33-values.json','all_requested164_root_values'))
    ell=dict(own.field(ROOT/'reports/AMBIENT_B18_C3_PURE_ELL_READBACK_VALUES_20261105_v63/all164-cached-pure-ell-readbacks-and-physical12-second-values.json','all164_actual_pure_ell_root_values3'))
    basis=own.field(ROOT/'reports/AMBIENT_B18_ELL2_RAW33_QUOTIENT_20261105_v63/actual-ell2-product33-and-raw45-actions.json','actual33_basis')
    idx=[basis.index(m) for m in [[0,0,0],[0,1,0],[0,0,1],[0,1,1]]]
    lower=[[i,[old[i].get(j,'0') for j in idx]+[ell[i][2]]] for i in roots]
    assert lower==c['expected6_lower_B6_coefficients_order_1_r_ell_rell_ell2']
    for meta in expected_meta:
        assert dict(lower)[meta['aA_node']]==dict(lower)[meta['aB_node']]
    p={'all7095_original_records_minimal_closure_and_six_full_root_labels_verified':True,
      'all630_given105_input_lower_overlaps_verified':True,'same_provided105_full10_jets_as_alpha4_SOURCE_peer':True,
      'all30_prior_lower_root_values_verified':True,'no_lambda_or_requested_positive_root_composed':True,
      'majorant_hash':hasher.hexdigest(),'coefficient_height':env.node_coarse,
      'lastuse_peak':peak,'DAG_before_operation_bits':env.pre,'inverse_diagnostic_bits':inverse_pre,
      'delta_difference_bits':diff.pre,'lambda_bounds':lambda_bounds,'literal_n_d_fixtures':fixtures,
      'formal_RIGHT2_source_peer_pinned':True,'source_static_numeric_upper':125000,'hard_source_numeric':150000,
      'prepared_code_postroot_height_and_wholecustody_pending':True}
    save(out/'six-root-source-bound-checkpoint.json',p)
    save(out/'phase.json',{'phase':'SOURCE_CORE_PASS_NO_GO_pending_combined_worker_custody','no_lambda_values':True})
    print(json.dumps({'core_passed':True,'DAG_pre':env.pre,'inverse_pre':inverse_pre,'delta_pre':diff.pre,'lambda_bounds':lambda_bounds,'seconds':time.monotonic()-start}));signal.alarm(0)
if __name__=='__main__':main()
