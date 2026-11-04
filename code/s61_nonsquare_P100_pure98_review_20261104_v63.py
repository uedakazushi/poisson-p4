#!/usr/bin/env python3
"""Independent source-coefficient recurrence and one fixed QQ[z] pure-row check."""
import argparse,gzip,hashlib,itertools,json,resource,signal,subprocess,time
from fractions import Fraction as F
from functools import lru_cache
from pathlib import Path
from s61_nonsquare_P100_H1_review_20261104_v63 import native_matrix
from s61_nonsquare_raw_H_B_curve_components_review_20261104_v63 import poly,mul
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_P100_Q5_PURE98_ROW_20261104_v63_REVIEW'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main(out):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
    out=out.resolve();out.mkdir(parents=True,exist_ok=False);work=ROOT/'work/P100-pure98-independent-review';work.mkdir(parents=True,exist_ok=True);hashes={};base=ROOT/'reports'
    def read(p):hashes[str(p.relative_to(ROOT))]=sha(p);return json.loads(p.read_text())
    def audit(p,passed=False):
        r=read(p)
        if passed:assert r['all_passed']
        for name,h in r.get('input_sha256',{}).items():assert sha(ROOT/name)==h,name;hashes[name]=h
        for field in ['evidence_sha256','native_evidence_sha256']:
            for name,h in r.get(field,{}).items():assert sha(p.parent/name)==h,name;hashes[str((p.parent/name).relative_to(ROOT))]=h
        if not passed:
            for key,name in [('native_input_sha256','native-input.txt'),('native_output_sha256','native.json'),('certificate_sha256','certificate.json'),('manifest_sha256','manifest.json'),('cubic_input_sha256','cubic-coefficient-input.json')]:
                if key in r:assert sha(p.parent/name)==r[key];hashes[str((p.parent/name).relative_to(ROOT))]=r[key]
        return r
    b=base/'S61_NONSQUARE_P100_Q5_PURE98_ROW_20261104_v63';r=audit(b/'report.json');assert r['complete']
    gb=base/'S61_NONSQUARE_P100_PURE98_ROW_GCD_20261104_v63';gr=audit(gb/'report.json');assert gr['row_unimodular_on_whole_P100']
    for stem in ['S61_NONSQUARE_P100_H1_20261104_v63_REVIEW','S61_NONSQUARE_ORACLE_FRAME_ADAPTER_20261103_v63_REVIEW',
        'S61_NONSQUARE_JOINT_LELL_FREE164_20261104_v63_REVIEW','S61_NONSQUARE_H_SQUARE_H1_PREFLIGHT_20261104_v63_REVIEW',
        'S61_NONSQUARE_JOINT_LELL_QPAIR_UNITS_20261104_v63_REVIEW']:audit(base/stem/'report.json',True)
    pf=base/'S61_NONSQUARE_P100_PURE98_SOURCE_PREFLIGHT_20261104_v63';audit(pf/'report.json')
    ci=read(pf/'cubic-coefficient-input.json');manifest=read(pf/'manifest.json')
    fd=read(base/'S61_NONSQUARE_JOINT_LELL_FREE164_20261104_v63/certificate.json')
    h=read(base/'S61_NONSQUARE_H_SQUARE_H1_PREFLIGHT_20261104_v63/certificate.json')
    ri=read(base/'S61_NONSQUARE_P100_H1_RIGHT_INVERSE_20261104_v63/native.json')
    binding=read(b/'source-binding.json');cube=read(b/'cubic504.json');cof=read(b/'cofactor42.json');jc=read(b/'joint-cubic504.json');j4=read(b/'joint210.json');pure=read(b/'pure98.json');gcd=read(gb/'native.json')
    qs=read(base/'S61_NONSQUARE_JOINT_LELL_QPAIR_UNITS_20261104_v63/checkpoint.json')
    cp=base/'S61_NONSQUARE_H1_KNOWN_UNIT_CONTINUATION_20261104_v63/checkpoint.json.gz';hashes[str(cp.relative_to(ROOT))]=sha(cp)
    with gzip.open(cp,'rt')as f:cache=json.load(f)
    q5=next(row for i,row in cache['materialized_raw_rows']if i==0);del cache
    classes=h['target122_Hdegree_classes'];remaining=qs['remaining_columns'];assert len(q5)==23
    for i in classes['0']+classes['2']:assert not qs['matrix58x136'][57][remaining[i]][0]
    for j,i in enumerate(classes['1']):
        n,d=qs['matrix58x136'][57][remaining[i]];nn,dd=q5[j]
        assert mul(poly(n),poly(dd))==mul(poly(nn),poly(d))
    overlap=read(base/'S61_NONSQUARE_JOINT_LELL_OVERLAP_UNIT_20261104_v63/remaining-pair-block.json')
    assert overlap['source_rows'][57]==[5,[0,0,0,0,0,0,0,1]]
    assert binding['actual_q5H_cache_row']==0 and binding['actual_qpair_source_row']==57
    assert cube['normalized_Tn_FrameAdapter_ROWUNIT_DIVISION']
    active=[i for i,row in enumerate(ri['right_inverse112x23'])if any(v['numerator']for v in row)]
    assert active==[7*s+j for s in[1,6,8,9,11,14]for j in range(7)]==binding['active_rightinverse_source_rows42']
    # Independent literal coefficient recurrence; only the accepted generic
    # atom/interner operations are reused, never Oracle.coefficient/core_value.
    from s61_nonsquare_lazy_quartic_oracle_20261103_v63 import Oracle
    import sympy as sy
    adapter=read(base/'S61_NONSQUARE_ORACLE_FRAME_ADAPTER_20261103_v63/adapter.json')
    oraclefile=ROOT/adapter['oracle_path'];original=read(oraclefile);assert sha(oraclefile)==adapter['oracle_sha256']
    o,_=Oracle.load(oraclefile);oldnodes=original['arithmetic']['nodes'];weights=original['arithmetic']['node_weights'];vw=original['arithmetic']['variable_weights']
    def forbidden(*args,**kwargs):raise AssertionError('author coefficient/core-value method is not used')
    o.coefficient=forbidden;o.core_value=forbidden
    @lru_cache(None)
    def coefficient(i,m):
        if not i or sum(w*k for w,k in zip(vw,m))!=weights[i]:return 0
        op,*a=oldnodes[i]
        if op=='v':return o.const(1)if m==tuple(int(j==a[0])for j in range(8))else 0
        if op=='z':return 0
        if op=='l':return o.lin([(o.atoms[c],coefficient(j,m))for c,j in a[0]])
        assert op=='m';terms=[]
        for left in itertools.product(*(range(k+1)for k in m)):
            if sum(w*k for w,k in zip(vw,left))!=weights[a[0]]:continue
            right=tuple(k-j for k,j in zip(m,left));x=coefficient(a[0],left);y=coefficient(a[1],right)
            if x and y:terms.append((1,o.mul(x,y)))
        return o.lin(terms)
    D=[o.K.from_expr(sy.sympify(c))for c in adapter['Tden']]
    units=[o.K.from_expr(sy.sympify(c))for c in adapter['basis57_encoded_source_units']]
    assert len(ci['sources6'])==6 and len(ci['old_physical_monomials84'])==84
    for source,expected_i in zip(ci['sources6'],[1,6,8,9,11,14]):
        pos=41+expected_i;record=adapter['basis57_source_definitions'][pos]
        assert source['basis57_position']==pos and source['genuine_index']==expected_i and record['kind']=='genuine_cubic'
        assert source['actual135_row']==record['original_equation']==h['all112_cubic_source319_definitions'][7*expected_i][2]
        rebuilt=[]
        for mon in ci['old_physical_monomials84']:
            assert sum(mon)==3 and mon[5]==0
            raw=coefficient(original['residual_weight_le4_roots'][record['original_equation']],tuple(mon))
            sc=o.K.one
            for v,k in zip(D,mon):sc*=v**k
            # Keep the two mathematically distinct operations visible.
            scaled=o.lin([(sc,raw)]);normalized=o.lin([(1/units[pos],scaled)]);rebuilt.append(normalized)
        assert rebuilt==source['roots84']
    for i,node in ci['needed_nodes']:assert o.nodes[i]==node
    def kpoly(p):return o.R.from_dict({tuple(m):o.R.domain.convert(F(c))for m,c in p})
    for i,(n,d)in ci['needed_atoms']:assert o.atoms[i]==o.K(kpoly(n))/o.K(kpoly(d))
    assert len(ci['needed_nodes'])==19960 and len(ci['needed_atoms'])==604
    del o;coefficient.cache_clear()
    # Rebuild both NF dependency closures from the accepted full projection.
    NF=fd['full_target_projection239_to122'];ln=fd['linear_nodes'];sn=fd['scalar_circuit']['nodes'];atoms=fd['scalar_circuit']['atoms_num_den']
    purecols=classes['0'];roots=[NF['output122_linear_roots'][j]for j in purecols]
    Lneed=set();Sroots=set();todo=roots[:]
    while todo:
        i=todo.pop()
        if i in Lneed:continue
        Lneed.add(i);op,*a=ln[i]
        if op=='l':
            for c,j in a[0]:Sroots.add(c);todo.append(j)
        else:assert op in['z','v']
    Sneed=set();SA=set();todo=list(Sroots)
    while todo:
        i=todo.pop()
        if i in Sneed:continue
        Sneed.add(i);op,*a=sn[i]
        if op=='c':SA.add(a[0])
        elif op in['a','s','m']:todo+=a
        else:assert op=='z'
    assert sorted(Lneed)==manifest['NF_needed_linear_ids']and sorted(Sneed)==manifest['NF_needed_scalar_ids']and sorted(SA)==manifest['NF_needed_atom_ids']
    assert roots==manifest['NF_pure98_linear_roots']and purecols==pure['target122_columns']
    assert [len(Lneed),len(Sneed),len(SA)]==[547,244,122]
    pre=dict(pinned_coefficient_queries=504,coefficient_DAG_nodes=19960,coefficient_atoms=604,
        source_only_coefficient_recurrence_checked=True,ROWUNIT_operation='DIVIDE',NF_nodes=[244,547],
        source_static_coefficient_lastuse_bound=manifest['bounds']['oracle_coefficient_live'],
        actual_RF_guard_includes_inputs_outputs_temporaries_unit_powers=150000,
        seconds_guard=60,memory_bytes_guard=2_000_000_000,output_bytes_guard=32*1024*1024,
        no_new_source_monomial_or_quartic_query=True)
    (out/'preflight-and-source-check.json').write_text(json.dumps(pre,indent=2)+'\n')
    invfile=base/'S61_NONSQUARE_JOINT_LELL_PREFLIGHT_20261104_v63/coordinate-inverse-numerators7x7.txt';forwardfile=invfile.with_name('coordinate-forward7x7.txt')
    n,m,inverse=native_matrix(invfile);assert(n,m)==(7,7);n,m,forward=native_matrix(forwardfile);assert(n,m)==(7,7)
    for p in[invfile,forwardfile]:hashes[str(p.relative_to(ROOT))]=sha(p)
    raw=read(base/'S61_NONSQUARE_RAW_H_LINEAR_UNITS_20261104_v63/preflight.json')
    oldmap={x['old112_index']:x['value']for x in cof['old_source42']};newmap={x['old112_index']:x['value']for x in cof['new_source42']};assert sorted(oldmap)==sorted(newmap)==active
    assert gcd['gcd_all98_numerators']==[[[0],'1']]and[g['position']for g in gcd['literal_Bezout_contributors']]==[0,1]
    assert[g['position']for g in gcd['global_regular_row_dual']]==[0,1]
    cpp=Path(__file__).with_suffix('.cpp')
    for p in[Path(__file__),cpp,ROOT/'code/s61_nonsquare_P100_H1_review_20261104_v63.cpp',ROOT/'code/s61_nonsquare_P100_H1_review_20261104_v63.py',ROOT/'code/s61_nonsquare_lazy_quartic_oracle_20261103_v63.py',ROOT/'code/s61_nonsquare_raw_H_B_curve_components_review_20261104_v63.py',ROOT/f'research/notes/{ID}.md',ROOT/'research/notes/S61_NONSQUARE_P100_Q5_PURE98_ROW_20261104_v63.md',ROOT/'research/notes/S61_NONSQUARE_P100_PURE98_ROW_GCD_20261104_v63.md']:hashes[str(p.relative_to(ROOT))]=sha(p)
    inp=work/'fixed-source-chain-input.txt'
    with inp.open('w')as f:
        def pw(p):
            p={int(m[0]):F(c)for m,c in p}if isinstance(p,list)else p
            f.write(str(len(p))+''.join(f' {i} {c}'for i,c in sorted(p.items()))+'\n')
        def sw(p):
            p={tuple(m):F(c)for m,c in p}if isinstance(p,list)else p
            f.write(str(len(p))+''.join(f' {a} {e} {d} {c}'for(a,e,d),c in sorted(p.items()))+'\n')
        def rw(r):pw(r['numerator']);pw(r['denominator'])
        def nw(items):
            f.write(str(len(items))+'\n')
            for i,node in items:
                if node[0]=='l':f.write(f'{i} l {len(node[1])}'+''.join(f' {c} {j}'for c,j in node[1])+'\n')
                else:f.write(str(i)+' '+' '.join(map(str,node))+'\n')
        def aw(items):
            f.write(str(len(items))+'\n')
            for i,(n,d)in items:f.write(str(i)+'\n');sw(n);sw(d)
        for c in[F(0),F(20),F(60),F(2700,203),F(12)]:pw({0:c,1:F(1)}if c else{1:F(1)})
        aw(ci['needed_atoms']);nw(ci['needed_nodes'])
        for row in ci['sources6']:f.write(' '.join(map(str,row['roots84']))+'\n')
        for mon in ci['old_physical_monomials84']:f.write(' '.join(str(mon[i])for i in[0,1,2,3,4,6,7])+'\n')
        for row in cube['values6x84']:
            for v in row:rw(v)
        for (n,d),expected in zip(q5,cof['q5H_actualROW0_values23']):sw(n);sw(d);rw(expected)
        f.write('42\n')
        for i in active:
            f.write(str(i)+'\n')
            for v in ri['right_inverse112x23'][i]:rw(v)
            rw(newmap[i])
        for i in active:f.write(str(i)+'\n');rw(oldmap[i])
        for p in inverse:sw(p)
        for p in forward:sw(p)
        for row in raw['constant_inverse7x7']:f.write(' '.join(map(str,row))+'\n')
        for row in jc['values6x84']:
            for v in row:rw(v)
        for mon in j4['monomials']:f.write(' '.join(map(str,mon))+'\n')
        for v in j4['values210']:rw(v)
        aw([(i,atoms[i])for i in sorted(SA)]);nw([(i,sn[i])for i in sorted(Sneed)]);nw([(i,ln[i])for i in sorted(Lneed)])
        for values in[roots,purecols,NF['input_variables']]:f.write(' '.join(map(str,values))+'\n')
        for mon in h['parent_joint_target239_monomials']:f.write(' '.join(map(str,mon))+'\n')
        for v in pure['values98']:rw(v)
        pw(gcd['gcd_all98_numerators'])
        for p in gcd['all98_quotients']:pw(p)
        f.write('2\n')
        for a,d in zip(gcd['literal_Bezout_contributors'],gcd['global_regular_row_dual']):f.write(str(a['position'])+'\n');pw(a['coefficient']);rw(d)
    assert inp.stat().st_size<16*1024*1024;exe=work/'verify-native'
    with(out/'compile.log').open('w')as log:subprocess.run(['c++','-O2','-std=c++17',str(cpp),'-o',str(exe),*subprocess.check_output(['pkg-config','--cflags','--libs','flint'],text=True).split()],stdout=log,stderr=subprocess.STDOUT,check=True,timeout=max(1,55-(time.monotonic()-start)))
    with inp.open()as f:proc=subprocess.run([str(exe)],stdin=f,text=True,capture_output=True,timeout=max(1,55-(time.monotonic()-start)))
    (out/'native.log').write_text(proc.stdout);(out/'native.stderr').write_text(proc.stderr)
    if proc.returncode:
        (out/'failure.json').write_text(json.dumps(dict(all_passed=False,returncode=proc.returncode,no_automatic_retry=True,reason=proc.stderr),indent=2)+'\n');raise SystemExit(proc.returncode)
    result=json.loads(proc.stdout);assert result['all_passed']and result['live_polynomial_slots_peak']<=150000
    proof=dict(actual_q5H_cache0_equals_QPAIR57_all23_polynomial_identities=True,q5H_pure_H2_zero_and_H2deletion_unchanged=True,
        independent504_source_coefficient_recurrence_against_original135=True,ROWUNIT_normalization_DIVISION_verified=True,
        accepted_generic_Oracle_atom_interner_only_reused_no_author_coefficient_or_core_value_calls=True,
        all504_normalized_and_all504_joint_values_independently_checked=True,all42_cofactor_sums_and_delta_undo=True,
        both49_physical_joint_inverse_products=True,all210_actual_negative_puretail_combinations=True,
        NF244_scalar547_linear_source_closures_independently_rebuilt=True,all98_actual_NF_values_checked=True,
        all98_nonzero_individual_named_units_none=True,literal_polynomial_Bezout_only_positions0_1_equals1=True,
        regular_rowdual_AiDi_and_pairing_equals1=True,actual_pure_row_unimodular_on_whole_P100=True,
        accepted141block_plus_unimodular_pure_row_gives_pointwise_original_quartic_rank_lower142=True,
        no_fixed_original142minor_columnGL_free_source_complement_rank98_full239_or_exhaustion_claim=True,
        no_other_source_monomial_newcorevalue_quartic134family_POWER_ranksearch_point_prime_GB_or_GCD=True,checker=result)
    (out/'independent-proof.json').write_text(json.dumps(proof,indent=2)+'\n')
    report=dict(id=ID,all_passed=True,actual_source_pure98_chain_independently_verified=True,actual_pure_row_unimodular_whole_P100=True,
        pointwise_original_quartic_rank_lower142=True,input_sha256=hashes,evidence_sha256={p.name:sha(p)for p in out.iterdir()},
        ignored_verifier_input_sha256=sha(inp),ignored_executable_sha256=sha(exe),seconds=time.monotonic()-start,max_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    p=out/'report.json';p.write_text(json.dumps(report,indent=2)+'\n');assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
    print(json.dumps(dict(all_passed=True,seconds=report['seconds'],rss=report['max_rss_kib'],report_sha256=sha(p))))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args();main(a.output)
