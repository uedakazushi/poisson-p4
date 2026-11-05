#!/usr/bin/env python3
"""Independent fixed supplied-input jet2 replay; no source/gauge solve/import."""
import hashlib,json,pathlib,resource,signal,time
from fractions import Fraction as Q
ROOT=pathlib.Path(__file__).resolve().parents[1]
ID='AMBIENT_B18_C3_PURE_ELL_READBACK_VALUES_20261105_v63_REVIEW'
SRC='AMBIENT_B18_C3_PURE_ELL_READBACK_SOURCE_20261105_v63'
PRI='AMBIENT_B18_C3_PURE_ELL_READBACK_VALUES_20261105_v63'
PRIMARY_SHA='341251643c4936ca11d818ee69f7c8a336898cbea86b3bebd4b6dc8ad825d208'
SOURCE_PEER_SHA='58abe6866cb8959d2b73479128c80f3dba79a06f8383a2fe3d1eec1419a32d4b'
def sha(p):
    h=hashlib.sha256()
    with pathlib.Path(p).open('rb') as f:
        for b in iter(lambda:f.read(65536),b''):h.update(b)
    return h.hexdigest()
def read(p):return json.loads(pathlib.Path(p).read_text())
def save(p,x):pathlib.Path(p).write_text(json.dumps(x,indent=2)+'\n')
def ncount(x):
    if isinstance(x,dict):return sum(ncount(v) for v in x.values())
    if isinstance(x,(list,tuple)):return sum(map(ncount,x))
    if isinstance(x,int):return 1
    if isinstance(x,str):
        try:Q(x);return 1
        except (ValueError,ZeroDivisionError):return 0
    return 0
def main():
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30))
    out=ROOT/'reports'/ID;out.mkdir(exist_ok=False)
    try:
        pp=ROOT/'reports'/PRI/'report.json';assert sha(pp)==PRIMARY_SHA;pr=read(pp);assert pr['all_passed']
        sp=ROOT/'reports'/(SRC+'_REVIEW')/'report.json';assert sha(sp)==SOURCE_PEER_SHA;spr=read(sp);assert spr['all_passed']
        pins={str(pp.relative_to(ROOT)):PRIMARY_SHA,str(sp.relative_to(ROOT)):SOURCE_PEER_SHA}
        for p,h in pr['input_sha256'].items():assert sha(ROOT/p)==h;pins[p]=h
        assert pr['code_sha256']=='29308014c69da649c66596454c95a749c14d8724f129c4ac97f9cfb939da1ca4'
        codepath=ROOT/'code/ambient_B18_C3_pure_ell_readback_values_v2_20261105_v63.py';assert sha(codepath)==pr['code_sha256'];pins[str(codepath.relative_to(ROOT))]=sha(codepath)
        for p,h in pr['evidence_sha256'].items():
            ep=ROOT/'reports'/PRI/p;assert sha(ep)==h;pins[str(ep.relative_to(ROOT))]=h
        actual=read(ROOT/'reports'/PRI/'all164-cached-pure-ell-readbacks-and-physical12-second-values.json')
        source=read(ROOT/'reports'/SRC/'compiled-original164-DAG-and-provided-pure-ell-inputs.json')
        static=ncount(source)+3*len(source['operator_records_originalID_op_args'])+ncount(actual)+ncount(spr)+15000
        assert static<130000
        opbits=0;prebound=0;calls=0;rootbits=0
        def height(x):return max(abs(x.numerator).bit_length(),x.denominator.bit_length())
        def rat(s,d=None):
            nonlocal opbits
            if d is None:
                v=str(s).split('/');n=int(v[0]);d=int(v[1]) if len(v)==2 else 1
            else:n=s
            assert d and max(abs(n).bit_length(),abs(d).bit_length())<=2048
            z=Q(n,d);opbits=max(opbits,height(z));return z
        def a(x,y):
            nonlocal opbits,prebound,calls
            if not x:return y
            if not y:return x
            nx,dx=abs(x.numerator).bit_length(),x.denominator.bit_length();ny,dy=abs(y.numerator).bit_length(),y.denominator.bit_length()
            b=max(max(nx+dy,ny+dx)+1,dx+dy);prebound=max(prebound,b);assert b<=2048
            z=x+y;assert height(z)<=2048;opbits=max(opbits,height(z));calls+=1;return z
        def m(x,y):
            nonlocal opbits,prebound,calls
            if not x or not y:return rat(0)
            b=max(abs(x.numerator).bit_length()+abs(y.numerator).bit_length(),x.denominator.bit_length()+y.denominator.bit_length())
            prebound=max(prebound,b);assert b<=2048
            z=x*y;assert height(z)<=2048;opbits=max(opbits,height(z));calls+=1;return z
        zero=(rat(0),rat(0),rat(0));one=(rat(1),rat(0),rat(0))
        def plus(x,y):return tuple(a(p,q) for p,q in zip(x,y))
        def minus(x,y):return plus(x,tuple(-p for p in y))
        def product(x,y):
            # Reverse coefficient and summand order, independent of the source worker.
            ans=[rat(0),rat(0),rat(0)]
            for k in [2,1,0]:
                for i in reversed(range(k+1)):ans[k]=a(ans[k],m(x[i],y[k-i]))
            return tuple(ans)
        phi={i:tuple(rat(x) for x in v) for i,v in source['provided105_originalID_pure_ell_coefficients']}
        records=source['operator_records_originalID_op_args'];vals=[None]*(1+max(r[0] for r in records));inverses=[];stream=hashlib.sha256()
        for rec in records:
            i,op,*arg=rec
            if op=='input':v=phi[i]
            elif op=='c':v=(rat(arg[0],arg[1]),rat(0),rat(0))
            elif op=='a':v=plus(vals[arg[0]],vals[arg[1]])
            elif op=='n':v=tuple(-p for p in vals[arg[0]])
            elif op=='m':v=product(vals[arg[0]],vals[arg[1]])
            else:
                assert op=='i';x=vals[arg[0]];assert x[0]==1
                u=minus(x,one);v=plus(minus(one,u),product(u,u))
                assert product(v,x)==one
                inverses.append({'inverse_node':i,'input_node':arg[0],'value3':list(map(str,v)),'reversed_unit_product3':list(map(str,product(v,x)))})
            assert all(height(p)<=2048 for p in v)
            rootbits=max(rootbits,*(height(p) for p in v));vals[i]=v
            stream.update((json.dumps([i,list(map(str,v))],separators=(',',':'))+'\n').encode())
        assert [[r['inverse_node'],r['input_node'],'1'] for r in inverses]==source['unit_inverse_centers6']
        targets={p:(rat(0),rat(1,5),rat(0)) if j==8 else zero for j,p in enumerate(source['CA_target_positions'])}
        assert all(vals[i]==targets.get(j,zero) for j,i in enumerate(source['equation105_order']))
        assert all(vals[i][:2]==tuple(rat(x) for x in v) for i,v in source['expected164_constant_and_ell_coefficients_from_C2'])
        roots=[[i,list(map(str,vals[i]))] for i in source['requested_root_nodes164']]
        assert roots==actual['all164_actual_pure_ell_root_values3']
        nu=[list(map(str,vals[i])) for i in source['nu_roots_CA_CB']]
        assert nu==actual['actual_nuA_B_pure_ell3']==source['expected_nu_pure_ell3_from_Q20']
        assert 2*rootbits+2<=2048
        physical=[]
        for row in source['physical12_coefficient_roots']:
            delta=minus(vals[row['aA_node']],vals[row['aB_node']])
            physical.append({'physical_row':row['physical_row'],'source_roots':row,'delta_jet3':list(map(str,delta)),
              'bA_jet3':list(map(str,vals[row['bA_node']])),'bB_jet3':list(map(str,vals[row['bB_node']]))})
        assert physical==actual['all12_actual_physical_coefficients_pure_ell3']
        assert all(r['delta_jet3'][2]=='0' and r['bA_jet3'][2]=='0' for r in physical)
        nonzero=[[r['physical_row'],r['bB_jet3'][2]] for r in physical if Q(r['bB_jet3'][2])]
        assert nonzero==[[18,'-1/15'],[60,'-29/120']]
        result={'all164_values3':roots,'all12_physical_values3':physical,'actual_nu2_values3':nu,
          'all18584_own_node_value_stream_sha256':stream.hexdigest(),'all6_inverse_values_and_reversed_products':inverses,
          'all105_original_equations_and164_C2_first_overlaps_verified':True,
          'all164_full_second_values_and36_physical_jets_match_primary':True,
          'delta_ell2_all12_zero':True,'bA_ell2_all12_zero':True,'nonzero_bB_ell2':nonzero,
          'own_max_saved_DAG_coefficient_bits':rootbits,'own_max_temporary_QQ_result_bits':opbits,
          'own_max_before_operation_bound_bits':prebound,'guarded_nontrivial_QQ_operations':calls,
          'strict2048_checks_before_all_Fraction_add_multiply':True,'static_numeric_upper':static,
          'source_frame':'SAME original105/F9 and zero-based physical150 array; source peer58abe686, z7=rA=0,sA=ell/5.',
          'source_or_gauge_or_jet_solve_count':0,'C3_rank_or_torsion_exponents_claimed':False,
          'scope':'Only existing-cache pure-ell second readback, not an all-plane second-zero or cubic primitive claim.'}
        save(out/'independent164-jet2-and-physical12-fixed-value-proof.json',result)
        raw=sum((ROOT/p).stat().st_size for p in pins)+sum(p.stat().st_size for p in out.iterdir() if p.is_file())
        assert raw<2<<20 and time.monotonic()-start<60
        report={'id':ID,'all_passed':True,'frozen':True,'primary_report_sha256':PRIMARY_SHA,
          'source_peer_sha256':SOURCE_PEER_SHA,'input_sha256':pins,'code_sha256':sha(__file__),
          'evidence_sha256':{'independent164-jet2-and-physical12-fixed-value-proof.json':sha(out/'independent164-jet2-and-physical12-fixed-value-proof.json')},
          'all_original_source_frame_and164_values_accepted':True,'no_author_producer_imports':True,
          'no_new_source_gauge_higherjet_solve':True,'no_C3_rank_or_exponent_claim':True,
          'raw_input_plus_output_bytes':raw,'static_numeric_upper':static,'elapsed_seconds':time.monotonic()-start}
        save(out/'report.json',report);signal.alarm(0)
        print(json.dumps({'all_passed':True,'report':str(out/'report.json'),'sha256':sha(out/'report.json'),'seconds':report['elapsed_seconds'],'nonzero_bB_ell2':nonzero}))
    except Exception as e:
        save(out/'failure.json',{'error':type(e).__name__+': '+str(e),'no_positive_claim':True});raise
if __name__=='__main__':main()
