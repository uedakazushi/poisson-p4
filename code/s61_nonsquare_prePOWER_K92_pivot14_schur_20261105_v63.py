#!/usr/bin/env python3
"""Source-only PIVOT14: finite determinant operand and faithful 1x8 operators."""
import argparse, hashlib, json, resource, signal, time
from fractions import Fraction as Q
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
ID='S61_NONSQUARE_PREPOWER_K92_PIVOT14_SCHUR_20261105_v63'

def sha(p):
    with p.open('rb') as f: return hashlib.file_digest(f,'sha256').hexdigest()
def read(p): return json.loads(p.read_text())
def dump(p,x): p.write_text(json.dumps(x,indent=2)+'\n')

def universal_check():
    """Exact free Laurent arithmetic, only s is invertible; no source values."""
    nv=18; z={}; one={(0,)*nv:Q(1)}
    def var(i):
        m=[0]*nv; m[i]=1; return {tuple(m):Q(1)}
    def add(a,b):
        c=a.copy()
        for m,v in b.items():
            c[m]=c.get(m,Q(0))+v
            if not c[m]: del c[m]
        return c
    def neg(a): return {m:-v for m,v in a.items()}
    def mul(a,b):
        c={}
        for m,v in a.items():
            for n,w in b.items():
                e=tuple(x+y for x,y in zip(m,n)); c[e]=c.get(e,Q(0))+v*w
                if not c[e]: del c[e]
        return c
    def mm(a,b):
        assert len(a[0])==len(b)
        c=[[{} for _ in b[0]] for _ in a]
        for i in range(len(a)):
            for j in range(len(b[0])):
                for k in range(len(b)): c[i][j]=add(c[i][j],mul(a[i][k],b[k][j]))
        return c
    def eye(n): return [[one if i==j else z for j in range(n)] for i in range(n)]
    s=var(0); im=[0]*nv; im[0]=-1; si={tuple(im):Q(1)}
    b=[var(1+i) for i in range(8)]; c=var(9); d=[var(10+i) for i in range(8)]
    A=[[s]+b,[c]+d]
    L=[[si,z],[neg(mul(c,si)),one]]; Li=[[s,z],[c,one]]
    R=eye(9); Ri=eye(9)
    for j in range(8): R[0][j+1]=neg(mul(si,b[j])); Ri[0][j+1]=mul(si,b[j])
    checks=[mm(L,Li)==eye(2),mm(Li,L)==eye(2),mm(R,Ri)==eye(9),mm(Ri,R)==eye(9)]
    expected=[[one]+[z]*8,[z]+[add(d[j],neg(mul(mul(c,si),b[j]))) for j in range(8)]]
    checks.append(mm(mm(L,A),R)==expected)
    assert all(checks)
    return dict(ring='QQ[s,s^-1,b0..b7,c,d0..d7]',all_five_matrix_identities=True,
                no_parameter_coefficients_evaluated=True,new_inverse_variables=['s'])

def main(out):
    start=time.monotonic(); resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000)); signal.alarm(60)
    out.mkdir(parents=True,exist_ok=False); pins={}
    def pin(rel):
        p=ROOT/rel; pins[rel]=sha(p); return p
    def load(st,name='report.json'): return read(pin(f'reports/{st}/{name}'))
    prev='S61_NONSQUARE_PREPOWER_K120_PIVOT13_SCHUR_20261104_v63'
    po=load(prev,'operator.json'); assert load(prev)['all_passed']; assert load(prev+'_REVIEW')['all_passed']
    en='S61_NONSQUARE_PREPOWER_K92_WHOLE_AXIS_BAREISS_ENTRY_20261105_v63'
    er=load(en); ax=load(en,'axis-rank-and-determinant-normalization.json')
    assert er['all_passed'] and er['rank_over_QQ_a']==11 and er['axis_determinant_NONZERO']
    assert ax['exact_identity']=='detM(a,0,0)=a^5*detCaxis(a)/prod(axis_rowLCD)'
    co=dict(ax['axis_detM_e0d0_QQ_coefficients']); assert Q(co[12])!=0
    raw='S61_NONSQUARE_PREPOWER_RAWCROSS18_ENTRY_20261104_v63'
    assert load(raw+'_REVIEW')['all_passed']
    p12='S61_NONSQUARE_PREPOWER_G115_PIVOT12_SCHUR_20261104_v63'
    assert load(p12+'_REVIEW')['all_passed']
    lift='S61_NONSQUARE_PREPOWER_FULL57_SOURCE_LIFT_20261104_v63'
    assert load(lift+'_REVIEW')['all_passed']; load(lift,'operator.json')
    src='S61_NONSQUARE_PREPOWER_K92_T18_SOURCE_PREFLIGHT_20261104_v63'
    recipe=load(src,'fixed-K92-target18-recipe.json')
    assert recipe['conditional_actual_next_entry']=='c8*g*det11/(Delta407*a²*N7) oncurrentchartifraw8prerequisiteholds'
    ts='S61_NONSQUARE_PREPOWER_K92_TROPICAL_SOURCE_PREFLIGHT_V2_20261105_v63'
    binding=load(ts,'source-binding.json')
    assert binding['det_normalization'].startswith('original fixed det11=a²')
    state_rel='reports/S61_NONSQUARE_PREPOWER_K92_DET11_GUARDED_ATTEMPT_20261104_v63/state.bin'
    state=pin(state_rel); assert pins[state_rel]==binding['state_sha256']=='cd5fb2431c2f65649153d82890d52b52ce026f9d0e478bc30592687f17b20f02'
    decoder=pin('code/s61_nonsquare_prePOWER_K92_tropical_source_preflight_v2_20261105_v63.py')
    original_binding='reports/S61_NONSQUARE_PREPOWER_K92_DET11_GUARDED_ATTEMPT_20261104_v63/source-binding.json'
    pin(original_binding)
    hc='S61_NONSQUARE_CUBIC_GENUINE14_H5_PREFLIGHT_20261104_v63'
    h=load(hc,'certificate.json'); c8=Q(h['known_row_units_rho14'][8]['constant']); assert c8==Q(-88000000,13851)
    assert len(h['raw_compact14x5'][6][1])==26
    N=[]
    def n(op,shape,args=(),**kw):
        assert all(0<=a<len(N) for a in args)
        i=len(N);N.append(dict(id=i,op=op,shape=shape,args=list(args),**kw));return i
    def ref(st,name,ptr,shape):
        rel=f'reports/{st}/{name}';p=pin(rel)
        return n('external_frozen_reference',shape,reference=dict(path=rel,sha256=sha(p),json_pointer=ptr))
    def root(st,name,shape):
        op=load(st,'operator.json');idx=op['exports'][name];assert op['nodes'][idx]['shape']==shape
        return ref(st,'operator.json','/exports/'+name,shape)
    def eye(m):return n('identity',[m,m])
    def zero(m,k):return n('zero',[m,k])
    def prod(a,b):
        assert N[a]['shape'][1]==N[b]['shape'][0]
        return n('matrix_product',[N[a]['shape'][0],N[b]['shape'][1]],[a,b])
    def neg(a):return n('negative',N[a]['shape'],[a])
    def sub(a,b):
        assert N[a]['shape']==N[b]['shape'];return n('difference',N[a]['shape'],[a,b])
    def scale(s,a):
        assert N[s]['shape']==[1,1];return n('central_scalar_multiple',N[a]['shape'],[s,a])
    def sel(a,rows,cols):return n('submatrix',[len(rows),len(cols)],[a],rows=rows,columns=cols)
    def block(rows):
        hs=[N[r[0]]['shape'][0] for r in rows];ws=[N[x]['shape'][1] for x in rows[0]]
        for i,r in enumerate(rows):
            for j,x in enumerate(r):assert N[x]['shape']==[hs[i],ws[j]]
        return n('block_matrix',[sum(hs),sum(ws)],[x for r in rows for x in r],blocks=rows)
    def diag(a,b):return block([[a,zero(N[a]['shape'][0],N[b]['shape'][1])],[zero(N[b]['shape'][0],N[a]['shape'][1]),b]])
    Y=root(prev,'residual2x9',[2,9]);Ls=root(prev,'fullsource15',[15,15]);Lsi=root(prev,'fullsource15_inverse',[15,15])
    Us=root(prev,'fulltarget22',[22,22]);Usi=root(prev,'fulltarget22_inverse',[22,22])
    mat=n('frozen_saved_integer_polynomial_matrix',[11,11],state=dict(path=state_rel,sha256=sha(state),version='K92SP001'),
          decoder=dict(path=str(decoder.relative_to(ROOT)),sha256=sha(decoder),function='selected_inputs',argument_roots_reference=original_binding+'#/actual_pure_matrix_roots'),
          matrix_state_ids=binding['matrix11x11_state_ids'],literal_zero_state_id=-1,
          entry_semantics='QQcommonfactor times saved integer polynomial in nonnegative (a,e,d)',
          normalization=binding['det_normalization'],physical_targets=binding['physical_targets'])
    P=n('universal_polynomial_determinant',[1,1],[mat],expanded=False,order=11,variable_order=['a','e','d'])
    iP=n('conditional_chart_inverse',[1,1],[P],chart='D(P_K92)',primary_nonzero=en,independent_nonzero_peer_pending=True)
    Delta=ref('S61_NONSQUARE_RAW_H_GENERIC_RANK_20261104_v63','certificate.json','/nonzero_selected_minor_Delta',[1,1])
    ng=ref(hc,'certificate.json','/raw_compact14x5/6/1',[1,1])
    ig=n('inherited_parameter_unit_inverse',[1,1],[ng],gate=p12+'_REVIEW',proof='g*D5*r=c10*a²*N6; N6 and c10*a² are units on old chart, hence each product factor is a unit; no new g chart')
    c8n=n('QQ_constant',[1,1],value=str(c8));ic8=n('QQ_constant',[1,1],value=str(1/c8))
    N7=ref('S61_NONSQUARE_PREPOWER_K120_T14_NESTED_ENTRY_20261104_v63','native.json','/normalized_N7',[1,1])
    ir=prod(prod(prod(Delta,N7),ic8),prod(ig,iP))
    cp=[8,0,1,2,3,4,5,6,7];Pc=n('column_permutation',[9,9],order=cp);Pci=n('transpose',[9,9],[Pc])
    YY=prod(Y,Pc);s=sel(YY,[0],[0]);b=sel(YY,[0],list(range(1,9)));c=sel(YY,[1],[0]);d=sel(YY,[1],list(range(1,9)))
    ci=prod(c,ir);ib=scale(ir,b);S=sub(d,prod(ci,b));cleared=sub(scale(s,d),prod(c,b))
    L=block([[ir,zero(1,1)],[neg(ci),eye(1)]]);Li=block([[s,zero(1,1)],[c,eye(1)]])
    U=block([[eye(1),neg(ib)],[zero(8,1),eye(8)]]);Ui=block([[eye(1),ib],[zero(8,1),eye(8)]])
    W=sel(L,[1],[0,1]);src15=prod(diag(eye(13),L),Ls);src15i=prod(Lsi,diag(eye(13),Li))
    tgt22=prod(Us,diag(eye(13),prod(Pc,U)));tgt22i=prod(diag(eye(13),prod(Ui,Pci)),Usi)
    ex=dict(residual1x8=S,cleared_residual1x8=cleared,pivot_scalar=s,pivot_scalar_inverse=ir,P_K92=P,
            fullsource15=src15,fullsource15_inverse=src15i,fulltarget22=tgt22,fulltarget22_inverse=tgt22i)
    for suffix,k in [('source_original25',25),('corrected_full57',57),('corrected_original63',63),('actual_pure84',84),('raw7',7)]:
        ex['surviving1_'+suffix]=prod(W,root(prev,'surviving2_'+suffix,[2,k]))
    op=dict(nodes=N,exports=ex,chart='FULLU D(Delta407*Saux*T*N6*N7bar*P_K92)',
            source_pivot_order_current15=po['source_pivot_order_current15']+[9],source_remaining_current15=[14],
            target_pivot_original32=po['target_pivot_original32']+[18],target_remaining_original32=[1,3,4,6,7,8,9,17],
            all_tail_effects_inherited=True,no_new_expanded_parameter_coefficients=True,
            only_new_chart_factor='P_K92',independent_nonzero_peer_pending=True)
    dump(out/'operator.json',op)
    proof=dict(exact_pivot_identity='s*Delta407*N7=c8*g*P_K92',old_det11_identity='det11_old=a²*P_K92',
               c8=str(c8),g_pointer='/raw_compact14x5/6/1',g_is_NOT_G20_anchor=True,
               old_chart='FULLU D(Delta407*Saux*T*N6*N7bar)',only_new_inverse='P_K92',
               no_new_g_D5_p_R_h22_cleanup_deltaP_deltaQ_N0_chart=True,
               nonzero_primary=dict(axis_coefficient_monomial=[12,0,0],QQ_coefficient=co[12],
                 proof='restriction of nonnegative-exponent actual generic polynomial matrix; nonzero coefficient implies P_K92 not zero',independent_peer_pending=True),
               column_permutation_sign=1,block_identity='L*(Y*Pc)*U=diag(1,residual1x8)',
               reverse_source_identity='all old source vectors multiplied by [-c/s,1]; raw7 stays zero by exact linearity',
               universal_free_Laurent_check=universal_check(),
               full57_scope='1x8 full row rank SUFFICES for actual57 on this chart via accepted FULL57_SOURCE_LIFT; not equality with actual bad Fitting ideal',
               no_residual_rank_global57_global239_or_exhaustion_claim=True)
    dump(out/'proof.json',proof)
    pin(str(Path(__file__).relative_to(ROOT)));pin(f'research/notes/{ID}.md')
    total=sum(p.stat().st_size for p in out.iterdir());assert total<1024*1024 and len(N)<200
    report=dict(id=ID,all_passed=True,frozen=True,source_only=True,conditional_on_independent_nonzero_peer=True,
                typed_nodes=len(N),same_cokernel_shape=[1,8],pivot_count=14,new_chart_factor='P_K92',
                no_new_scalar_source_or_determinant_calculation=True,actual_full57_source_lift_retained=True,
                no_global57_or_actual_bad_ideal_equality=True,total_new_RAW_before_report=total,
                input_sha256=pins,evidence_sha256={p.name:sha(p) for p in out.iterdir()},seconds=time.monotonic()-start)
    dump(out/'report.json',report);signal.alarm(0)
    print(json.dumps({k:v for k,v in report.items() if k not in ['input_sha256','evidence_sha256']},indent=2))

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    main(parser.parse_args().output.resolve())
