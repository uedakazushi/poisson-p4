#!/usr/bin/env python3
"""Source-only corrected A4+Schur10 last-row formula, no parameter expansion."""
import argparse,functools,hashlib,json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
ID='S61_NONSQUARE_PREPOWER_LASTROW17_CORRECTED_PAIR_SOURCE_20261105_v63'
def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def dump(p,x):p.write_text(json.dumps(x,indent=2)+'\n')

def main(out):
    started=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60)
    out.mkdir(parents=True,exist_ok=False);pins={};files={}
    def load(rel):
        if rel not in files:
            p=ROOT/rel;pins[rel]=sha(p);files[rel]=json.loads(p.read_text())
        return files[rel]
    def st(s,fn='operator.json'):return f'reports/{s}/{fn}'
    p14='S61_NONSQUARE_PREPOWER_K92_PIVOT14_SCHUR_20261105_v63';pp=st(p14);op=load(pp)
    gate='S61_NONSQUARE_PREPOWER_K92_PIVOT14_NONZERO_GATE_20261105_v63'
    assert load(st(gate,'report.json'))['all_passed']
    sf='S61_NONSQUARE_PREPOWER_PIVOT14_LASTROW_SOURCE_PROFILE_20261105_v63'
    profiles=load(st(sf,'eight-tail-source-map.json'));pr=next(x for x in profiles['profiles'] if x['original32_target']==17)
    staged='S61_NONSQUARE_PREPOWER_STAGED443_PREFLIGHT_20261104_v63';sop=load(st(staged))
    k11='S61_NONSQUARE_PREPOWER_KAPPA11_SCHUR_20261104_v63';ko=load(st(k11))
    qr='reports/S61_NONSQUARE_CUBIC_Q2_L_UNIT_COMPLETION_20261104_v63/certificate.json';q=load(qr)
    allcols=q['new22_target_original32_columns'];pivrows=op['source_pivot_order_current15'];pivcols=op['target_pivot_original32']
    freecols=op['target_remaining_original32'];freeidx=[allcols.index(x) for x in freecols];pividx=[allcols.index(x) for x in pivcols]
    assert pivrows[:4]==[2,4,5,8] and pivcols[:4]==[27,31,22,2]
    assert pivrows[4:]==[0,1,3,6,7,13,10,11,12,9]
    assert set(pivrows)==set(range(14)) and freecols==[1,3,4,6,7,8,9,17]
    # Structural abstract arithmetic: known rational constants or UNKNOWN.
    # UNKNOWN never cancels or becomes zero; thus every certified zero/one is
    # a source-independent identity, not a coefficient specialization.
    count=[0];cache={}
    def plus(a,b):return None if a is None or b is None else a+b
    def times(a,b):
        if a==0 or b==0:return Q(0)
        return None if a is None or b is None else a*b
    def mm(A,B):
        assert len(A[0])==len(B)
        C=[]
        for row in A:
            rr=[]
            for j in range(len(B[0])):
                v=Q(0)
                for k,a in enumerate(row):v=plus(v,times(a,B[k][j]))
                rr.append(v)
            C.append(rr)
        return C
    def eye(n):return [[Q(i==j) for j in range(n)] for i in range(n)]
    def ptr(obj,p):
        for a in p.strip('/').split('/'):obj=obj[int(a)] if isinstance(obj,list) else obj[a]
        return obj
    def ev(path,idx):
        key=(path,idx)
        if key in cache:return cache[key]
        ob=load(path);n=ob['nodes'][idx];kind=n['op'];m,k=n['shape']
        unknown=lambda:[[None for _ in range(k)] for _ in range(m)]
        aa=n.get('args',[])
        if kind.startswith('external_frozen'):
            ref=n['reference'];rp=ref['path'];target=load(rp)
            assert pins[rp]==ref['sha256']
            if rp.endswith('operator.json'):
                j=ref['node_index'] if 'node_index' in ref else ptr(target,ref['json_pointer'])
                ans=ev(rp,j)
            else:ans=unknown()
        elif kind=='identity':ans=eye(m)
        elif kind=='zero':ans=[[Q(0) for _ in range(k)] for _ in range(m)]
        elif kind=='one':ans=[[Q(1)]]
        elif kind=='QQ_constant':ans=[[Q(n['value'])]]
        elif kind=='row_permutation':ans=[[Q(n['order'][i]==j) for j in range(k)] for i in range(m)]
        elif kind=='column_permutation':ans=[[Q(n['order'][j]==i) for j in range(k)] for i in range(m)]
        elif kind=='transpose':ans=[list(r) for r in zip(*ev(path,aa[0]))]
        elif kind=='matrix_product':ans=mm(ev(path,aa[0]),ev(path,aa[1]))
        elif kind in ['difference','sum']:
            A=ev(path,aa[0]);B=ev(path,aa[1]);sgn=-1 if kind=='difference' else 1
            ans=[[plus(a,times(Q(sgn),b)) for a,b in zip(ar,br)] for ar,br in zip(A,B)]
        elif kind=='negative':ans=[[times(Q(-1),a) for a in row] for row in ev(path,aa[0])]
        elif kind=='central_scalar_multiple':
            v=ev(path,aa[0])[0][0];ans=[[times(v,a) for a in row] for row in ev(path,aa[1])]
        elif kind=='submatrix':
            A=ev(path,aa[0]);ans=[[A[i][j] for j in n['columns']] for i in n['rows']]
        elif kind=='block_matrix':
            blocks=n.get('blocks',n.get('block_rows'));ans=[]
            for row in blocks:
                mats=[ev(path,j) for j in row]
                for ii in range(len(mats[0])):ans.append([v for A in mats for v in A[ii]])
        elif kind=='diagonal':
            vals=[ev(path,i)[0][0] for i in n['diagonal']];ans=[[vals[i] if i==j else Q(0) for j in range(k)] for i in range(m)]
        elif kind=='diagonal_known_FULLU_unit_inverse':
            A=ev(path,aa[0]);assert all(A[i][j]==0 for i in range(m) for j in range(k) if i!=j)
            ans=[[((None if A[i][i] is None else 1/A[i][i]) if i==j else Q(0)) for j in range(k)] for i in range(m)]
        elif kind=='horizontal_concatenation':
            mats=[ev(path,i) for i in aa];ans=[[v for A in mats for v in A[i]] for i in range(m)]
        else:ans=unknown() # Every other source scalar/adjugate is opaque.
        assert len(ans)==m and all(len(r)==k for r in ans)
        count[0]+=m*k;assert count[0]+15000<150000
        cache[key]=ans;return ans
    L=ev(pp,op['exports']['fullsource15']);U=ev(pp,op['exports']['fulltarget22'])
    assert all(L[i][14]==0 for i in range(14)) and L[14][14]==1
    assert all(U[i][j]==0 for i in freeidx for j in range(14))
    assert [[U[i][j] for j in range(14,22)] for i in freeidx]==eye(8)
    # Column normalization is diagonal and equals one on all eight free rows.
    Qmat=ev(st(k11),ko['exports']['normalized_column_units'])
    assert all(Qmat[i][j]==0 for i in range(22) for j in range(22) if i!=j)
    assert all(Qmat[i][i]==1 for i in freeidx)
    structural=dict(unknown_never_cancels=True,cached_structural_slots=count[0],
      L_top14_source14_zero=True,L_last_source14_one=True,U_free_rows_first14_zero=True,
      U_free_rows_last8_identity=True,column_normalization_free8_identity=True,
      extracted_inverse_proof='L00*K14*(Qdiag^-1*U)00=I14; square matrices over the commutative parameter ring give a two-sided inverse',
      exact_lastrow_entry_proof='last L coefficient at source14 is1 and free quotient of U isI8; hence scalar Schur complement equals existing residual1x8, not just a unit multiple')
    dump(out/'source-target-flag-proof.json',structural)
    N=[]
    def n(kind,shape,args=(),**kw):
        assert all(0<=i<len(N) for i in args);i=len(N);N.append(dict(id=i,op=kind,shape=shape,args=list(args),**kw));return i
    def ext(path,idx,shape):
        x=load(path);assert x['nodes'][idx]['shape']==shape
        return n('external_frozen_operator_root',shape,reference=dict(path=path,sha256=pins[path],node_index=idx))
    def root(s,name,shape):
        path=st(s);x=load(path);return ext(path,x['exports'][name],shape)
    def subm(a,rr,cc):return n('submatrix',[len(rr),len(cc)],[a],rows=rr,columns=cc)
    def mul(a,b):
        assert N[a]['shape'][1]==N[b]['shape'][0]
        return n('matrix_product',[N[a]['shape'][0],N[b]['shape'][1]],[a,b])
    def diff(a,b):assert N[a]['shape']==N[b]['shape'];return n('difference',N[a]['shape'],[a,b])
    inp=root(k11,'input15x22',[15,22]);diag=root(k11,'normalized_column_units',[22,22]);dinv=n('diagonal_known_FULLU_unit_inverse',[22,22],[diag])
    norm=mul(inp,diag);lroot=root(p14,'fullsource15',[15,15]);uroot=root(p14,'fulltarget22',[22,22]);V=mul(dinv,uroot)
    L00=subm(lroot,list(range(14)),pivrows);V00=subm(V,pividx,list(range(14)));Kinv=mul(V00,L00)
    Sinv=subm(Kinv,list(range(4,14)),list(range(4,14)))
    pair=pivrows[:4];other=pivrows[4:];first=pividx[:4];rest=pividx[4:];j=allcols.index(17)
    A=root(staged,'A4',[4,4]);Ai=root(staged,'A4inverse',[4,4]);B=subm(norm,pair,rest)
    C=subm(norm,other,first);D=subm(norm,other,rest);c=subm(norm,[14],first);d=subm(norm,[14],rest)
    v4=subm(norm,pair,[j]);v10=subm(norm,other,[j]);v1=subm(norm,[14],[j])
    x=mul(Ai,v4);corrected10=diff(v10,mul(C,x));corrected1=diff(v1,mul(c,x));row10=diff(d,mul(mul(c,Ai),B))
    S=diff(D,mul(mul(C,Ai),B));final=diff(corrected1,mul(mul(row10,Sinv),corrected10))
    original=root(p14,'residual1x8',[1,8]);original17=subm(original,[0],[freecols.index(17)])
    operator=dict(nodes=N,exports=dict(pair_A4=A,pair_A4_inverse=Ai,pair_target17=v4,pair_correction_vector=x,
      remaining10_corrected_target17=corrected10,last_source_corrected_target17=corrected1,
      Schur10=S,Schur10_inverse_from_existing_pivots=Sinv,last_source_corrected_pivot10=row10,
      corrected_final_target17=final,original_PIVOT14_target17=original17),
      chart=op['chart'],new_inverse_factors=[],all_four_pair_components_retained=True,
      exact_equality_to_original_lastrow=True,no_naive_M11_G134_12border=True,
      source_pair_rows=pair,source_pivot10=other,source_last=14,pair_targets_original32=pivcols[:4],
      pivot10_targets_original32=pivcols[4:],final_target_original32=17)
    dump(out/'operator.json',operator)
    coeffgate=dict(actual_scalar_job_authorized=False,
      first_finite_target='x=A4_inverse*v17 using four actual pair components; detA4=-T and cached last adjugate row=-kappa',
      literal_pair_source_counts=pr['pair4'],literal_all15_source_counts=pr['all15'],
      no_parameter_coefficient_expansion_in_this_job=True,
      largest_remaining_coefficient_gate='expanded S10_inverse is NOT computed; source/operator access and bounded type construction do not bound expanded rational coefficients',
      required_before_scalar_worker=['choose a literal source-bound contraction schedule','propagate all denominator/content/degree/Newton supports','bound ALL live source/intermediate/output/scratch integers and serialized state'],
      no_empty_bad_locus_or_global57_claim=True)
    dump(out/'future-scalar-gate.json',coeffgate)
    code=Path(__file__);pins[str(code.relative_to(ROOT))]=sha(code);note=ROOT/f'research/notes/{ID}.md';pins[str(note.relative_to(ROOT))]=sha(note)
    total=sum(p.stat().st_size for p in out.iterdir());assert total<1024*1024
    rep=dict(id=ID,all_passed=True,frozen=True,source_only=True,corrected_pair4_plus_pivot10_decomposition=True,
      no_new_parameter_inverse=True,all_three_additional_pair_tail_corrections_retained=True,
      exact_source_target_flag_slots=count[0],typed_nodes=len(N),no_new_scalar_or_determinant_value=True,
      lastrow_scalar_worker_gate_pending=True,total_new_RAW_before_report=total,input_sha256=pins,
      evidence_sha256={p.name:sha(p) for p in out.iterdir()},seconds=time.monotonic()-started)
    dump(out/'report.json',rep);signal.alarm(0);print(json.dumps({k:v for k,v in rep.items() if k not in ['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
