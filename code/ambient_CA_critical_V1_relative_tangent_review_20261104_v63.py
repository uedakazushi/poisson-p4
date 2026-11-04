#!/usr/bin/env python3
"""Independent exact physical150 source and fixed-rank-certificate audit at V1.
No author imports, new point, source jet, implicit solve, rank search or GB.
"""
import argparse,gzip,hashlib,itertools,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
PRIMARY='AMBIENT_CA_CRITICAL_V1_RELATIVE_TANGENT_20261104_v63'
STEM=PRIMARY+'_REVIEW'
PAIRS=list(itertools.combinations(range(5),2));MONS=list(itertools.combinations_with_replacement(range(5),2))
PI={p:i for i,p in enumerate(PAIRS)};MI={p:i for i,p in enumerate(MONS)}
def read(p):
    b=p.read_bytes();return json.loads(gzip.decompress(b) if p.suffix=='.gz' else b)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def mat(M):return [[F(v) for v in row] for row in M]
def transpose(M):return [list(c) for c in zip(*M)]
def col(M,j):return [r[j] for r in M]
def cols(cs):return transpose(cs)
def mv(M,v):return [sum((c*v[j] for j,c in enumerate(row) if c),F(0)) for row in M]
def eye(n):return [[F(i==j) for j in range(n)] for i in range(n)]
def determinant(M):
    a=[row[:] for row in M];answer=F(1)
    for k in range(len(a)):
        p=next((i for i in range(k,len(a)) if a[i][k]),None);assert p is not None
        if p!=k:a[k],a[p]=a[p],a[k];answer=-answer
        x=a[k][k];answer*=x
        for i in range(k+1,len(a)):
            if not a[i][k]:continue
            q=a[i][k]/x
            for j in range(k+1,len(a)):a[i][j]-=q*a[k][j]
            a[i][k]=F(0)
    return answer
def push(v,G,Gi):
    """Literal coefficient substitution in G P(G^-1 X) G^T, no CAS."""
    result=[F(0)]*150
    for k,coef in enumerate(v):
        if not coef:continue
        a,b=PAIRS[k//15];u,w=MONS[k%15];poly={}
        for i in range(5):
            for j in range(5):
                if Gi[u][i] and Gi[w][j]:
                    m=tuple(sorted((i,j)));poly[m]=poly.get(m,F(0))+coef*Gi[u][i]*Gi[w][j]
        for i,j in PAIRS:
            c=G[i][a]*G[j][b]-G[i][b]*G[j][a]
            if c:
                for m,x in poly.items():result[15*PI[i,j]+MI[m]]+=c*x
    return result
def entry(v,a,b):
    if a==b:return [F(0)]*15
    sign=1 if a<b else -1;a,b=sorted((a,b))
    return [sign*x for x in v[15*PI[a,b]:15*PI[a,b]+15]]
def action(v,i,j):
    """Coefficient of t in (I+t Eij).P, retaining both tensor and variable terms."""
    ans=[F(0)]*150
    for a,b in PAIRS:
        offset=15*PI[a,b]
        if a==i:
            for k,x in enumerate(entry(v,j,b)):ans[offset+k]+=x
        if b==i:
            for k,x in enumerate(entry(v,a,j)):ans[offset+k]+=x
        for k,x in enumerate(v[offset:offset+15]):
            if not x:continue
            m=list(MONS[k]);count=m.count(i)
            if count:
                m.remove(i);m.append(j);ans[offset+MI[tuple(sorted(m))]]-=count*x
    return ans
def fixed_rank(M,rec):
    r=rec['rank'];ci=rec['independent_columns'];ri=rec['minor_rows']
    assert len(ci)==len(ri)==r and len(set(ci))==len(set(ri))==r
    minor=[[M[i][j] for j in ci] for i in ri];det=determinant(minor)
    assert det==F(rec['minor_determinant']) and det
    factor=mat(rec['all_column_factorization']);assert len(factor)==r
    B=[[row[j] for j in ci] for row in M]
    for i,row in enumerate(B):
        value=[F(0)]*len(M[0])
        for k,c in enumerate(row):
            if c:
                for j,x in enumerate(factor[k]):
                    if x:value[j]+=c*x
        assert value==M[i]
    return {'rank':r,'fixed_minor_determinant':str(det),'all_column_factorization_verified':True}

def main(out):
    start=time.monotonic();signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False)
    folder=ROOT/'reports'/PRIMARY;primary=folder/'report.json';A=read(primary)
    assert sha(primary)=='6e6245a6b2115c0b938d466ce0f00c1adac532f501c3e384ad41d352d22c9558' and A['all_passed']
    for p,h in A['input_sha256'].items():assert sha(ROOT/p)==h,p
    for p,h in A['evidence_sha256'].items():assert sha(folder/p)==h,p
    E=read(folder/'actual-V1-source-relative-tangent.json')
    (out/'phase.json').write_text(json.dumps({'phase':'all_frozen_inputs_hashed','primary_sha256':sha(primary)})+'\n')
    C=read(ROOT/'reports/AMBIENT_CA_V1_ODD_ORDER_CHECKPOINT_20261103_v63/checkpoint.json.gz')
    S=read(ROOT/'reports/AMBIENT_CA_V1_ODD_ORDER_V2_20261103_v63/setup.json.gz')
    B=read(ROOT/'reports/AMBIENT_CA_V1_ACTUAL_ODD_RANK9_20261103_v63/certificate.json.gz')
    P=read(ROOT/'reports/AMBIENT_B18_PROPER_POLYNOMIAL_MODEL_PREFLIGHT_20261104_v63/literal-proper-source14-preflight.json')
    CB=read(ROOT/'reports/AMBIENT_B18_CB_CRITICAL_SCHEME_20261103_v63/certificate.json.gz')
    point=list(map(F,C['point150']));assert len(point)==150
    G=eye(5);Gi=eye(5);G[1][0]=F(-1,10);Gi[1][0]=F(1,10)
    assert determinant(G)==determinant(Gi)==1
    assert G==mat(E['proper_Npair_frame5']) and Gi==mat(E['inverse_frame'])
    prime=push(point,G,Gi);assert push(prime,Gi,G)==point
    fixed=mat(P['all150_CB_fixed_flag19']);rows=P['zero_coordinate_rows36'];br=P['CB_constant_left_inverse_rows19'];inverse=mat(P['CB_left_inverse19'])
    selected=[[fixed[rows[i]][j] for j in range(19)] for i in br]
    assert all(mv(inverse,col(selected,j))==[F(i==j) for i in range(19)] for j in range(19))
    coeff=mv(inverse,[prime[rows[i]] for i in br])
    assert mv(fixed,coeff)==prime and coeff==list(map(F,E['CB_fixed19_coefficients']))
    assert coeff[:3]==list(map(F,E['CB_binary_f_coefficients']))==[F(-1),F(-1,5),F(-1,100)]
    # The original CB chart is f(h,x); fixed19 is AFTER h/x swap. The ratio
    # derivative g=-f has only the z-w component and zero divergence.
    q=[F(0)]*150
    for (u,v),c in zip([(1,1),(0,1),(0,0)],coeff[:3]):q[15*PI[3,4]+MI[u,v]]=-c
    fiber=cols([*transpose(fixed)[:3],q,*transpose(fixed)[3:]])
    assert fiber==mat(E['CB_full20fiber_columns_prime_frame'])
    flags=CB['flag_positions'];assert flags==E['CB_flag9_positions'] and len(flags)==9
    Jprime=transpose(fiber)+[action(prime,i,j) for i,j in flags]
    JCB=cols([push(v,Gi,G) for v in Jprime])
    assert JCB==mat(E['CB_normalizer_J150_by29_original_frame'])
    # All three f coefficients form the same square; this is not a copied
    # CA ratio tangent or a newly chosen point.
    assert coeff[1]**2-4*coeff[0]*coeff[2]==0
    JCA=mat(C['full150_by29_target_Jacobian']);indices=C['complement28_columns']
    wc=list(map(F,S['corrected_quadratic150']));odd=mat(B['all150_actual_odd_linear']);oc=B['nine_columns']
    RCA=S['actual_target_Y28']+[S['actual_branch29th_row']]
    assert len(indices)==28 and len(oc)==9 and len(RCA)==29
    assert all(not x for i in RCA for x in odd[i])
    ca=cols([*[col(JCA,j) for j in indices],wc,*[col(odd,j) for j in oc]])
    assert ca==mat(E['CA_complete_target_tangent150_by38_original_frame'])
    assert E['CA_tangent_columns_rule']=={'source28_indices':indices,'even_branch_raw_row':22,'actual_odd9_base_columns':oc}
    assert E['CA_even29_actual_raw_rows']==RCA
    proj=[JCB[i] for i in RCA];assert proj==mat(E['CB_to_CA_even_projection29'])
    union=[x+y for x,y in zip(ca,JCB)]
    (out/'phase.json').write_text(json.dumps({'phase':'all150_source_frame_and_complete_tangent_matrices_verified'})+'\n')
    matrices={'CA38_certificate':ca,'CB29_certificate':JCB,'actual_union_certificate':union,'CB_to_CA_even_projection_certificate':proj}
    cert={k:fixed_rank(M,E[k]) for k,M in matrices.items()}
    assert [cert[k]['rank'] for k in matrices]==[38,29,38,28]
    # The Ac8 complete-image theorem supplies the target tangent completion.
    # The normalizer rank28 alone does not suffice; Wc and odd9 are essential.
    gates=[ROOT/p for p in A['input_sha256'] if p.startswith('reports/') and ('REVIEW' in p or 'GENERIC_V_FULL_COMPLETION' in p)]
    for p in gates:assert read(p)['all_passed'],p
    proper_gate=ROOT/'reports/AMBIENT_B18_PROPER_POLYNOMIAL_MODEL_20261104_v63_REVIEW/report.json'
    assert read(proper_gate)['all_passed']
    assert A['actual_D_tangent_dimension']==E['actual_tangent_intersection_dimension']==38+29-38==29
    assert not A['generic_CAcritical_kappa_zero_proved'] and not A['actual_D_and_N_smooth28_proved']
    assert not A['CB_to_CA_even_projection_unit']
    termcount=sum(bool(c) for M in matrices.values() for row in M for c in row)
    termcount+=sum(bool(F(c)) for k in matrices for row in E[k]['all_column_factorization'] for c in row)
    assert termcount<150000
    ev={'own_all150_CB_fixed19_reconstruction':True,'own_inverse_shear_returns_old_V1':True,
      'frame_determinant':1,'CB_f3':list(map(str,coeff[:3])),'CB_f_discriminant':'0',
      'own_actual_CB_ratio_derivative_minus_f_G':True,'own_all9_physical_Lie_flag_directions':True,
      'own_all150_by29_source_differential_and_frame_pullback':True,
      'CA_complete_tangent_source28_corrected_p2_fold_odd9':True,
      'own_odd_zero_on29_even_rows':True,'fixed_rank_receipts':cert,
      'actual_D_tangent_dimension':29,'CB_even_projection_rank':28,
      'generic_kappa_or_K_rank_D0_depth_not_determined':True,
      'initial_tuple_reader_error_preserved_before_actual_source_rank_work':True,
      'stored_nonzero_source_and_rank_factor_coefficients':termcount}
    ep=out/'independent-V1-source-rank-checks.json';ep.write_text(json.dumps(ev,indent=2)+'\n')
    note=ROOT/'research/notes'/f'{STEM}.md';files=[primary,Path(__file__),note,proper_gate]
    files += [ROOT/p for p in A['input_sha256']]+[folder/p for p in A['evidence_sha256']]
    report={'id':STEM,'all_passed':True,'primary_report_sha256':sha(primary),
      'independent_physical150_coefficient_source_frame_and_fixed_rank_certificates':True,
      'CA_complete_target_tangent_rank':38,'CB_complete_normalizer_tangent_rank':29,
      'union_rank':38,'actual_D_tangent_dimension':29,'CB_to_CA_even29_projection_rank':28,
      'generic_kappa_K_rank_D0_depth_or_higher_order_claimed':False,
      'no_author_producer_new_point_root_jet_IFT_GB_or_rank_search':True,
      'input_sha256':{str(p.relative_to(ROOT)):sha(p) for p in files},
      'evidence_sha256':{ep.name:sha(ep),'phase.json':sha(out/'phase.json')},
      'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:report[k] for k in ['id','all_passed','seconds','max_rss_kib']}))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
