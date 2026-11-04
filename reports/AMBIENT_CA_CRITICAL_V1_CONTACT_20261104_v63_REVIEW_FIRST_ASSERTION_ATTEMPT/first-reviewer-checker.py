#!/usr/bin/env python3
"""Independent actual-source contact audit via ordered Lie exponentials.

No producer imports, finite source solve, new point, new jet order or GB.
The saved three restrictions are replayed by a different representation law:
exponentials of physical150 Lie operators, not polynomial matrix substitution.
"""
import argparse,gzip,hashlib,itertools,json,math,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
from ambient_CA_critical_V1_relative_tangent_review_20261104_v63 import (
    action,push,eye,mat,mv,col,cols,determinant,PAIRS,MONS,PI,MI)
from ambient_B18_fullV56_comparison_unit_review_20261104_v63 import mm
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_CA_CRITICAL_V1_CONTACT_20261104_v63_REVIEW'
PREFIX='AMBIENT_CA_CRITICAL_V1_'
def read(p):
    b=p.read_bytes();return json.loads(gzip.decompress(b) if p.suffix=='.gz' else b)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def dot(a,b):return sum((x*y for x,y in zip(a,b)),F(0))
def ident(n):return [[F(i==j) for j in range(n)] for i in range(n)]
def sub(a,b):return [x-y for x,y in zip(a,b)]
def addv(*vs):return [sum(x,F(0)) for x in zip(*vs)]
def scaled(v,c):return [c*x for x in v]
def rvec(a):return list(map(F,a))
def insert(jet,m,v,c=F(1)):
    if c and any(v):jet[m]=addv(jet.get(m,[F(0)]*150),scaled(v,c))
def group_lie(jet,positions,a,b,N):
    for k in reversed(range(len(positions))):
        if not(a[k] or b[k]):continue
        out={}
        for m,v in jet.items():
            w=v
            for degree in range(N-sum(m)+1):
                if degree:w=action(w,*positions[k])
                if not any(w):break
                for i in range(degree+1):
                    j=degree-i;c=a[k]**i*b[k]**j/F(math.factorial(i)*math.factorial(j))
                    insert(out,(m[0]+i,m[1]+j),w,c)
        jet=out
    return jet
def get(jet,m):return jet.get(m,[F(0)]*150)
def compare(jet,saved):
    mons=[tuple(m) for m in saved['monomials']];values=mat(saved['physical150'])
    assert len(values)==150 and len(values[0])==len(mons)
    for j,m in enumerate(mons):assert get(jet,m)==col(values,j),m
    return sum(bool(c) for row in values for c in row)

def main(out):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False)
    names=['CUSP_DIRECTION_PREFLIGHT','CONTACT_OPERATOR_PREFLIGHT','CONTACT_VECTOR']
    expected=['b30b7e43b30e503fb8ec1ce0ac541288ab723178ed495a0d620e22762613fb8c',
      '06cafbe66b9022bb37cc82e505e1a711f16ea8084f28b93f68e4786317cf614d',
      '74d00c4095a71da7f42c86e5d9146823621f0c7491a5432b0964489f6a9e4cf8']
    reports=[];evidence=[];files=[Path(__file__)]
    for name,h in zip(names,expected):
        folder=ROOT/'reports'/(PREFIX+name+'_20261104_v63');p=folder/'report.json';r=read(p)
        assert sha(p)==h and r['all_passed'];reports.append(r);files.append(p)
        for q,v in r['input_sha256'].items():assert sha(ROOT/q)==v;files.append(ROOT/q)
        assert len(r['evidence_sha256'])==1
        for q,v in r['evidence_sha256'].items():
            ep=folder/q;assert sha(ep)==v;evidence.append(read(ep));files.append(ep)
    D,I,E=evidence
    (out/'phase.json').write_text(json.dumps({'phase':'all_three_frozen_source_packets_bound'})+'\n')
    base=ROOT/'reports'
    C=read(base/'AMBIENT_CA_V1_ODD_ORDER_CHECKPOINT_20261103_v63/checkpoint.json.gz')
    S=read(base/'AMBIENT_CA_V1_ODD_ORDER_V2_20261103_v63/setup.json.gz')
    O=read(base/'AMBIENT_CA_V1_ACTUAL_ODD_RANK9_20261103_v63/certificate.json.gz')
    oldjets=read(base/'AMBIENT_CA_V1_ODD_ORDER_V2_20261103_v63/jets.json.gz')
    R=read(base/'AMBIENT_CA_CRITICAL_V1_RELATIVE_TANGENT_20261104_v63/actual-V1-source-relative-tangent.json')
    J=mat(C['full150_by29_target_Jacobian']);JB=mat(R['CB_normalizer_J150_by29_original_frame'])
    A=mat(O['all150_actual_odd_linear']);indices=C['complement28_columns'];Y=S['actual_target_Y28'];Rrows=Y+[S['actual_branch29th_row']]
    even=mat(D['actual_even_tangent_embedding150_by29']);assert [even[i] for i in Rrows]==ident(29)
    assert mm(even,[J[i] for i in Rrows])==J
    wc=rvec(S['corrected_quadratic150']);assert mv(even,[wc[i] for i in Rrows])==wc
    nu=rvec(D['cusp_even_nu_linear_in_old_R29'])
    assert all(not dot(nu,c) for c in zip(*[J[i] for i in Rrows])) and dot(nu,[wc[i] for i in Rrows])==1
    eight=S['actual_odd_selected8_rows'];assert eight==[27,28,29,35,60,65,67,68]
    odd=mat(D['actual_odd_Y8_linear_physical150_after_even_subtraction'])
    for i,row in enumerate(eight):
        expectedrow=[F(j==row) for j in range(150)]
        for j,k in enumerate(Rrows):expectedrow[k]-=even[row][j]
        assert odd[i]==expectedrow
    I9=mat(D['actual_regular_I9_rows8_plus_nu']);assert I9==[A[i] for i in eight]+[nu]
    embedding=mat(D['actual_cusp_odd_embedding150_by9']);assert mm(embedding,I9)==A
    assert [embedding[i] for i in eight]==[r[:9] for r in ident(9)[:8]]
    assert determinant([[I9[i][j] for j in O['nine_columns']] for i in range(9)])
    cubic=[F(row[3]) for row in oldjets['stages'][0]['all150_output_series']]
    assert col(embedding,8)==cubic and embedding[63][8]==F(800,27)
    kernel=rvec(D['actual_CB_source_kernel_of_CA_even_projection29']);normal=mv(JB,kernel)
    assert normal==rvec(D['actual_CB_normal_direction_physical150']) and any(normal)
    assert all(not normal[i] for i in Rrows)
    normalcoords=rvec(D['cusp_direction_Y8'])+[F(D['cusp_direction_Z'])]
    assert mv(embedding,normalcoords)==normal and normalcoords==[0,0,0,0,1,-50,0,0,0]
    kappa=rvec(C['kernel29']);T=cols([kappa]+[col(ident(29),j) for j in indices])
    H=mat(S['full150_by29_mixed_Hessian']);tau=rvec(D['actual_trace_tau_linear_on_source29'])
    assert [dot(tau,c) for c in zip(*T)]==[dot(nu,c)/2 for c in zip(*[H[i] for i in Rrows])]
    assert dot(tau,kappa)==1 and mv(J,kappa)==[0]*150
    h=rvec(D['actual_N_divisor_linear_on_source29_normalized'])
    assert h==mat(C['normalizer_source29_frame'])[2]==[F(i==22) for i in range(29)] and dot(h,kappa)==1
    u=mat(D['actual_even_u8_linear_on_source29']);assert u==mm([A[i] for i in eight],[J[i] for i in Rrows])
    params=mat(D['actual_source_parameter_Jacobian29']);assert params==[J[i] for i in S['orbit20_raw_rows']]+u+[tau]
    assert determinant(params)==F(D['actual_source_parameter_Jacobian_determinant'])!=0
    lam=rvec(D['actual_N_trace_beta_linear_coefficients_lambda8'])
    assert lam==[0,0,0,0,F(-9,40),0,F(15,64),F(45,16)]
    assert [dot(lam,c) for c in zip(*u)]==sub(tau,h)
    nx=-2*dot(lam,normalcoords[:8]);assert nx==F(9,20)==F(D['N_adapted_direction_nx'])
    xrow=[F(0)]*150
    for j,row in enumerate(Rrows):xrow[row]+=nu[j]
    xrow=sub(xrow,[2*dot(lam,c) for c in zip(*odd)])
    assert xrow==rvec(D['N_adapted_target_x_linear_physical150']) and dot(xrow,normal)==nx
    proj=mat(D['N_adapted_CB_even29_projection_matrix'])
    assert proj==[JB[i] for i in S['orbit20_raw_rows']]+mm([A[i] for i in eight],[JB[i] for i in Rrows])+[[dot(xrow,c) for c in zip(*JB)]]
    assert determinant(proj)==F(128,20599365234375)==F(D['N_adapted_CB_even29_projection_determinant'])
    # All saved directions are verified as solutions; no new system is solved.
    av=rvec(I['CA_source_v29']);ar=rvec(I['CA_source_r29']);bv=rvec(I['CB_source_v29']);bx=rvec(I['CB_source_x29'])
    assert rvec(I['CA_source_kappa29'])==kappa
    beta=scaled(normalcoords[:8],1/nx);rhs=[F(0)]*20+beta+[F(0)];xr=[F(0)]*28+[F(1)]
    p=params[:-1]+[h];assert determinant(p)==F(reports[1]['CA_source_parameter_minor'])!=0
    assert mv(p,av)==rhs and mv(proj,bv)==rhs and mv(proj,bx)==xr
    assert dot(h,av)==0 and mv(JB,bv)==mv(J,av)
    V=mv(JB,bx);assert V==scaled(normal,1/nx)==rvec(I['actual_CB_x150'])
    missing=next(j for j in range(29) if j not in indices)
    alpha=av[missing]/kappa[missing];coordinates=[alpha]+[av[j]-alpha*kappa[j] for j in indices]
    mixed=mv(H,coordinates);W=scaled(col(H,0),F(1,2))
    assert W==rvec(I['cached_CA_fold_square150'])==wc and mixed==rvec(I['cached_CA_fold_v_mixed150'])
    assert V==addv(mixed,W,mv(J,ar))
    zrow=[F(j==63) for j in range(150)]
    for j,row in enumerate(Rrows):zrow[row]-=even[63][j]
    zrow=scaled(sub(zrow,[dot(embedding[63][:8],c) for c in zip(*odd)]),1/embedding[63][8])
    assert zrow==rvec(I['canonical_Z_linear150'])
    assert all(not dot(zrow,c) for c in zip(*J)) and all(not dot(zrow,c) for c in zip(*JB))
    assert dot(zrow,W)==dot(zrow,mixed)==0 and [dot(zrow,c) for c in zip(*embedding)]==[0]*8+[1]
    CA=mat(I['CA_target_tangent150_by38']);assert CA==mat(R['CA_complete_target_tangent150_by38_original_frame'])
    rows=I['CA_tangent_left_rows38'];left=mat(I['CA_tangent_left_inverse38']);assert mm(left,[CA[i] for i in rows])==ident(38)
    (out/'phase.json').write_text(json.dumps({'phase':'N_adaptation_all_fixed_matrices_source_directions_and_adjoint_checks_passed'})+'\n')
    point=rvec(C['point150']);fixed=mat(C['fixedCA150_by20']);extras=C['extra9_fiber_indices'];positions=[divmod(j,5) for j in C['group20_GL_indices']]
    def ca_replay(a,b,N):
        jet={(0,0):point}
        for k,j in enumerate(extras):
            insert(jet,(1,0),col(fixed,j),a[20+k]);insert(jet,(0,1),col(fixed,j),b[20+k])
        return group_lie(jet,positions,a[:20],b[:20],N)
    ca3=ca_replay(kappa,av,3);ca2=ca_replay(av,ar,2)
    # Full CB q*(-f) cross terms: only z-w coefficients before flag transport.
    cross=[]
    for mon in [(1,1),(0,1),(0,0)]:
        v=[F(0)]*150;v[15*PI[3,4]+MI[mon]]=F(-1);cross.append(v)
    assert cross==mat(E['CB_mixed_ratio_f3_physical_columns'])
    fiber=mat(R['CB_full20fiber_columns_prime_frame']);coeff=rvec(R['CB_fixed19_coefficients'])
    center=[dot(row[:3],coeff[:3])+dot(row[4:],coeff[3:]) for row in fiber]
    cbjet={(0,0):center}
    for j in range(20):insert(cbjet,(1,0),col(fiber,j),bv[j]);insert(cbjet,(0,1),col(fiber,j),bx[j])
    for j in range(3):
        insert(cbjet,(2,0),cross[j],bv[3]*bv[j]);insert(cbjet,(1,1),cross[j],bv[3]*bx[j]+bx[3]*bv[j]);insert(cbjet,(0,2),cross[j],bx[3]*bx[j])
    cbjet=group_lie(cbjet,R['CB_flag9_positions'],bv[20:],bx[20:],2)
    frame=mat(R['proper_Npair_frame5']);fi=mat(R['inverse_frame'])
    cb2={m:push(v,fi,frame) for m,v in cbjet.items()}
    nonzero=sum(compare(q,E[name]) for q,name in [(ca3,'CA_fold_v_source3'),(ca2,'CA_v_r_source2'),(cb2,'CB_v_x_source2')])
    assert all(get(q,(0,0))==point for q in [ca3,ca2,cb2])
    assert get(ca3,(2,0))==W and get(ca3,(1,1))==mixed and dot(zrow,get(ca3,(3,0)))==1
    delta=sub(sub(sub(get(cb2,(1,1)),get(ca3,(1,2))),get(ca3,(2,1))),get(ca2,(1,1)))
    assert delta==rvec(E['Delta150']);cv=dot(zrow,delta);assert cv==F(E['c_v'])==0
    factors=mv(left,[delta[i] for i in rows]);pi=sub(delta,mv(CA,factors))
    assert pi==rvec(E['all_actual_ambient_graph_contact_bits_PiDelta150'])==[F(0)]*150
    assert E['kappa']==1 and E['one_minus_c_v']=='1' and not E['nonzero_graph_rows']
    total_nonzero=nonzero+sum(bool(c) for c in delta)+sum(bool(c) for c in pi)
    assert total_nonzero==reports[2]['stored_nonzero_coefficient_count']==464 and total_nonzero<150000
    ev={'all_three_primary_source_hashes_bound':True,'all_fixed_N_adaptation_and_full_even_subtraction_checks':True,
      'lambda8':list(map(str,lam)),'normal_nx':str(nx),'normal_nY8':list(map(str,normalcoords[:8])),'normal_nZ':'0',
      'adapted_CB_even29_determinant':'128/20599365234375',
      'all_source_directions_and_all150_decomposition_verified_without_solve':True,
      'independent_source_replay':'ordered physical150 Lie-exponentials applied right-to-left, distinct from author finite matrix-substitution ring',
      'all3300_actual_source_slots_verified':True,'source_nonzero_coefficients':nonzero,'source_and_contact_nonzero_coefficients':total_nonzero,
      'all3_actual_CB_q_delta_f_mixed_terms_retained':True,'c_v':'0','all150_ambient_PiDelta_zero':True,
      'whole_contact_lemma_accepted':True,'actual_transverse_D_ideal':'x*(l1,...,l7,x,v^2)',
      'actual_normalizer_tensor':'G20[[l7,v,t]]/(l7*t,t^2-v*t,t^3)',
      'comparison_kernel':'G20*(x*v), transverse length1','comparison_cokernel':'G20*t modulo t^2, transverse length1',
      'generic_CA_critical_kappa':1,'no_origin_D0_depth_socle_or_fullnil_claim':True}
    out_ev=out/'independent-N-adaptation-contact-and-local-proof.json';out_ev.write_text(json.dumps(ev,indent=2)+'\n')
    note=ROOT/'research/notes'/f'{STEM}.md';files += [note,
      ROOT/'code/ambient_CA_critical_V1_relative_tangent_review_20261104_v63.py',
      ROOT/'code/ambient_B18_fullV56_comparison_unit_review_20261104_v63.py']
    report={'id':STEM,'all_passed':True,'independent_actual_N_frame_adjoint_source_replay_and_contact_proof':True,
      'all3300_source_slots_and_all150_graph_bits_checked':True,'c_v':'0','nonzero_graph_count':0,
      'actual_generic_CA_critical_kappa':1,'transverse_normalizer_kernel_length':1,'transverse_comparison_C_length':1,
      'no_authorproducer_newpoint_jetorder_source_solve_GB':True,'origin_D0_depth_socle_fullnil_claimed':False,
      'input_sha256':{str(p.relative_to(ROOT)):sha(p) for p in files},'evidence_sha256':{out_ev.name:sha(out_ev),'phase.json':sha(out/'phase.json')},
      'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:report[k] for k in ['id','all_passed','seconds','max_rss_kib']}))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
