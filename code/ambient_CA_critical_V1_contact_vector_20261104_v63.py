#!/usr/bin/env python3
"""ONE bounded exact relative contact job at the old accepted V1 point."""
import argparse,gzip,hashlib,itertools,json,math,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
import sympy as s
from ambient_CA_global_flag_model_20261103_v63 import X,canonical,vector
from ambient_CB_global_flag_model_20261103_v63 import raw as rawB
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_CA_CRITICAL_V1_CONTACT_VECTOR_20261104_v63'
PAIRS=list(itertools.combinations(range(5),2));MONS=list(itertools.combinations_with_replacement(range(5),2))
def read(p):return json.loads(gzip.decompress(p.read_bytes()) if p.suffix=='.gz' else p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def q(v):return Q(str(v))
def jsonwrite(p,x):p.write_text(json.dumps(x,indent=2)+'\n')
class Ring:
    def __init__(self,N):
        self.N=N;self.mons=[(i,d-i) for d in range(N+1) for i in range(d,-1,-1)];self.idx={m:i for i,m in enumerate(self.mons)};self.L=len(self.mons)
        self.pairs=[(i,j,self.idx[(a+c,b+d)]) for i,(a,b) in enumerate(self.mons) for j,(c,d) in enumerate(self.mons) if a+b+c+d<=N]
        self.mul_calls=0;self.scalar_products=0
    def zero(self):return [Q(0)]*self.L
    def const(self,a):return [q(a)]+[Q(0)]*(self.L-1)
    def linear(self,a,b):
        r=self.zero();r[self.idx[1,0]]=q(a);r[self.idx[0,1]]=q(b);return r
    def add(self,a,b):return [x+y for x,y in zip(a,b)]
    def neg(self,a):return [-x for x in a]
    def scale(self,a,c):return [x*c for x in a]
    def mul(self,a,b):
        self.mul_calls+=1;r=self.zero()
        for i,j,k in self.pairs:
            if a[i] and b[j]:r[k]+=a[i]*b[j];self.scalar_products+=1
        return r
    def exp(self,a):
        assert not a[0];r=self.const(1);power=self.const(1)
        for k in range(1,self.N+1):power=self.mul(power,a);r=self.add(r,self.scale(power,Q(1,math.factorial(k))))
        return r
    def ident(self):return [[self.const(int(i==j)) for j in range(5)] for i in range(5)]
    def mm(self,A,B):return [[self.sum(self.mul(A[i][k],B[k][j]) for k in range(5)) for j in range(5)]for i in range(5)]
    def sum(self,values):
        r=self.zero()
        for a in values:r=self.add(r,a)
        return r
def group(r,parameters,positions,inverse=False):
    G=r.ident();order=list(range(len(positions)));order=order[::-1] if inverse else order
    for k in order:
        i,j=positions[k];v=parameters[k]
        if inverse:v=r.neg(v)
        if i==j:v=r.add(r.exp(v),r.const(-1))
        col=[r.mul(G[t][i],v) for t in range(5)]
        for t in range(5):G[t][j]=r.add(G[t][j],col[t])
    return G
def push(r,P,G,Gi):
    forms=[]
    for a,b in MONS:
        row=[]
        for c,d in MONS:
            z=r.mul(Gi[a][c],Gi[b][d])
            if c!=d:z=r.add(z,r.mul(Gi[a][d],Gi[b][c]))
            row.append(z)
        forms.append(row)
    middle=[]
    for k in range(10):
        row=[r.zero() for _ in range(15)]
        for m in range(15):
            if any(P[15*k+m]):
                for n in range(15):
                    if any(forms[m][n]):row[n]=r.add(row[n],r.mul(P[15*k+m],forms[m][n]))
        middle.append(row)
    out=[r.zero() for _ in range(150)]
    for k,(i,j) in enumerate(PAIRS):
        for l,(a,b) in enumerate(PAIRS):
            c=r.add(r.mul(G[i][a],G[j][b]),r.neg(r.mul(G[i][b],G[j][a])))
            if any(c):
                for m in range(15):
                    if any(middle[l][m]):out[15*k+m]=r.add(out[15*k+m],r.mul(c,middle[l][m]))
    return out
def ca_eval(r,C,d1,d2):
    t=[r.linear(a,b) for a,b in zip(d1,d2)];P=[r.const(x) for x in C['point150']];basis=C['fixedCA150_by20'];extras=C['extra9_fiber_indices']
    for k,bj in enumerate(extras):
        for i in range(150):
            z=q(basis[i][bj])
            if z:P[i]=r.add(P[i],r.scale(t[20+k],z))
    pos=[divmod(i,5) for i in C['group20_GL_indices']]
    return push(r,P,group(r,t[:20],pos),group(r,t[:20],pos,True))
def cb_eval(r,R,d1,d2,cross):
    t=[r.linear(a,b) for a,b in zip(d1,d2)];coeff=[q(x) for x in R['CB_fixed19_coefficients']]
    F=[[q(x) for x in row]for row in R['CB_full20fiber_columns_prime_frame']]
    # Center in prime frame = fixed19*f19; q center is zero.
    p0=[sum((row[i]*coeff[i] for i in range(3)),Q(0))+sum((row[4+i]*coeff[3+i]for i in range(16)),Q(0))for row in F]
    P=[r.const(x)for x in p0]
    for i in range(150):
        for j in range(20):
            if F[i][j]:P[i]=r.add(P[i],r.scale(t[j],F[i][j]))
        for j in range(3):
            if cross[j][i]:P[i]=r.add(P[i],r.scale(r.mul(t[3],t[j]),cross[j][i]))
    G=group(r,t[20:],R['CB_flag9_positions']);Gi=group(r,t[20:],R['CB_flag9_positions'],True)
    frame=s.Matrix([[s.Rational(x)for x in row]for row in R['proper_Npair_frame5']]);fi=frame.inv()
    A=[[r.const(fi[i,j])for j in range(5)]for i in range(5)];Ai=[[r.const(frame[i,j])for j in range(5)]for i in range(5)]
    return push(r,P,r.mm(A,G),r.mm(Gi,Ai))
def col(out,r,mon):return [a[r.idx[mon]] for a in out]
def dot(a,b):return sum((q(x)*q(y)for x,y in zip(a,b)),Q(0))
def serial(out):return [[str(v)for v in row]for row in out]
def main(out):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));out.mkdir(parents=True,exist_ok=False)
    files={
      'preflight':'reports/AMBIENT_CA_CRITICAL_V1_CONTACT_OPERATOR_PREFLIGHT_20261104_v63/literal-V1-contact-vector-input.json',
      'checkpoint':'reports/AMBIENT_CA_V1_ODD_ORDER_CHECKPOINT_20261103_v63/checkpoint.json.gz',
      'relative':'reports/AMBIENT_CA_CRITICAL_V1_RELATIVE_TANGENT_20261104_v63/actual-V1-source-relative-tangent.json'}
    I,C,R=map(read,[ROOT/p for p in files.values()]);assert I['source_coefficient_evaluation_executed'] is False
    phase={'id':STEM,'status':'early_phase_inputs_bound','input_sha256':{v:sha(ROOT/v)for v in files.values()},'no_source_solve_or_new_point':True}
    jsonwrite(out/'phase.json',phase)
    # These are the full actual q*delta_f columns, including canonicalization
    # and the SAME fixed h/x swap used by the frozen proper CB source.
    h,x,y,z,w=X;sw=s.eye(5);sw[0,0]=sw[1,1]=0;sw[0,1]=sw[1,0]=1
    cross=[]
    for f in [h*h,h*x,x*x]:
        M=canonical(rawB(0,-f,0));M=sw*M.subs({h:x,x:h},simultaneous=True)*sw.T
        cross.append([q(v)for v in vector(M)])
    F=R['CB_full20fiber_columns_prime_frame'];fc=list(map(q,R['CB_fixed19_coefficients'][:3]))
    assert all(sum((fc[j]*cross[j][i]for j in range(3)),Q(0))==q(F[i][3])for i in range(150))
    d=lambda k:list(map(q,I[k]))
    r3=Ring(3);r2=Ring(2);rb=Ring(2)
    ca3=ca_eval(r3,C,d('CA_source_kappa29'),d('CA_source_v29'))
    jsonwrite(out/'CA-fold-v-source3.json',{'monomials':r3.mons,'physical150':serial(ca3)})
    phase['status']='CA_10_monomials_complete';jsonwrite(out/'phase.json',phase)
    ca2=ca_eval(r2,C,d('CA_source_v29'),d('CA_source_r29'))
    cb2=cb_eval(rb,R,d('CB_source_v29'),d('CB_source_x29'),cross)
    point=list(map(q,C['point150']));J=C['full150_by29_target_Jacobian'];JB=R['CB_normalizer_J150_by29_original_frame']
    assert col(ca3,r3,(0,0))==col(ca2,r2,(0,0))==col(cb2,rb,(0,0))==point
    for values,ring,jac,d1,d2 in [(ca3,r3,J,d('CA_source_kappa29'),d('CA_source_v29')),(ca2,r2,J,d('CA_source_v29'),d('CA_source_r29')),(cb2,rb,JB,d('CB_source_v29'),d('CB_source_x29'))]:
        assert col(values,ring,(1,0))==[dot(row,d1)for row in jac]
        assert col(values,ring,(0,1))==[dot(row,d2)for row in jac]
    assert col(ca3,r3,(2,0))==list(map(q,I['cached_CA_fold_square150']))
    assert col(ca3,r3,(1,1))==list(map(q,I['cached_CA_fold_v_mixed150']))
    zlinear=list(map(q,I['canonical_Z_linear150']));assert dot(zlinear,col(ca3,r3,(3,0)))==1
    C12=col(ca3,r3,(1,2));C21=col(ca3,r3,(2,1));H=col(ca2,r2,(1,1));HB=col(cb2,rb,(1,1))
    delta=[b-c-d-e for b,c,d,e in zip(HB,C12,C21,H)]
    cv=dot(zlinear,delta);rows=I['CA_tangent_left_rows38'];left=[[q(x)for x in row]for row in I['CA_tangent_left_inverse38']]
    coeff=[dot(row,[delta[i]for i in rows])for row in left]
    pi=[v-dot(row,coeff)for v,row in zip(delta,I['CA_target_tangent150_by38'])]
    nonzero=[i for i,v in enumerate(pi)if v]
    result={'CA_fold_v_source3':{'monomials':r3.mons,'physical150':serial(ca3)},
      'CA_v_r_source2':{'monomials':r2.mons,'physical150':serial(ca2)},
      'CB_v_x_source2':{'monomials':rb.mons,'physical150':serial(cb2)},
      'CB_mixed_ratio_f3_physical_columns':serial(cross),
      'Delta150':list(map(str,delta)),'c_v':str(cv),'one_minus_c_v':str(1-cv),
      'all_actual_ambient_graph_contact_bits_PiDelta150':list(map(str,pi)),
      'nonzero_graph_rows':nonzero,'first_nonzero_graph_bit':None if not nonzero else [nonzero[0],str(pi[nonzero[0]])],
      'all_source0_1_and_CA_cached_mixed2_verified':True,
      'source_monomial_slots':3300,'scalar_multiplications':sum(r.scalar_products for r in [r3,r2,rb]),
      'ring_product_calls':sum(r.mul_calls for r in [r3,r2,rb]),
      'conditional_lemma_gate_one_minus_c_v_unit':bool(1-cv),
      'kappa':None if cv==1 else int(not nonzero),
      'stops_if_degenerate_no_extra_order':True}
    ev=out/'actual-V1-contact-vector-and-source-jets.json';jsonwrite(ev,result)
    totalterms=sum(bool(c)for arr in [ca3,ca2,cb2]for row in arr for c in row)+sum(bool(x)for x in delta)+sum(bool(x)for x in pi)
    assert totalterms<150000 and sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
    deps=[ROOT/'reports/AMBIENT_CA_CRITICAL_V1_CONTACT_OPERATOR_PREFLIGHT_20261104_v63/report.json']
    assert all(read(p)['all_passed']for p in deps)
    paths=[Path(__file__),ROOT/'research/notes'/f'{STEM}.md']+[ROOT/p for p in files.values()]+deps+[ROOT/'code/ambient_CA_global_flag_model_20261103_v63.py',ROOT/'code/ambient_CB_global_flag_model_20261103_v63.py']
    report={'id':STEM,'all_passed':True,'status':'coefficient_gate_degenerate_inconclusive' if cv==1 else 'conditional_contact_coefficients_nondegenerate',
      'c_v':str(cv),'one_minus_c_v':str(1-cv),'nonzero_graph_count':len(nonzero),'first_nonzero_graph_bit':result['first_nonzero_graph_bit'],
      'kappa_by_preflight_lemma':result['kappa'],'independent_review_required':True,
      'no_source_solve_new_point_or_extra_jet_order':True,'D0_origin_socle_or_depth_claimed':False,
      'stored_nonzero_coefficient_count':totalterms,'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},
      'evidence_sha256':{ev.name:sha(ev)},'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    jsonwrite(out/'report.json',report);phase['status']='FINAL';jsonwrite(out/'phase.json',phase)
    print(json.dumps({k:report[k]for k in ['id','all_passed','status','c_v','nonzero_graph_count','kappa_by_preflight_lemma','seconds','max_rss_kib']}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
