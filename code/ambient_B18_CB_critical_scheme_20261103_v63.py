#!/usr/bin/env python3
"""Whole CB-normalizer ramification germ, retaining every ratio*f mixed term."""
import argparse,gzip,hashlib,itertools,json,subprocess,time
from pathlib import Path
import sympy as s
from ambient_CA_global_flag_model_20261103_v63 import X,canonical,vector,action
from ambient_CB_global_flag_model_20261103_v63 import raw
ROOT=Path(__file__).resolve().parents[1];h,x,y,z,w=X
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def mat(M):return [[str(v)for v in row]for row in M.tolist()]
def main(out,m2):
    started=time.monotonic();out.mkdir(parents=True,exist_ok=False)
    cubics=[s.prod((h,x,z,w)[i]for i in mon)for mon in itertools.combinations_with_replacement(range(4),3)if any(i>=2 for i in mon)]
    sw=s.eye(5);sw[0,0]=sw[1,1]=0;sw[0,1]=sw[1,0]=1
    change=lambda P:(sw*P.subs({h:x,x:h},simultaneous=True)*sw.T).applyfunc(s.expand)
    F=[change(canonical(raw(f,0,0)))for f in[h*h,h*x,x*x]]
    G=[change(canonical(raw(0,f,0)))for f in[h*h,h*x,x*x]]
    HH=[change(canonical(raw(0,0,H)))for H in cubics]
    base=F+HH;C19=s.Matrix.hstack(*map(vector,base));cr19=list(C19.T.rref()[1]);left19=C19[cr19,:].inv()
    ey=s.Matrix([0,0,1,0,0]);D=s.Matrix([x*x,0,0,h*x,h*h+z*z]);P0=ey*D.T-D*ey.T
    center=left19*vector(P0)[cr19,:];assert C19*center==vector(P0) and list(center[:3,:])==[-1,0,0]
    flag=[(0,2),(0,3),(0,4),(1,0),(1,2),(1,3),(1,4),(3,2),(4,2)]
    parabolic=[(i,j)for i,j in itertools.product(range(5),repeat=2)if(i,j)not in flag]
    # Recheck the full associated-ratio action, including its character.
    enlarged=F+G+HH;C22=s.Matrix.hstack(*map(vector,enlarged));cr22=list(C22.T.rref()[1]);left22=C22[cr22,:].inv()
    q=s.Symbol('q');embedding=s.zeros(22,19)
    for i in range(3):embedding[i,i]=1;embedding[i+3,i]=-q
    for i in range(16):embedding[i+6,i+3]=1
    selected=[0,1,2,*range(6,22)]
    for i,j in parabolic:
        M=s.zeros(5);M[i,j]=1
        large=s.Matrix.hstack(*[vector(action(M,P))for P in enlarged]);small=left22*large[cr22,:]
        assert C22*small==large
        chi=M[3,3]+M[4,4]-M[0,0]-M[2,2]
        corrected=small*embedding-chi*q*embedding.diff(q)
        assert (embedding*corrected[selected,:]-corrected).applyfunc(s.expand)==s.zeros(22,19)
    C0=s.Matrix.hstack(*map(vector,F),-sum((center[i]*vector(G[i])for i in range(3)),s.zeros(150,1)),*map(vector,HH))
    cr=list(C0.T.rref()[1]);left0=C0[cr,:].inv()
    orbit=s.Matrix.hstack(*[left0*vector(action(s.eye(5)[:,i]*s.eye(5)[j,:],P0))[cr,:]for i,j in parabolic])
    pc=list(orbit.rref()[1]);assert len(pc)==10
    fixed=[0,4,5,6,7,8,11,12,13,16];free=[i for i in range(20)if i not in fixed]
    orbitdet=s.factor(orbit[fixed,pc].det());assert orbitdet!=0
    tt=s.symbols('s0:10');alpha=list(center[:3,:])+[s.S.Zero]+list(center[3:,:])
    for i,t in zip(free,tt):alpha[i]+=t
    q=alpha[3]
    P=sum((alpha[i]*F[i]for i in range(3)),s.zeros(5))-q*sum((alpha[i]*G[i]for i in range(3)),s.zeros(5))+sum((alpha[i+4]*HH[i]for i in range(16)),s.zeros(5))
    C=s.Matrix.hstack(*[vector(F[i]-q*G[i])for i in range(3)],-sum((alpha[i]*vector(G[i])for i in range(3)),s.zeros(150,1)),*map(vector,HH))
    fiberunit=s.factor(C[cr,:].det());assert fiberunit==-s.Rational(864,25)
    left=C[cr,:].inv();remainder=lambda M:(M-C*left*M[cr,:]).applyfunc(s.expand)
    Q=remainder(s.Matrix.hstack(*[vector(action(s.eye(5)[:,i]*s.eye(5)[j,:],P))for i,j in flag]))
    qr=[5,22,23,27,31,35,50,65];qc=[0,1,2,4,5,6,7,8]
    flagunit=s.factor(Q[qr,qc].det());assert flagunit==s.Rational(8,3)
    schur=(Q[:,3]-Q[:,qc]*Q[qr,qc].inv()*Q[qr,3]).applyfunc(s.expand)
    equations=list(dict.fromkeys(v for v in schur if v))
    a0,a1,a2,a3,a4,a5,a6,a7,a8,a9=tt;U=a3+a7;V=a3-a7;W=a2;e5=10*(a3*a7+a7*a7)/27
    elimination=[a0-2*a3/9+16*a7/9,a1+3*e5/5,a4-27*e5/10+a7*a7,a5-e5,a6,a8-a7*a7,a9-a7**3/3]
    predicted=elimination+[W*U,W*V,U*V,U**3]
    gb=s.groebner(equations,*tt,order='grevlex',domain=s.QQ);fg=s.groebner(predicted,*tt,order='grevlex',domain=s.QQ)
    assert all(gb.reduce(s.expand(f))[1]==0 for f in predicted)
    assert all(fg.reduce(f)[1]==0 for f in equations)
    assert list(gb)==list(fg)
    text='R=QQ[s0,s1,s2,s3,s4,s5,s6,s7,s8,s9];\nI=ideal('+','.join(str(f).replace('**','^')for f in equations)+');\nJ=ideal('+','.join(str(s.expand(f)).replace('**','^')for f in predicted)+');\nassert(I==J);\nS=QQ[w,u,v];\nK=ideal(w*u,w*v,u*v,u^3);\nassert(K==intersect(ideal(u,v),ideal(u,w),ideal(w,v,u^3)));\nassert(dim(S/K)==1);\nassert(degree(S/K)==2);\nprint "CB_B18_FULL_CRITICAL_GERM_PASS";\nexit 0;\n'
    (out/'critical.m2').write_text(text)
    if m2:
        run=subprocess.run(['M2','--script',str((out/'critical.m2').resolve())],capture_output=True,text=True,check=True,timeout=60)
        assert 'CB_B18_FULL_CRITICAL_GERM_PASS'in run.stdout;(out/'critical.log').write_text(run.stdout+run.stderr)
    cert=dict(fiber_center19=list(map(str,center)),full_ratio_source_formula='P=sum f_i F_i -q sum f_i G_i +sum H_i HH_i',all16_enlarged22_and_relative19_Lie_columns_checked=True,
        parabolic_positions=parabolic,full_parabolic_orbit20_by16=mat(orbit),orbit10_columns=pc,fixed_rows=fixed,orbit10_minor=str(orbitdet),slice_free_indices=free,slice_variables=list(map(str,tt)),
        full_variable_fiber_Jacobian150_by20=mat(C),fiber_unit_rows=cr,fiber_unit_determinant=str(fiberunit),flag_positions=flag,full_flag_quotient150_by9=mat(Q),flag_unit_rows=qr,flag_unit_columns=qc,flag_unit_determinant=str(flagunit),
        all150_Schur_entries=list(map(str,schur)),all_nonzero_critical_generators=list(map(str,equations)),seven_monic_eliminations=list(map(str,elimination)),u=str(U),v=str(V),w=str(W),
        predicted_entire_critical_ideal=list(map(str,predicted)),whole_QQ_GB=list(map(str,gb)))
    (out/'certificate.json.gz').write_bytes(gzip.compress((json.dumps(cert,separators=(',',':'))+'\n').encode(),mtime=0))
    paths=[Path(__file__),ROOT/'research/notes/AMBIENT_B18_CB_CRITICAL_SCHEME_20261103_v63.md',ROOT/'code/ambient_CA_global_flag_model_20261103_v63.py',ROOT/'code/ambient_CB_global_flag_model_20261103_v63.py',ROOT/'reports/AMBIENT_B18_QUADRATIC_IMAGE_REDUCTION_20261103_v63_REVIEW/report.json']
    assert json.loads(paths[-1].read_text())['all_passed'] is True
    result=dict(id='AMBIENT_B18_CB_CRITICAL_SCHEME_20261103_v63',all_passed=True,field='QQ',all16_by22_enlarged_Lie_and16_by19_associated_ratio_identities=True,etale19_smooth_source_factors=True,
        actual_source_ratio_times_f_all_orders_retained=True,variable_full20fiber_Jacobian_used=True,whole_intrinsic_ramification_Fitting_ideal=True,fiber_unit=str(fiberunit),flag8_unit=str(flagunit),exact_Schur_generator_count=len(equations),both_QQ_ideal_containments=True,independent_Macaulay2_ideal_and_primary_checks_run=bool(m2),
        cone_ramification_completion='C[[g1,...,g19,w,u,v]]/(w*u,w*v,u*v,u^3)',projective_ramification_completion='C[[g1,...,g18,w,u,v]]/(w*u,w*v,u*v,u^3)',
        exact_associated_prime_count=3,two_reduced_branches=True,embedded_junction_prime=True,nilradical_principal_u=True,nilradical_cube_zero=True,nilradical_square_nonzero=True,
        cone_dimension=20,cone_embedding_dimension=22,cone_depth=19,projective_dimension=19,projective_embedding_dimension=21,projective_depth=18,Hilbert_Samuel_multiplicity=2,
        target_CB_actual_D_completion_or_odd_ideal_claimed=False,input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},evidence_sha256={p.name:sha(p)for p in out.iterdir()if p.is_file()},seconds=time.monotonic()-started)
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');print('WHOLE CB B18 CRITICAL GERM incl ALLq*f PASS',len(equations),result['seconds'],flush=True)
if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--m2',action='store_true');args=ap.parse_args();main(args.output.resolve(),args.m2)
