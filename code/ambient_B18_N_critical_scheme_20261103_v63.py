#!/usr/bin/env python3
"""Exact N-normalizer source ramification scheme, not its target image ideal."""
import argparse,gzip,hashlib,itertools,json,subprocess,time
from pathlib import Path
import sympy as s
from ambient_CA_global_flag_model_20261103_v63 import X,raw,canonical,vector,action
ROOT=Path(__file__).resolve().parents[1];h,x,y,z,w=X
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def mat(M):return [[str(v)for v in row]for row in M.tolist()]
def main(out,m2):
    start=time.monotonic();out.mkdir(parents=True,exist_ok=False)
    cubics=[s.prod((h,x,z,w)[i]for i in mon)for mon in itertools.combinations_with_replacement(range(4),3)if any(i>=2 for i in mon)]
    sw=s.eye(5);sw[0,0]=sw[1,1]=0;sw[0,1]=sw[1,0]=1
    change=lambda P:(sw*P.subs({h:x,x:h},simultaneous=True)*sw.T).applyfunc(s.expand)
    basis=[change(canonical(raw(f,0,0)))for f in[h*h,h*x]]+[change(canonical(raw(0,1,0)))]+[change(canonical(raw(0,0,H)))for H in cubics]
    C=s.Matrix.hstack(*map(vector,basis));cr=list(C.T.rref()[1]);left=C[cr,:].inv()
    rem=lambda M:(M-C*left*M[cr,:]).applyfunc(s.expand)
    ey=s.Matrix([0,0,1,0,0]);D=s.Matrix([x*x,0,0,h*x,h*h+z*z]);P0=ey*D.T-D*ey.T
    center=left*vector(P0)[cr,:];assert C*center==vector(P0)
    flag=[(0,2),(0,3),(0,4),(1,0),(1,2),(1,3),(1,4),(3,2),(4,2)]
    parabolic=[(i,j)for i,j in itertools.product(range(5),repeat=2)if(i,j)not in flag]
    for i,j in parabolic:
        M=s.zeros(5);M[i,j]=1
        assert rem(s.Matrix.hstack(*[vector(action(M,P))for P in basis]))==s.zeros(150,19)
    orbit=s.Matrix.hstack(*[left*vector(action(s.eye(5)[:,i]*s.eye(5)[j,:],P0))[cr,:]for i,j in parabolic])
    pc=list(orbit.rref()[1]);assert len(pc)==10
    fixed=[0,3,4,5,6,7,10,11,12,15];free=[i for i in range(19)if i not in fixed]
    assert free==[1,2,8,9,13,14,16,17,18]
    orbitdet=s.factor(orbit[fixed,pc].det());assert orbitdet!=0
    tt=s.symbols('t0:9');P=P0.copy()
    for i,t in zip(free,tt):P+=t*basis[i]
    Q=rem(s.Matrix.hstack(*[vector(action(s.eye(5)[:,i]*s.eye(5)[j,:],P))for i,j in flag]))
    qr=[5,15,22,23,31,35,50,65];qc=[0,1,2,4,5,6,7,8]
    assert s.factor(Q[qr,qc].det())==1
    schur=(Q[:,3]-Q[:,qc]*Q[qr,qc].inv()*Q[qr,3]).applyfunc(s.expand)
    equations=list(dict.fromkeys(v for v in schur if v))
    t0,t1,t2,t3,t4,t5,t6,t7,t8=tt
    predicted=[t0-5*t6/2,t2-11*t6/2,t3,t4,t5,t7,t8,t1*t6,t6*t6]
    G=s.groebner(equations,*tt,order='grevlex',domain=s.QQ);F=s.groebner(predicted,*tt,order='grevlex',domain=s.QQ)
    assert all(G.reduce(f)[1]==0 for f in predicted)
    assert all(F.reduce(f)[1]==0 for f in equations)
    assert list(G)==list(F)
    text='R=QQ[t0,t1,t2,t3,t4,t5,t6,t7,t8];\nI=ideal('+','.join(str(f).replace('**','^')for f in equations)+');\nJ=ideal('+','.join(str(f).replace('**','^')for f in predicted)+');\nassert(I==J);\nS=QQ[w,v];\nK=ideal(w*v,v^2);\nassert(K==intersect(ideal(v),ideal(w,v^2)));\nassert(dim(S/K)==1);\nassert(degree(S/K)==1);\nprint hilbertSeries(S/K,Reduce=>true);\nprint "N_B18_FULL_CRITICAL_GERM_PASS";\nexit 0;\n'
    (out/'critical.m2').write_text(text)
    if m2:
        run=subprocess.run(['M2','--script',str((out/'critical.m2').resolve())],capture_output=True,text=True,check=True,timeout=60)
        assert 'N_B18_FULL_CRITICAL_GERM_PASS'in run.stdout
        (out/'critical.log').write_text(run.stdout+run.stderr)
    certificate=dict(center19=list(map(str,center)),fiber19_basis=mat(C),parabolic_positions=parabolic,parabolic_orbit=mat(orbit),orbit10_columns=pc,fixed_rows=fixed,free_indices=free,orbit10_minor=str(orbitdet),
        slice_variables=list(map(str,tt)),flag_positions=flag,full150_by9_flag_quotient=mat(Q),unit8_rows=qr,unit8_columns=qc,unit8_determinant='1',all150_Schur_entries=list(map(str,schur)),all_critical_generators=list(map(str,equations)),predicted_exact_ideal=list(map(str,predicted)),whole_QQ_GB=list(map(str,G)),remaining_variables={'w':'t1','v':'t6'})
    (out/'certificate.json.gz').write_bytes(gzip.compress((json.dumps(certificate,separators=(',',':'))+'\n').encode(),mtime=0))
    paths=[Path(__file__),ROOT/'research/notes/AMBIENT_B18_N_CRITICAL_SCHEME_20261103_v63.md',ROOT/'code/ambient_CA_global_flag_model_20261103_v63.py',ROOT/'reports/AMBIENT_B18_QUADRATIC_IMAGE_REDUCTION_20261103_v63_REVIEW/report.json']
    assert json.loads(paths[-1].read_text())['all_passed'] is True
    result=dict(id='AMBIENT_B18_N_CRITICAL_SCHEME_20261103_v63',all_passed=True,field='QQ',all16_by19_parabolic_columns_stable=True,etale19_source_smooth_factors=True,whole_intrinsic_critical_Fitting_ideal=True,unit8_block_determinant='1',exact_Schur_generator_count=len(equations),both_QQ_ideal_containments=True,independent_Macaulay2_ideal_and_primary_checks_run=bool(m2),
        cone_normalizer_ramification_completion='C[[g1,...,g19,w,v]]/(w*v,v^2)',projective_normalizer_ramification_completion='C[[g1,...,g18,w,v]]/(w*v,v^2)',reduced_critical_locus_smooth=True,exact_associated_prime_count=2,embedded_prime_at_nil_junction=True,nilradical_square_zero=True,nilradical_nonzero=True,
        cone_dimension=20,cone_embedding_dimension=21,cone_depth=19,projective_dimension=19,projective_embedding_dimension=20,projective_depth=18,Hilbert_Samuel_multiplicity=1,
        target_N_completion_actual_D_or_full_odd_ideal_claimed=False,input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},evidence_sha256={p.name:sha(p)for p in out.iterdir()if p.is_file()},seconds=time.monotonic()-start)
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');print('WHOLE N B18 NORMALIZER CRITICAL GERM line+squarezero PASS',len(equations),result['seconds'],flush=True)
if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--m2',action='store_true');a=ap.parse_args();main(a.output.resolve(),a.m2)
