#!/usr/bin/env python3
"""Independent full CB source-critical audit using coefficient Lie operators."""
import gzip,hashlib,itertools,json,time
from pathlib import Path
import sympy as s
from ambient_B18_CA_N_critical_review_20261103_v63 import ROOT,X,h,x,y,z,w,vec,canon,coefficient_action
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def raw(f=0,g=0,H=0):
 P=s.zeros(5);P[0,2]=f;P[3,4]=g;P[2,3]=-s.diff(f,h)*z/2-s.diff(H,w);P[2,4]=-s.diff(f,h)*w/2+s.diff(H,z)
 return canon(P-P.T)
def main():
 start=time.monotonic();folder=ROOT/'reports/AMBIENT_B18_CB_CRITICAL_SCHEME_20261103_v63';p=json.loads((folder/'report.json').read_text());assert p['all_passed']
 for f,want in p['input_sha256'].items():assert sha(ROOT/f)==want,f
 for f,want in p['evidence_sha256'].items():assert sha(folder/f)==want,f
 d=json.loads(gzip.decompress((folder/'certificate.json.gz').read_bytes()))
 FF=[vec(raw(f=f))for f in[x*x,h*x,h*h]];GG=[vec(raw(g=f))for f in[x*x,h*x,h*h]]
 HH=[vec(raw(H=s.prod((x,h,z,w)[i]for i in mon)))for mon in itertools.combinations_with_replacement(range(4),3)if any(i>=2 for i in mon)]
 C19=s.Matrix.hstack(*FF,*HH);r19=list(C19.T.rref()[1]);iv19=C19[r19,:].inv()
 P0=s.zeros(5)
 for j,v in[(0,x*x),(3,h*x),(4,h*h+z*z)]:P0[2,j]=v;P0[j,2]=-v
 vp=vec(P0);center=iv19*vp[r19,:];assert C19*center==vp and center==s.Matrix(d['fiber_center19'])
 E=s.Matrix.hstack(*FF,*GG,*HH);rr=list(E.T.rref()[1]);iv=E[rr,:].inv();q=s.Symbol('q');embed=s.zeros(22,19)
 for i in range(3):embed[i,i]=1;embed[i+3,i]=-q
 for i in range(16):embed[i+6,i+3]=1
 flag=[tuple(v)for v in d['flag_positions']];para=[tuple(v)for v in d['parabolic_positions']];acts={ab:coefficient_action(*ab)for ab in flag+para}
 select=[0,1,2,*range(6,22)]
 for a,b in para:
  M=acts[a,b]*E;rep=iv*M[rr,:];assert E*rep==M
  chi=int(a==b==3)+int(a==b==4)-int(a==b==0)-int(a==b==2)
  corrected=rep*embed-chi*q*embed.diff(q);assert(embed*corrected[select,:]-corrected).applyfunc(s.expand)==s.zeros(22,19)
 C0=s.Matrix.hstack(*FF,-sum((center[i]*GG[i]for i in range(3)),s.zeros(150,1)),*HH);cr=d['fiber_unit_rows'];iv0=C0[cr,:].inv()
 O=s.Matrix.hstack(*(iv0*(acts[ab]*vp)[cr,:]for ab in para));assert O==s.Matrix(d['full_parabolic_orbit20_by16'])
 assert O[d['fixed_rows'],d['orbit10_columns']].det()==s.Rational(d['orbit10_minor'])!=0
 vv=s.symbols(' '.join(d['slice_variables']));alpha=list(center[:3,:])+[s.S.Zero]+list(center[3:,:])
 for i,t in zip(d['slice_free_indices'],vv):alpha[i]+=t
 q=alpha[3];vp=sum((alpha[i]*(FF[i]-q*GG[i])for i in range(3)),s.zeros(150,1))+sum((alpha[i+4]*HH[i]for i in range(16)),s.zeros(150,1))
 C=s.Matrix.hstack(*[FF[i]-q*GG[i]for i in range(3)],-sum((alpha[i]*GG[i]for i in range(3)),s.zeros(150,1)),*HH)
 assert C==s.Matrix(d['full_variable_fiber_Jacobian150_by20'])
 assert s.expand(C[cr,:].det())==-s.Rational(864,25);iv=C[cr,:].inv()
 Q=s.Matrix.hstack(*(acts[ab]*vp for ab in flag));Q=(Q-C*iv*Q[cr,:]).applyfunc(s.expand);assert Q==s.Matrix(d['full_flag_quotient150_by9'])
 qr=d['flag_unit_rows'];qc=d['flag_unit_columns'];extra=next(iter(set(range(9))-set(qc)));B=Q[qr,qc];assert s.expand(B.det())==s.Rational(8,3)
 schur=(Q[:,extra]-Q[:,qc]*B.inv()*Q[qr,extra]).applyfunc(s.expand);assert schur==s.Matrix(d['all150_Schur_entries'])
 expected=list(map(s.sympify,d['predicted_entire_critical_ideal']));G=s.groebner(list(schur),*vv,order='lex',domain=s.QQ);F=s.groebner(expected,*vv,order='lex',domain=s.QQ)
 assert all(G.reduce(s.expand(f))[1]==0 for f in expected)and all(F.reduce(f)[1]==0 for f in schur)
 note=ROOT/'research/notes/AMBIENT_B18_CB_CRITICAL_SCHEME_20261103_v63_REVIEW.md';files=[Path(__file__),note,ROOT/'code/ambient_B18_CA_N_critical_review_20261103_v63.py',ROOT/'reports/AMBIENT_B18_CA_N_CRITICAL_SCHEMES_20261103_v63_REVIEW.json',folder/'report.json',folder/'certificate.json.gz']
 result=dict(id='AMBIENT_B18_CB_CRITICAL_SCHEME_20261103_v63_REVIEW',all_passed=True,own_physical22_basis_and_actual_center=True,all16_enlarged22_and_relative_ratio19_stability_identities=True,all_ratio_times_f_mixed_terms_included=True,orbit10_etale_minor=True,variable_fiber20_and_flag9_reconstructed_all150_coefficients=True,both_global_constant_unit_blocks_verified=True,all150_schur_polynomials_and_both_whole_QQ_ideal_containments=True,source_ramification_completion='C[[g19,w,u,v]]/(wu,wv,uv,u^3)',target_CB_or_actual_D_or_odd_ideal_claimed=False,input_sha256={str(f.relative_to(ROOT)):sha(f)for f in files},seconds=time.monotonic()-start)
 (ROOT/'reports/AMBIENT_B18_CB_CRITICAL_SCHEME_20261103_v63_REVIEW.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
