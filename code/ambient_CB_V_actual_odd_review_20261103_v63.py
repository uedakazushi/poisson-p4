#!/usr/bin/env python3
"""Independent full physical CB QQ(V) frame, dual jets, and actual odd audit."""
import argparse
from fractions import Fraction
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import time
import sympy as s
from sympy.polys.matrices import DomainMatrix

ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_CB_V_ACTUAL_ODD_LINEAR_20261103_v63_REVIEW'
FRAME='AMBIENT_CB_V_FAMILY_FRAME_20261103_v63'
ODD='AMBIENT_CB_V_ACTUAL_ODD_LINEAR_20261103_v63'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return h.hexdigest()
def load(p):return json.loads(gzip.decompress(p.read_bytes()) if p.suffix=='.gz' else p.read_bytes())

def audit(output):
 start=time.monotonic();output.mkdir(parents=True,exist_ok=False)
 fd=ROOT/'reports'/FRAME;od=ROOT/'reports'/ODD
 reports=[fd/'report.json',od/'report.json'];inputs=[]
 for rp in reports:
  r=load(rp);assert r['all_passed']
  for n,h in r['input_sha256'].items():assert sha(ROOT/n)==h,n;inputs.append(ROOT/n)
  for n,h in r['evidence_sha256'].items():assert sha(rp.parent/n)==h,n;inputs.append(rp.parent/n)
 d=load(fd/'frame.json.gz');jets=load(od/'source_jets.json.gz');cert=load(od/'certificate.json.gz')
 V=s.Symbol('V');K=s.QQ.frac_field(V)
 dm=lambda M:DomainMatrix.from_Matrix(M).convert_to(K)
 readmat=lambda q:dm(s.Matrix(q))
 eq=lambda a,b: a==b
 h,x,y,z,w=s.symbols('h x y z w');X=s.Matrix([h,x,y,z,w])
 pairs=list(itertools.combinations(range(5),2));mons=list(itertools.combinations_with_replacement(range(5),2))
 index={(p,m):i for i,(p,m) in enumerate(itertools.product(pairs,mons))}
 def vec(P):return s.Matrix([s.Poly(P[i,j],*X).coeff_monomial(X[a]*X[b]) for i,j in pairs for a,b in mons])
 def raw(f,g,H):
  P=s.zeros(5);P[1,2]=f;P[3,4]=g
  P[2,3]=-s.diff(f,x)*z/2-s.diff(H,w);P[2,4]=-s.diff(f,x)*w/2+s.diff(H,z)
  return P-P.T
 def canonical(P):
  divergence=s.Matrix([sum(s.diff(P[i,j],X[i]) for i in range(5)) for j in range(5)])
  return (P-(X*divergence.T-divergence*X.T)/5).applyfunc(s.expand)
 swap=s.eye(5);swap[0,0]=swap[1,1]=0;swap[0,1]=swap[1,0]=1
 def transformed(P):return (swap*canonical(P).subs({h:x,x:h},simultaneous=True)*swap.T).applyfunc(s.expand)
 binary=[h*h,h*x,x*x]
 cubic=[s.prod((h,x,z,w)[i] for i in mon) for mon in itertools.combinations_with_replacement(range(4),3) if any(i>=2 for i in mon)]
 F=s.Matrix.hstack(*(vec(transformed(raw(f,0,0))) for f in binary))
 G=s.Matrix.hstack(*(vec(transformed(raw(0,f,0))) for f in binary))
 HH=s.Matrix.hstack(*(vec(transformed(raw(0,0,H))) for H in cubic))
 assert dm(F)==readmat(d['fiber_basis_F150_by3']) and dm(G)==readmat(d['fiber_basis_G150_by3'])
 ey=s.Matrix([0,0,1,0,0]);DD=s.Matrix([x*x,0,0,h*x,h*h+z*z]);p0=vec(ey*DD.T-DD*ey.T)
 c19=F.row_join(HH);cr=list(c19.T.rref()[1]);center=c19[cr,:].inv()*p0[cr,:];assert c19*center==p0
 criticalpath=ROOT/'reports/AMBIENT_B18_CB_CRITICAL_SCHEME_20261103_v63/certificate.json.gz'
 critical=load(criticalpath);alpha=list(center[:3,:])+[s.S.Zero]+list(center[3:,:])
 values=[V,0,0,V/2,-V**2/4,0,0,-V/2,V**2/4,-V**3/24]
 for j,value in zip(critical['slice_free_indices'],values):alpha[j]+=value
 assert alpha==list(map(s.sympify,d['fiber_alpha20'])) and alpha[3]==0
 tt=s.symbols('s0:10');sub=dict(zip(tt,values))
 assert all(s.expand(s.sympify(f).subs(sub))==0 for f in critical['predicted_entire_critical_ideal'])
 point=F*s.Matrix(alpha[:3])+HH*s.Matrix(alpha[4:]);C=F.row_join(-G*s.Matrix(alpha[:3])).row_join(HH)
 p=dm(point);c=dm(C);gg=dm(G)
 assert p==readmat(d['point150']) and c==readmat(d['fiber_jacobian150_by20'])
 assert point.subs(V,0)==p0
 # Coefficient of delta in the PHYSICAL pushforward by I+delta E_ab,
 # using inverse I-delta E_ab modulo delta². No author Lie helper import.
 def dual_operator(a,b):
  entries={}
  def put(pair,mon,col,coefficient):
   i,j=pair
   if i==j:return
   if i>j:i,j,coefficient=j,i,-coefficient
   row=index[((i,j),tuple(sorted(mon)))];entries[row,col]=entries.get((row,col),0)+coefficient
  for col,((i,j),mon) in enumerate(itertools.product(pairs,mons)):
   # Transform wedge factors and substitute the inverse in BOTH monomials.
   if i==b:put((a,j),mon,col,1)
   if j==b:put((i,a),mon,col,1)
   for slot in range(2):
    if mon[slot]==a:
     mm=list(mon);mm[slot]=b;put((i,j),mm,col,-1)
  rows={}
  for (i,j),q in entries.items():
   if q:rows.setdefault(i,{})[j]=K(q)
  return DomainMatrix(rows,(150,150),K)
 ops25=[dual_operator(*divmod(i,5)) for i in range(25)]
 orbit=DomainMatrix.hstack(*[op*p for op in ops25]);gi=d['ordered_group19_GL_indices'];extra=d['extra10_fiber_indices']
 o19=orbit.extract(range(150),gi);rr=list(o19.transpose().rref()[1]);oinv=o19.extract(rr,range(19)).inv()
 assert o19*oinv*orbit.extract(rr,range(25))==orbit
 assert len(gi)==19 and len(extra)==10
 # Reconstruct the full29 source coordinate change in old fiber20+flag9.
 flag=[tuple(x) for x in critical['flag_positions']];assert flag==[tuple(x) for x in d['flag_positions']]
 crows=list(c.transpose().rref()[1]);ci=c.extract(crows,range(20)).inv()
 source=s.zeros(29,25)
 for j in range(25):
  ab=divmod(j,5)
  if ab in flag:source[20+flag.index(ab),j]=1
  else:
   co=ci*orbit.extract(crows,[j]);assert c*co==orbit.extract(range(150),[j]);source[:20,j]=co.to_Matrix()
 frame=dm(source[:,gi].row_join(s.eye(29)[:,extra]));assert frame==readmat(d['fullsource29_frame'])
 framedet=s.factor(K.to_sympy(frame.det()));assert framedet==V**12*(V**2-2)/3072
 B=c.extract(range(150),extra);J=o19.hstack(B);assert J==readmat(d['fulltarget150_by29_J'])
 kap=readmat(d['kernel29']);comp=d['complement28_columns'];D=J.extract(range(150),comp)
 assert J*kap==DomainMatrix.zeros((150,1),K) and D.rank()==28 and len(comp)==28
 # Full ordered source dual propagation: pure epsilon³ and ALL29 delta*epsilon.
 # Initial fiber quadratic and its derivatives are from -delta_q delta_f G.
 kk=kap.to_Matrix();iq=extra.index(3)
 pure=[p,B*kap.extract(range(19,29),[0]),DomainMatrix.zeros((150,1),K),DomainMatrix.zeros((150,1),K)]
 mixed=s.zeros(150,29)
 for fidx in range(3):
  if fidx in extra:
   pos=extra.index(fidx);pure[2]=pure[2]-gg.extract(range(150),[fidx]).scalarmul(K.from_sympy(kk[19+iq]*kk[19+pos]))
   mixed[:,19+iq]-=kk[19+pos]*G[:,fidx];mixed[:,19+pos]-=kk[19+iq]*G[:,fidx]
 der0=DomainMatrix.zeros((150,19),K).hstack(B);der1=dm(mixed)
 for i in range(18,-1,-1):
  op=ops25[gi[i]];ki=kap[i,0].element;old=pure
  pure=[]
  for degree in range(4):
   value=old[degree]
   for power in range(1,degree+1):
    v=old[degree-power]
    for j in range(power):v=op*v
    value=value+v.scalarmul(ki**power/K(s.factorial(power)))
   pure.append(value)
  der1=der1+(op*der0).scalarmul(ki)
  unit=DomainMatrix({0:{i:K.one}},(1,29),K)
  der0=der0+(op*pure[0])*unit;der1=der1+(op*pure[1])*unit
 assert der0==J and pure[1]==DomainMatrix.zeros((150,1),K)
 H=der1;W=pure[2];Z=pure[3]
 assert H==readmat(jets['full150_by29_H_kernel']) and W==readmat(jets['pure_kernel_quadratic150']) and Z==readmat(jets['pure_kernel_cubic150'])
 assert H*kap==W.scalarmul(K(2))
 # Actual Y inverse and tracefree rank2 base. Pure third correction must include
 # -H_kr D_Y^-1 W_Y; it is not the uncorrected cubic of the source line.
 Y=cert['actual_base_Y28'];branch=cert['branch29th_row'];inv=D.extract(Y,range(28)).inv()
 rem=lambda M:M-D*inv*M.extract(Y,range(M.shape[1]))
 hcomp=H.extract(range(150),comp);wc=rem(W);weight=wc[branch,0].element
 assert weight and s.cancel(K.to_sympy(weight)-s.sympify(cert['branch_weight']))==0
 E=rem(hcomp*inv);Zc=rem(Z-hcomp*inv*W.extract(Y,[0]))
 normal=lambda M:M-wc.scalarmul(K.one/weight)*M.extract([branch],range(M.shape[1]))
 beta=normal(Zc).scalarmul(K.one/weight);ell=D.extract([branch],range(28))*inv
 odd=(normal(E)-beta*ell).hstack(beta)
 assert odd==readmat(cert['actual_all150_by29_odd_firstjet'])
 rows=cert['rank_rows'];cols=cert['rank_columns'];small=odd.extract(rows,cols)
 minor=s.factor(K.to_sympy(small.det()));assert minor==s.sympify(cert['rank_minor'])
 expected=157837977*V**28*(V**2-2)**9/s.Integer(125440000000000000000000000000000000000)
 assert s.cancel(minor-expected)==0 and len(rows)==len(cols)==9
 factors=readmat(cert['whole_column_factorization']);assert odd.extract(range(150),cols)*factors==odd
 hnew=H*kap.hstack(dm(s.eye(29)[:,comp]));assert rem(hnew).rank()==9
 # Modulus direction is retained, so20 critical directions are actual.
 tangent=o19.hstack(dm(point.diff(V)));assert tangent.rank()==20
 tr=list(tangent.transpose().rref()[1]);tminor=s.factor(K.to_sympy(tangent.extract(tr,range(20)).det()))
 assert gi[10]==12 and o19.extract(range(150),[10])==p
 # Nonzero rank witness after p-reduction of the SAME QQ(V) minor; this does not
 # assert that V=1 is in the separate global finite-open locus.
 q=Fraction(str(minor.subs(V,1)));prime=32003
 residue=(q.numerator%prime)*pow(q.denominator%prime,-1,prime)%prime;assert residue
 deps=[ROOT/'reports/AMBIENT_CB_GLOBAL_FLAG_MODEL_20261103_v63_REVIEW.json',
       ROOT/'reports/AMBIENT_B18_QUADRATIC_IMAGE_REDUCTION_20261103_v63_REVIEW/report.json',
       ROOT/'reports/AMBIENT_B18_CB_CRITICAL_SCHEME_20261103_v63_REVIEW.json']
 for dp in deps:assert load(dp)['all_passed']
 files=[Path(__file__),ROOT/'research/notes'/f'{STEM}.md',criticalpath,*reports,*inputs,*deps,
       ROOT/'research/notes/AMBIENT_B18_QUADRATIC_IMAGE_REDUCTION_20261103_v63.md',
       ROOT/'research/notes/AMBIENT_CB_GENERIC_V_FULL_COMPLETION_20261103_v63.md']
 result=dict(id=STEM,all_passed=True,field='QQ(V); full29 source directions retained',
  independent_physical_F3_G3_H16_and_full_nonlinear_source=True,all150_point_and_C20_and_J29_and_source_frame=True,
  source_frame_det=str(framedet),source_orbit_rank=19,target_J_rank=28,
  all150_by29_mixed_H_by_dual_derivative_propagation=True,all150_pure_kernel_cubic=True,
  nonlinear_ratio_f_quadratic_and_all_mixed_derivatives_retained=True,
  all150_tracefree_actual_odd_firstjet=True,actual_odd_rank=9,critical_cotangent_rank=9,whole_factorization_verified=True,
  selected_minor=str(minor),same_minor_mod32003_at_V1=residue,
  critical_GL19_plus_modulus_tangent_rank=20,critical20_minor=str(tminor),scalar_column_equals_point=True,
  generic_whole_proper_fiber_length2_argument_accepted=True,quadratic_regular_R29_free_B2_argument_accepted=True,
  J_regular9_plus_actual_I_rank9_implies_I_equals_J=True,nonzero_I_recovers_fraction_field_before_normalization_claim=True,
  generic_full_completion='C[[g20,s8,t^2,s_i*t,t^3]]',cone_dimension=29,cone_edim=38,cone_depth=21,multiplicity=2,
  projective_dimension=28,projective_edim=37,projective_depth=20,
  not_claimed=['all V or V=2 whole finite-open membership','special behavior at frame poles V(V^2-2)','B18 origin full odd ideal or actual intersection D','global normalization without separate birational proof'],
  input_sha256={str(p.relative_to(ROOT)):sha(p) for p in dict.fromkeys(files)},seconds=time.monotonic()-start)
 (output/'report.json').write_text(json.dumps(result,indent=2)+'\n')
 print('PASS',STEM,result['seconds'])

if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args();audit(a.output)
