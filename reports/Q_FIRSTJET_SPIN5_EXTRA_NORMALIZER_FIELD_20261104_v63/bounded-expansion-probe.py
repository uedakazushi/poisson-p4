import sympy as S,json,time,signal,resource
from sympy.polys.rings import ring
from pathlib import Path
resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60)
start=time.monotonic();out=Path('work/spin4-extra-normalizer-v63/action-probe.json');result={'complete':False,'stage':'start'}
def save(stage):
 result['stage']=stage;result['seconds']=time.monotonic()-start;out.write_text(json.dumps(result,indent=2)+'\n');print(stage,flush=True)
def alarm(*_):save('TIME_LIMIT_AT_'+result['stage']);raise SystemExit(124)
signal.signal(signal.SIGALRM,alarm)
RR,p,q,r,s,t=ring('p,q,r,s,t',S.QQ);K=RR.to_field();p,q,r,s,t=map(K,[p,q,r,s,t]);R=(p*q-s*t)**2-(p*r-s)*(q-r*t)
J0=2*(p*r-s);J1=4*(p*q-s*t)
z=(K.zero,K.zero);one=(K.one,K.zero)
def const(a):return(K(a),K.zero)
def add(*args):return tuple(sum((x[i]for x in args),K.zero)for i in range(2))
def neg(a):return(-a[0],-a[1])
def mul(a,b):return(a[0]*b[0]+R*a[1]*b[1],a[0]*b[1]+a[1]*b[0])
def scale(a,c):return(a[0]*c,a[1]*c)
def pw(a,n):
 ans=one
 for i in range(n):ans=mul(ans,a)
 return ans
A00=(K(S.Rational(1,2)),-J1/(8*R));A10=(K.zero,J0/(4*R));A01=(-J1/(2*J0),-2/J0);A11=one
assert add(A00,neg(mul(A01,A10)))==one
M00=neg(add(scale(pw(A00,2),p),mul(A00,A10),scale(pw(A10,2),t)))
M01=neg(add(scale(pw(A00,2),s),scale(mul(A00,A10),r),scale(pw(A10,2),q)))
M10=neg(add(scale(pw(A01,2),p),A01,const(t)))
M11=neg(add(scale(pw(A01,2),s),scale(A01,r),const(q)))
assert add(mul(M00,M11),neg(mul(M01,M10)))==(K.zero,-K.one)
assert add(scale(mul(A00,A01),2*p),A00,mul(A10,A01),scale(A10,2*t))==z
assert add(scale(mul(A00,A01),2*s),scale(A00,r),scale(mul(A10,A01),r),scale(A10,2*q))==z
save('exact_frame_and_pure_square_M')
fac=(K.zero,-1/R)
b00=mul(fac,M11);b01=neg(mul(fac,M01));b10=neg(mul(fac,M10));b11=mul(fac,M00)
w=[add(mul(A00,pw(b00,2)),mul(A10,pw(b10,2))),scale(add(mul(A00,mul(b00,b01)),mul(A10,mul(b10,b11))),2),add(mul(A00,pw(b01,2)),mul(A10,pw(b11,2))),add(mul(A01,pw(b00,2)),pw(b10,2)),scale(add(mul(A01,mul(b00,b01)),mul(b10,b11)),2),add(mul(A01,pw(b01,2)),pw(b11,2))]
save('six_normalized_W_coefficients')
images={}
expressions=[mul((K.zero,-K.one),mul(w[1],w[4])),scale(add(mul(w[0],pw(w[1],3)),mul(w[5],pw(w[4],3))),R**2),mul((K.zero,-R),add(mul(w[0],pw(w[1],3)),neg(mul(w[5],pw(w[4],3))))),scale(add(mul(w[3],pw(w[1],5)),neg(mul(w[2],pw(w[4],5)))),R**3),mul((K.zero,-R**2),add(mul(w[3],pw(w[1],5)),mul(w[2],pw(w[4],5))))]
for i,a in enumerate(expressions,1):
 assert not a[1],i
 images[str(i)]=str(a[0]);result['images_before_old_z_scalars']=images;save('image_'+str(i))
result['complete']=True;result['R']=str(R);result['old_z_scalar_powers_for_images1_to5']=[1,2,2,3,3];save('COMPLETE')
