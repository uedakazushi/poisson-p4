#!/usr/bin/env python3
"""Whole Anil3 rank-three-quadric sources: reduced support and orbit identities.
The original complete linear sources and their Jacobi lists are inputs certified
independently in CLASSIFICATION_DEGENERATE_QUADRIC_SOURCES_REVIEW_20261107.json.
No arbitrary source slice replaces the algebraic branch decomposition below.
"""
import hashlib,json,resource,time
from pathlib import Path
import sympy as s
ROOT=Path(__file__).resolve().parents[1]
xx=s.symbols('x0:5'); nn=s.symbols('n0:12'); a,b,c,r,z=xx
ENV={str(v):v for v in xx+nn}; q=2*a*c+b*b
TYPES=['CNIL2_UPPER','CNIL2_LOWER','CNIL2_ZERO','CZERO_NONZERO','CZERO_ZERO']

def read(typ):
    path=ROOT/('reports/CLASSIFICATION_DEGENERATE_QUADRIC_Q3_ANIL3_'+typ+'_20261107.json')
    raw=json.loads(path.read_text());P=s.zeros(5)
    for key,value in raw['brackets'].items():
        i,j=map(int,key);P[i,j]=s.sympify(value,locals=ENV)
    P-=P.T
    return path,raw,P,[s.sympify(v,locals=ENV) for v in raw['residual_jacobi_coefficients']],s.cancel(s.sympify(raw['euler_contracted_pfaffian_not_yet_actual_divisor'],locals=ENV)/q)

def brackets(P):return {str(i)+str(j):str(s.factor(P[i,j])) for i in range(5) for j in range(i+1,5) if P[i,j]!=0}
def sub(P,at):return P.subs(at,simultaneous=True).applyfunc(s.expand)
def branch(data,at):
    _,_,P,J,C=data
    assert all(s.expand(v.subs(at,simultaneous=True))==0 for v in J)
    return sub(P,at),s.factor(C.subs(at,simultaneous=True))

def identity(P,target,coords):
    T=s.Matrix(coords).jacobian(xx)
    assert all(s.expand(v)==0 for v in T*P*T.T-target.subs(dict(zip(xx,coords)),simultaneous=True))
    return s.factor(T.det())

def main():
    resource.setrlimit(resource.RLIMIT_CPU,(80,81));resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2))
    start=time.monotonic();data={k:read(k) for k in TYPES}; results={}; families={}
    # Necessary upper constraints give n6 in {2,11/8}, n5=0, then these complete branches.
    up={nn[6]:2,nn[5]:0,nn[2]:0,nn[1]:-4*nn[3],nn[0]:s.Rational(3,5)*nn[4]**2}
    upper,Cup=branch(data['CNIL2_UPPER'],up)
    _,Cempty=branch(data['CNIL2_UPPER'],{nn[6]:s.Rational(11,8),nn[5]:0,nn[2]:-nn[4],nn[3]:0,nn[0]:0});assert Cempty==0
    baseup=sub(upper,{nn[3]:0,nn[4]:0})
    coordup=[a,b,c,r-s.Rational(2,5)*nn[4]*a,z-s.Rational(2,5)*nn[4]*b+s.Rational(6,5)*nn[3]*a]
    assert identity(baseup,upper,coordup)==1
    assert s.expand(s.diff(Cup,z)-s.Rational(5,12)*q)==0
    assert Cup.subs({a:0,b:0})==s.Rational(5,36)*c**3
    families['upper']={'source':up,'free':['n3','n4'],'open':'all','P':upper,'C':Cup,'base':baseup,'coords':coordup,'det':s.Integer(1)}
    # Lower constraints (2n2-1)(4n8+3), (n8+2)(n2+3n8+3), n7(2n2+6n8+1).
    lower_ats=[{nn[2]:s.Rational(1,2),nn[8]:-2,nn[7]:0,nn[4]:-4*nn[6],nn[5]:0,nn[3]:s.Rational(27,5)*nn[6]**2},
      {nn[2]:s.Rational(1,2),nn[8]:s.Rational(-7,6),nn[7]:0,nn[6]:0,nn[4]:3*nn[5],nn[3]:s.Rational(9,5)*nn[5]**2},
      {nn[2]:s.Rational(-3,4),nn[8]:s.Rational(-3,4),nn[7]:0,nn[6]:0,nn[4]:nn[5],nn[3]:0}]
    lower=[]
    for at in lower_ats:
        _,C=branch(data['CNIL2_LOWER'],at)
        lower.append({'substitution':{str(k):str(v) for k,v in at.items()},'cubic':str(C)})
    assert s.expand(s.sympify(lower[0]['cubic'],locals=ENV)+(6*nn[6]*a+5*r)*(18*nn[6]*a*a+10*q+15*a*r)/120)==0
    assert s.rem(s.sympify(lower[1]['cubic'],locals=ENV),a,a)==0 and lower[2]['cubic']=='0'
    results['CNIL2_LOWER']={'all_three_branches':lower,'irreducible_cubic_support_empty':True}
    # Zero Cnil2: exact Jacobi coefficients imply n10=n5=0 on reduced support.
    _,_,_,J,C=data['CNIL2_ZERO'];assert any(s.expand(v-s.Rational(5,4)*nn[10])==0 for v in J)
    assert any(s.expand(v-nn[5]**2)==0 for v in J)
    Cz=s.factor(C.subs({nn[10]:0,nn[5]:0}));assert s.rem(Cz,a,a)==0
    results['CNIL2_ZERO']={'forced_reduced_equations':['n10=0','n5=0'],'remaining_cubic':str(Cz),'irreducible_cubic_support_empty':True}
    # Czero_zero: only two necessary exact Jacobi coefficients are needed.
    _,_,_,J,C=data['CZERO_ZERO'];assert nn[10]*nn[9] in J and nn[11]*nn[9] in J
    assert s.expand(C.subs(a,0)-b*b*(nn[10]*nn[9]*r+nn[11]*nn[9]*z)/5)==0
    results['CZERO_ZERO']={'necessary_jacobi':['n10*n9','n11*n9'],'cubic_mod_x0':str(s.factor(C.subs(a,0))),'irreducible_cubic_support_empty':True}
    # Remaining Czero_nonzero branches: n7!=0 and the exceptional n7=0 resonance.
    first={nn[8]:2,nn[2]:s.Rational(4,3),nn[4]:2*nn[6],nn[5]:-s.Rational(2,5)*nn[1]*nn[7],nn[3]:-s.Rational(3,5)*nn[0]*nn[7]}
    P1,C1=branch(data['CZERO_NONZERO'],first);base1=sub(P1,{nn[0]:0,nn[1]:0,nn[6]:0,nn[7]:1})
    A=-nn[0]/5+s.Rational(6,25)*nn[1]*nn[6];B=-nn[1]*nn[7]/5;C=-s.Rational(6,5)*nn[6]
    coords1=[a,b,c,r+A*a+B*z,nn[7]*z+C*a]
    assert identity(base1,P1,coords1)==nn[7]
    assert s.expand(s.diff(C1,r)-nn[7]*a*a/2)==0
    assert C1.subs(a,0)!=0
    families['Czero_nonzero_I']={'source':first,'free':['n0','n1','n6','n7'],'open':'n7 !=0','P':P1,'C':C1,'base':base1,'coords':coords1,'det':nn[7]}
    second={nn[8]:2,nn[7]:0,nn[2]:s.Rational(7,9),nn[5]:0,nn[4]:8*nn[6]}
    P2,C2=branch(data['CZERO_NONZERO'],second);base2=sub(P2,{nn[0]:0,nn[1]:0,nn[6]:0,nn[3]:1})
    A=-s.Rational(3,5)*nn[0]+s.Rational(54,25)*nn[1]*nn[6];B=-s.Rational(3,10)*nn[3]*nn[1];C=-s.Rational(18,5)*nn[6]
    coords2=[a,b,c,r+A*a+B*z,nn[3]*z+C*a]
    assert identity(base2,P2,coords2)==nn[3]
    assert s.expand(C2+(6*nn[3]*a**3+(18*nn[6]*a+5*z)*q)/36)==0
    families['Czero_nonzero_II']={'source':second,'free':['n0','n1','n3','n6'],'open':'n3 !=0','P':P2,'C':C2,'base':base2,'coords':coords2,'det':nn[3]}
    _,_,_,J,C=data['CZERO_NONZERO']
    # n7=0, n8!=2 forces n2=n8-3/2: cubic always has a factor x0.
    Cother=s.factor(C.subs({nn[7]:0,nn[2]:nn[8]-s.Rational(3,2)},simultaneous=True));assert s.rem(Cother,a,a)==0
    # n7=0,n8=2,n2=1/2: full branch is reducible.
    _,Chalf=branch(data['CZERO_NONZERO'],{nn[8]:2,nn[7]:0,nn[2]:s.Rational(1,2),nn[6]:0,nn[3]:s.Rational(3,5)*nn[4]*nn[5]})
    assert s.expand(Chalf+nn[5]*a*a*(2*nn[4]*a+5*z)/20)==0
    # On the complementary n2!=1/2 branch, n5=0; n3=0 is q-divisible.
    Cres=s.factor(C.subs({nn[8]:2,nn[7]:0,nn[5]:0,nn[3]:0}));assert s.div(Cres,q,*xx)[1]==0
    results['CZERO_NONZERO']={'n7_zero_n8_not_two_cubic':str(Cother),'half_resonance_cubic':str(Chalf),'n3_zero_remaining_cubic':str(Cres),'remaining_full_families':['Czero_nonzero_I','Czero_nonzero_II']}
    output={'task_id':'E_reduced_nonnormal','all_passed':True,'scope':'All five Anil3 pure nilpotent rank-three quadric cofactor sources reduce to exactly three full irreducible-cubic families, each a single projective orbit. Membership of those orbits is a separate obligation.',
      'original_source_verification_dependency':'reports/CLASSIFICATION_DEGENERATE_QUADRIC_SOURCES_REVIEW_20261107.json',
      'inputs_sha256':{str(v[0].relative_to(ROOT)):hashlib.sha256(v[0].read_bytes()).hexdigest() for v in data.values()},
      'producer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'empty_and_complement_checks':results,
      'families':{name:{'source_substitution':{str(k):str(v) for k,v in f['source'].items()},'free_parameters':f['free'],'open_condition':f['open'],'universal_brackets':brackets(f['P']),'actual_cubic':str(f['C']),'canonical_brackets':brackets(f['base']),'whole_family_coordinate_identity':[str(v) for v in f['coords']],'coordinate_determinant':str(f['det']),'all_original_jacobi_coefficients_zero':True,'whole_source_one_orbit':True} for name,f in families.items()},'seconds':time.monotonic()-start}
    path=ROOT/'reports/CLASSIFICATION_DEGENERATE_QUADRIC_Q3_ANIL3_SUPPORT_20261107.json'
    if path.exists():raise FileExistsError(path)
    path.write_text(json.dumps(output,indent=2)+'\n');print(json.dumps({'all_passed':True,'seconds':output['seconds'],'remaining_whole_orbits':3}))
if __name__=='__main__':main()
