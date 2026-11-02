"""Explicit rational squarezero graph and six-jet lifts of111 actual kernel generators."""
import json,os,subprocess,tempfile,time
from fractions import Fraction
from pathlib import Path
import sympy as s
from sympy.polys.rings import ring
from sympy.polys.domains import QQ
from s61_model import S61,X,J as J5
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from weight1_structure import COORDS,sp_basis,J
ROOT=Path(__file__).resolve().parents[1]

def source_geometry(initial):
    start=time.monotonic();m=S61();z,bu,beta=primitive_beta(m)
    point=list(map(s.Rational,initial['point61']))
    one=[i for i,w in enumerate(m.reduced['parameter_weights'])if w==1]
    tensor=s.Matrix([s.diff(beta[i,j],z[k])for (i,j),k in COORDS]);B=tensor.jacobian(bu)
    rr=list(B.T.to_DM().rref()[1]);left=B[rr,:].inv()
    b0=beta.subs(dict(zip(bu,[point[i]for i in one])))
    a,b,c,w=s.symbols('a b c w');sp=sp_basis();g=s.eye(4)
    for v,j in [(a,4),(b,6),(c,9),(w,5)]:
        assert sp[j]*sp[j]==s.zeros(4)
        g=g*(s.eye(4)+v*sp[j])
    assert g.T*J*g==J
    # The source jet for pushforward by g is g^* beta.
    source=(g.T*b0.subs(dict(zip(z,g*s.Matrix(z))),simultaneous=True)*g).applyfunc(s.expand)
    tv=s.Matrix([s.diff(source[i,j],z[k])for (i,j),k in COORDS])
    up=(left*tv[rr,:]).applyfunc(s.expand)
    assert B*up==tv
    u4=up[list(map(str,bu)).index('u4')]
    u5=up[list(map(str,bu)).index('u5')]
    assert s.expand(s.diff(u4,w)+u5)==0 and not u5.has(w)
    record=dict(group= [[str(v)for v in row]for row in g.tolist()],parameters=['a','b','c','w'],source_jet_variables=list(map(str,bu)),source_jet=list(map(str,up)),u4=s.factor(u4).__str__(),degree_w=int(s.degree(u4,w)),seconds=time.monotonic()-start)
    return record

def quadratic_representation():
    start=time.monotonic()
    m=S61();two=[i for i,w in enumerate(m.reduced['parameter_weights'])if w==2]
    q=s.symbols('q:45');p=[s.S.Zero]*61
    for i,v in zip(two,q):p[i]=v
    P=m.reconstruct(p)-X[0]**2*J5
    records={};exps={};a,b,c,w=s.symbols('a b c w')
    for j in (4,6,9,5):
        A=s.zeros(5);A[1:,1:]=sp_basis()[j]
        AQ=A*s.Matrix(X)
        dP=A*P+P*A.T-sum((AQ[i]*P.diff(X[i])for i in range(5)),s.zeros(5))
        out=m.normalize(X[0]**2*J5+dP)
        D=s.Matrix([out['parameters'][i]for i in two]).jacobian(q)
        assert not D.free_symbols
        powers=[s.eye(45)];k=1
        while True:
            v=powers[-1]*D
            if v==s.zeros(45):break
            assert k<=5;powers.append(v);k+=1
        records[str(j)]=dict(matrix=[[str(v)for v in row]for row in D.tolist()],nilpotence_index=k)
        z={4:a,6:b,9:c,5:w}[j]
        exps[j]=sum((z**i/s.factorial(i)*M for i,M in enumerate(powers)),s.zeros(45))
        print('generator',j,'nilindex',k,'seconds',time.monotonic()-start,flush=True)
    free=[19,26,27,28,33,34,35,37,38];rows=free+[25]
    R=s.eye(45)[rows,:]
    for j in (4,6,9,5): R=(R*exps[j]).applyfunc(s.expand)
    aa,bb,cc,ww=a,b,c,w
    expected=[-cc,2*aa-2*bb*ww,2*ww,2*bb-2*cc*ww,-ww*(2*bb-2*cc*ww),-(2*bb-2*cc*ww)**2/4,-bb,2*bb*ww,bb*(2*bb-2*cc*ww),-1]
    assert all(s.expand(-R[i,25]-expected[i])==0 for i in range(10))
    result=dict(generators=records,selected_rows=rows,selected_finite_action=[[str(v)for v in row]for row in R.tolist()],seconds=time.monotonic()-start)
    return result

def kernel_candidates(initial):
    started=time.monotonic();saved=initial
    variables=[2,3,4,5,6,7,8];idx={v:i for i,v in enumerate(variables)}
    polys=[]
    for row in saved['initial_relations']:
        terms=[]
        for mon,c in row:
            if any(i in (0,1)for i in mon):continue
            terms.append('('+c+')*'+'*'.join('x_'+str(idx[i])for i in mon))
        if terms:polys.append('+'.join(terms))
    header='R=QQ[x_0..x_6];J=ideal('+','.join(polys)+');G=gb J;\n'
    forms=['x_0','x_1+2*x_4','x_3','x_6']
    generators=['x_2','x_5']+['('+a+')*('+b+')'for i,a in enumerate(forms)for b in forms[i:]]
    gen_degrees=[1,1]+[2]*10
    env=os.environ.copy();env['DYLD_LIBRARY_PATH']=str(ROOT/'work/m2-runtime/lib')+':/opt/homebrew/lib:/opt/homebrew/opt/givaro/lib'
    chosen=[];by_degree=[];logs=[]
    with tempfile.TemporaryDirectory(prefix='zero12-kernel-generators-')as temp:
        td=Path(temp)
        for degree in range(1,8):
            cand=[(i,())for i,d in enumerate(gen_degrees)if d==degree]
            if by_degree:
                cand += [(g,tuple(sorted(m+(j,))))for g,m in by_degree[-1]for j in range(7)]
            cand=list(dict.fromkeys(cand));expr=[]
            for g,m in cand:
                expr.append('('+generators[g]+')'+''.join('*x_'+str(j)for j in m))
            output=td/'normalforms.txt'
            program=header+'f=openOut "'+str(output)+'";\n'
            for i,e in enumerate(expr):
                program+='a=listForm(('+e+')%G); f << toString(#a) << " "; scan(a,z->{f << toString(last z) << " "; scan(first z,v->f << toString(v) << " ");}); f << endl;\n'
            program+='close f;\n'
            source=td/'verify.m2';source.write_text(program)
            run=subprocess.run(['/opt/homebrew/Cellar/macaulay2/1.26.05/bin/M2-binary','--script',str(source)],env=env,capture_output=True,text=True,timeout=30)
            assert run.returncode==0,run.stdout+run.stderr
            lines=output.read_text().splitlines();assert len(lines)==len(cand)
            echelon={};selected=[]
            for candidate,line in zip(cand,lines):
                a=line.split();count=int(a[0]);assert len(a)==1+8*count
                v={tuple(map(int,a[2+8*i:9+8*i])):Fraction(a[1+8*i])for i in range(count)}
                for lead,row in sorted(echelon.items()):
                    if lead in v:
                        c=v[lead]
                        for m,z in row.items():
                            v[m]=v.get(m,Fraction(0))-c*z
                            if not v[m]:del v[m]
                if v:
                    lead=min(v);c=v[lead];echelon[lead]={m:z/c for m,z in v.items()};selected.append(candidate)
            by_degree.append(selected);chosen+=selected
            logs.append(dict(degree=degree,candidates=len(cand),rank=len(selected)))
            print(logs[-1],flush=True)
        assert [len(r)for r in by_degree]==[2,15,28,38,18,10,0]
        assert len(chosen)==111
        records=[dict(B_generator=g,multiplier=[variables[j]for j in m],degree=len(m)+gen_degrees[g])for g,m in chosen]
        report=dict(scope='Polynomial generators of initial E; actual lifts require the formal graph equations of D.',
                    base_variables=[0,1],B_generator_formulas=['v4','v7']+['('+a+')*('+b+')'for i,a in enumerate(['v2','v3+2*v6','v5','v8'])for b in ['v2','v3+2*v6','v5','v8'][i:]],
                    generators=records,degree_counts=[len(r)for r in by_degree[:-1]],checks=logs,seconds=time.monotonic()-started)
        return report

def actual_lifts(source,action,normal,candidates):
    started=time.monotonic()
    a,b,c,w=s.symbols('a b c w');v=s.symbols('v:9')
    R,A,B,C=ring('a,b,c',QQ);degree=5

    def trunc(p,n=degree):return R.from_dict({m:z for m,z in p.items()if sum(m)<=n})
    def mul(p,q):return trunc(p*q)
    def power(p,n):
        if n<0:return power(inverse(p),-n)
        ans=R.one
        while n:
            if n&1:ans=mul(ans,p)
            p=mul(p,p);n//=2
        return ans
    def inverse(p):
        c=p.get((0,0,0),QQ.zero);assert c
        q=R.one-p/c;term=R.one;out=term
        for _ in range(degree):term=mul(term,q);out+=term
        assert mul(p,out/c)==R.one
        return out/c
    def evaluate(expr,values):
        expr=s.sympify(expr)
        if expr.is_Number:return R.ground_new(QQ.from_sympy(expr))
        if expr.is_Symbol:return values[expr]
        if expr.is_Add:return sum((evaluate(z,values)for z in expr.args),R.zero)
        if expr.is_Mul:
            out=R.one
            for z in expr.args:out=mul(out,evaluate(z,values))
            return out
        if expr.is_Pow:return power(evaluate(expr.base,values),int(expr.exp))
        raise ValueError(expr)
    up=list(map(s.sympify,source['source_jet']));ui=list(map(s.Symbol,source['source_jet_variables']))
    fsource=s.sympify(source['u4']);den=-s.diff(fsource,w)
    assert not den.has(w)
    values={a:A,b:B,c:C}
    W=mul(evaluate(fsource.subs(w,0),values),inverse(evaluate(den,values)))
    values[w]=W
    uv={x:evaluate(z,values)for x,z in zip(ui,up)}
    assert uv[s.Symbol('u4')]==0
    print('source',time.monotonic()-started,flush=True)
    N=[[evaluate(z,uv)for z in row]for row in normal['squarezero_normal_directions']]
    Rm=[[evaluate(z,values)for z in row]for row in action['selected_finite_action']]
    Np=[[sum((mul(Rm[i][k],N[k][j])for k in range(45)),R.zero)for j in range(4)]for i in range(10)]
    vc=[-C,2*A-2*mul(B,W),2*W,2*B-2*mul(C,W),R.zero,R.zero,-B,2*mul(B,W),R.zero]
    vc[4]=-mul(W,vc[3]);vc[5]=-mul(vc[3],vc[3])/4;vc[8]=-mul(vc[6],vc[3])
    assert all(vc[i]==-Rm[i][25] for i in range(9))
    Nh=[[Np[i][j]+mul(vc[i],Np[9][j])for j in range(4)]for i in range(9)]
    ai=(v[1]-v[6]*v[2])/2;bi=-v[6];ci=-v[0];wi=v[2]/2
    fs=[v[3]+2*v[6]-v[0]*v[2],2*v[4]+v[2]*v[3],4*v[5]+v[3]**2,v[7]+v[6]*v[2],v[8]+v[6]*v[3],s.expand(fsource.subs({a:ai,b:bi,c:ci,w:wi},simultaneous=True))]
    vv=dict(zip(v,vc));assert all(evaluate(f,vv)==0 for f in fs)
    DM=[[sum((mul(evaluate(s.diff(f,v[k]),vv),Nh[k][j])for k in range(9)),R.zero)for j in range(4)]for f in fs]
    sel=[0,2,4,5];mat=[DM[k]for k in sel]
    M0=s.Matrix([[s.Rational(z.get((0,0,0),0))for z in row]for row in mat]);M0i=M0.inv()
    MI0=[[R.ground_new(QQ.from_sympy(z))for z in row]for row in M0i.tolist()]
    def mm(X,Y):return [[sum((mul(X[i][k],Y[k][j])for k in range(len(Y))),R.zero)for j in range(len(Y[0]))]for i in range(len(X))]
    I=[[R.one if i==j else R.zero for j in range(4)]for i in range(4)]
    NM=mm(mat,MI0);Z=[[I[i][j]-NM[i][j]for j in range(4)]for i in range(4)]
    term=I;geom=[row[:]for row in I]
    for _ in range(degree):
        term=mm(term,Z)
        geom=[[geom[i][j]+term[i][j]for j in range(4)]for i in range(4)]
    MI=mm(MI0,geom);assert mm(mat,MI)==I and mm(MI,mat)==I
    alpha=mm([DM[1]],MI)[0];gamma=mm([DM[3]],MI)[0]
    assert mm([alpha],mat)==[DM[1]] and mm([gamma],mat)==[DM[3]]
    assert all(p.get((0,0,0),0)==0 for p in alpha+gamma)
    print('matrix inverse',M0.det(),time.monotonic()-started,flush=True)
    # Compose coefficient jets with the flag big-cell inverse. All discarded terms
    # have order>=6, and multiply one f_i of order>=1.
    S,*V=ring(','.join('v'+str(i)for i in range(9)),QQ)
    def tr9(p,n=6):return S.from_dict({m:z for m,z in p.items()if sum(m)<=n})
    def mu9(p,q):return tr9(p*q)
    sub=[(V[1]-V[6]*V[2])/2,-V[6],-V[0]]
    def compose(p):
        out=S.zero
        for mon,coef in p.items():
            z=S.ground_new(coef)
            for t,n in zip(sub,mon):
                for _ in range(n):z=mu9(z,t)
            out+=z
        return out
    F=[S.from_expr(f)for f in fs]
    G4=tr9((F[1]-sum((mu9(compose(alpha[j]),F[sel[j]])for j in range(4)),S.zero))/2)
    G7=tr9(F[3]-sum((mu9(compose(gamma[j]),F[sel[j]])for j in range(4)),S.zero))
    H=[-2*F[5],F[0],F[2]/4,F[4]]
    coordinate_expressions=[ai,bi,ci]+[z.as_expr()for z in H]+[G4.as_expr(),G7.as_expr()]
    coordinate_jacobian=s.Matrix(coordinate_expressions).jacobian(v).subs(dict.fromkeys(v,s.S.Zero))
    assert coordinate_jacobian.det()!=0
    gens=[G4,G7]+[mu9(x,y)for i,x in enumerate(H)for y in H[i:]]
    initial=[V[4],V[7]]+[x*y for i,x in enumerate([V[2],V[3]+2*V[6],V[5],V[8]])for y in [V[2],V[3]+2*V[6],V[5],V[8]][i:]]
    for i,(g,h)in enumerate(zip(gens,initial)):
        n=1 if i<2 else 2
        assert S.from_dict({m:z for m,z in g.items()if sum(m)==n})==h
    # Verify entire value and first-order normal derivative of each jet modulo
    # degree6/value and degree5/derivative in the source parameters.
    for g in gens:
        assert evaluate(g.as_expr(),vv)==0
        for j in range(4):assert sum((mul(evaluate(g.diff(V[k]).as_expr(),vv),Nh[k][j])for k in range(9)),R.zero)==0
    cand=candidates['generators'];lifts=[]
    for rec in cand:
        p=gens[rec['B_generator']]
        for j in rec['multiplier']:p=mu9(p,V[j])
        lifts.append(p)
    def enc(p):return [[list(m),str(z)]for m,z in sorted(p.items())]
    report=dict(actual_lifts_defined_by_exact_rational_conormal_formula=True,reduced_curve_equations=list(map(str,fs)),bigcell_inverse=list(map(str,[ai,bi,ci,wi])),coordinate_jacobian_determinant=str(coordinate_jacobian.det()),jet_order=6,conormal_selected_rows=sel,constant_conormal_determinant=str(M0.det()),alpha_sixjet_coefficient_series=[enc(p)for p in alpha],gamma_sixjet_coefficient_series=[enc(p)for p in gamma],source_parameters=['a','b','c'],ambient_variables=list(map(str,v)),kernel_generators_sixjet=[enc(p)for p in gens],kernel_111_lifts_sixjet=[enc(p)for p in lifts],generator_metadata=cand,all_values_and_normal_derivatives_verified=True,seconds=time.monotonic()-started)
    return report
