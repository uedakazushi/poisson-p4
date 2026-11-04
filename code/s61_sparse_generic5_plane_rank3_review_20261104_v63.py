#!/usr/bin/env python3
"""Independent sparse5 source and fixed-operation certificate audit, no search."""
import argparse,gzip,hashlib,itertools,json,math,resource,signal,subprocess,time
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from pathlib import Path
import sympy as s

ROOT=Path(__file__).resolve().parents[1]
ID='S61_SPARSE_GENERIC5_PLANE_RANK3_CIRCUIT_20261104_v63_REVIEW'
STEMS=['S61_SPARSE_GENERIC5_INPUT_PREFLIGHT_20261104_v63','S61_SPARSE_GENERIC5_CONSTANT_UNIT_PREFIX_20261104_v63','S61_SPARSE_GENERIC5_SMALL_TAIL_FIELD_20261104_v63','S61_SPARSE_GENERIC5_PLANE_RANK3_CIRCUIT_20261104_v63']
ZERO=(0,)*5
def read(p):
    with gzip.open(p,'rt') if str(p).endswith('.gz') else open(p) as f:return json.load(f)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def add(a,b):
    c=dict(a)
    for m,v in b.items():c[m]=c.get(m,F(0))+v
    return {m:v for m,v in c.items() if v}
def scale(a,c):return {m:v*c for m,v in a.items() if v*c}
def mul(a,b):
    c={}
    for m,v in a.items():
        for n,w in b.items():
            e=tuple(x+y for x,y in zip(m,n));c[e]=c.get(e,F(0))+v*w
    return {m:v for m,v in c.items() if v}
def dec(a):
    if isinstance(a,dict):return {tuple(t[1:]):F(t[0])*F(a['content']) for t in a['terms'] if F(t[0])*F(a['content'])}
    if not a:return {}
    if isinstance(a[0][0],list):return {tuple(e):F(c) for e,c in a if F(c)}
    return {tuple(t[1:]):F(t[0]) for t in a if F(t[0])}
def eqmatrix(a,b):assert a==b

class Program:
    """Literal saved field operations, plus independent source/kernel identities."""
    def __init__(self):self.lines=[];self.n=0;self.memo={};self.equalities=0
    def register(self,line):
        i=self.n;self.n+=1;self.lines.append(line.format(i=i));return i
    def rf(self,n,d):
        key=(tuple(sorted(n.items())),tuple(sorted(d.items())))
        if key not in self.memo:
            def txt(p):return ' '.join([str(len(p)),*(' '.join(map(str,(c,*e))) for e,c in sorted(p.items()))])
            self.memo[key]=self.register('L {i} '+txt(n)+' '+txt(d))
        return self.memo[key]
    def p(self,a):
        den=math.lcm(*(v.denominator for v in a.values())) if a else 1
        return self.rf({e:int(v*den) for e,v in a.items()},{ZERO:den})
    def saved(self,a):return self.rf({tuple(t[1:]):int(t[0]) for t in a['numerator']},{tuple(t[1:]):int(t[0]) for t in a['denominator']})
    def op(self,o,a,b):
        if o in ['A','M'] and b<a:a,b=b,a
        key=(o,a,b)
        if key not in self.memo:self.memo[key]=self.register(f'{o} {{i}} {a} {b}')
        return self.memo[key]
    def total(self,a):
        z=self.p({})
        for v in a:z=self.op('A',z,v)
        return z
    def eq(self,a,b):self.lines.append(f'E {a} {b}');self.equalities+=1
    def nz(self,a):self.lines.append(f'N {a}')
    def write(self,p):
        # Drop registers after last use, keeping the checker far below its cap.
        last={}
        for k,line in enumerate(self.lines):
            a=line.split();o=a[0]
            if o in ['A','S','M','D']:refs=list(map(int,a[2:]))
            elif o=='E':refs=list(map(int,a[1:]))
            elif o=='N':refs=[int(a[1])]
            else:refs=[]
            for r in refs:last[r]=k
        kill={}
        for r,k in last.items():kill.setdefault(k,[]).append(r)
        with gzip.open(p,'wt',compresslevel=6) as f:
            for k,line in enumerate(self.lines):
                f.write(line+'\n')
                for r in kill.get(k,[]):f.write(f'X {r}\n')

def source_audit(dirs):
    from s61_model import S61
    model=S61();cert=read(dirs[0]/'certificate.json');known=read(dirs[0]/'known-five-input.json.gz');pulled=read(dirs[0]/'pulled134-input.json.gz')
    weights=model.reduced['parameter_weights'];one=[i for i,w in enumerate(weights) if w==1];two=[i for i,w in enumerate(weights) if w==2];assert one==cert['historical_weight1_indices'] and two==cert['historical_weight2_indices']
    pairs=list(itertools.combinations(range(5),2));mons=list(itertools.combinations_with_replacement(range(5),2));z=s.symbols('z0:4');U=s.symbols('u0:16');J=s.zeros(4);J[0,1]=J[2,3]=1;J[1,0]=J[3,2]=-1
    raw=s.zeros(150,16)
    for c,i in enumerate(one):raw[model.certificate['raw_free_indices'][i],c]=1
    for row in model.data['eliminated_raw_variables']:
        for mon,n,d in row['polynomial']:
            if len(mon)==1 and mon[0] in one:raw[row['raw_index'],one.index(mon[0])]+=s.Rational(n,d)
    assert raw==s.Matrix(cert['original_raw_weight1_columns150x16'])
    P=s.zeros(5)
    for k,(ij,mm) in enumerate(itertools.product(pairs,mons)):
        i,j=ij;mon=s.prod(1 if x==0 else z[x-1] for x in mm);P[i,j]+=sum(raw[k,c]*U[c] for c in range(16))*mon
    for i,j in pairs:P[j,i]=-P[i,j]
    affine=s.Matrix(4,4,lambda i,j:s.expand(P[i+1,j+1]-z[i]*P[0,j+1]+z[j]*P[0,i+1]))
    linear=affine.applyfunc(lambda x:s.Poly(x,*z).homogeneous_component(1) if False else sum(s.diff(x,y).subs(dict.fromkeys(z,0))*y for y in z))
    alpha=s.Matrix([s.diff(linear[0,1]+linear[2,3],y)/5 for y in z]);assert alpha.jacobian(U)==s.Matrix(cert['Lee_alpha4x16'])
    v=J*alpha;zz=s.Matrix(z);beta=(-J*(linear-2*alpha.dot(zz)*J+zz*v.T-v*zz.T)*J).applyfunc(s.expand)
    coordinates=[((i,j),k) for i,j in itertools.combinations(range(4),2) for k in range(4)]
    M=s.Matrix([s.diff(beta[i,j],z[k]) for (i,j),k in coordinates]).jacobian(U);assert M==s.Matrix(cert['primitive_projection24x16'])
    K=s.Matrix(cert['primitive24_basis']);A=s.Matrix(cert['adapter_b16_equals_A_u16']);Ai=s.Matrix(cert['adapter_u16_equals_Ainv_b16']);assert K.rank()==16 and K*A==M and A*Ai==s.eye(16) and A.det()==s.Rational(36,625)
    # Basis independently satisfies all four closed equations and four traces.
    for c in range(16):
        B=s.zeros(4)
        for r,((i,j),k) in enumerate(coordinates):B[i,j]+=K[r,c]*z[k];B[j,i]-=K[r,c]*z[k]
        assert B[0,1]+B[2,3]==0
        for i,j,k in itertools.combinations(range(4),3):assert s.diff(B[j,k],z[i])-s.diff(B[i,k],z[j])+s.diff(B[i,j],z[k])==0
    H,L,ss,qq,tt=s.symbols('H L s0 q0 t0');pars=[H,L,ss,qq,tt];b=s.zeros(16,1)
    for i,x in {0:-45*ss/16,1:-30*L,2:2400,4:tt,5:H*(5*qq+ss+1)/240,6:ss/54000,8:-32*H/3,9:-s.Rational(32,675),14:H*(ss+1)/120,15:-64*H/3}.items():b[i]=x
    assert (b-s.Matrix(cert['primitive_beta16'])).applyfunc(s.expand)==s.zeros(16,1)
    us=Ai*b;up=[{e:F(str(c)) for e,c in s.Poly(x,*pars).terms() if c} for x in us]
    assert sum(bool(x) for x in up)==9 and sum(map(len,up))==10
    assert all(s.expand(x-s.sympify(y))==0 for x,y in zip(us,cert['source_u16']))
    u={i:{} for i in range(61)};u.update(dict(zip(one,up)));twoidx={i:j for j,i in enumerate(two)}
    @lru_cache(None)
    def mono(mon):
        return {ZERO:F(1)} if not mon else mul(u[mon[0]],mono(mon[1:]))
    count=0
    for row,eq in enumerate(model.reduced['equations']):
        f={};ww={sum(weights[i] for i in mon) for mon,n,d in eq};assert len(ww)==1 and next(iter(ww)) in [4,5]
        for mon,n,d in eq:
            qm=tuple(twoidx[i] for i in mon if i in twoidx);bm=tuple(i for i in mon if i not in twoidx)
            for e,c in mono(bm).items():key=(qm,e);f[key]=f.get(key,F(0))+F(n,d)*c
        f={k:v for k,v in f.items() if v};saved={(tuple(qm),tuple(e)):F(c) for qm,e,c in pulled['equations'][row]};assert f==saved;count+=len(f)
        lin=[{} for _ in range(45)]
        for (qm,e),c in f.items():
            if len(qm)==1:lin[qm[0]][e]=c
        assert lin==[dec(x) for x in pulled['linear_q_matrix134x45'][row]]
    assert Counter(next(iter({sum(weights[i] for i in mon) for mon,n,d in eq})) for eq in model.reduced['equations'])=={4:70,5:64}
    # All150 reconstruction terms have the physical Gm weights.
    rawweights=[]
    for (i,j),mm in itertools.product(pairs,mons):rawweights.append(2-int(i>0)-int(j>0)+sum(x>0 for x in mm))
    for k,i in zip(model.certificate['raw_free_indices'],range(61)):assert rawweights[k]==weights[i]
    for row in model.data['eliminated_raw_variables']:assert all(sum(weights[i] for i in mon)==rawweights[row['raw_index']] for mon,n,d in row['polynomial'])
    historical=s.symbols(' '.join(cert['historical_weight1_names']))
    def substitute(text):
        out={}
        for ex,c in s.Poly(s.sympify(text),*historical).terms():
            mon=tuple(i for i,p in zip(one,ex) for _ in range(p));out=add(out,scale(mono(mon),F(str(c))))
        return out
    p51=read(ROOT/'reports/Q51_UNIVERSAL_SECTION_20261101.json');p52=read(ROOT/'reports/Q52_UNIVERSAL_SECTION_20261101.json');plane=read(ROOT/'reports/TORIC_UNIVERSAL_THREE_PLANE_20261101.json')
    for key,source in [('Q51',p51['quadratic_section_q']),('Q52',p52['quadratic_section_q']),('q_center',plane['q_center'])]:assert [substitute(x) for x in source]==[dec(x) for x in known[key]]
    for key,source in [('L51x20',plane['all51_homogeneous_linear_consequences']),('N45x20',plane['q_linear_matrix_N'])]:assert [[substitute(x) for x in row] for row in source]==[[dec(x) for x in row] for row in known[key]]
    for a,b,c,d in zip(known['Q51'],known['Q52'],known['q_center'],known['Q52_minus_center']):assert scale(dec(c),3)==add(scale(dec(a),2),dec(b)) and dec(d)==add(dec(b),scale(dec(c),-1))
    return known,dict(adapter='36/625',raw150_weight1_columns=2400,primitive24_by16=384,pulled134_terms=count,known_plane_L_N_entries=51*20+45*20,weighted134=[70,64])

def field_audit(dirs,known,P):
    z=P.p({});one=P.p({ZERO:F(1)});minus=P.p({ZERO:F(-1)})
    eq=P.eq;op=P.op
    L=[[P.p(dec(x)) for x in r] for r in known['L51x20']];A=[r[:] for r in L];rr=list(range(51));cc=list(range(20));prefix=read(dirs[1]/'constant-prefix.json')
    for k,step in enumerate(prefix['elementary_operations']):
        assert step['stage']==k;pi=step['swap_row_with'];pj=step['swap_column_with'];A[pi],A[k]=A[k],A[pi];rr[pi],rr[k]=rr[k],rr[pi]
        for r in A:r[pj],r[k]=r[k],r[pj]
        cc[pj],cc[k]=cc[k],cc[pj];pivot=P.p({ZERO:F(step['pivot_QQ'])});P.nz(pivot);eq(A[k][k],pivot)
        for j in range(k,20):A[k][j]=op('D',A[k][j],pivot)
        for i in range(51):
            m=P.p(dec(step['elimination_multipliers51'][i]));eq(m,z if i==k else A[i][k])
            if i!=k:
                for j in range(k+1,20):A[i][j]=op('S',A[i][j],op('M',m,A[k][j]))
                A[i][k]=z
    assert cc==prefix['column_order'] and rr==prefix['row_order']
    for i in range(51):
        for j in range(20):eq(A[i][j],P.p(dec(prefix['reduced51x20'][i][j])))
    primitive=read(dirs[2]/'row-primitive.json');assert not primitive['literal_exact_division_and_remultiplication_count']
    tail=[[P.p(dec(x)) for x in r] for r in primitive['primitive40x9']]
    for i in range(40):
        for j in range(9):eq(tail[i][j],A[11+i][11+j])
    tc=list(range(9));tr=list(range(40))
    for k in range(7):
        step=read(dirs[2]/f'field-op-{k}.json');pi=step['swap_row_with'];pj=step['swap_column_with'];tail[pi],tail[k]=tail[k],tail[pi];tr[pi],tr[k]=tr[k],tr[pi]
        for r in tail:r[pj],r[k]=r[k],r[pj]
        tc[pj],tc[k]=tc[k],tc[pj];pivot=P.saved(step['pivot']);P.nz(pivot);eq(tail[k][k],pivot)
        normal=[P.saved(x) for x in step['normalized_row9']]
        for j in range(9):eq(tail[k][j],op('M',pivot,normal[j]) if j>=k else normal[j])
        tail[k]=normal
        for i in range(40):
            m=P.saved(step['elimination_multipliers40'][i]);eq(m,z if i==k else tail[i][k])
            if i!=k:
                for j in range(k+1,9):tail[i][j]=op('S',tail[i][j],op('M',m,normal[j]))
                tail[i][k]=z
    checkpoint=read(dirs[2]/'field-prefix-checkpoint.json');assert tc==checkpoint['column_order'] and tr==checkpoint['row_order']
    for i in range(40):
        for j in range(9):eq(tail[i][j],P.saved(checkpoint['rref40x9'][i][j]));eq(tail[i][j],z) if i>=7 else None
    kernel=read(dirs[2]/'generic-field-kernel.json');T=[[P.saved(x) for x in r] for r in kernel['T20x2_rational']];free=kernel['free_original_columns'];assert free==[cc[11+tc[7]],cc[11+tc[8]]]
    for i in range(51):
        for j in range(2):eq(P.total(op('M',L[i][h],T[h][j]) for h in range(20)),z)
    for i in range(2):
        for j in range(2):eq(T[free[i]][j],one if i==j else z)
    old=read(ROOT/'reports/S61_SPARSE_GENERIC5_PLANE_ADAPTER_20261104_v63/kernel20x2.json');oldT=[[P.p(dec(x)) for x in r] for r in old['T_numerators20x2']];oldD=P.p(dec(old['Delta']));P.nz(oldD)
    for i in range(20):
        for j in range(2):
            eq(P.total(op('M',T[i][h],oldT[free[h]][j]) for h in range(2)),oldT[i][j])
            eq(P.total(op('M',oldT[i][h],T[old['free_columns'][h]][j]) for h in range(2)),op('M',oldD,T[i][j]))
    N=[[P.p(dec(x)) for x in r] for r in known['N45x20']];V=[[P.saved(x) for x in r] for r in read(dirs[2]/'plane-rational45x3.json')['V45x3']]
    shared=read(dirs[2]/'shared-denominator-kernel-plane.json.gz');den=P.p(dec(shared['common_denominator_P']));P.nz(den)
    for i in range(20):
        for j in range(2):eq(op('M',T[i][j],den),P.p(dec(shared['kernel20x2_numerators'][i][j])))
    W0=read(dirs[3]/'plane-W-original.json.gz');eq(den,P.p(dec(W0['common_P'])));W=[[P.p(dec(x)) for x in r] for r in W0['W45x3']]
    for i in range(45):
        for j in range(2):eq(P.total(op('M',N[i][h],T[h][j]) for h in range(20)),V[i][j])
        eq(V[i][2],P.p(dec(known['Q52_minus_center'][i])));eq(P.p(dec(W0['center45'][i])),P.p(dec(known['q_center'][i])))
        for j in range(3):
            eq(W[i][j],op('M',den,V[i][j]));mixed=P.p(dec(shared['plane45x3_mixed_numerators'][i][j]));eq(W[i][j],mixed if j<2 else op('M',den,mixed))
    prim=read(dirs[3]/'row-primitive45.json.gz');B=[[P.p(dec(x)) for x in r] for r in prim['quotients45x3']];factors=[P.p(dec(x)) for x in prim['factors45']]
    for i in range(45):
        P.nz(factors[i])
        for j in range(3):eq(W[i][j],op('M',factors[i],B[i][j]))
    br=list(range(45));bc=list(range(3));piv=[]
    for k in range(3):
        step=read(dirs[3]/f'field-op-{k}.json.gz');pi=step['swap_row_with'];pj=step['swap_column_with'];B[pi],B[k]=B[k],B[pi];br[pi],br[k]=br[k],br[pi]
        for r in B:r[pj],r[k]=r[k],r[pj]
        bc[pj],bc[k]=bc[k],bc[pj];pivot=P.saved(step['pivot']);piv.append(pivot);P.nz(pivot);eq(B[k][k],pivot);normal=[P.saved(x) for x in step['normalized_row3']]
        for j in range(3):eq(B[k][j],op('M',pivot,normal[j]) if j>=k else normal[j])
        B[k]=normal
        for i in range(45):
            m=P.saved(step['multipliers45'][i]);eq(m,B[i][k] if i>k else z)
            if i>k:
                for j in range(k+1,3):B[i][j]=op('S',B[i][j],op('M',m,normal[j]))
                B[i][k]=z
    rank=read(dirs[3]/'rank-certificate.json.gz');assert br[:3]==rank['selected_original_q_rows']==[3,6,13] and bc==rank['column_order']==[2,1,0]
    assert rank['column_permutation_sign']==-1
    for i in range(3):eq(piv[i],P.saved(rank['nonzero_pivot_factors'][i]));eq(factors[br[i]],P.p(dec(rank['selected_row_factors'][i])))
    return dict(constant_pivots=11,tail_pivots=7,rank_L=18,kernel_dimension=2,old_span_both_directions=80,plane_rows=br[:3],column_order=bc,minor_sign=-1)

def slp_audit(path):
    d=read(path);nodes=d['nodes'];vals=[]
    # Independent commutative sparse polynomial evaluator on abstract leaves.
    def a(x,y):return add(x,y)
    for i,node in enumerate(nodes):
        o,*v=node
        if o=='qq':x={():F(v[0])} if F(v[0]) else {}
        elif o in ['original_W_poly','center','z','n','selected_row_factor_poly','rank_pivot_rf']:x={((o,*v),):F(1)}
        elif o=='inv':
            assert nodes[v[0]][0] in ['selected_row_factor_poly','rank_pivot_rf'];x={(('inv',*nodes[v[0]]),):F(1)}
        else:
            assert all(k<i for k in v)
            if o in ['add','sub']:x=a(vals[v[0]],scale(vals[v[1]],-1 if o=='sub' else 1))
            else:
                assert o=='mul';x={}
                for m,c in vals[v[0]].items():
                    for n,e in vals[v[1]].items():key=tuple(sorted(m+n));x[key]=x.get(key,F(0))+c*e
                x={m:c for m,c in x.items() if c}
        vals.append(x)
    def leaf(o,*v):return {((o,*v),):F(1)}
    def times(x,y):
        out={}
        for m,c in x.items():
            for n,e in y.items():key=tuple(sorted(m+n));out[key]=out.get(key,F(0))+c*e
        return {m:c for m,c in out.items() if c}
    rows=d['pivot_q_rows'];B=[[leaf('original_W_poly',i,j) for j in range(3)] for i in rows];g={}
    for perm in itertools.permutations(range(3)):
        t={():F((-1)**sum(perm[i]>perm[j] for i in range(3) for j in range(i+1,3)))}
        for i in range(3):t=times(t,B[i][perm[i]])
        g=a(g,t)
    assert vals[d['Gamma3_node']]==g
    inv={():F(-1)}
    for o in ['selected_row_factor_poly','rank_pivot_rf']:
        for j in range(3):inv=times(inv,{(('inv',o,j),):F(1)})
    assert vals[d['inverse_Gamma3_node']]==inv
    adj=[[None]*3 for _ in range(3)]
    for i in range(3):
        for j in range(3):
            r=[k for k in range(3) if k!=j];c=[k for k in range(3) if k!=i];adj[i][j]=scale(a(times(B[r[0]][c[0]],B[r[1]][c[1]]),scale(times(B[r[0]][c[1]],B[r[1]][c[0]]),-1)),(-1)**(i+j));assert vals[d['adjugate3x3_nodes'][i][j]]==adj[i][j]
    for i in range(3):
        for j in range(3):
            x={}
            for k in range(3):x=a(x,times(B[i][k],adj[k][j]))
            assert x==(g if i==j else {})
            x={}
            for k in range(3):x=a(x,times(adj[i][k],B[k][j]))
            assert x==(g if i==j else {})
    normal=[i for i in range(45) if i not in rows];assert d['normal_q_rows42']==normal
    for i in range(45):
        E=[]
        for j in range(3):
            x={}
            for k in range(3):x=a(x,times(leaf('original_W_poly',i,k),adj[k][j]))
            E.append(times(inv,x));assert E[-1]==vals[d['E45x3_nodes'][i][j]]
        q=leaf('center',i)
        for j in range(3):q=a(q,times(E[j],leaf('z',j)))
        if i in normal:q=a(q,leaf('n',normal.index(i)))
        assert q==vals[d['q45_nodes'][i]]
    from s61_model import S61
    two=[i for i,w in enumerate(S61().reduced['parameter_weights']) if w==2]
    assert d['q_to_u_indices_zero_based']==two
    assert d['pivot_u_numbers_one_based']==[two[i]+1 for i in rows]==[6,9,19] and len(nodes)==1608
    return dict(nodes=len(nodes),all45_forward_coordinate_identities=True,abstract_B_adjugate_identity=True,inverse_affine_map_proved=True,no_unknown_inverse=True)

def main(out,work,verified_checkpoint=None):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False);work.mkdir(parents=True,exist_ok=True)
    dirs=[ROOT/'reports'/x for x in STEMS];hashes={}
    for d in dirs:
        report=read(d/'report.json');assert report['all_passed']
        for p,h in report['input_sha256'].items():assert sha(ROOT/p)==h;hashes[p]=h
        for p,h in report['evidence_sha256'].items():assert sha(d/p)==h;hashes[str((d/p).relative_to(ROOT))]=h
        hashes[str((d/'report.json').relative_to(ROOT))]=sha(d/'report.json')
    known,source=source_audit(dirs);P=Program();field=field_audit(dirs,known,P);slp=slp_audit(dirs[3]/'normal42-adapter-SLP.json');inp=out/'independent-fixed-operation-checks.txt.gz';P.write(inp)
    assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
    native=work/'verify-field-certificates';cpp=ROOT/'code/s61_sparse_generic5_plane_rank3_review_20261104_v63.cpp'
    if verified_checkpoint is None:
        subprocess.run(['bash','-c','source /workspace/cloud-setup/activate.sh; c++ -std=c++17 -O2 "$1" -o "$2" $(pkg-config --cflags --libs flint) -lgmp','compile',str(cpp),str(native)],check=True,timeout=15)
        # No producer binary, search, pivot choice, determinant expansion, or point.
        with gzip.open(inp,'rb') as f:
            result=subprocess.run(['bash','-c','source /workspace/cloud-setup/activate.sh; exec "$1"','verify',str(native)],input=f.read(),capture_output=True,check=True,timeout=45)
        (out/'native.log').write_bytes(result.stdout+result.stderr);arith=json.loads(result.stdout)
    else:
        # Reuse only this review's already completed arithmetic on byte-identical
        # commands and checker source. Never reuse the producer's result.
        checkpoint=read(verified_checkpoint)
        with gzip.open(inp,'rb') as f:assert hashlib.sha256(f.read()).hexdigest()==checkpoint['uncompressed_input_sha256']
        assert sha(cpp)==checkpoint['independent_checker_cpp_sha256']
        assert hashlib.sha256(checkpoint['native_log'].encode()).hexdigest()==checkpoint['native_log_sha256']
        (out/'native.log').write_text(checkpoint['native_log']);arith=json.loads(checkpoint['native_log']);assert arith['all_passed']
        (out/'completed-arithmetic-checkpoint.json').write_text(json.dumps(checkpoint,indent=2)+'\n')
    proof=dict(source=source,field=field,SLP=slp,arithmetic=arith,independent_fixed_operation_commands=len(P.lines),no_elimination_search_or_points=True,scope='QQ(H,L,s0,q0,t0), with Sp4 and physical Gm dominance; no exhaustion or all-parameter chart claim')
    (out/'independent-proof.json').write_text(json.dumps(proof,indent=2)+'\n')
    paths=[Path(__file__),cpp,ROOT/f'research/notes/{ID}.md'];hashes.update({str(p.relative_to(ROOT)):sha(p) for p in paths})
    report=dict(id=ID,all_passed=True,source_and_all_QQ_fixed_certificates_independently_verified=True,rank_L18_plane3_normal42=True,no_author_solver_or_new_point_prime_GB=True,input_sha256=hashes,evidence_sha256={p.name:sha(p) for p in out.iterdir()},native_arithmetic_seconds=arith['seconds'],publication_reuses_exact_completed_check=verified_checkpoint is not None,seconds=time.monotonic()-start)
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--work',type=Path,required=True);p.add_argument('--verified-checkpoint',type=Path);a=p.parse_args();main(a.output,a.work,a.verified_checkpoint)
