#!/usr/bin/env python3
"""Independent sparse C36/C44 and fixed-rank witness verification of honest T_nil.
No producer imports, source solve, new jet, rank search, or quotient derivative.
"""
import hashlib,itertools,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
PRIMARY='AMBIENT_B18_HONEST_NIL_T_FIRST_VALUES_20261104_v63';ID=PRIMARY+'_REVIEW';OUT=ROOT/'reports'/ID
def sha(p):
    with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def put(d,k,c):
    c=d.get(k,F(0))+c
    if c:d[k]=c
    else:d.pop(k,None)
def plus(*ps):
    d={}
    for p in ps:
        for k,c in p.items():put(d,k,c)
    return d
def scaled(p,c):return {k:a*c for k,a in p.items()if a*c}
def sparse(v):return {i:F(c)for i,c in enumerate(v)if F(c)}
def pairs(v):return {int(i):F(c)for i,c in v if F(c)}
def save(n,d):(OUT/n).write_text(json.dumps(d,indent=2)+'\n')
def mvec(A,v):
    return {i:c for i,r in enumerate(A)if(c:=sum((a*v.get(k,F(0))for k,a in r.items()),F(0)))}
def mprod(A,B):
    return [plus(*(scaled(B[k],a)for k,a in row.items()))for row in A]
def eye(n):return [{i:F(1)}for i in range(n)]
def det_fixed(A):
    A=[dict(r)for r in A];D=F(1)
    for j in range(len(A)):
        i=next(i for i in range(j,len(A))if A[i].get(j,F(0)))
        if i!=j:A[i],A[j]=A[j],A[i];D=-D
        p=A[j][j];D*=p
        for i in range(j+1,len(A)):
            q=A[i].get(j,F(0))
            if q:
                for k,c in A[j].items():put(A[i],k,-q*c/p)
    return D

class FiniteC:
    """Generic monomial quotient B[tA,tB]/(tA²-nuA,tB²-nuB)."""
    def __init__(self,mons,nu):
        self.mons=[tuple(m)for m in mons];idx={m:i for i,m in enumerate(self.mons)}
        self.base={(i,j):idx.get(tuple(a+b for a,b in zip(m,n)))for i,m in enumerate(self.mons)for j,n in enumerate(self.mons)}
        self.nu=[sparse(v)for v in nu];self.dim=4*len(mons);self.table={}
        assert self.mons[0]==(0,)*8 and all(not p.get(0,F(0))for p in self.nu)
        for i,j in itertools.product(range(self.dim),repeat=2):
            a,u=divmod(i,4);b,v=divmod(j,4);k=self.base[a,b]
            if k is None:continue
            p={k:F(1)}
            for bit in range(2):
                if (u&v)&(1<<bit):p=self.bmul(p,self.nu[bit])
            if p:self.table[i,j]={4*k+(u^v):c for k,c in p.items()}
    def bmul(self,p,q):
        r={}
        for i,a in p.items():
            for j,b in q.items():
                k=self.base[i,j]
                if k is not None:put(r,k,a*b)
        return r
    def c(self,c):return {0:F(c)}if c else{}
    def mul(self,p,q):
        r={}
        for i,a in p.items():
            for j,b in q.items():
                for k,c in self.table.get((i,j),{}).items():put(r,k,a*b*c)
        return r
    def power(self,p,n):
        q=self.c(1)
        for _ in range(n):q=self.mul(q,p)
        return q
    def polynomial(self,poly,variables):
        q={}
        for e,c in poly:
            t=self.c(F(c))
            for p,n in zip(variables,e):
                if n:t=self.mul(t,self.power(p,n))
            q=plus(q,t)
        return q
    def component(self,p,i):return {k//4:c for k,c in p.items()if k%4==i}

def main():
    start=time.monotonic();signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));OUT.mkdir(exist_ok=True);pins={}
    def load(path):
        p=ROOT/path;pins[path]=sha(p);return json.loads(p.read_text())
    prefix='reports/'+PRIMARY+'/'
    report=load(prefix+'report.json');assert pins[prefix+'report.json']=='6d171143810a747bf83ad76af12fa1aec219db98be2243bb17bc2f37c7c3a76f'
    for p,h in report['input_sha256'].items():assert sha(ROOT/p)==h
    A1=load(prefix+'phase1-actual-m0-and-section.json');A2=load(prefix+'phase2-actual-radical-X-firsts.json');A3=load(prefix+'phase3-actual-W84-rank-and-mu.json')
    for n,h in report['evidence_sha256'].items():assert pins[prefix+n]==h
    def artifact(st,name):return load('reports/'+st+'_20261104_v63/'+name)
    H=artifact('AMBIENT_B18_HONEST_NIL_T_OPERATOR','actual-honest-Tnil128-972-formal-operator.json')
    honest_peer=load('reports/AMBIENT_B18_HONEST_NIL_T_OPERATOR_20261104_v63_REVIEW/report.json');assert honest_peer['all_passed']
    R=load(H['source_artifacts']['radical']['path']);K=load(H['source_artifacts']['new_resolution']['path'])
    N1=artifact('AMBIENT_B18_NATIVE128_FIRST','phase1-source-firsts.json')
    G=artifact('AMBIENT_B18_NATIVE128_FIRST','phase5-native-adapter-firsts.json')
    S1=artifact('AMBIENT_B18_NATIVE128_FIRST','phase4-native-source-and-R12.json')
    SB=artifact('AMBIENT_B18_T_NATIVE_S2_B11','phase3-actual-proper-G8-B11.json')
    NB=artifact('AMBIENT_B18_T_NATIVE_S2_B11','phase1-source-B11-values.json')
    BB=artifact('AMBIENT_B18_T_NATIVE_B11_SOURCE_PREFLIGHT','actual-B11-source-and-gauge-preflight.json')
    NIL=artifact('AMBIENT_B18_NORMALIZER_NIL_GENERATORS','actual-E4-nil-two-generators-and-relations.json')
    Bfirst=artifact('AMBIENT_B18_COMPARISON_X_FIRST640_LINEAR','actual-first640-X-mod-m2.json')
    for st in ['AMBIENT_B18_NATIVE128_FIRST','AMBIENT_B18_T_NATIVE_S2_B11','AMBIENT_B18_COMPARISON_X_FIRST640_LINEAR']:
        assert load('reports/'+st+'_20261104_v63_REVIEW/report.json')['all_passed']
    assert NIL['forward_GL4']==R['actual_C_variable_rules']['actual_core_GL4']
    mons9=[[0]*8]+[[int(i==j)for i in range(8)]for j in range(8)]
    C36=FiniteC(mons9,N1['nu_first9']);C44=FiniteC(BB['B11_basis_monomials8'],NB['nu_A_B_B11'])
    assert BB['B11_basis_monomials8'][2]==[0,0,0,0,0,0,0,2]
    def variables(C,kappa,u):
        core=[sparse(kappa[i])for i in [0,1,2,6]]
        return [sparse(v)for v in u[:5]]+[C.polynomial(NIL['forward_GL4'][n],core)for n in ['h','c','r','l']]
    va=variables(C36,S1['kappa9'],G['G8']);vb=variables(C44,SB['kappa9_C_B44'],SB['actualG8_C_B44'])
    def new_q(C,vs):
        d=C.polynomial(R['d_polynomial'],vs)
        q=[C.mul(d,C.polynomial(p,vs))for p in R['R9_generators_polynomials9']]
        return q,[C.mul(p,{b:F(1)})for p in q for b in range(4)]
    qa,qa36=new_q(C36,va);qb,qb44=new_q(C44,vb)
    # V=[Qrad,-U8]. Seven U_i=zi*tauA; last gA order>=3 is zero in both caches.
    for i in range(8):
        qa36.append({4*(i+1)+1:F(-1)}if i<7 else{})
        z=[0]*8
        if i<7:z[i]=1
        idx=BB['B11_basis_monomials8'].index(z)if i<7 and z in BB['B11_basis_monomials8']else None
        qb44.append({4*idx+1:F(-1)}if idx is not None else{})
    def gamma(C,q):
        x,y,z=[C.component(q,i)for i in [1,3,2]]
        return [x,y,plus(z,scaled(x,F(-11,12)))]
    gv=[gamma(C36,q)for q in qa36];gb=[gamma(C44,q)for q in qb44]
    assert all(not v.get(0,F(0))for col in gv for v in col)
    assert gv==[[pairs(v)for v in col]for col in A1['Vgamma_firsts_all44']]
    assert [c[2].get(2,F(0))for c in gb]==list(map(F,A1['third_gamma_r_squared_all44']))
    m0=[[F(0)]*44 for _ in range(24)]
    for j in range(44):
        for s in range(3):
            for i in range(7):m0[8*s+i][j]=gv[j][s].get(i+1,F(0))
            if s<2:m0[8*s+7][j]=gv[j][s].get(8,F(0))
            else:assert not gv[j][s].get(8,F(0));m0[23][j]=gb[j][s].get(2,F(0))
    assert m0==[list(map(F,r))for r in A1['m0_24_by44']]
    selected=A1['selected_source_columns24'];small=[sparse([r[j]for j in selected])for r in m0];smi=[sparse(r)for r in A1['selected_minor_inverse']]
    assert len(selected)==24 and mprod(small,smi)==mprod(smi,small)==eye(24)
    assert det_fixed(small)==F(A1['selected_minor_determinant'])
    # Validate the SAVED A by N*A=[D;0] and F=A*(dR), not by solving/inverting anew.
    savedA=[[sparse(v)for v in row]for row in A2['new_A12_by9_C36']]
    Nc=[[sparse(v)for v in row]for row in G['N']]
    D=[[C36.polynomial(p,va)for p in row]for row in R['G8_polynomial_factor_D8_by9']]
    detN=det_fixed([{j:v[0]for j,v in enumerate(row)if 0 in v}for row in Nc]);assert detN
    for i,j in itertools.product(range(12),range(9)):
        value=plus(*(C36.mul(Nc[i][k],savedA[k][j])for k in range(12)))
        assert value==(D[i][j]if i<8 else{})
    for i in range(12):assert plus(*(C36.mul(savedA[i][j],qa[j])for j in range(9)))==sparse(S1['physicalF12'][i])
    # Reconstruct all640 Res4(A^t)B/P and all140 restricted radical syzygy columns.
    old=Bfirst['source_columns640_rows56_sparse_dual9'];assert len(old)==640
    computed=[{}for _ in range(780)]
    def acc(col,row,p):
        q=plus(computed[col].get(row,{}),p)
        if q:computed[col][row]=q
        else:computed[col].pop(row,None)
    byrow={}
    for j,col in enumerate(old):
        for row,sp in col:
            if row<48:byrow.setdefault(row,[]).append((j,pairs(sp)))
            else:acc(j,row-12,pairs(sp))
    for row,entries in sorted(byrow.items()):
        oldgen,trace=divmod(row,4)
        for j in range(9):
            v=C36.mul(savedA[oldgen][j],{trace:F(1)})
            for component in range(4):
                p=C36.component(v,component)
                if p:
                    for col,c in entries:acc(col,4*j+component,C36.bmul(p,c))
    native=R['native_tensor_Koszul6_HB2_resolution'][1]
    assert native['shape']==[9,35] and len(native['columns'])==35
    for j,col in enumerate(native['columns']):
        for i,p in col:
            value=C36.polynomial(p,va)
            for trace in range(4):
                v=C36.mul(value,{trace:F(1)})
                for comp in range(4):
                    p=C36.component(v,comp)
                    if p:acc(640+4*j+trace,4*i+comp,p)
    actual=[{int(i):pairs(v)for i,v in col}for col in A2['all780_source_columns44_sparse_dual9']]
    assert computed==actual and sum(len(v)for c in actual for v in c.values())==3217
    # Fixed X0 certificate: inverse gives lower rank; all780 factorizations give upper.
    x0=[{i:v[0]for i,v in col.items()if v.get(0,F(0))}for col in actual]
    mrows=[sparse(r)for r in m0];mcolumns=[{i:m0[i][j]for i in range(24)if m0[i][j]}for j in range(44)]
    assert all(not mvec(mrows,v)for v in x0)
    xp=A2['X0_pivot_columns'];xr=A2['X0_pivot_rows'];xi=[sparse(r)for r in A2['X0_fixed_minor_inverse']]
    assert len(xp)==len(xr)==18
    xminor=[{j:x0[c].get(r,F(0))for j,c in enumerate(xp)if x0[c].get(r,F(0))}for r in xr]
    assert mprod(xminor,xi)==mprod(xi,xminor)==eye(18) and det_fixed(xminor)==F(A2['X0_fixed_minor_determinant'])
    xc=[]
    for j,col in enumerate(x0):
        co=mvec(xi,{i:col.get(r,F(0))for i,r in enumerate(xr)if col.get(r,F(0))});xc.append(co)
        assert [co.get(i,F(0))for i in range(18)]==list(map(F,A2['all780_X0_factorization_coefficients'][j]))
        assert plus(*(scaled(x0[xp[i]],c)for i,c in co.items()))==col
    # Zero-origin constant combinations, so the unknown m1 term cancels EXACTLY.
    residual=[i for i in range(780)if i not in xp];assert residual==A3['residual_original_X_columns'] and len(residual)==762
    wedges=K['basis_sector_wedge'][2];W=[];equations=0
    for j in residual:
        res={i:{t-1:c for t,c in v.items()if t>0}for i,v in actual[j].items()}
        for p,c in xc[j].items():
            for i,v in actual[xp[p]].items():
                res[i]=plus(res.get(i,{}),{t-1:-c*a for t,a in v.items()if t>0})
        first={}
        for i,v in res.items():
            for s,b in mcolumns[i].items():
                for t,a in v.items():put(first,(s,t),b*a)
        w={k:first.get((8*s+J[1],J[0]),F(0))for k,(s,J)in enumerate(wedges)};w={k:c for k,c in w.items()if c}
        replay={}
        for k,c in w.items():
            for r,mon,sgn in K['all8_differentials_sparse_columns'][1][k]:
                if sum(mon)==1:put(replay,(r,mon.index(1)),c*sgn)
        assert replay==first and not first.get((23,7),F(0));equations+=192;W.append(w)
    assert W==[pairs(c)for c in A3['all_W84_columns_sparse']]
    assert equations==A3['all_Koszul_linear_coefficients_verified']==146304
    # Fixed W81 minor and ALL762 column factorizations, no rank/pivot search.
    wp=A3['W_pivot_residual_columns'];wr=A3['W_pivot_rows'];wi=[sparse(r)for r in A3['W_fixed_minor_inverse']]
    assert len(wp)==len(wr)==81
    wminor=[{j:W[c].get(r,F(0))for j,c in enumerate(wp)if W[c].get(r,F(0))}for r in wr]
    assert mprod(wminor,wi)==mprod(wi,wminor)==eye(81) and det_fixed(wminor)==F(A3['W_minor_determinant'])
    factorhash=hashlib.sha256()
    for col in W:
        co=mvec(wi,{i:col[r]for i,r in enumerate(wr)if r in col})
        assert plus(*(scaled(W[wp[k]],a)for k,a in co.items()))==col
        text=json.dumps([[i,str(c)]for i,c in sorted(co.items())],separators=(',',':'))+'\n';factorhash.update(text.encode())
    assert factorhash.hexdigest()==A3['W_factorization_stream_sha256']
    assert 128-24-18-81==5==report['mu_R0_Tnil']==A3['mu_S_Tnil']
    save('independent-source-and-fixed-rank-proof.json',{
        'independent_generic_monomial_quotient_C36_C44_arithmetic':True,
        'tau_squares_nu_firsts_and_C44_r_squared_preserved':True,
        'all44_Vgamma_firsts_third_r_squared_and_m0_values_verified':True,
        'm0_rank24_fixed_minor_both_inverse_products_and_determinant':True,
        'actual_N12_origin_unit_determinant':str(detN),
        'all_N_A_E8D_C36_and_all12_F12_A_dR_identities':True,
        'all640_actual_Res4_Atranspose_B_and_unchanged_positive_P_sources':True,
        'all140_complete_radical_native_syzygies_not_old128':True,
        'all780_X_dual9_QQ_nonzero_slots':3217,
        'X0_rank18_fixed_minor_and_all780_factorizations':True,
        'constant_origin_cancellation_removes_m1_without_zero_assumption':True,
        'all762_residuals_full192_linear_Koszul_identities':equations,
        'third_r_squared_sector_component7_r_coefficient_explicitly_zero':True,
        'W0_rank81_fixed_minor_both_inverse_products_determinant_and_all762_factorizations':True,
        'QQ_linearity_of_original_column_Euler_lift_suffices_not_S_linearity':True,
        'mu_R0_Tnil':5,'unit_pair_count':123,
        'abstract_minimal_target_shape_only':[5,849],
        'new_scalar_source_jet_or_quotient_derivative_not_evaluated':True,
        'no_length_nil_generator_list_or_cube_claim':True})
    rep={'id':ID,'all_passed':True,'frozen':True,'m0_rank':24,'actual_X0_rank':18,'actual_W0_rank':81,'mu_R0_Tnil':5,
        'all780_actual_source_firsts_all_residual192_identities_and_fixed_rank_certificates_accepted':True,
        'complete128x972_source_unchanged_minimal5x849_not_built':True,
        'generator_count_only_not_length_nilgens_or_cube':True,
        'no_producer_import_rerun_source_solve_newjet_point_GB_or_rank_search':True,
        'input_sha256':pins,'code_sha256':sha(Path(__file__)),
        'note_sha256':sha(ROOT/'research/notes'/f'{ID}.md'),
        'evidence_sha256':{'independent-source-and-fixed-rank-proof.json':sha(OUT/'independent-source-and-fixed-rank-proof.json')},
        'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    assert sum(p.stat().st_size for p in OUT.iterdir())<32*1024*1024
    save('report.json',rep);signal.alarm(0);print(json.dumps({k:v for k,v in rep.items()if k!='input_sha256'},indent=2))

if __name__=='__main__':main()
