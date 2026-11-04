#!/usr/bin/env python3
"""Independent ambient16 gradients/Poisson coefficient basis; pure Python arithmetic."""
import argparse
from collections import Counter,defaultdict
from fractions import Fraction as F
import gzip,hashlib,itertools,json,resource,signal,time
from pathlib import Path
import sympy as s

ROOT=Path(__file__).resolve().parents[1]
STEM='Q_NEW_FAMILY_DEGREE26_POISSON_BASIS_20261104_v63_REVIEW'
P=32003;BASE=64
def read(p):
    with(gzip.open(p,'rt')if p.suffix=='.gz'else p.open())as f:return json.load(f)
def sha(p):
    h=hashlib.sha256()
    with p.open('rb')as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()
def mod(q):
    q=F(str(q));assert q.denominator%P;return q.numerator*pow(q.denominator,-1,P)%P
def add(a,b,c=1):
    z=dict(a)
    for m,v in b.items():
        z[m]=(z.get(m,0)+c*v)%P
        if not z[m]:del z[m]
    return z
def mul(a,b):
    z=defaultdict(int)
    for m,x in a.items():
        for n,y in b.items():z[m+n]+=x*y
    return {m:v%P for m,v in z.items()if v%P}
def total(values):
    z={}
    for p in values:z=add(z,p)
    return z
def scale(a,c):return {m:v*c%P for m,v in a.items()if v*c%P}
def mm(a,b):return [[total(mul(a[i][k],b[k][j])for k in range(4))for j in range(4)]for i in range(4)]
def madd(a,b):return [[add(x,y)for x,y in zip(r,t)]for r,t in zip(a,b)]
def pair(a,b):return total(mul(a[i][j],b[j][i])for i in range(4)for j in range(4))
def coef(p,d):return [p.get(BASE*i+j,0)for i in range(d+1)for j in range(d-i+1)]
def fromcoef(a,d):return {BASE*i+j:c for (i,j),c in zip(((i,j)for i in range(d+1)for j in range(d-i+1)),a)if c}
def det(a):
    a=[r.copy()for r in a];d=1
    for j in range(len(a)):
        i=next(i for i in range(j,len(a))if a[i][j])
        if i!=j:a[i],a[j]=a[j],a[i];d=-d%P
        d=d*a[j][j]%P;iv=pow(a[j][j],-1,P)
        for i in range(j+1,len(a)):
            c=a[i][j]*iv%P
            if c:a[i]=[(x-c*y)%P for x,y in zip(a[i],a[j])]
    return d
def main(out):
    started=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
    out.mkdir(parents=True,exist_ok=False)
    folder=ROOT/'reports/Q_NEW_FAMILY_DEGREE26_POISSON_BASIS_20261104_v63'
    ar,cert,witness,definitions=[read(folder/p)for p in ['report.json','native-certificate.json.gz','minor-witnesses.json','definitions.json']]
    assert ar['all_passed']and ar['complete_QQ_Inv26_basis']
    info_path=ROOT/'reports/Q_NEW_FAMILY_HYPERSURFACE_RECOGNITION_20261103_v63/invariant-certificates.json'
    spin_path=ROOT/'reports/Q_FIRSTJET_SPIN5_STRUCTURE_20261103_v63/certificate.json.gz'
    plane_path=ROOT/'reports/Q_NEW_FAMILY_DEGREE16_INVARIANT_PROBE_20261103_v63/report.json'
    oldgate=ROOT/'reports/Q_NEW_FAMILY_FULL_DEGREE24_HESSIAN_BASIS_20261104_v63_REVIEW/report.json'
    oldcolumns=oldgate.parent/'own-columns.json'
    spinegate=ROOT/'reports/Q_FIRSTJET_SPIN5_STRUCTURE_20261103_v63_REVIEW/report.json'
    info,spin,plane=map(read,[info_path,spin_path,plane_path]);og=read(oldgate)
    assert og['all_passed']and read(spinegate)['all_passed']
    assert sha(oldcolumns)==og['evidence_sha256']['own-columns.json']
    b=plane['ambient_beta16'][:3]
    beta=[{0:b[0][i]%P,BASE:b[1][i]%P,1:b[2][i]%P}for i in range(16)]
    mon=[tuple(x)for x in spin['quadratic_monomials']]
    # Derivatives are taken in all16 independent beta coordinates, then the
    # already fixed formal plane is substituted into their degree1 formulas.
    matrices=[];derivatives=[]
    for key in ['spin_adjoint_C_coefficients16x136','vector_adjoint_D_coefficients16x136']:
        values=[];ds=[[{}for _ in range(16)]for _ in range(16)]
        for pos,row in enumerate(spin[key]):
            value={}
            for (i,j),c in zip(mon,row):
                c=mod(c)
                if not c:continue
                value=add(value,mul(beta[i],beta[j]),c)
                ds[i][pos]=add(ds[i][pos],beta[j],c)
                ds[j][pos]=add(ds[j][pos],beta[i],c)
            values.append(value)
            assert total(mul(beta[k],ds[k][pos])for k in range(16))==scale(value,2)
        matrices.append([values[4*i:4*i+4]for i in range(4)])
        derivatives.append([[row[4*i:4*i+4]for i in range(4)]for row in ds])
    eye=[[{0:1}if i==j else{}for j in range(4)]for i in range(4)]
    cache={'':eye,'0':matrices[0],'1':matrices[1]}
    def word(w):
        if w not in cache:cache[w]=mm(word(w[:-1]),matrices[int(w[-1])])
        return cache[w]
    def valgrad(w):
        value=total(word(w)[i][i]for i in range(4))
        partial=[[[{}for _ in range(4)]for _ in range(4)]for _ in range(2)]
        for k,letter in enumerate(w):partial[int(letter)]=madd(partial[int(letter)],word(w[k+1:]+w[:k]))
        grad=[total(pair(partial[l],derivatives[l][k])for l in range(2))for k in range(16)]
        assert total(mul(beta[k],grad[k])for k in range(16))==scale(value,2*len(w))
        return value,grad
    w4=['0000','0001','0011','0101','0111'];w6=['000101','000111','011011']
    traces,grads=map(list,zip(*(valgrad(w)for w in w4+w6)))
    assert [coef(p,8)for p in traces[:5]]==cert['Inv8_trace4_columns']
    assert [coef(p,12)for p in traces[5:]]==cert['trace6_columns']
    ownGrad=[[coef(p,7 if i<5 else 11)for p in row]for i,row in enumerate(grads)]
    assert ownGrad==cert['ambient16_gradient_columns']
    O=s.Matrix(info['symplectic_form16']);Oi=O.inv();assert O.T==-O and O*Oi==s.eye(16)
    inverse=[(i,j,mod(Oi[i,j]))for i in range(16)for j in range(16)if Oi[i,j]]
    def bracket(f,g):return total(scale(mul(grads[f][i],grads[g][j]),c)for i,j,c in inverse)
    qC=valgrad('00')[0];qD=valgrad('11')[0];qCD=valgrad('01')[0]
    F4=add(add(qC,qCD,-4),qD,4)
    pairs=list(itertools.combinations(range(5),2));degree14=[bracket(i,j)for i,j in pairs]
    assert [coef(p,14)for p in degree14]==cert['degree14_columns']
    cols14=cert['degree14_selected_columns'];assert cols14==[0,1,4]
    basis14=[degree14[i]for i in cols14]
    degree18=[mul(F4,p)for p in basis14]+[bracket(i,5+j)for i in range(5)for j in range(3)]
    assert [coef(p,18)for p in degree18]==cert['degree18_columns']
    primitive=cert['primitive18_selected_candidate_columns'];assert primitive==list(range(3,9))
    old=read(oldcolumns);selected12=[0,1,2,3,4,5,9,10,21]
    I12=[fromcoef(old['degree12_columns'][i],12)for i in selected12]
    degree26=[mul(p,q)for p in basis14 for q in I12]+[mul(degree18[i],q)for i in primitive for q in traces[:5]]
    assert len(degree26)==57
    assert [coef(p,26)for p in degree26]==cert['degree26_columns']
    semantics=[dict(kind='Inv14_times_Inv12',Inv14_basis_index=i,Inv12_basis_index=j)for i in range(3)for j in range(9)]+[dict(kind='primitive18_times_Inv8',primitive18_basis_index=i,Inv8_index=j)for i in range(6)for j in range(5)]
    assert semantics==definitions['degree26']
    assert [semantics[i]for i in cert['degree26_selected_columns']]==definitions['Inv26_basis']
    minors={}
    for degree,columns in [(14,degree14),(18,degree18),(26,degree26)]:
        ww=witness[str(degree)];values=[coef(p,degree)for p in columns]
        minors[str(degree)]=det([[values[j][i]for j in ww['columns']]for i in ww['coefficient_rows']])
        assert minors[str(degree)]==ww['minor_mod32003']
    assert minors=={'14':7584,'18':26750,'26':7998}
    # Independent complete symmetric-power weight DP and all eight C2 Weyl
    # denominator terms, including the new degree26 upper bound.
    weights=[tuple(w)for w in info['Cartan_weights']];dp=[Counter()for _ in range(27)];dp[0][0,0]=1
    for w in weights:
        for degree in range(1,27):
            for (x,y),c in dp[degree-1].items():dp[degree][x+w[0],y+w[1]]+=c
    rho=(2,1);weyl=[]
    for perm in [(0,1),(1,0)]:
        for signs in itertools.product([-1,1],repeat=2):
            wr=[signs[i]*rho[perm[i]]for i in range(2)];sgn=(1 if perm==(0,1)else-1)*signs[0]*signs[1]
            weyl.append(((rho[0]-wr[0],rho[1]-wr[1]),sgn))
    dims=[sum(sign*dp[d][w]for w,sign in weyl)for d in range(27)]
    assert [dims[d]for d in [8,12,14,18,24,25,26]]==[5,9,3,9,67,0,45]
    saved=read(folder/'character-records.json')
    assert dims==[r['invariant_dimension']for r in saved]
    for r in saved:
        assert {tuple(t['weight']):(t['sign'],t['multiplicity'])for t in r['Weyl_terms']}=={w:(sign,dp[r['degree']][w])for w,sign in weyl}
    evidence=dict(ambient16_gradient_columns=ownGrad,degree14_columns=[coef(p,14)for p in degree14],
        degree18_columns=[coef(p,18)for p in degree18],degree26_columns=[coef(p,26)for p in degree26],
        exact_minors=minors,independent_character_dimensions_through26=dims,
        semantic57_columns=semantics,all_derivatives_ambient16_before_plane=True,
        polynomial_engine='independent Python integer dictionaries mod32003; packed exponent base64, maximal degree26')
    with gzip.open(out/'independent-coefficients.json.gz','wt')as f:json.dump(evidence,f);f.write('\n')
    paths=[Path(__file__),ROOT/f'research/notes/{STEM}.md',info_path,spin_path,plane_path,oldgate,oldcolumns,spinegate,
        *[folder/p for p in ['report.json','native-certificate.json.gz','minor-witnesses.json','definitions.json','character-records.json']],
        ROOT/'code/q_new_family_degree26_poisson_basis_20261104_v63.py',ROOT/'code/q_new_family_degree26_poisson_basis_20261104_v63.cpp',
        ROOT/'research/notes/Q_NEW_FAMILY_DEGREE26_POISSON_BASIS_20261104_v63.md']
    report=dict(id=STEM,all_passed=True,own_all32_CD_quadratic_tensors_and16_ambient_directional_derivatives=True,
        all8x16_word_gradients_and_Euler_identities_independently_verified=True,
        own_full_Omega_inverse_Poisson_contractions_all10_and18_columns=True,
        all57_degree26_columns_and45_semantic_basis_definitions_match=True,
        exact_coefficient_minors=minors,independent_full_Weyl_DP_through26=True,
        complete_QQ_Inv14_3_Inv18_9_Inv26_45_bases_certified=True,
        no_NEW_point_source26_expansion_producer_native_replay_or_higher_probe=True,
        NEW25_degree_bound_unchanged=26,NEW26_equation_claimed=False,
        input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},
        evidence_sha256={'independent-coefficients.json.gz':sha(out/'independent-coefficients.json.gz')},seconds=time.monotonic()-started)
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');signal.alarm(0)
    print(json.dumps({k:report[k]for k in ['all_passed','exact_coefficient_minors','seconds']},indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
