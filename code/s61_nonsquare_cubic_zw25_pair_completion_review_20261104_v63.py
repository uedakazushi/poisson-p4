#!/usr/bin/env python3
"""Independent exact pair-dot and free-polynomial operator audit. No solver."""
import argparse,hashlib,itertools,json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
import sympy as s
from sympy.polys.rings import ring
ROOT=Path(__file__).resolve().parents[1]
BASE='S61_NONSQUARE_CUBIC_ZW25_PAIR_UNIT_COMPLETION_20261104_v63'
ID=BASE+'_REVIEW'
def read(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

# A separate polynomial dictionary engine, with arbitrary formal leaf symbols.
# No actual core rational functions are substituted into the other sixteen rows.
def plus(a,b):
    c=dict(a)
    for m,v in b.items():c[m]=c.get(m,Q(0))+v
    return {m:v for m,v in c.items() if v}
def times(a,b):
    c={}
    for m,v in a.items():
        for n,w in b.items():
            k=tuple(sorted(m+n));c[k]=c.get(k,Q(0))+v*w
    return {m:v for m,v in c.items() if v}
def neg(a):return {m:-v for m,v in a.items()}
def sub(a,b):return plus(a,neg(b))
ZERO={};ONE={():Q(1)}
def var(*v):return {(v,):Q(1)}
def eye(n):return [[ONE if i==j else ZERO for j in range(n)] for i in range(n)]
def total(values):
    out={}
    for v in values:out=plus(out,v)
    return out
def mm(A,B):return [[total(times(A[i][k],B[k][j]) for k in range(len(B))) for j in range(len(B[0]))] for i in range(len(A))]

def universal():
    b=[var('b',i) for i in range(5)];v=[var('v',i) for i in range(5)];k=var('k');t=total(times(x,y) for x,y in zip(b,v));tm=sub(t,ONE)
    U=[[sub(ONE if i==j else ZERO,times(v[i],b[j])) for j in range(5)]+[v[i]] for i in range(5)]+[b+[ZERO]]
    E=eye(6);Ei=eye(6)
    for j in range(5):E[5][j]=neg(times(k,b[j]));Ei[5][j]=times(k,b[j])
    W=[[times(v[i],b[j]) for j in range(5)]+[neg(v[i])] for i in range(5)]+[[neg(x) for x in b]+[ONE]]
    U2=mm(U,U);right=mm(mm(U,E),mm(Ei,U));left=mm(mm(Ei,U),mm(U,E));leftW=mm(mm(Ei,W),E)
    assert mm(E,Ei)==mm(Ei,E)==eye(6)
    for i in range(6):
        for j in range(6):
            ident=ONE if i==j else ZERO
            assert sub(U2[i][j],ident)==times(tm,W[i][j])
            assert sub(right[i][j],ident)==times(tm,W[i][j])
            assert sub(left[i][j],ident)==times(tm,leftW[i][j])
    output=mm([b+[k]],mm(U,E))[0]
    assert output==[times(times(sub(ONE,t),plus(ONE,k)),x) for x in b]+[t]
    determinant={}
    for p in itertools.permutations(range(6)):
        term=ONE
        for i in range(6):term=times(term,U[i][p[i]])
        determinant=plus(determinant,neg(term) if sum(p[i]>p[j] for i in range(6) for j in range(i+1,6))%2 else term)
    assert determinant==neg(t)
    return dict(independent_integer_polynomial_engine=True,both_U_inverse_products36=True,U0_square36=True,E_inverse72=True,row_output6=True,determinant720='-sum(b_i*v_i)')

def main(out):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False)
    folder=ROOT/'reports'/BASE;report=read(folder/'report.json');cert=read(folder/'certificate.json');assert report['frozen']
    hashes={}
    for p,h in report['input_sha256'].items():assert sha(ROOT/p)==h;hashes[p]=h
    assert sha(folder/'certificate.json')==report['certificate_sha256']
    # Independently bind the previously unreviewed pair-preflight coefficient
    # transport to the accepted P22 extractor, not merely its metadata.
    pairfolder=ROOT/'reports/S61_NONSQUARE_CUBIC_ZW25_PAIR_DUAL_PREFLIGHT_20261104_v63';pr=read(pairfolder/'report.json');pc=read(pairfolder/'certificate.json')
    for p,h in pr['input_sha256'].items():assert sha(ROOT/p)==h;hashes[p]=h
    assert sha(pairfolder/'certificate.json')==pr['certificate_sha256']
    cp=read(ROOT/'reports/S61_NONSQUARE_CUBIC_ZW25_UNITS_20261104_v63/checkpoint.json');inp=read(ROOT/'reports/S61_NONSQUARE_CUBIC_ZW25_AUGMENTED_INPUT_20261104_v63/input.json')
    assert read(ROOT/'reports/S61_NONSQUARE_CUBIC_ZW25_UNITS_20261104_v63_REVIEW/report.json')['all_passed']
    assert read(ROOT/'reports/S61_NONSQUARE_KERNEL22_EXTRACTOR_20261104_v63_REVIEW/report.json')['all_passed']
    R,a,e,d=ring('a,e,d',s.QQ);K=R.to_field()
    def poly(ts):return R.from_dict({tuple(m):s.QQ(c) for m,c in ts})
    def dec(q):return K(poly(q[0]))/K(poly(q[1]))
    old=read(ROOT/'reports/S61_NONSQUARE_KERNEL22_EXTRACTOR_20261104_v63/global-P22.json');roots=old['pair_dual6x2_scalar_roots'];cir=old['circuit'];need=set()
    def visit(i):
        if i in need:return
        need.add(i);n=cir['scalar_nodes'][i]
        if n[0] in ['a','s','m']:visit(n[1]);visit(n[2])
        else:assert n[0] in ['z','c']
    for row in roots:
        for i in row:visit(i)
    vals={}
    for i in sorted(need):
        n=cir['scalar_nodes'][i];o=n[0]
        if o=='z':v=K.zero
        elif o=='c':v=dec(cir['atoms_num_den'][n[1]])
        elif o=='a':v=vals[n[1]]+vals[n[2]]
        elif o=='s':v=vals[n[1]]-vals[n[2]]
        else:v=vals[n[1]]*vals[n[2]]
        vals[i]=v
    dual=[[vals[roots[5-i][0]],vals[roots[5-i][1]]/1600] for i in range(6)]
    assert dual==[[dec(x) for x in r] for r in pc['Lframe_dual6x2']]
    lines=(ROOT/'reports/S61_NONSQUARE_ELL_LFRAME_PREFLIGHT_20261104_v63/pair2x6.txt').read_text().splitlines();assert lines[0].split()==['2','6']
    pair=[]
    for i in range(2):
        row=[]
        for j in range(6):
            t=lines[1+i*6+j].split();assert len(t)==1+4*int(t[0]);row.append(K(R.from_dict({tuple(map(int,t[k+1:k+4])):s.QQ(t[k]) for k in range(1,len(t),4)})))
        pair.append(row)
    assert pair==[[dec(x) for x in r] for r in pc['Lframe_pair2x6']]
    assert [[sum((pair[i][k]*dual[k][j] for k in range(6)),K.zero) for j in range(2)] for i in range(2)]==[[K.one,K.zero],[K.zero,K.one]]
    assert [[K(poly(inp['matrix25x32'][i][j])) for j in range(20,26)] for i in [8,9]]==pair
    steps=[(step,co) for step in cp['all_completed_unit_steps'] for row,co in step['eliminations'] if row==4]
    assert len(steps)==1 and steps[0][0]['pivot_row']==0 and steps[0][0]['pivot_column']==0
    coefficient=dec(steps[0][1]);av=[dec(cp['matrix25x32'][4][j]) for j in range(4,10)];v=[x[1] for x in dual];assert v[1]==0
    for k,j in enumerate(range(4,10)):
        initial=K(poly(inp['matrix25x32'][4][j]));pivotrow=dec(cp['matrix25x32'][0][j]);assert initial==pair[1][k] and pivotrow==pair[0][k]/K(25*a/4)
        assert av[k]==initial-coefficient*pivotrow
    assert sum((x*y for x,y in zip(av,v)),K.zero)==1
    literal=cert['actual_six_dot_and_binding'];assert av==list(map(dec,literal['row4_coefficients_on_x0_Sym2_y'])) and v==list(map(dec,literal['transported_pair_second_column'])) and dec(literal['six_term_dot'])==1
    named=[a,d,500*a+90*e**2+3*e-24,(50*a+18*e**2-9*e-2)*d+60*e-40,100*a+18*e**2-9*e-2,5*a*d-6*e-4]
    for q,proof in zip(v+av,literal['all_12_new_atom_denominator_proofs']):
        denominator=R.from_expr(s.sympify(proof['constant']))
        for f,n in zip(named,proof['a_d_F_B_C_J_powers']):assert n>=0;denominator*=f**n
        assert denominator==q.denom
        assert all(n==0 for n in proof['a_d_F_B_C_J_powers'][1:])
    univ=universal()
    # Replay every new DAG node over independent formal checkpoint/dual leaves.
    nodes=cert['scalar_DAG']['nodes'];values=[]
    for i,n in enumerate(nodes):
        o,*t=n
        if o=='z':x=ZERO
        elif o=='one':x=ONE
        elif o=='checkpoint':
            r,c=t;assert cp['matrix25x32'][r][c][0];x=var('checkpoint',r,c)
        elif o=='accepted_pair_second':
            assert v[t[0]];x=var('dual',t[0])
        else:
            assert o in ['a','s','m'] and all(0<=j<i for j in t)
            x=plus(values[t[0]],values[t[1]]) if o=='a' else sub(values[t[0]],values[t[1]]) if o=='s' else times(values[t[0]],values[t[1]])
        values.append(x)
    def cpv(i,j):return var('checkpoint',i,j) if cp['matrix25x32'][i][j][0] else ZERO
    def get(A):return [[values[x] for x in row] for row in A]
    group=[4,6,7,8,9,5];duals=[0,2,3,4,5];b=[cpv(4,j) for j in group[:5]];vb=[var('dual',j) for j in duals];k=cpv(4,5)
    U0=[[sub(ONE if i==j else ZERO,times(vb[i],b[j])) for j in range(5)]+[vb[i]] for i in range(5)]+[b+[ZERO]]
    E=eye(6);Ei=eye(6)
    for j in range(5):E[5][j]=neg(times(k,b[j]));Ei[5][j]=times(k,b[j])
    U6=mm(U0,E);Ui6=mm(Ei,U0)
    assert get(cert['U0_6x6_roots'])==U0 and get(cert['E_6x6_roots'])==E and get(cert['E_inverse6x6_roots'])==Ei
    assert get(cert['U6_roots'])==U6 and get(cert['U6_inverse_roots'])==Ui6
    rows=cp['remaining_rows'];cols=cp['remaining_columns'];assert len(rows)==17 and len(cols)==24;gp=[cols.index(j) for j in group];pivot=cols.index(5)
    T=eye(24);Ti=eye(24)
    for i in range(6):
        for j in range(6):T[gp[i]][gp[j]]=U6[i][j];Ti[gp[i]][gp[j]]=Ui6[i][j]
    outside=[j for j in cols if j not in range(4,10) and cp['matrix25x32'][4][j][0]];assert outside==[1,2,3]
    tail=eye(24);taili=eye(24)
    for j in outside:tail[pivot][cols.index(j)]=neg(cpv(4,j));taili[pivot][cols.index(j)]=cpv(4,j)
    assert mm(tail,taili)==mm(taili,tail)==eye(24)
    T=mm(T,tail);Ti=mm(taili,Ti);assert get(cert['target24_roots'])==T and get(cert['target24_inverse_roots'])==Ti
    A=[[cpv(i,j) for j in cols] for i in rows];AT=mm(A,T);rp=rows.index(4);assert [values[x] for x in cert['pivot_row_transformed_before_simplification_roots']]==AT[rp]
    S=eye(17)
    for i in range(17):
        if i!=rp:S[i][rp]=neg(AT[i][pivot])
    assert get(cert['source17_roots'])==S
    newrows=[i for i in rows if i!=4];newcols=[j for j in cols if j!=5]
    expected=[[AT[rows.index(i)][cols.index(j)] for j in newcols] for i in newrows]
    assert get(cert['remaining16x23_scalar_roots'])==expected
    assert cert['old17_source_original25_rows']==rows and cert['old24_target_original32_columns']==cols
    assert cert['new16_source_original25_rows']==newrows and cert['new23_target_original32_columns']==newcols
    assert cert['new16_source_labels']==[inp['source25_labels'][i] for i in newrows]
    assert cert['pivot_row_exact_simplified']==[int(i==pivot) for i in range(24)]
    assert cert['source_pivot']==dict(original25_row=4,original32_column=5,value='1')
    assert cert['tail_clear_original_columns']==outside and cert['permuted_six_target_original_columns']==group
    for j in [27,28,29,30,31]:
        pos=cols.index(j)
        for i in range(24):assert T[i][pos]==Ti[i][pos]==(ONE if i==pos else ZERO)
    assert len(nodes)==2138 and not any(n[0] in ['i','inv'] for n in nodes)
    proof=dict(actual_P22_dual_ancestors=len(need),reverse6_and_q5_scale1600=True,actual_pair_four_identities=True,all_six_row4_source_bindings=True,actual_dot_one=True,all12_denominators_only_QQ_and_a=True,universal=univ,all2138_DAG_nodes_independently_replayed=True,all368_reduced_entries_bound=True,source_and_target_operators_bound=True,heavy5_unchanged=True,actual_same_cokernel_shapes=[[17,24],[16,23]],scope='One additional global unit; no original9minor, whole-FULLU57, quartic239, or new kernel count claim')
    (out/'independent-proof.json').write_text(json.dumps(proof,indent=2)+'\n')
    for p in [folder/'certificate.json',folder/'report.json',pairfolder/'certificate.json',pairfolder/'report.json',Path(__file__),ROOT/f'research/notes/{ID}.md']:hashes[str(p.relative_to(ROOT))]=sha(p)
    result=dict(id=ID,all_passed=True,actual_one_new_unit_and_same_cokernel_verified=True,no_author_solver_Gaussian_search_point_prime_or_new_inverse=True,input_sha256=hashes,evidence_sha256={'independent-proof.json':sha(out/'independent-proof.json')},seconds=time.monotonic()-start)
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(dict(all_passed=True,seconds=result['seconds'],report_sha256=sha(out/'report.json'))))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
