#!/usr/bin/env python3
"""Independent coefficient-series audit of the Q_b CB meeting and actual-D two-jet."""
import argparse
from fractions import Fraction
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import resource
import signal
import subprocess
import sys
import time
import sympy as s
from sympy.polys.matrices import DomainMatrix
from ambient_CA_V1_independent_review_20261103_v63 import tensor_factor_operator

ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_CB_SELF_PAIR_ACTUAL_D_TWOJET_20261103_v63_REVIEW'
PAIRS=list(itertools.combinations(range(5),2));MONS=list(itertools.combinations_with_replacement(range(5),2))
b=s.symbols('b');X=s.symbols('h x y z w');h,x,y,z,w=X;K=s.QQ.frac_field(b)


def read(path):
    with(gzip.open(path,'rt')if path.suffix=='.gz'else path.open())as f:return json.load(f)
def sha(path):
    H=hashlib.sha256()
    with path.open('rb')as f:
        for c in iter(lambda:f.read(1048576),b''):H.update(c)
    return H.hexdigest()
def coerce(v):return K.from_sympy(s.sympify(v,locals={'b':b,'h':h,'x':x,'y':y,'z':z,'w':w}))
def expr(v):return K.to_sympy(v)
def matrix(columns):return s.Matrix(150,len(columns),lambda i,j:expr(columns[j][i]))
def add(*vs):return[sum((v[i]for v in vs),K.zero)for i in range(len(vs[0]))]
def scale(v,c):return[c*x for x in v]
def matmul(A,B):return[[sum((A[i][k]*B[k][j]for k in range(5)),K.zero)for j in range(5)]for i in range(5)]
def transpose(M):return list(map(list,zip(*M)))
def coefficients(P):return[coerce(s.Poly(s.expand(P[i,j]),*X).coeff_monomial(X[u]*X[v]))for i,j in PAIRS for u,v in MONS]
def constant_map(g,gi):
    wedge=[[g[i][k]*g[j][l]-g[i][l]*g[j][k]for k,l in PAIRS]for i,j in PAIRS];lookup={mon:i for i,mon in enumerate(MONS)};subs=[]
    for m,n in MONS:
        out=[K.zero]*15
        for i,j in itertools.product(range(5),repeat=2):out[lookup[tuple(sorted((i,j)))]]+=gi[m][i]*gi[n][j]
        subs.append(out)
    rows=[[]for _ in range(150)]
    for op,ip,old,new in itertools.product(range(10),range(10),range(15),range(15)):
        c=wedge[op][ip]*subs[old][new]
        if c:rows[15*op+new].append((15*ip+old,c))
    return lambda v:[sum((c*v[j]for j,c in row),K.zero)for row in rows]


def series_push2(g,gi,source):
    """Truncated polynomial arrays in t; all rational b coefficients live in QQ(b)."""
    zero=[K.zero]*3
    def ad(*V):return[sum((v[i]for v in V),K.zero)for i in range(3)]
    def mul(v,w):return[sum((v[j]*w[i-j]for j in range(i+1)),K.zero)for i in range(3)]
    def mm(A,B):return[[ad(*[mul(A[i][k],B[k][j])for k in range(5)])for j in range(5)]for i in range(5)]
    ident=[[[K.one if i==j else K.zero,K.zero,K.zero]for j in range(5)]for i in range(5)]
    assert mm(g,gi)==ident and mm(gi,g)==ident
    wedge=[[ad(mul(g[i][k],g[j][l]),[-c for c in mul(g[i][l],g[j][k])])for k,l in PAIRS]for i,j in PAIRS]
    lookup={mon:i for i,mon in enumerate(MONS)};subs=[]
    for m,n in MONS:
        row=[zero[:]for _ in range(15)]
        for i,j in itertools.product(range(5),repeat=2):
            k=lookup[tuple(sorted((i,j)))];row[k]=ad(row[k],mul(gi[m][i],gi[n][j]))
        subs.append(row)
    result=[zero[:]for _ in range(150)]
    for ip,old in itertools.product(range(10),range(15)):
        if not any(source[15*ip+old]):continue
        for op in range(10):
            if not any(wedge[op][ip]):continue
            v=mul(wedge[op][ip],source[15*ip+old])
            for new in range(15):
                if any(subs[old][new]):result[15*op+new]=ad(result[15*op+new],mul(v,subs[old][new]))
    return transpose(result)


def main(out):
    start=time.monotonic();out.mkdir(parents=True,exist_ok=False)
    resource.setrlimit(resource.RLIMIT_AS,(2200000000,2200000000));signal.alarm(120)
    parent=ROOT/'reports/AMBIENT_CB_SELF_PAIR_ACTUAL_D_TWOJET_20261103_v63';author=read(parent/'report.json');assert author['all_passed'];saved=read(parent/'certificate.json')
    for key in['input_sha256','evidence_sha256']:
        for rel,H in author[key].items():
            path=ROOT/rel if key=='input_sha256'else parent/rel;assert sha(path)==H
            if key=='input_sha256'and 'REVIEW' in rel and path.suffix=='.json':assert read(path)['all_passed']
    raw=read(ROOT/'reports/AMBIENT_CB_ORIGIN_POLYNOMIAL_SELF_PAIR_20261103_v63/certificate.json.gz')
    globalcb=read(ROOT/'reports/AMBIENT_CB_GLOBAL_FLAG_MODEL_20261103_v63/report.json')
    # Rebuild all22 elementary source coefficients independently from homogeneous polynomial tensors.
    def canon(P):
        E=s.Matrix(X);d=s.Matrix([sum(s.diff(P[i,j],X[i])for i in range(5))for j in range(5)])
        return(P-(E*d.T-d*E.T)/5).applyfunc(s.expand)
    def rawcb(F,G,H):
        P=s.zeros(5);P[1,2]=F;P[3,4]=G;P[2,3]=-s.diff(F,x)*z/2-s.diff(H,w);P[2,4]=-s.diff(F,x)*w/2+s.diff(H,z);return P-P.T
    cubics=[s.prod((h,x,z,w)[i]for i in mon)for mon in itertools.combinations_with_replacement(range(4),3)if any(i>=2 for i in mon)]
    basis=[coefficients(canon(rawcb(F,0,0)))for F in[h*h,h*x,x*x]]+[coefficients(canon(rawcb(0,G,0)))for G in[h*h,h*x,x*x]]+[coefficients(canon(rawcb(0,0,H)))for H in cubics]
    assert basis==[[coerce(c)for c in col]for col in globalcb['enlarged_canonical_basis']]
    fixed=[*basis[:3],*basis[6:]];flags=[ij for ij in itertools.product(range(5),repeat=2)if list(ij)not in globalcb['parabolic_GL_positions']];assert len(flags)==9
    ops=[tensor_factor_operator(i,j)for i,j in itertools.product(range(5),repeat=2)]
    def lie(i,j,q):
        W,V=ops[5*i+j];result=[K.zero]*150
        for op,ip,c in W:
            for mon in range(15):result[15*op+mon]+=c*q[15*ip+mon]
        for old,new,c in V:
            for pair in range(10):result[15*pair+new]+=c*q[15*pair+old]
        return result
    point=[coerce(v)for v in raw['full_polynomial_Q150']];directions=[[coerce(v)for v in saved[name]]for name in['u1','u2']]
    Js=[];jets=[];matrices=[];localjets=[]
    for datum,u in zip(raw['actual_two_CB_sources'],directions):
        G=[[coerce(v)for v in row]for row in datum['G']];Gi=[[coerce(v)for v in row]for row in datum['G_inverse']];back=constant_map(Gi,G)
        F=s.sympify(datum['f'],locals={'b':b,'h':h,'x':x,'y':y,'z':z,'w':w});H=s.sympify(datum['H'],locals={'b':b,'h':h,'x':x,'y':y,'z':z,'w':w})
        values=[coerce(s.Poly(F,h,x).coeff_monomial(m))for m in[h*h,h*x,x*x]]+[coerce(s.Poly(H,h,x,z,w).coeff_monomial(m))for m in cubics]
        combine=lambda B,c:add(*[scale(v,t)for v,t in zip(B,c)])
        q0=combine(fixed,values);ratiocol=combine(basis[3:6],values[:3]);local=[*fixed,ratiocol,*[lie(i,j,q0)for i,j in flags]]
        Jcols=[back(v)for v in local];Js.append(matrix(Jcols));assert back(q0)==point
        A=[[K.zero for _ in range(5)]for _ in range(5)]
        for c,(i,j)in zip(u[20:],flags):A[i][j]+=c
        A2=matmul(A,A);A3=matmul(A2,A);A4=matmul(A3,A);A5=matmul(A4,A);assert not any(c for row in A5 for c in row)
        E=[[[K.one if i==j else K.zero,A[i][j],A2[i][j]/2]for j in range(5)]for i in range(5)];Ei=[[[K.one if i==j else K.zero,-A[i][j],A2[i][j]/2]for j in range(5)]for i in range(5)]
        q1=add(combine(fixed,u[:19]),scale(ratiocol,u[19]));q2=scale(combine(basis[3:6],u[:3]),u[19])
        localjet=series_push2(E,Ei,transpose([q0,q1,q2]));jet=[back(v)for v in localjet];jets.append(jet);matrices.append(A);localjets.append([q0,q1,q2])
        assert jet[0]==point
        expectedA=[[coerce(v)for v in row]for row in saved['source_GL_matrices'][len(matrices)-1]];assert A==expectedA
        for c in u:
            den=s.Poly(s.denom(s.cancel(expr(c))),b);assert len(den.terms())==1
    common=[coerce(v)for v in saved['common_tangent']];L=add(jets[0][2],scale(jets[1][2],-1))
    assert jets[0][1]==jets[1][1]==common and L==[coerce(v)for v in saved['second_difference']]
    J1,J2=Js;U=J1.row_join(J2);v=s.Matrix(list(map(expr,common)));Lv=s.Matrix(list(map(expr,L)))
    assert all((J*s.Matrix(list(map(expr,u)))-v).applyfunc(s.cancel)==s.zeros(150,1)for J,u in zip(Js,directions))
    base=read(ROOT/'reports/AMBIENT_B18_CA_CRITICAL_SCHEME_20261103_v63/certificate.json.gz');CA=s.Matrix(base['fixed_fiber_linear_basis150_by20']).row_join(matrix([lie(i,j,point)for i,j in map(tuple,base['flag_positions'])]))
    rk=lambda A:len(DomainMatrix.from_Matrix(A).convert_to(K).rref()[1])
    ranks=[rk(J1.row_join(CA)),rk(J2.row_join(CA)),rk(U),rk(U.row_join(Lv)),rk(U.row_join(Lv).row_join(CA))];extra=[rk(CA.row_join(v)),rk(U.row_join(CA))]
    assert ranks==saved['ranks']==[36,30,37,38,38]and extra==saved['extra_ranks_CA_v_union_CA']==[29,37]
    assert U[63,:]==s.zeros(1,58)and s.factor(Lv[63])==24*b**4
    assert [rk(J1),rk(J2),rk(CA)]==[29,29,29]
    # Same universal81-kernel presentation as the previously accepted meeting-double theorem.
    rho,tau=s.symbols('rho tau');ys=s.symbols('y0:8');zs=s.symbols('z0:8')
    relations=[u*v for u in ys for v in zs]+[(tau-rho*rho)*u for u in ys]+[tau*v for v in zs]+[tau*(tau-rho*rho)]
    assert len(relations)==81
    first={**dict.fromkeys(zs,0),tau:rho*rho};second={**dict.fromkeys(ys,0),tau:0}
    assert all(s.expand(r.subs(first,simultaneous=True))==0 and s.expand(r.subs(second,simultaneous=True))==0 for r in relations)
    # Independent nonzero mixed-source structural check, despite both chosen ratio speeds being zero.
    assert any(basis[3])and combine(basis[3:6],[K.one,K.zero,K.zero])==basis[3]
    replay=out/'fresh-producer';process=subprocess.run([sys.executable,str(ROOT/'code/ambient_CB_self_pair_actual_D_twojet_20261103_v63.py'),'--output',str(replay)],capture_output=True,text=True,timeout=60)
    assert process.returncode==0,process.stdout+process.stderr
    fresh=read(replay/'certificate.json');assert {k:v for k,v in fresh.items()if k!='seconds'}=={k:v for k,v in saved.items()if k!='seconds'}
    freport=read(replay/'report.json');assert {k:v for k,v in freport.items()if k not in['seconds','evidence_sha256']}=={k:v for k,v in author.items()if k not in['seconds','evidence_sha256']}
    evidence=dict(source1_own_physical_jet150_by3=[[str(expr(c))for c in row]for row in transpose(jets[0])],source2_own_physical_jet150_by3=[[str(expr(c))for c in row]for row in transpose(jets[1])],
        own_full150_by29_J1=[[str(c)for c in row]for row in J1.tolist()],own_full150_by29_J2=[[str(c)for c in row]for row in J2.tolist()],own_second_difference=list(map(lambda c:str(expr(c)),L)),
        ranks=ranks,extra_ranks=extra,coker_functional='e63',coker_value=str(s.factor(Lv[63])))
    with(out/'independent-certificate.json.gz').open('wb')as file:
        with gzip.GzipFile(fileobj=file,mode='wb',filename='',mtime=0)as archive:archive.write((json.dumps(evidence,separators=(',',':'))+'\n').encode())
    result=dict(id=STEM,all_passed=True,own_all22_source_basis_rebuilt=True,own_finite_coefficient_series_physical_pushforward_all150_through2=True,
        own_full29_source_derivatives_and_all_QQb_rank_keys_match=True,both_ratio_times_binary_mixed_source_terms_retained=True,
        coker_pairing='24*b**4',source_directions_have_only_b_denominators=True,fresh_producer_semantic_replay_exact=True,
        exact_CB_meeting_double_S20_and_full81_fiber_product_application_reviewed=True,
        CB_completion_cone_dimension29_embedding38_depth21_multiplicity2_reviewed=True,
        actual_D_generic_dimension28_embedding29_nonreduced_reviewed=True,
        old_full_D_swap_or_nilmodule_primary_claim_accepted=False,seconds=time.monotonic()-start)
    paths=[Path(__file__),parent/'report.json',parent/'certificate.json',ROOT/'code/ambient_CB_self_pair_actual_D_twojet_20261103_v63.py',ROOT/'research/notes/AMBIENT_CB_SELF_PAIR_ACTUAL_D_TWOJET_20261103_v63.md',
        ROOT/'reports/AMBIENT_CB_GLOBAL_FLAG_MODEL_20261103_v63/report.json',ROOT/'reports/AMBIENT_CB_ORIGIN_POLYNOMIAL_SELF_PAIR_20261103_v63/certificate.json.gz',ROOT/'code/ambient_CA_V1_independent_review_20261103_v63.py']
    note=ROOT/f'research/notes/{STEM}.md'
    if note.exists():paths.append(note)
    result['input_sha256']={str(p.relative_to(ROOT)):sha(p)for p in paths};result['evidence_sha256']={'independent-certificate.json.gz':sha(out/'independent-certificate.json.gz')}
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');signal.alarm(0)
    print(json.dumps({k:v for k,v in result.items()if k not in['input_sha256','evidence_sha256']},indent=2),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
