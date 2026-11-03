#!/usr/bin/env python3
"""Full actual CB29-source rank at the explicit CA self-pair's CB lift."""
import argparse
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import resource
import signal
import time
import sympy as s
from sympy.polys.matrices import DomainMatrix

ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_CA_SELF_PAIR_CB_SOURCE_RANK_20261103_v63'
X=s.symbols('h x y z w');h,x,y,z,w=X;a=s.symbols('a',nonzero=True)
PAIRS=list(itertools.combinations(range(5),2));MONS=list(itertools.combinations_with_replacement(range(5),2))


def read(path):
    with(gzip.open(path,'rt')if path.suffix=='.gz'else path.open())as f:return json.load(f)
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def vector(P):return s.Matrix([s.Poly(s.expand(P[i,j]),*X).coeff_monomial(X[m]*X[n])for i,j in PAIRS for m,n in MONS])
def canonical(P):
    d=s.Matrix([sum(s.diff(P[i,j],X[i])for i in range(5))for j in range(5)]);E=s.Matrix(X)
    return(P-(E*d.T-d*E.T)/5).applyfunc(s.expand)
def rawCB(f,g,H):
    P=s.zeros(5);P[1,2]=f;P[3,4]=g;P[2,3]=-s.diff(f,x)*z/2-s.diff(H,w);P[2,4]=-s.diff(f,x)*w/2+s.diff(H,z)
    return P-P.T
def action(i,j,P):
    M=s.zeros(5);M[i,j]=1
    return(M*P+P*M.T-X[j]*P.diff(X[i])).applyfunc(s.expand)
def main(out):
    start=time.monotonic();out.mkdir(parents=True,exist_ok=False)
    resource.setrlimit(resource.RLIMIT_AS,(2200000000,2200000000));signal.alarm(60)
    parent=ROOT/'reports/AMBIENT_CA_ORIGIN_SELF_PAIR_N_CB_INCIDENCE_20261103_v63';assert read(parent/'report.json')['all_passed']
    data=read(parent/'certificate.json.gz');globalpath=ROOT/'reports/AMBIENT_CB_GLOBAL_FLAG_MODEL_20261103_v63/report.json';globaldata=read(globalpath)
    reviewed=ROOT/'reports/AMBIENT_CB_GLOBAL_FLAG_MODEL_20261103_v63_REVIEW.json';assert read(reviewed)['all_passed']
    cubics=[s.prod((h,x,z,w)[i]for i in mon)for mon in itertools.combinations_with_replacement(range(4),3)if any(i>=2 for i in mon)];assert len(cubics)==16
    variables=s.symbols('u0:19');mu=s.symbols('mu');f=sum(c*b for c,b in zip(variables[:3],[h*h,h*x,x*x]));H=sum(c*b for c,b in zip(variables[3:],cubics))
    # Keep actual mu*f before differentiating; mu is not an independent free binary g.
    universal=vector(canonical(rawCB(f,mu*f,H)));sourcevars=[*variables,mu]
    actualH=s.sympify(data['actual_CB_cubic_H'],locals=dict(zip(['a',*map(str,X)],[a,*X])));actualF=3*a*x*x
    actual=[0,0,3*a,*[s.Poly(actualH,h,x,z,w).coeff_monomial(mon)for mon in cubics],0];sub=dict(zip(sourcevars,actual))
    J20=universal.jacobian(sourcevars).subs(sub).applyfunc(s.cancel);q=universal.subs(sub).applyfunc(s.cancel)
    expected=s.Matrix([s.sympify(p,locals={'a':a})for p in data['all150_CB_source']]);assert q==expected
    ratio_expected=vector(canonical(rawCB(0,actualF,0)));assert J20[:,19]==ratio_expected
    # Full finite nonlinear mu*f model: all three mixed dmu/du_i derivatives are retained and checked.
    mixed=[]
    for i,b in enumerate([h*h,h*x,x*x]):
        v=universal.diff(mu).diff(variables[i]);assert v==vector(canonical(rawCB(0,b,0)));mixed.append(v)
    P=s.zeros(5)
    for k,(i,j)in enumerate(PAIRS):P[i,j]=sum(q[15*k+m]*X[u]*X[v]for m,(u,v)in enumerate(MONS));P[j,i]=-P[i,j]
    parabolic=list(map(tuple,globaldata['parabolic_GL_positions']));flags=[pos for pos in itertools.product(range(5),repeat=2)if pos not in parabolic];assert len(flags)==9
    Jflag=s.Matrix.hstack(*[vector(action(i,j,P))for i,j in flags]);full=J20.row_join(Jflag)
    fixed=J20[:,:19];assert not any(c.free_symbols for c in fixed)
    rows=list(DomainMatrix.from_Matrix(fixed.T).convert_to(s.QQ).rref()[1]);assert len(rows)==19
    left=fixed[rows,:].inv();quot=(full[:,19:]-fixed*left*full[rows,19:]).applyfunc(s.cancel)
    domain=s.QQ.frac_field(a);dm=DomainMatrix.from_Matrix(quot).convert_to(domain);pc=list(dm.rref()[1]);pr=list(DomainMatrix.from_Matrix(quot[:,pc].T).convert_to(domain).rref()[1])
    sq=DomainMatrix.from_Matrix(quot[pr,pc]).convert_to(domain);det=s.factor(domain.to_sympy(sq.det()));assert det
    payload=dict(variable='a',actual_CB_source_parameter_variables=list(map(str,sourcevars)),actual_CB_source_parameter_values=list(map(str,actual)),
        universal_actual_CB150=list(map(str,universal)),all150_by20_actual_source_derivative=[[str(c)for c in row]for row in J20.tolist()],
        all_three_mu_times_f_mixed_derivatives150_by3=[[str(c)for c in row]for row in s.Matrix.hstack(*mixed).tolist()],
        actual_ratio_derivative150=list(map(str,ratio_expected)),full_flag9_positions=flags,
        full150_by29_actual_CB_source_derivative=[[str(c)for c in row]for row in full.tolist()],
        fixed19_selected_rows=rows,fixed19_inverse=[[str(c)for c in row]for row in left.tolist()],
        all150_by10_quotient_after_fixed19=[[str(c)for c in row]for row in quot.tolist()],
        quotient_rank=len(pc),quotient_minor_rows=pr,quotient_minor_columns=pc,quotient_minor_determinant=str(det))
    with(out/'certificate.json.gz').open('wb')as file:
        with gzip.GzipFile(fileobj=file,mode='wb',filename='',mtime=0)as archive:archive.write((json.dumps(payload,separators=(',',':'))+'\n').encode())
    result=dict(id=STEM,all_passed=True,full_CB_source_dimension=29,fixed_fiber19_rank=19,
        actual_binary_ratio_product_mu_times_f_retained=True,all_three_mixed_mu_f_derivatives_checked=True,
        actual_ratio_column_is_g_equals_3a_xsquared=True,all9_flag_directions_retained=True,
        exact_generic_CB_source_differential_rank=19+len(pc),quotient10_minor=str(det),
        CB_whole_fiber_type_or_global_normality_claimed=False,seconds=time.monotonic()-start)
    paths=[Path(__file__),parent/'report.json',parent/'certificate.json.gz',globalpath,reviewed];note=ROOT/f'research/notes/{STEM}.md'
    if note.exists():paths.append(note)
    result['input_sha256']={str(p.relative_to(ROOT)):sha(p)for p in paths};result['evidence_sha256']={'certificate.json.gz':sha(out/'certificate.json.gz')}
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');signal.alarm(0)
    print(json.dumps({k:v for k,v in result.items()if k not in ['input_sha256','evidence_sha256']},indent=2),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
