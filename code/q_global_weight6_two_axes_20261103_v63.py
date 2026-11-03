#!/usr/bin/env python3
"""Exact pure three-weight6 closed slice over QQ[a,a^-1], with two axes."""
import argparse
from collections import defaultdict
from fractions import Fraction
import gzip
import hashlib
import json
from pathlib import Path
import time
import sympy as s
from sympy.polys.rings import ring
from sympy.polys.matrices import DomainMatrix
from sympy.polys.matrices.normalforms import smith_normal_decomp

ROOT=Path(__file__).resolve().parents[1]
MONS=[(2,0,0),(1,1,0),(1,0,1),(0,2,0),(0,1,1),(0,0,2)]
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def evaluate(nodes,point,R):
    values=[]
    for op,*data in nodes:
        if op=='q':value=R(s.Rational(data[0]))
        elif op=='v':value=point[data[0]]
        elif op=='+':value=sum((values[i]for i in data[0]),R.zero)
        elif op=='*':
            value=R.one
            for i in data[0]:value*=values[i]
        else:raise ValueError(op)
        values.append(value)
    return values


def encode(M,R):return [[str(R.to_sympy(v))for v in row]for row in M.to_list()]


def build(output):
    started=time.monotonic();folder=ROOT/'reports/Q_SMITH_GLOBAL_FINITE_MODEL_20261103_v63'
    parent_path=folder/'model37.json.gz';parent=json.loads(gzip.decompress(parent_path.read_bytes()))
    kept=[parent['input_names'].index(name)for name in ['z31','z32','z33']]
    assert all(parent['input_weights'][i]==6 for i in kept)
    R,a,x,y,z=ring('a,x,y,z',s.QQ);point=[a]+[R.zero]*36
    for i,value in zip(kept,[x,y,z]):point[i]=value
    values=evaluate(parent['nodes'],point,R)
    print('weight6 source pullback ready',time.monotonic()-started,flush=True)
    residuals=[values[i]for i in parent['all350_residual_nodes']];raw=[values[i]for i in parent['raw150_reconstruction_nodes']]
    assert all(sum(mon[1:])in[2,3]for f in residuals for mon in f)
    quad_rows=[i for i,f in enumerate(residuals)if f and all(sum(mon[1:])==2 for mon in f)]
    cube_rows=[i for i,f in enumerate(residuals)if f and all(sum(mon[1:])==3 for mon in f)]
    assert len(quad_rows)==29 and cube_rows==[33,69]
    aa=s.Symbol('a');base=s.QQ.poly_ring(aa);A=base.gens[0]
    def coefficient(poly,mon):
        return sum((base.convert(c)*A**power[0]for power,c in poly.items()if power[1:]==mon),base.zero)
    M=DomainMatrix.from_list([[coefficient(residuals[i],mon)for mon in MONS]for i in quad_rows],base).to_dense()
    D,U,V=smith_normal_decomp(M);assert U*M*V==D
    vdet=base.to_sympy(V.det());assert vdet and not vdet.has(aa)
    print('Smith identity ready',time.monotonic()-started,flush=True)
    diagonal=[s.factor(base.to_sympy(D.to_list()[i][i]))for i in range(6)]
    assert diagonal[:5]==[1,1,1,aa,aa]
    assert s.cancel(diagonal[5]/(aa*(9*aa-8)*(9*aa+8))).is_Rational
    alpha=base.from_sympy(s.Rational(90,41)*(9*aa-8));beta=base.from_sympy(-s.Rational(360,41)*(9*aa-8))
    p=base.from_sympy((9*aa-8)*(9*aa+8));one,zero=base.one,base.zero
    C=DomainMatrix.from_list([[one,zero,2*alpha,zero,zero,alpha**2],
          [zero,one,beta,zero,alpha,alpha*beta],[zero,zero,zero,one,2*beta,beta**2],
          [zero,zero,one,zero,zero,alpha],[zero,zero,zero,zero,one,beta],[zero,zero,zero,zero,zero,p]],base).to_dense()
    CV=(C*V).to_list();diag=D.to_list();L=[[zero]*len(quad_rows)for _ in range(6)]
    for i in range(6):
        for j in range(6):L[i][j]=base.exquo(A*CV[i][j],diag[j][j])
    reverse_numerator=DomainMatrix.from_list(L,base).to_dense()*U
    assert reverse_numerator*M==C.scalarmul(A)
    print('reverse quadratic coefficients ready',time.monotonic()-started,flush=True)
    S,b,W,T,e=ring('a,W1,W2,eta',s.QQ)
    al=S.from_expr(base.to_sympy(alpha));be=S.from_expr(base.to_sympy(beta));pp=(9*b-8)*(9*b+8)
    generators=[W*W,W*T,T*T,W*e,T*e,pp*e*e]
    def change(poly):
        return sum((S(c)*b**mon[0]*(W-al*e)**mon[1]*(T-be*e)**mon[2]*e**mon[3]for mon,c in poly.items()),S.zero)
    transformed=[change(f)for f in residuals];forward=[]
    for f in transformed:
        grouped=defaultdict(lambda:base.zero)
        for mon,c in f.items():grouped[mon[1:]]+=base.convert(c)*A**mon[0]
        multipliers=[S.zero]*6
        for mon,c in grouped.items():
            i,j,k=mon
            if i>=2:which,rest=0,(i-2,j,k)
            elif i and j:which,rest=1,(i-1,j-1,k)
            elif j>=2:which,rest=2,(i,j-2,k)
            elif i and k:which,rest=3,(i-1,j,k-1)
            elif j and k:which,rest=4,(i,j-1,k-1)
            else:
                assert not i and not j and k>=2
                which,rest=5,(0,0,k-2);c=base.exquo(c,p)
            multipliers[which]+=S.from_expr(base.to_sympy(c))*W**rest[0]*T**rest[1]*e**rest[2]
        assert f==sum((q*g for q,g in zip(multipliers,generators)),S.zero)
        forward.append(multipliers)
    reverse=[[S.zero]*350 for _ in range(6)]
    reverse_num=reverse_numerator.to_list()
    for i in range(6):
        for j,row in enumerate(quad_rows):reverse[i][row]=S.from_expr(base.to_sympy(reverse_num[i][j]))
        assert sum((q*f for q,f in zip(reverse[i],transformed)),S.zero)==b*generators[i]
    print('all350 two-sided ideal identities ready',time.monotonic()-started,flush=True)
    # Verify the negative axis in the frozen old-frame T40 parametrization.
    eta_ring,eta=ring('eta',s.QQ)
    def along(poly,a_value,vector):
        answer=eta_ring.zero
        for mon,c in poly.items():
            term=eta_ring(c)*eta_ring(a_value)**mon[0]
            for coefficient,power in zip(vector,mon[1:]):
                if power:term*= (eta_ring(coefficient)*eta)**power
            answer+=term
        return answer
    minus=[along(f,s.Rational(-8,9),[s.Rational(1440,41),-s.Rational(5760,41),1])for f in raw]
    plus=[along(f,s.Rational(8,9),[0,0,1])for f in raw]
    assert all(not along(f,s.Rational(-8,9),[s.Rational(1440,41),-s.Rational(5760,41),1])for f in residuals)
    assert all(not along(f,s.Rational(8,9),[0,0,1])for f in residuals)
    t40_path=ROOT/'reports/Q_T40_FINITE_QUOTIENT_20261103_v63/T40-nine-parameter-model.json.gz'
    t40=json.loads(gzip.decompress(t40_path.read_bytes()));old_point=[eta_ring.zero]*9
    old_point[t40['input_weights'].index(6)]=eta
    old_values=evaluate(t40['nodes'],old_point,eta_ring);old_raw=[old_values[i]for i in t40['raw150_nodes']]
    ratio=next(f.get((1,),s.QQ.zero)/g.get((1,),s.QQ.zero)for f,g in zip(minus,old_raw)if g.get((1,),s.QQ.zero))
    assert ratio
    assert all(f==sum((eta_ring(c*ratio**mon[0])*eta**mon[0]for mon,c in g.items()),eta_ring.zero)for f,g in zip(minus,old_raw))
    plus_nonzero=[i for i,f in enumerate(plus)if f]
    assert plus_nonzero==[64,83,90,105,116,121,137,140]
    assert plus[64]==-s.Rational(256,3)*eta and plus[83]==-64*eta and plus[116]==64*eta
    output.mkdir(parents=True,exist_ok=False)
    certificate=dict(parent_kept_input_indices=kept,closed_variables=['a','x','y','z'],
       quadratic_monomial_basis=MONS,quadratic_rows=quad_rows,quadratic_matrix=encode(M,base),
       Smith_left=encode(U,base),Smith_right=encode(V,base),Smith_diagonal=list(map(str,diagonal)),
       normalized_coordinates={'W1':'x+90*(9*a-8)*z/41','W2':'y-360*(9*a-8)*z/41','eta':'z'},
       canonical_generators=list(map(str,generators)),source_all350=list(map(str,residuals)),
       transformed_all350=list(map(str,transformed)),raw150=list(map(str,raw)),
       forward350_by6=[[str(v)for v in row]for row in forward],
       reverse6_by350_numerator_over_a=[[str(v)for v in row]for row in reverse],
       negative_axis_T40_parameter_ratio=str(ratio),plus_raw150=list(map(str,plus)),minus_raw150=list(map(str,minus)))
    archive=output/'certificate.json.gz'
    with archive.open('wb')as stream:
        with gzip.GzipFile(fileobj=stream,mode='wb',filename='',mtime=0)as zipped:zipped.write((json.dumps(certificate,separators=(',',':'))+'\n').encode())
    note=ROOT/'research/notes/Q_GLOBAL_WEIGHT6_TWO_AXES_20261103_v63.md'
    inputs=[Path(__file__),note,parent_path,t40_path,
            ROOT/'reports/Q_SMITH_GLOBAL_FINITE_MODEL_20261103_v63_REVIEW.json',
            ROOT/'reports/Q_T40_GRAPH_OPEN_20261103_v63_REVIEW.json']
    report=dict(id='Q_GLOBAL_WEIGHT6_TWO_AXES_20261103_V63',all_passed=True,
       exact_base='QQ[a,a^-1]',closed_positive_inputs='global37 z31=x,z32=y,z33=z; all other positive inputs zero',
       exact_all350_ideal='(W1,W2)*(W1,W2,eta) + ((9a-8)*(9a+8)*eta^2)',
       all350_forward_and6_reverse_over_a_certificates_verified=True,
       quadratic_Smith_diagonal=list(map(str,diagonal)),cubic_rows=cube_rows,cubics_redundant_globally_after_a_unit=True,
       exact_primary_decomposition='(W1,W2,eta)^2 intersection (9a-8,W1,W2) intersection (9a+8,W1,W2)',
       no_embedded_prime_in_closed_slice=True,generic_zero_section_length=4,
       nilradical_square_zero_and_free_rank3_over_QQ_a_unit=True,
       nilradical_generators=['W1','W2','(9a-8)*(9a+8)*eta'],
       actual_integrable_axis_parameters=['a=8/9, (x,y,z)=(0,0,1)*eta','a=-8/9, (x,y,z)=(1440/41,-5760/41,1)*eta'],
       negative_axis_identified_with_existing_T40_by_all150_QQ_eta_comparison=True,
       negative_axis_parameter_ratio=str(ratio),positive_axis_component_not_identified=True,
       other_full_R37_nonflatness_or_components_inferred=False,
       full_weight6_four_variable_slice_at_a_two_over_nine_exhausted=False,
       input_sha256={str(p.relative_to(ROOT)):sha(p)for p in inputs},
       evidence_sha256={'certificate.json.gz':sha(archive)},seconds=time.monotonic()-started)
    (output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print('PASS exact closed weight6 primary decomposition and two axes;',report['seconds'],'seconds')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();build(args.output.resolve())
