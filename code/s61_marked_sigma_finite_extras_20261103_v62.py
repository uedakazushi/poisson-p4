#!/usr/bin/env python3
"""Full-incidence transverse certificates at the two finite-h extra points.

Candidate number-field coordinates were found by p-adic exploration; this
producer verifies the exact original135 identities and the actual ambient
incidence jets independently of that reconstruction procedure.
"""
import argparse
from collections import Counter
from functools import lru_cache
import hashlib
import json
from pathlib import Path
import time
import sympy as s
from sympy.polys.matrices import DomainMatrix
from s61_marked_sigma_point_lift_20261103_v62 import load_model,coerce
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from toric_closed_jet_coordinates_20261102 import coords,pairs,m2
from weight1_structure import J,COORDS,sp_basis,infinitesimal_action
from s61_model import X

ROOT=Path(__file__).resolve().parents[1]
PRIME=32003


def verify_case(lifter,data,entry):
    started=time.monotonic();K=s.QQ.algebraic_field(s.sqrt(s.Rational(entry['h_square'])),alias='h');h=K.unit
    conv=lambda value:coerce(K,value)
    expr=lambda value:s.factor(K.to_sympy(value))
    def residue(value,root):
        result=0
        for c in value.rep:
            assert int(c.denominator)%PRIME
            result=(result*root+int(c.numerator)*pow(int(c.denominator),-1,PRIME))%PRIME
        return result
    roots=[entry['h_mod_p'],(-entry['h_mod_p'])%PRIME]
    def integral(value):return all(int(c.denominator)%PRIME for c in value.rep)
    def unit(value):return integral(value) and all(residue(value,root) for root in roots)
    def dm(matrix):
        if isinstance(matrix,DomainMatrix):return matrix.convert_to(K)
        if isinstance(matrix,s.MatrixBase):return matrix.to_DM().convert_to(K)
        return DomainMatrix(matrix,(len(matrix),len(matrix[0])),K)
    def enc(matrix):
        if isinstance(matrix,DomainMatrix):matrix=matrix.to_Matrix()
        return [[str(c)for c in row]for row in matrix.tolist()]
    def minor(matrix,size=None):
        matrix=dm(matrix);columns=list(matrix.rref()[1]);rows=list(matrix.extract(list(range(matrix.shape[0])),columns).transpose().rref()[1])
        if size is not None:assert len(columns)==size
        determinant=matrix.extract(rows,columns).det();assert unit(determinant)
        return dict(rows=rows,columns=columns,determinant=str(expr(determinant)))

    normals=[conv(s.Rational(c))*(h if e else K.one)for c,e in zip(entry['normal_coefficients'],entry['normal_h_exponents'])]
    point=lifter.lift(h,normals,rho=1,domain=K);assert point['all135_zero'] and not point['t']
    q,beta=point['q'],point['beta'];assert all(integral(c)for c in q+beta)
    one,two,source=lifter.one,lifter.two,lifter.original;fixed=25;assert q[fixed]==-K.one
    qexpr=s.Matrix(list(map(K.to_sympy,q)));space=s.Matrix(X[1:]);x,y,z,w=space;physical=data['E']*qexpr;Q=s.zeros(4)
    for k,(i,j)in enumerate(pairs):Q[i,j]=sum(physical[k*10+a]*space[v]*space[t]for a,(v,t)in enumerate(m2))
    Q-=Q.T;Q=Q.applyfunc(s.expand)
    first=entry['h_square']=='-52/125';orbit_label='T4'if first else'T7';length=2 if first else 4
    r=conv(s.Rational(25,2)if first else s.Rational(25,4))*h;re=K.to_sympy(r)
    aa=s.Rational(5,8)if first else s.Rational(5,11);bb=-re/5 if first else 4*re/5
    darboux=s.Matrix([[0,aa,1,0],[0,bb,0,1],[1,0,-bb,aa],[0,1,0,0]])
    assert (darboux.T*J*darboux-J).applyfunc(s.expand)==s.zeros(4)
    shear=s.eye(4);shear[2,3]=-5*re/24 if first else 10*re/33
    swap=s.Matrix([[1,0,0,0],[0,1,0,0],[0,0,0,1],[0,0,-1,0]])
    frame=swap*shear*darboux.inv()
    assert (frame*J*frame.T-J).applyfunc(s.expand)==s.zeros(4)
    transformed=(frame*Q.subs(dict(zip(space,frame.inv()*space)),simultaneous=True)*frame.T).applyfunc(s.expand)
    A=-s.Rational(125,156)if first else -s.Rational(125,2112)*re;B=s.Rational(1,3)
    G=A*x*z*w+B*w**3 if first else A*x*x*z+B*w**3
    expected=s.zeros(4);expected[1,2]=s.diff(G,w);expected[1,3]=-s.diff(G,z);expected-=expected.T
    assert transformed==expected and unit(conv(A)) and unit(conv(B))
    assert all(integral(conv(c))for c in frame)

    M=[[K.zero]*16 for _ in range(135)];quadratic=[[K.zero]*16 for _ in range(16)]
    for row,equation in enumerate(source):
        for mon,n,d in equation:
            assert d%PRIME
            if sum(i in two for i in mon)!=2:continue
            bs=[one.index(i)for i in mon if i in one];value=conv(s.Rational(n,d))
            for i in mon:
                if i in two:value*=q[two.index(i)]
            if len(bs)==1:M[row][bs[0]]+=value
            if len(bs)==2:quadratic[bs[0]][bs[1]]+=value
    MD=dm(M);H=MD.nullspace().transpose();hrows=list(H.transpose().rref()[1]);H=H*H.extract(hrows,list(range(H.shape[1]))).inv()
    hdim=9 if first else 11;odim=8 if first else 7;assert H.shape==(16,hdim)
    assert (MD*H).is_zero_matrix and (MD*dm([[c]for c in beta])).is_zero_matrix
    QQ=dm(quadratic);assert (H.transpose()*(QQ+QQ.transpose())*H).is_zero_matrix
    assert all(integral(c)for row in H.to_list()for c in row)
    kernel_minor=minor(MD,16-hdim)
    motions=[]
    for Z in sp_basis():
        dQ=(sum(((Z*space)[i]*Q.diff(space[i])for i in range(4)),s.zeros(4))-Z*Q-Q*Z.T).applyfunc(s.expand)
        coefficients=s.Matrix([s.Poly(dQ[i,j],*space).coeff_monomial(space[v]*space[t])for i,j in pairs for v,t in m2])
        dq=data['Pi']*coefficients;motions.append(dq-qexpr*dq[fixed]/qexpr[fixed])
    O=s.Matrix.hstack(*motions);orbit_minor=minor(O,odim);qcuts=orbit_minor['rows']
    zz,bu,form=primitive_beta(lifter.model);vector=s.Matrix([form[i,j].diff(zz[k])for(i,j),k in COORDS]);embedding=vector.jacobian(bu)
    erows=list(embedding.T.to_DM().rref()[1]);left=embedding[erows,:].inv();sample=H.to_Matrix()*s.Matrix(range(2,hdim+2));sv=vector.subs(dict(zip(bu,sample)))
    beta_actions=s.Matrix.hstack(*[(left*infinitesimal_action(sv,Z)[erows,:]).applyfunc(s.expand)for Z in sp_basis()])
    dominance=minor(H.to_Matrix().row_join(beta_actions),16)

    frozen={two[fixed],*[two[i]for i in qcuts],*[one[i]for i in hrows]};variables=[i for i in range(62)if i not in frozen]
    N=len(variables);location=dict(zip(variables,range(N)));values=dict(zip(one,beta))|dict(zip(two,q))|{61:K.zero}
    linear=[[K.zero]*N for _ in range(135)];encoded=[]
    for row,equation in enumerate(source):
        terms=[]
        for mon,n,d in equation:
            ids=tuple(mon)+(61,)*(2-sum(i in two for i in mon));coefficient=conv(s.Rational(n,d));terms.append((ids,coefficient))
            for variable,power in Counter(ids).items():
                if variable not in location:continue
                rest=list(ids);rest.remove(variable);value=coefficient*power
                for i in rest:value*=values[i]
                linear[row][location[variable]]+=value
        encoded.append(terms)
    LD=dm(linear);pivots=list(LD.rref()[1]);assert len(pivots)==N-1
    selected=list(LD.extract(list(range(135)),pivots).transpose().rref()[1]);block=LD.extract(selected,pivots);jacobian_det=block.det();assert unit(jacobian_det)
    inverse=block.inv();free=next(i for i in range(N)if i not in pivots)
    direction=[K.zero]*N;direction[free]=K.one;tail=inverse*LD.extract(selected,[free])
    for index,column in enumerate(pivots):direction[column]=-tail.to_list()[index][0]
    assert (LD*dm([[c]for c in direction])).is_zero_matrix
    assert all(integral(c)for c in direction)
    lifts={1:direction};jet_records=[]
    for degree in range(2,length+1):
        @lru_cache(None)
        def product(ids,d):
            if not ids:return K.one if d==0 else K.zero
            i=ids[0];total=values[i]*product(ids[1:],d)
            if i in location:
                total+=sum((lifts[k][location[i]]*product(ids[1:],d-k)for k in range(1,min(d,degree-1)+1)),K.zero)
            return total
        raw=[[sum((coefficient*product(ids,degree)for ids,coefficient in row),K.zero)]for row in encoded]
        correction=inverse*dm([raw[i]for i in selected]);lift=[K.zero]*N
        for i,column in enumerate(pivots):lift[column]=-correction.to_list()[i][0]
        lifts[degree]=lift
        residual=dm(raw)-LD.extract(list(range(135)),pivots)*correction;residuals=[row[0]for row in residual.to_list()]
        assert all(residuals[i]==K.zero for i in selected)
        if degree<length:assert not any(residuals)
        else:
            witness=next(i for i,c in enumerate(residuals)if unit(c))
        jet_records.append(dict(degree=degree,residuals=[str(expr(c))for c in residuals],
                                implicit_lift=[str(expr(c))for c in lift]))
    print(orbit_label,'length',length,'unitrow',witness,'coefficient',expr(residuals[witness]),flush=True)
    return dict(h_square=entry['h_square'],residue_h_roots=roots,all135_original_equations_zero=True,
                normal_coordinates=[str(expr(c))for c in normals],
                normal_residues_by_h={str(root):[residue(c,root)for c in normals]for root in roots},
                beta=[str(expr(c))for c in beta],q=[str(expr(c))for c in q],t='0',Q=enc(Q),
                symplectic_frame=enc(frame),transformed_Q=enc(transformed),Hamiltonian=str(G),
                orbit=orbit_label,H_dimension=hdim,orbit_dimension=odim,incidence_dimension=hdim+odim,
                highest_kernel=enc(H),kernel_rank_minor=kernel_minor,orbit_minor=orbit_minor,
                beta_dominance_minor=dominance,beta_cut_indices=hrows,q_cut_indices=qcuts,
                variables=variables,Jacobian_rows=selected,Jacobian_columns=pivots,
                Jacobian_determinant=str(expr(jacobian_det)),normal_direction=[str(expr(c))for c in direction],
                jets=jet_records,unit_residual_row=witness,unit_residual=str(expr(residuals[witness])),
                transverse_length_upper=length,known_generic_same_incidence_length=length,
                ambient_free_rank=length,ambient_t_zero=True,seconds=time.monotonic()-started)


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--case',choices=['all','52','16'],default='all');args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=False)
    lifter=load_model(ROOT/'reports/S61_MARKED_SIGMA_POINT_LIFT_20261103_v62_FINAL/model');data=coords()
    path=ROOT/'code/certificates/s61_marked_sigma_extra_points_20261103_v62.json';entries=json.loads(path.read_text())
    if args.case!='all':entries=[entry for entry in entries if entry['h_square']==f'-{args.case}/125']
    results=[verify_case(lifter,data,entry)for entry in entries]
    paths=[Path(__file__),path,ROOT/'code/s61_marked_sigma_point_lift_20261103_v62.py',
           ROOT/'code/s61_marked_sigma_pencil_model_20261102.py',ROOT/'code/s61_model.py',
           ROOT/'code/primitive_jet_elliptic_curve_20261101.py',ROOT/'code/toric_closed_jet_coordinates_20261102.py',
           ROOT/'code/weight1_structure.py',ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json',
           ROOT/'legacy/v02/data/symplectic_slice_reduced_equations.json',
           ROOT/'research/notes/S61_R0_MIXED_CURVE_EXHAUSTION_20261101.md',
           ROOT/'research/notes/S61_R0_MIXED_SURFACE_T7_EXHAUSTION_20261101.md',
           ROOT/'research/notes/S61_BINARY_A_ZERO_LIFT_20261103.md',
           ROOT/'research/notes/S61_MARKED_SIGMA_FINITE_EXTRAS_20261103_v62.md']
    report=dict(id='S61_MARKED_SIGMA_FINITE_EXTRAS_20261103_v62',all_passed=True,cases=results,
                input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in paths},
                prepared_model_input_sha256=hashlib.sha256(lifter.input_path.read_bytes()).hexdigest(),
                scope='Two finite-h nonzero rho points; full marked pencil exhaustion is separate.')
    (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__':main()
