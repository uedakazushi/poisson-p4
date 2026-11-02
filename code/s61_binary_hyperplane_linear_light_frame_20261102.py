#!/usr/bin/env python3
"""A core-linear light frame with no additional localization."""
import argparse,hashlib,json,time,os,subprocess,tempfile
from pathlib import Path
from itertools import combinations,combinations_with_replacement
import sympy as s
from sympy.polys.rings import ring
from s61_model import S61
from s61_flag_transverse_obstruction_20261101 import specialize
root=Path(__file__).resolve().parents[1]
NOTE=root/'research/notes/S61_BINARY_HYPERPLANE_LINEAR_LIGHT_FRAME_20261102.md'
def verify():
    begin=time.time()
    data=json.loads((root/'reports/S61_BINARY_HYPERPLANE_WEIGHTED_CHART_20261102.json').read_text())
    old=json.loads((root/'reports/S61_BINARY_DEEPER_TRIANGULAR_CHART_20261101.json').read_text())['linear']
    R,a,e,d=ring('a,e,d',s.QQ);K=R.to_field();pol=lambda x:R.from_expr(s.sympify(x));frac=lambda x:K.from_expr(s.sympify(x))
    T=[[frac(row[j]) for j in data['normal34']] for row in data['T34']]
    U=[[frac(v) for v in row] for row in data['U34']]
    td=[]
    for j in range(8):
        den=R.one
        for row in T:den=den.lcm(row[j].denom)
        td.append(den)
    ud=R.one
    for row in U:
        for v in row:ud=ud.lcm(v.denom)
    Tn=[[v.numer*td[j].exquo(v.denom) for j,v in enumerate(row)] for row in T]
    Un=[[v.numer*ud.exquo(v.denom) for v in row] for row in U]
    m=S61();beta={int(k[1:])-1:s.sympify(v) for k,v in data['beta'].items()};two,F=specialize(m,beta,True)
    fixed=old['fixed'];cols=old['cols'];J=[[R.zero for j in range(45)] for i in range(135)];H=[]
    for i,f in enumerate(F):
        row=[]
        for mon,c in f.items():
            cc=pol(c)
            if len(mon)==2:
                x,y=mon
                if x!=fixed and y==fixed:J[i][cols.index(x)]+=cc
                if y!=fixed and x==fixed:J[i][cols.index(y)]+=cc
            elif len(mon)==1 and mon[0]==fixed:J[i][44]+=cc
            if fixed not in mon:
                ij=[cols.index(j) for j in mon]+[44]*(2-len(mon))
                if cc:row.append((ij,cc))
        H.append(row)


    candidate=root/'code/certificates/s61_binary_hyperplane_linear_light_frame_20261102.json'
    cert=json.loads(candidate.read_text());assert cert['variables']==['a','e','d']
    assert cert['light_positions']==[0,1,2,3,4,6,7]
    cw=data['local_coordinate_weights'];nw=data['normal_weights34'];zero=R.zero
    assert cert['raw_coordinate_indices']==[i for i,w in enumerate(cw) if w==1]
    new=[[zero for j in range(8)] for i in range(45)]
    for i in range(45):new[i][5]=Tn[i][5]
    for j,vec in zip(cert['light_positions'],cert['light_columns']):
        assert len(vec)==12
        for i,terms in zip(cert['raw_coordinate_indices'],vec):
            new[i][j]=R.from_dict({tuple(mon):s.QQ.convert(s.Rational(c)) for mon,c in terms})
    assert all(sum((v*new[j][k] for j,v in enumerate(row) if v),zero)==0
               for row in J for k in range(8))
    selected=data['selected34_rows']
    assert all(sum((v*Un[j][k] for j,v in enumerate(J[i]) if v),zero)==(ud if r==k else zero)
               for r,i in enumerate(selected) for k in range(34))
    W=s.Matrix(old['chart_basis']).applyfunc(s.sympify);Wi=W.inv()
    oldT=s.Matrix([[v.as_expr() for v in row] for row in T])
    newT=s.Matrix([[v.as_expr() for v in row] for row in new])
    normal_rows=[data['free34'][j] for j in data['normal34']]
    pi=Wi[normal_rows,:]
    assert (pi.to_DM()*oldT.to_DM()).to_Matrix()==s.eye(8)
    transition=(pi.to_DM()*newT.to_DM()).to_Matrix()
    tr=[[frac(v) for v in row] for row in transition.tolist()]
    assert all(sum((T[i][k]*tr[k][j] for k in range(8)),K.zero)==K(new[i][j])
               for i in range(45) for j in range(8))
    cut_rows=[data['free34'][j] for j in data['family_parameters']]
    assert (Wi[cut_rows,:].to_DM()*newT.to_DM()).is_zero_matrix
    assert all(not v or cw[i]==nw[j] for i,row in enumerate(new) for j,v in enumerate(row))
    assert all(not transition[i,j] or nw[i]==nw[j] for i in range(8) for j in range(8))
    det=s.factor(transition.to_DM().to_field().det().as_expr())
    B=50*a*d+18*d*e**2-9*d*e-2*d+60*e-40
    assert pol(det)==51200000*a*B
    light=cert['light_positions'];ldet=s.factor(transition[light,light].to_DM().to_field().det().as_expr())
    assert pol(ldet)==2560000000*a
    light_degree=max(sum(mon) for row in new for j in light for mon in row[j])
    heavy_degree=max(sum(mon) for row in new for mon in row[5])
    assert light_degree==1 and heavy_degree==3
    coeffs=[s.Rational(c) for row in new for j in light for c in row[j].values()]
    assert all(c.q==1 for c in coeffs) and max(map(abs,coeffs))<=240
    paths=[Path(__file__),NOTE,candidate,
           root/'reports/S61_BINARY_HYPERPLANE_WEIGHTED_CHART_20261102.json',
           root/'reports/S61_BINARY_DEEPER_TRIANGULAR_CHART_20261101.json',
           root/'code/s61_model.py',root/'code/s61_flag_transverse_obstruction_20261101.py']
    return dict(id='S61_BINARY_HYPERPLANE_LINEAR_LIGHT_FRAME_20261102',all_passed=True,
        core_variables=['a','e','d'],candidate_generation_replayed=False,
        original_equations=135,all_kernel_identities_checked=True,
        old_selected_right_inverse_rechecked=True,old_cut_conditions_checked=True,
        old_frame_times_transition_equals_new_frame=True,
        T_new=[[str(v.as_expr()) for v in row] for row in new],
        transition=[[str(v) for v in row] for row in transition.tolist()],
        transition_determinant=str(det),light_determinant=str(ldet),
        light_core_degree=light_degree,heavy_core_degree=heavy_degree,
        maximum_light_integer_coefficient=int(max(map(abs,coeffs))),
        normal_weights=nw,additional_localization_factors=[],
        full_implicit_models_identified_by_linear_substitution=True,
        proof_note=str(NOTE.relative_to(root)),
        limitations=['No new uniform rank coverage, length65, or ambient-removal assertion.'],
        input_sha256={str(path.relative_to(root)):hashlib.sha256(path.read_bytes()).hexdigest() for path in paths},
        seconds=time.time()-begin)

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();assert not args.output.exists()
    result=verify()
    with args.output.open('x') as stream:json.dump(result,stream,indent=2);stream.write('\n')
    print(json.dumps({k:result[k] for k in ['all_passed','transition_determinant','light_core_degree','seconds']}))
