#!/usr/bin/env python3
"""Exact5-parameter first-jet slice and secant/Wronski generic projections.

Verifies the S61 coordinate map, dominance, weighted scaling, and universal
Wronskian identity. Generates full exact solver inputs, but does not run a GB
or assert a numerical generic-rank bound.
"""
import argparse
from collections import Counter
import gzip
import hashlib
import json
from pathlib import Path
import time

import sympy as sp
from sympy.polys.domains import QQ
from sympy.polys.rings import ring

from s61_model import S61
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from weight1_structure import COORDS,infinitesimal_action,sp_basis

ROOT=Path(__file__).resolve().parents[1]


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def run(output):
    output.mkdir(parents=True,exist_ok=False);started=time.monotonic()
    model=S61();z,u,beta=primitive_beta(model)
    old=sp.Matrix([[sp.expand(beta[i,j]).coeff(z[k]).coeff(v)for v in u]for(i,j),k in COORDS])
    rows=list(old.T.rref()[1]);assert len(rows)==16;left=old[rows,:].inv()
    slice_path=ROOT/'reports/PRIMITIVE_CURVE_DOMINANT_SLICE_20261101.json'
    saved=json.loads(slice_path.read_text());assert saved['all_passed']
    wb=sp.Matrix(saved['weight_basis_24_by16']);sm=sp.Matrix(saved['slice_matrix_b']);target=wb*sm
    section=left*target[rows,:];assert old*section==target
    point=target[:,0]
    tangent24=sp.Matrix.hstack(*[infinitesimal_action(point,A)for A in sp_basis()],point,target[:,1:])
    tangent=left*tangent24[rows,:];det=sp.factor(tangent.det());assert det==sp.Rational(438994543750,3)
    source=ROOT/'legacy/v02/data/symplectic_slice_reduced_equations.json';data=json.loads(source.read_text())
    weights=data['parameter_weights'];first=[i for i,w in enumerate(weights)if w==1];second=[i for i,w in enumerate(weights)if w==2]
    assert len(first)==16 and len(second)==45
    assert list(map(str,u))==['u'+str(i+1)for i in first]
    equation_weights=[set(sum(weights[i]for i in mon)for mon,_,_ in equation)for equation in data['equations']]
    assert len(equation_weights)==134 and all(len(v)==1 for v in equation_weights)
    eqweights=[next(iter(v))for v in equation_weights];assert set(eqweights)=={4,5}
    R,*ss=ring('ss1,ss2,ss3,ss4,ss5',QQ)
    symbols=sp.symbols('ss1:6')
    forms=[R.from_expr(section[i,0]+sum(section[i,j+1]*symbols[j]for j in range(5)))for i in range(16)]
    fi={i:j for j,i in enumerate(first)};qi={i:j for j,i in enumerate(second)};cache={():R.one}
    polynomials=[];encoded=[]
    for equation in data['equations']:
        row={}
        for mon,n,d in equation:
            betamon=tuple(i for i in mon if i in fi);qmon=tuple(qi[i]for i in mon if i in qi)
            if betamon not in cache:
                value=R.one
                for i in betamon:value*=forms[fi[i]]
                cache[betamon]=value
            row[qmon]=row.get(qmon,R.zero)+QQ(n,d)*cache[betamon]
        row={m:c for m,c in row.items()if c}
        assert all(len(m)<=2 for m in row)
        polynomials.append('+'.join('('+str(c).replace('**','^')+')'+''.join('*x'+str(i)for i in mon)for mon,c in row.items())or'0')
        encoded.append([[list(mon),[[list(ex),str(c)]for ex,c in poly.items()]]for mon,poly in row.items()])
    with gzip.open(output/'restricted-original134.json.gz','wt')as stream:json.dump(encoded,stream,separators=(',',':'))
    matrices=dict(primitive_old_map=[[str(c)for c in row]for row in old.tolist()],
                  normalized_first_jet_map=[[str(c)for c in row]for row in section.tolist()],
                  normalized_parameters=['ss1','ss2','ss3','ss4','ss5'],
                  old_weight_one_indices=first,old_weight_two_indices=second,
                  dominance_matrix=[[str(c)for c in row]for row in tangent.tolist()],dominance_determinant=str(det))
    (output/'slice-matrices.json').write_text(json.dumps(matrices,indent=2)+'\n')
    a=sp.symbols('A0:45');b=sp.symbols('B0:45');zz=sp.Symbol('zz')
    P=sum(a[i]*zz**i for i in range(45));Q=sum(b[i]*zz**i for i in range(45))
    wronski=sp.expand(P*sp.diff(Q,zz)-sp.diff(P,zz)*Q)
    explicit=sum((j-i)*(a[i]*b[j]-a[j]*b[i])*zz**(i+j-1)for i in range(45)for j in range(i+1,45))
    assert sp.expand(wronski-explicit)==0
    variables=','.join(f'x{i}'for i in range(45));script_info={}
    for mode in ['secant','wronski']:
        parameters=[*(f'ss{i}'for i in range(1,6)),'zz1']+(['zz2']if mode=='secant'else[])+['cc1','cc2']
        L1='+'.join('zz1^'+str(i)+'*x'+str(i)for i in range(45))+'-cc1'
        L2=L1.replace('zz1','zz2').replace('cc1','cc2')if mode=='secant'else '+'.join(str(i)+'*zz1^'+str(i-1)+'*x'+str(i)for i in range(1,45))+'-cc2'
        text='ring rr=(0,'+','.join(parameters)+'),('+variables+'),dp;\noption(redSB);option(redTail);\nideal I='+','.join(polynomials+[L1,L2])+';\n'
        text+=f'print("EXACT_{mode.upper()}_SMALL_PARAMETER_START");ideal G=slimgb(I);print("EXACT_{mode.upper()}_SMALL_PARAMETER_DONE");\n'
        text+='print("NGENS="+string(size(G)));print("DIM="+string(dim(G)));print("VDIM="+string(vdim(G)));\n'
        text+='if(size(G)==1 && G[1]==1){matrix C=lift(I,G);matrix M=matrix(I)*C;if(M[1,1]!=1){ERROR("UNIT_IDENTITY_FAILED");quit;}print("EXACT_UNIT_COFACTOR_PASS");}\nquit;\n'
        path=output/(mode+'-exact.sing');path.write_text(text);script_info[mode]=dict(parameters=parameters,parameter_count=len(parameters),effective_unknown_count=43)
    proof_dependencies=[ROOT/'research/notes/S61_AFFINE_TWO_COORDINATE_PROJECTION_20261102.md',
                        ROOT/'research/notes/S61_AFFINE_SIGMA_FREE_DEGREE_BOUND_20261103_v62.md']
    inputs=[Path(__file__),source,slice_path,ROOT/'code/primitive_jet_elliptic_curve_20261101.py',
            ROOT/'code/s61_model.py',ROOT/'code/s61_infinity_flags.py',ROOT/'code/weight1_structure.py',*proof_dependencies]
    report=dict(id='S61_AFFINE_SMALL_PARAMETER_PROJECTION_20261103_v62',all_passed=True,field='QQ',
                normalized_slice_dimension=5,primitive_map_rank=16,dominance_group='Sp4 times Gm',
                dominance_determinant=str(det),original_equation_weight_counts=dict(Counter(eqweights)),
                all134_original_equations_retained=True,restricted_q_count=45,wronskian_universal_identity=True,
                secant_linear_elimination_determinant='zz2-zz1',wronski_linear_elimination_determinant='1',
                solver_models=script_info,
                secant_rank_criterion='With the known true-boundary dim<=1, center-free generic secant projection has rank equal to TOTAL surface degree; rank<=29 would exclude original generic surfaces.',
                wronski_unit_criterion='In characteristic0, exact generic UNIT for the Wronski pair is equivalent to support dimension<=1; no projective degree comparison is asserted for that pair.',
                field_extension_scope='Preserve total2D cycle degree under the dominant field extension; do not apply degree>=30 separately to slice-field components.',
                generic_rank_bound_computed=False,positive_dimensional_components_excluded=False,five_point_exhaustion_claimed=False,
                input_sha256={str(p.relative_to(ROOT)):sha(p)for p in inputs},
                evidence_sha256={p.name:sha(p)for p in output.iterdir()if p.is_file()},seconds=time.monotonic()-started)
    (output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print('PASS',det,script_info,report['seconds'])


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',required=True,type=Path)
    run(parser.parse_args().output)
