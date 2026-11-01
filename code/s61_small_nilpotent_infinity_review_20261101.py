#!/usr/bin/env python3
"""Independent leading-matrix audit of small nilpotent strata and vertex necessity."""
import argparse,hashlib,json
from pathlib import Path
import sympy as s
from s61_model import S61,X,J
ROOT=Path(__file__).resolve().parents[1]

def verify():
    model=S61();one=[i for i,w in enumerate(model.reduced['parameter_weights']) if w==1];two=[i for i,w in enumerate(model.reduced['parameter_weights']) if w==2]
    def matrix(Q):
        P=X[0]**2*J;P[1:,1:]+=Q;point=model.normalize(P)['parameters'];M=s.zeros(134,16)
        for row,terms in enumerate(model.reduced['equations']):
            constant=0
            for mon,n,d in terms:
                qi=[i for i in mon if i in two];bi=[i for i in mon if i in one]
                if len(qi)!=2:continue
                value=s.Rational(n,d)*s.prod(point[i] for i in qi)
                if bi:
                    assert len(bi)==1;M[row,one.index(bi[0])]+=value
                else:constant+=value
            assert s.cancel(constant)==0
        return M.applyfunc(s.expand)
    def parse(value):return s.Matrix(value).applyfunc(s.sympify)
    def det(M):
        domain=s.QQ.poly_ring(*sorted(M.free_symbols,key=str)) if M.free_symbols else s.QQ
        return domain.to_sympy(M.to_DM().convert_to(domain).det())
    small_path=ROOT/'reports/S61_OTHER_NILPOTENT_INFINITY_EXCLUSION_20261101.json'
    vertex_path=ROOT/'reports/S61_MINIMAL_INFINITY_VERTEX_NECESSITY_20261101.json'
    small=json.loads(small_path.read_text());vertex=json.loads(vertex_path.read_text());x,y,z,w=X[1:];V=s.Matrix([x,y,z,w])
    e=small['minimal_open'];Q=parse(e['Q_after_SL2_normalization']);Q0=parse(e['normalized_Q']);c,p,q,r,B=s.symbols('c p q r B')
    archive=json.loads((ROOT/'legacy/v05/data/minimal_nilpotent_certificate.json').read_text())
    oldlocals={f'x{i+1}':V[i] for i in range(4)};oldlocals.update({f't{i}':s.Symbol(f't{i}') for i in range(1,16)})
    replacement={oldlocals[f't{i}']:0 for i in range(1,16)}
    replacement.update({oldlocals['t1']:-3*c*c+3*c/4,oldlocals['t4']:p,oldlocals['t5']:q,oldlocals['t6']:4*c-s.Rational(1,2),oldlocals['t7']:r,oldlocals['t8']:5,oldlocals['t10']:c,oldlocals['t11']:1})
    from itertools import combinations
    for (i,j),formula in zip(combinations(range(4),2),archive['normal_form_brackets'].values()):
        assert s.expand(s.sympify(formula,locals=oldlocals).subs(replacement)-Q[i,j])==0
    T=parse(e['translation_matrix']);assert T*J[1:,1:]*T.T==J[1:,1:]
    pulled=(T.inv()*Q.subs(dict(zip(V,T*V)),simultaneous=True)*T.inv().T).applyfunc(s.expand)
    action={s.Symbol(k):s.sympify(v) for k,v in e['translation_parameter_action'].items()}
    assert (pulled-Q.subs(action,simultaneous=True)).applyfunc(s.expand)==s.zeros(4)
    M=matrix(Q0);assert M==parse(e['beta_matrix']);minor=det(M.extract(e['minor_rows'],e['minor_columns']))
    assert minor==-s.Rational(130977,327680)
    regular=[]
    for e in small['regular_nilpotent']:
        Q=parse(e['Q']);M=matrix(Q);assert M==parse(e['beta_matrix']);witness=e['rank_minor']
        # Recompute pivots independently rather than reuse the stored sample's rank.
        _,cols=M.to_DM().rref();cols=list(cols);_,rows=M[:,cols].T.to_DM().rref();rows=list(rows)
        assert len(cols)==12;minor=det(M[rows,cols]);assert minor!=0
        regular.append(dict(label=e['label'],rank=12,independent_minor=str(minor)))
    contradictions=[];g=s.symbols('g0:9');b=s.symbols('b0:16')
    assert vertex['beta_coordinate_order']==[f'u{i+1}' for i in one]
    for e in vertex['families']:
        Q=parse(e['Q']);M=matrix(Q)[70:,:];assert M==parse(e['beta_matrix'])
        K=parse(e['constant_row_combinations']);other=e['cancelled_beta_columns'];assert (K*M[:,other]).applyfunc(s.expand)==s.zeros(K.rows,len(other))
        relations=list(K*M[:,0]);gb=s.groebner(relations,*g,order='lex')
        offset=-s.Rational(6 if e['label']=='M0' else 1,5)
        expected=s.groebner([g[2],g[3]-offset,g[5],g[6],g[7]],*g,order='lex');assert list(gb)==list(expected)
        values={g[2]:0,g[3]:offset,g[5]:0,g[6]:0,g[7]:0,b[0]:1,b[2]:0};rows=(M*s.Matrix(b)).subs(values)
        assert s.expand(rows[16]-s.sympify(e['remaining_row86']))==0
        last=s.factor(rows[2].subs(b[3],0));assert last==s.sympify(e['row72_after_row86'])!=0
        contradictions.append(dict(label=e['label'],contradiction=str(last)))
    paths=[Path(__file__),small_path,vertex_path,ROOT/'code/s61_model.py',ROOT/'legacy/v05/data/minimal_nilpotent_certificate.json',ROOT/'legacy/v05/proofs/minimal_nilpotent_theorem.txt',ROOT/'legacy/v05/proofs/regular_nilpotent_theorem.txt',ROOT/'research/notes/S61_OTHER_NILPOTENT_INFINITY_EXCLUSION_20261101.md',ROOT/'research/notes/S61_MINIMAL_INFINITY_VERTEX_NECESSITY_20261101.md']
    return dict(id='S61_SMALL_NILPOTENT_INFINITY_REVIEW_20261101',all_passed=True,
        independent_original134_matrix_reconstructions=5,minimal_normal_form_verified_from_immutable_archive=True,
        uniform_minimal_open_rank=12,minimal_open_incidence_bound=13,regular_checks=regular,regular_incidence_bound=14,
        vertex_necessity_contradictions=contradictions,proof_note='research/notes/S61_SMALL_NILPOTENT_INFINITY_REVIEW_20261101.md',
        input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);a=p.parse_args();r=verify()
    with a.output.open('x') as f:json.dump(r,f,indent=2);f.write('\n')
    print(json.dumps(dict(all_passed=True,matrices=5)))
if __name__=='__main__':main()
