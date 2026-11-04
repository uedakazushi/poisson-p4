#!/usr/bin/env python3
"""Explicit T34 -> T34*diag(Tden) coefficient adapter; frozen inputs untouched."""
import argparse
from fractions import Fraction
import gzip
import json
from pathlib import Path
import resource
import signal
import time
import sympy as s
from s61_nonsquare_lazy_quartic_oracle_20261103_v63 import Oracle
from s61_nonsquare_structural_reductions_review_20261103_v63 import read,sha

ROOT=Path(__file__).resolve().parents[1]
STEM='S61_NONSQUARE_ORACLE_FRAME_ADAPTER_20261103_v63'
ORACLE=ROOT/'reports/S61_NONSQUARE_LAZY_QUARTIC_ORACLE_20261103_v63/oracle.json'


class FrameAdapter:
    def __init__(self,o,data,frame):
        self.oracle,self.data,self.frame=o,data,frame
        self.D=[o.K.from_expr(s.sympify(c))for c in frame['Tden']]
        self.q_units=[o.K.from_expr(s.sympify(c))for c in frame['quadratic_encoded_source_units']]
        self.basis_units=[o.K.from_expr(s.sympify(c))for c in frame['basis57_encoded_source_units']]
    @classmethod
    def load(cls,frame_path=None):
        if frame_path is None:frame_path=ROOT/f'reports/{STEM}/adapter.json'
        frame=read(Path(frame_path));o,data=Oracle.load(ROOT/frame['oracle_path'])
        assert sha(ROOT/frame['oracle_path'])==frame['oracle_sha256']
        return cls(o,data,frame)
    def scalar(self,mon):
        ans=self.oracle.K.one
        for c,k in zip(self.D,mon):ans*=c**k
        return ans
    def equation_coefficient(self,row,mon):
        """[v^mon] F_row(D*v), in the encoded source coordinates."""
        node=self.oracle.coefficient(self.data['residual_weight_le4_roots'][row],tuple(mon))
        return self.oracle.lin([(self.scalar(mon),node)])
    def quad_coefficient(self,i,mon):
        node=self.equation_coefficient(self.frame['quadric_original_rows'][i],mon)
        return self.oracle.lin([(1/self.q_units[i],node)])
    def basis57_coefficient(self,i,mon):
        record=self.frame['basis57_source_definitions'][i];mon=list(mon)
        if record['kind']=='product':
            mon[record['light_index']]-=1
            if min(mon)<0:return 0
        node=self.equation_coefficient(record['original_equation'],tuple(mon))
        return self.oracle.lin([(1/self.basis_units[i],node)])
    def quartic_coefficient(self,row,mon):
        """[v^mon] Qhat_row = F_row(D*v)/D_H^2, for the actual122 recipe."""
        assert row in self.data['genuine_quartic_original_rows']
        node=self.equation_coefficient(row,mon)
        return self.oracle.lin([(1/self.D[5]**2,node)])
    def core_value(self,node):return self.oracle.core_value(node)


def build(out):
    start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60)
    out.mkdir(parents=True,exist_ok=False)
    o,data=Oracle.load(ORACLE);n0=len(o.nodes);oracle_report=ORACLE.parent/'report.json';assert read(oracle_report)['all_exact_source_and_graph_checks_passed']
    source_path=ROOT/'reports/S61_BINARY_NONSQUARE_CUBIC_REPLAY_20261103_v62/source-report.json';source=read(source_path)
    chart_path=ROOT/'reports/S61_BINARY_HYPERPLANE_WEIGHTED_CHART_20261102.json';chart=read(chart_path)
    quad_path=ROOT/'reports/S61_BINARY_HYPERPLANE_QUADRATIC_CODIMENSION_20261102.json';quad=read(quad_path)
    basis_path=ROOT/'reports/S61_BINARY_NONSQUARE_CUBIC_ROW_BASIS_20261103_v62/report.json';basis=read(basis_path)
    review=ROOT/'reports/S61_BINARY_NONSQUARE_CUBIC_ROW_BASIS_20261103_v62_REVIEW/report.json';assert read(review)['all_passed']
    for p in[oracle_report,review,basis_path]:
        for rel,h in read(p).get('input_sha256',{}).items():assert sha(ROOT/rel)==h,rel
    D=[o.K.from_expr(s.sympify(c))for c in source['Tden']]
    a,e,d=o.a,o.e,o.d;F=500*a+90*e**2+3*e-24;B=(50*a+18*e**2-9*e-2)*d+60*e-40;C=100*a+18*e**2-9*e-2
    assert D==list(map(o.K,[a,a,a,1,a,B/50,a,a]))
    def unit(c):
        c=o.K(c);assert c;counts=[];ground=[]
        for p in[c.numer,c.denom]:
            powers=[]
            for f in o.known:
                k=0
                while p.rem(f)==0:p=p.exquo(f);k+=1
                powers.append(k)
            assert p.is_ground and p;counts.append(powers);ground.append(p)
        return dict(constant=str(o.K(ground[0])/o.K(ground[1])),a_F_B_C_powers=[x-y for x,y in zip(*counts)])
    for c in D:unit(c)
    def scalar(mon):
        ans=o.K.one
        for c,k in zip(D,mon):ans*=c**k
        return ans
    # The oracle's chart really is T34; the old native source starts with Tn.
    for i,row in enumerate(data['T_atom_ids']):
        for k,node in enumerate(row):assert o.atoms[node]==o.K.from_expr(s.sympify(chart['T34'][i][chart['normal34'][k]]))
    for i,row in enumerate(data['U_atom_ids']):
        for k,node in enumerate(row):assert o.atoms[node]==o.K.from_expr(s.sympify(chart['U34'][i][k]))
    pd=lambda terms:o.K(o.R.from_dict({tuple(ex):s.QQ(c)for ex,c in terms}))
    light=[0,1,2,3,4,6,7]
    from itertools import combinations_with_replacement
    mon2=[]
    for pair in combinations_with_replacement(light,2):mon2.append(tuple(pair.count(i)for i in range(8)))
    hm=[0]*8;hm[5]=1;mon2.append(tuple(hm))
    q_units=[];q_proofs=[];quad_checks=0
    for i,row in enumerate(quad['quadratic_rows_sparse']):
        saved={tuple(mon):pd(terms)for mon,terms in row};seed=min(saved,key=lambda mon:len(saved[mon].numer))
        original=quad['quadratic_rows'][i];root=data['residual_weight_le4_roots'][original]
        factor=o.core_value(o.coefficient(root,seed))*scalar(seed)/saved[seed];q_proofs.append(unit(factor));q_units.append(factor)
        for mon in mon2:
            actual=o.core_value(o.coefficient(root,mon))*scalar(mon)
            assert actual==factor*saved.get(mon,o.K.zero);quad_checks+=1
        (out/'progress.json').write_text(json.dumps({'quadrics_checked':i+1,'seconds':time.monotonic()-start})+'\n')
    assert quad_checks==174
    # Normalization factors in the old63 ->57 source refer to Tn coordinates.
    norm_path=ROOT/'reports/S61_BINARY_NONSQUARE_CUBIC_ROW_BASIS_20261103_v62/normalization.txt';nl=norm_path.read_text().splitlines();assert nl[0]=='63'
    def native_poly(line):
        ts=line.split();n=int(ts[0]);assert len(ts)==1+4*n
        return o.K(o.R.from_dict({tuple(map(int,ts[2+4*k:5+4*k])):s.QQ(ts[1+4*k])for k in range(n)}))
    norms=[native_poly(line)for line in nl[1:]];assert len(norms)==63
    for value in norms:unit(value)
    matrix_path=ROOT/'reports/S61_BINARY_NONSQUARE_CUBIC_ROW_BASIS_20261103_v62/remaining.txt.gz'
    raw=gzip.decompress(matrix_path.read_bytes());assert __import__('hashlib').sha256(raw).hexdigest()==basis['remaining_matrix_sha256']
    lines=raw.decode().splitlines();assert lines[0]=='57 91 57'and len(lines)==1+57*91
    b_units=[];b_defs=[];anchors=[];product_checks=0
    for pos,original in enumerate(basis['remaining_original_row_indices']):
        label=source['labels'][original]
        if label[0]=='product':
            _,qi,k=label;factor=q_units[qi]*norms[original];definition=dict(kind='product',original_equation=quad['quadratic_rows'][qi],light_index=k)
            q_saved={tuple(mon):pd(terms)for mon,terms in quad['quadratic_rows_sparse'][qi]}
            for j,mon in enumerate(source['normal_monomials']):
                lower=list(mon);lower[k]-=1;expected=q_saved.get(tuple(lower),o.K.zero)if min(lower)>=0 else o.K.zero
                assert native_poly(lines[1+pos*91+j])*norms[original]==expected;product_checks+=1
        else:
            original_equation=label[1];root=data['residual_weight_le4_roots'][original_equation]
            choices=[]
            for j,mon in enumerate(source['normal_monomials']):
                n=int(lines[1+pos*91+j].split()[0])
                if n:choices.append((0 if mon[5]else 1,n,j))
            _,_,j=min(choices);mon=tuple(source['normal_monomials'][j]);saved=native_poly(lines[1+pos*91+j])
            factor=o.core_value(o.coefficient(root,mon))*scalar(mon)/saved
            definition=dict(kind='genuine_cubic',original_equation=original_equation)
            anchors.append(dict(basis_position=pos,original_equation=original_equation,monomial=list(mon),
                known_unit_factor=unit(factor),expanded_seed_core_terms=len(saved.numer)))
        unit(factor);b_units.append(factor);b_defs.append(definition)
        (out/'progress.json').write_text(json.dumps({'basis57_units_fixed':pos+1,'seconds':time.monotonic()-start})+'\n')
    assert len(anchors)==16 and product_checks==41*91
    # Actual quartics need the common known unit D_H^-2 to retain old H² data.
    hm2=[0]*8;hm2[5]=2;hs_path=ROOT/'reports/S61_BINARY_HYPERPLANE_HEAVY_SQUARE_PIVOT_20261102.json';hs=read(hs_path)
    for row,saved in zip([77,131],hs['coefficient_formulas']):
        actual=o.core_value(o.coefficient(data['residual_weight_le4_roots'][row],tuple(hm2)))
        assert actual*scalar(hm2)/D[5]**2==o.K.from_expr(s.sympify(saved))
    frame=dict(id=STEM,oracle_path=str(ORACLE.relative_to(ROOT)),oracle_sha256=sha(ORACLE),Tden=list(map(str,D)),
        encoded_frame='y=D*v; old quadratic/cubic sources use T34*D*v, oracle uses T34*y',
        quadric_original_rows=quad['quadratic_rows'],quadratic_encoded_source_units=list(map(str,q_units)),quadratic_unit_factorizations=q_proofs,
        basis57_original_row_indices=basis['remaining_original_row_indices'],basis57_source_definitions=b_defs,
        basis57_encoded_source_units=list(map(str,b_units)),sixteen_genuine_anchor_checks=anchors,
        all57_identity_justification='Old exact full63 reconstruction uses Tn=T34D, Un=U*ud. Native A2/ud and A3/ud² are the same unique implicit graph after y=Dv; extraction uses common ud³ and known-unit row normalizations. The16 anchors fix those already-proved global row scalars.',
        actual_quartic_original_rows=data['genuine_quartic_original_rows'],
        actual_quartic_encoded_recipe='Qhat_k(v)=F_k(D*v)/D_H² for all33 original quartic rows',
        common_quartic_unit_rescale=str(1/D[5]**2),old_H2_coefficients_and_u_v_remain_literal=True,
        physical_oracle_H2_comparison_was_already_T34_and_correct=True)
    (out/'adapter.json').write_text(json.dumps(frame,indent=2)+'\n')
    result=dict(id=STEM,all_exact_coordinate_and_row_unit_checks_passed=True,
        all174_quadratic_coefficients_compared=True,all3731_product_basis_coefficients_compared=True,
        all16_genuine_global_row_units_fixed_with_exact_anchors_and_source_provenance=True,
        full57_row_units_known_a_F_B_C=True,all33_actual_quartics_have_explicit_common_unit_row_normalization=True,
        old_H2_actual_formulas_correct_in_original_T34_frame=True,
        actual_C122_mixed_frame_rank_probe_not_performed=True,no_generic239_rank_repeated=True,
        frozen_inputs_not_mutated=True,initial_nodes=n0,query_nodes_added=len(o.nodes)-n0,
        seconds=time.monotonic()-start)
    paths=[Path(__file__),ROOT/f'research/notes/{STEM}.md',ORACLE,oracle_report,source_path,chart_path,quad_path,basis_path,review,norm_path,matrix_path,hs_path,
        ROOT/'code/parametric_implicit_jets_20261101.cpp',ROOT/'code/s61_binary_hyperplane_cubic_extract_20261102.cpp',ROOT/'code/s61_binary_hyperplane_cubic_matrix_20261102.py']
    result['input_sha256']={str(p.relative_to(ROOT)):sha(p)for p in paths};result['evidence_sha256']={'adapter.json':sha(out/'adapter.json')}
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');signal.alarm(0)
    print(json.dumps({k:v for k,v in result.items()if k not in ['input_sha256','evidence_sha256']},indent=2),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);build(p.parse_args().output)
