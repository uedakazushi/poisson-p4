#!/usr/bin/env python3
"""Certify the eta!=0 chart of the true through-degree3 initial ideal.

Exact QQ constant-row cofactors and triangular polynomial division only.
No Groebner basis, modular inference, or claim about the full completed
Poisson ideal is used.
"""
import argparse
from fractions import Fraction as F
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import time

import sympy as sp
from sympy.polys.matrices import DomainMatrix

ROOT=Path(__file__).resolve().parents[1]


def add(p,q,c=F(1)):
    r=p.copy()
    for m,a in q.items():
        value=r.get(m,F(0))+c*a
        if value:r[m]=value
        else:r.pop(m,None)
    return r


def scale(p,c):return {m:a*c for m,a in p.items()if a*c}


def mul(p,q):
    r={}
    for m,a in p.items():
        for n,b in q.items():
            k=tuple(sorted(m+n));v=r.get(k,F(0))+a*b
            if v:r[k]=v
            else:r.pop(k,None)
    return r


def power(p,k):
    r={():F(1)}
    for _ in range(k):r=mul(r,p)
    return r


def substitute(p,i,rhs):
    out={};quotient={}
    for m,c in p.items():
        k=m.count(i)
        if not k:out=add(out,{m:c});continue
        rem=tuple(j for j in m if j!=i)
        out=add(out,mul({rem:c},power(rhs,k)))
        for e in range(k):
            quotient=add(quotient,mul({tuple(sorted(rem+(i,)*(k-1-e))):c},power(rhs,e)))
    relation=add({(i,):F(1)},rhs,F(-1))
    assert add(p,out,F(-1))==mul(relation,quotient)
    return out,quotient


def compose(p,rules):
    out={}
    for m,c in p.items():
        term={():c}
        for i in m:term=mul(term,rules[i])
        out=add(out,term)
    return out


def pack(p):return [[list(m),str(c)]for m,c in sorted(p.items())]


def unpack(data):return {tuple(m):F(c)for m,c in data}


def rank(M):return len(M.rref()[1])


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=False);started=time.monotonic()
    qpath=ROOT/'reports/Q_ORDINARY_QUADRATIC_CONE_20261103_v63/quadrics.json.gz'
    gpath=ROOT/'reports/Q_ORDINARY_CUBIC_MINIMAL_GENERATORS_20261103_v63/generators.json.gz'
    q=json.loads(gzip.decompress(qpath.read_bytes()));g=json.loads(gzip.decompress(gpath.read_bytes()))
    names=q['variables'];assert names==['u']+[f'z{i}'for i in range(34)]and g['variables']==names
    ix={name:i for i,name in enumerate(names)};eta=ix['z33'];assert eta==34
    allq=[{(i,j):F(c)for i,j,c in row}for row in q['all350_quadratic_rows']]
    qindices=g['quadratic_generators'];raw=[allq[k]for k in qindices]+[unpack(row['terms'])for row in g['cubic_generators']]
    assert len(qindices)==57 and len(raw)==74
    mons=list(itertools.combinations_with_replacement(range(35),2));mi={m:i for i,m in enumerate(mons)};Q=sp.QQ
    def qq(c):return Q(c.numerator,c.denominator)
    qrows=[]
    for p in raw[:57]:
        row=[Q.zero]*len(mons)
        for m,c in p.items():row[mi[m]]=qq(c)
        qrows.append(row)
    M=DomainMatrix(qrows,(57,len(mons)),Q);_,piv=M.rref();assert len(piv)==57
    inverse=M.extract(list(range(57)),list(piv)).inv()
    linears=[add({(ix[name],):F(1)},{(eta,):F(-c)})for name,c in [('z11',4),('z12',6),('z23',-4)]]
    pure_targets=[]
    pure=[ix[name]for name in ['z11','z12','z23','z33']];values=[4,6,-4,1]
    for i in range(4):
        for j in range(i,4):
            if (i,j)==(3,3):continue
            p=add({(pure[i],pure[j]):F(1)},{(eta,eta):F(-values[i]*values[j])})
            pure_targets.append((f'pure_{names[pure[i]]}_{names[pure[j]]}',p))
    square_targets=[(f'square_{i}',mul(p,p))for i,p in enumerate(linears)]
    cofactor=[]
    for label,p in pure_targets+square_targets:
        full=[Q.zero]*len(mons)
        for m,c in p.items():full[mi[m]]=qq(c)
        coeff=DomainMatrix([[full[j]for j in piv]],(1,57),Q)*inverse
        assert coeff*M==DomainMatrix([full],(1,len(mons)),Q)
        cofactor.append({'label':label,'polynomial':pack(p),'constant_cofactors':[
            [j,qindices[j],str(Q.to_sympy(c))]for j,c in enumerate(coeff.to_list()[0])if c]})
    # The three eta*l polynomials occur verbatim among the nine pure forms.
    assert all(mul({(eta,):F(1)},p)in[p for _,p in pure_targets]for p in linears)
    # Independent reconstruction of the complete pure-four-variable rowspace.
    outside=[j for j,m in enumerate(mons)if not set(m)<=set(pure)]
    pure_coefficients=M.extract(list(range(57)),outside).transpose().nullspace()
    assert pure_coefficients.shape==(9,57)
    assert rank(pure_coefficients*M)==9
    rows=[p.copy()for p in raw]
    initial_rules=[(eta,{():F(1)})]+[(ix[name],{():F(c)})for name,c in [('z11',4),('z12',6),('z23',-4)]]
    for i,rhs in initial_rules:rows=[substitute(p,i,rhs)[0]for p in rows]
    prepared=[p.copy()for p in rows];eliminated={i:rhs for i,rhs in initial_rules};chain=[]
    while any(rows):
        options=[]
        for j,p in enumerate(rows):
            for m,c in p.items():
                if len(m)!=1:continue
                i=m[0]
                if i in eliminated or any(i in k for k in p if k!=m):continue
                options.append((len(p),i,j,c))
        assert options,'No polynomial triangular pivot; this producer does not use GB.'
        _,i,j,c=min(options);chosen=rows[j]
        rhs={m:-a/c for m,a in chosen.items()if m!=(i,)}
        assert chosen==scale(add({(i,):F(1)},rhs,F(-1)),c)
        divisions=[];after=[]
        for k,p in enumerate(rows):
            v,d=substitute(p,i,rhs);after.append(v)
            if d:divisions.append([k,pack(d)])
        chain.append({'variable':names[i],'index':i,'source_generator':j,
                      'source_id':f'Q_raw{qindices[j]}'if j<57 else f'G_raw{g["cubic_generators"][j-57]["source"]}',
                      'pivot_coefficient':str(c),'selected_polynomial':pack(chosen),
                      'rhs':pack(rhs),'division_cofactors':divisions,
                      'nonzero_rows_after':sum(bool(p)for p in after),'terms_after':sum(map(len,after))})
        rows=after;eliminated[i]=rhs
        assert len(chain)<=31 and sum(map(len,rows))<100000
    free=[i for i in range(35)if i not in eliminated]
    expected=['z22','z24','z27','z28','z29','z30','z31','z32']
    assert [names[i]for i in free]==expected and len(chain)==23
    # Reverse substitution gives the entire exact graph on the final eight
    # variables, rather than leaving intermediate RHS expressions unresolved.
    graph={i:{(i,):F(1)}for i in free}
    for i,rhs in reversed(list(eliminated.items())):
        graph[i]=compose(rhs,graph)
    assert set(graph)==set(range(35)) and all(all(j in free for j in m)for p in graph.values()for m in p)
    assert all(len(m)<=1 for p in graph.values()for m in p)
    assert not graph[0]and graph[eta]=={():F(1)}
    assert all(not compose(p,graph)for p in raw+allq)
    # Homogenize the final graph.  The ordinary degree is distinct from the
    # older auxiliary weight.  Eta and the eight free coordinates all have
    # ordinary degree one.
    hgraph={i:{(eta,)if not m else m:c for m,c in p.items()}for i,p in graph.items()}
    assert hgraph[eta]=={(eta,):F(1)}
    assert all(not compose(p,hgraph)for p in raw+allq)
    equations=[add({(i,):F(1)},hgraph[i],F(-1))for i in range(35)if i!=eta and i not in free]
    assert len(equations)==26 and all(equations)
    # Compare with the already certified exact smooth9D T39 quotient.  The
    # earlier packet proves ALL350 actual equations vanish, with the nine
    # displayed parent coordinates as an exact left inverse.  Here we freshly
    # differentiate its DAG and compare every one of the35x9 entries.
    hpath=ROOT/'reports/Q_LOCAL_Q51_QUOTIENT_20261103_v63/certificate.json'
    tpath=ROOT/'reports/Q_LOCAL_Q51_QUOTIENT_20261103_v63/T39-nine-parameter-model.json.gz'
    t_report_path=ROOT/'reports/Q_LOCAL_Q51_QUOTIENT_20261103_v63/report.json'
    old=json.loads(hpath.read_text());tm=json.loads(gzip.decompress(tpath.read_bytes()));oldreport=json.loads(t_report_path.read_text())
    assert tm['exact_QQ9_all350_and45_zero']and tm['physical_free_coordinate_left_inverse']
    assert oldreport['all_passed']and oldreport['T39_all350_and45_exact_QQ9_zero']
    assert oldreport['T39_explicit_physical_coordinate_left_inverse']
    hfree=free+[eta];assert tm['input_names']==[names[i]for i in hfree]
    values=[]
    for op,arg in tm['nodes']:
        if op=='q':value=[F(arg)]+[F(0)]*9
        elif op=='v':value=[F(0)]+[F(int(i==arg))for i in range(9)]
        elif op=='+':value=[sum((values[k][i]for k in arg),F(0))for i in range(10)]
        elif op=='*':
            value=[F(1)]+[F(0)]*9
            for k in arg:
                b=values[k];value=[value[0]*b[0]]+[value[0]*b[i]+b[0]*value[i]for i in range(1,10)]
        else:raise AssertionError(op)
        values.append(value)
    parent=[values[k].copy()for k in tm['parent35_input_nodes']];parent[0][0]-=F(2,9)
    assert all(not row[0]for row in parent)
    tangent=[[row[j+1]for j in range(9)]for row in parent]
    graph_tangent=[[hgraph[i].get((j,),F(0))for j in hfree]for i in range(35)]
    assert tangent==graph_tangent
    D=DomainMatrix([[qq(c)for c in row]for row in tangent],(35,9),Q)
    H=DomainMatrix([[Q.from_sympy(sp.Rational(c))for c in row]for row in old['Q51_jacobian_45x35']],(45,35),Q)
    assert rank(D)==9 and rank(H)==26 and (H*D).is_zero_matrix
    assert (D.extract(hfree,list(range(9)))-DomainMatrix.eye(9,Q)).is_zero_matrix
    payload={'variables':names,'eta':'z33','raw_quadratic_generator_indices':qindices,
             'constant_row_certificates':cofactor,'pure_four_variable_quadratic_rowspace_dimension':9,
             'prepared_74_rows':list(map(pack,prepared)),'triangular_chain':chain,
             'free_variables':expected,'dehomogenized_graph':{names[i]:pack(graph[i])for i in range(35)},
             'homogeneous_graph':{names[i]:pack(hgraph[i])for i in range(35)},
             'homogeneous_graph_26_equations':list(map(pack,equations)),
             'T39_tangent_input_order':tm['input_names'],
             'T39_DAG_first_order_parent35':[[str(c)for c in row]for row in parent],
             'T39_exact_linear_map_35x9':[[str(c)for c in row]for row in tangent]}
    cp=args.output/'certificate.json.gz'
    with gzip.open(cp,'wt')as handle:json.dump(payload,handle,separators=(',',':'))
    inputs=[Path(__file__),qpath,gpath,ROOT/'code/q_ordinary_cubic_initial_20261103_v63.py',
            ROOT/'reports/Q_ORDINARY_CUBIC_INITIAL_20261103_v63/report.json',
            hpath,tpath,t_report_path,ROOT/'reports/Q_LOCAL_Q51_QUOTIENT_20261103_v63_REVIEW.json']
    report={'id':'Q_TRUE_INITIAL_ETA_OPEN_20261103_V63','all_passed':True,'coefficient_ring':'QQ',
            'homogeneous_ideal':'J=(57 exact quadratic initials,17 exact new cubic initials) in35 variables',
            'eta':'z33','constant_QQ_row_certificates_count':len(cofactor),
            'eta_times_three_linear_forms_in_J':True,'three_linear_squares_in_J':True,
            'pure_four_variable_quadratic_rowspace_dimension':9,
            'polynomial_triangular_steps':len(chain),'every_row_division_identity_verified':True,
            'reverse_all74_generators_zero':True,'reverse_all350_raw_quadratic_rows_zero':True,
            'final_graph_is_linear':True,'u_zero':True,'dehomogenized_free_variables':expected,
            'original_J_eta_equals_linear_prime_graph_eta':True,
            'original_J_eta_reduced_and_smooth':True,'original_J_eta_dimension':9,
            'dehomogenized_J_eta1_polynomial_ring_dimension':8,
            'scheme_isomorphism':'(QQ[u,z0,...,z33]/J)[1/eta] = QQ[eta,1/eta,z22,z24,z27,z28,z29,z30,z31,z32]',
            'fresh_T39_DAG_first_order_matches_all35x9_linear_graph_entries':True,
            'Q51_linear_relation_matrix_rank':26,'T39_linear_map_rank':9,
            'exact_full_T39_quotient_dependency_reused':True,
            'full_ordinary_initial_ideal_eta_equals_T39_linear_prime_eta':True,
            'full_initial_eta_proof':'J subset in_m(I) subset P_T39; localized J_eta=P_eta, and fresh exact T39 tangent map gives P_T39=P.',
            'all_other_ordinary_tangent_cone_components_lie_eta_zero':True,
            'full_actual_completed_ideal_still_distinguished_from_ordinary_initial':True,
            'full_actual_R35_eta_ideal_equality_claimed':False,'complete_initial_ideal_globally_claimed':False,
            'certificate_sha256':hashlib.sha256(cp.read_bytes()).hexdigest(),
            'input_sha256':{str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest()for path in inputs},
            'seconds':time.monotonic()-started}
    (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({key:report[key]for key in ['all_passed','constant_QQ_row_certificates_count',
          'polynomial_triangular_steps','dehomogenized_free_variables','original_J_eta_equals_linear_prime_graph_eta','seconds']},indent=2))


if __name__=='__main__':main()
