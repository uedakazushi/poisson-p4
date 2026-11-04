#!/usr/bin/env python3
"""One sequential unit-linear core elimination, retaining all seven residuals."""
import argparse,hashlib,itertools,json,math,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
import sympy as s
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_PROPER_CORE_ELIMINATION_20261104_v63'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def decode(terms,vars):return sum((s.Rational(n,d)*s.prod(v**k for v,k in zip(vars,m)) for m,n,d in terms),s.Integer(0))
def encode(poly):return [[list(m),int(s.numer(c)),int(s.denom(c))] for m,c in poly.terms() if c]
def main(out):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False)
    pp=ROOT/'reports/AMBIENT_B18_PROPER_RATIONAL_ELIMINATION_20261104_v63/actual12-flag-divided-polynomials.json';pr=pp.parent/'report.json'
    S=read(pp);assert read(pr)['all_passed']
    a,d,v2,v3,v4,v5,v6,v7,v8=vars=s.symbols('a d v2 v3 v4 v5 v6 v7 v8');R=(a,d,v2,v6)
    H=[decode(q['divided_primitive_polynomial'],vars) for q in S['all12_actual_remaining_and_divided_records']]
    HP=[s.Poly(q,*vars,domain=s.QQ) for q in H];D=s.Poly(1-d*d,*R,domain=s.QQ);P=s.Poly(d*d+5,*R,domain=s.QQ)
    origins=dict.fromkeys(vars,0);pivots=[v4,v3,v5,v7,v8];source_rows=[0,4,1,2,3]
    lin2=s.Matrix([[s.diff(H[i],x) for x in [v3,v4]] for i in [0,4]])
    assert all(s.Poly(H[i],v3,v4,v5,v7,v8).total_degree()<=1 for i in [0,4])
    assert s.factor(lin2.det()+16*(1-d*d)*(d*d+5))==0 and lin2.subs(origins).det()==-80
    tests=[];coefficients=[];constant_rows=[]
    # These are the only three degree/unit tests after the first fixed2x2.
    for k,(row,q) in enumerate(zip(source_rows,pivots)):
        f=s.Poly(H[row],q);degree=f.degree();coef=s.expand(f.nth(1)) if degree==1 else None
        if degree!=1 or not coef.subs(origins):
            blocked={'id':STEM,'all_passed':False,'status':'STOP_non_unit_linear_pivot','row':row,'pivot':str(q),'degree':degree,'tests':tests}
            (out/'report.json').write_text(json.dumps(blocked,indent=2)+'\n');print(json.dumps(blocked));return
        assert not(coef.free_symbols-set(R)),('coefficient depends on unsolved pivot',row)
        tail=s.expand(H[row]-coef*q);future=set(pivots[k+1:]);assert not(tail.free_symbols&future)
        coeff=s.Poly(coef,*R,domain=s.QQ);dc=pc=0
        while True:
            quotient,rem=coeff.div(D)
            if rem.is_zero:coeff=quotient;dc+=1
            else:break
        while True:
            quotient,rem=coeff.div(P)
            if rem.is_zero:coeff=quotient;pc+=1
            else:break
        assert coeff.total_degree()==0 and coeff.as_expr()!=0
        c0=coeff.as_expr();coefficients.append((c0,dc,pc));constant_rows.append(s.Poly(tail,*vars,domain=s.QQ))
        tests.append({'row':row,'pivot':str(q),'degree':degree,'coefficient':str(coef),
          'coefficient_at_origin':str(coef.subs(origins)),'constant_factor':str(c0),'D_power':dc,'P_power':pc,
          'future_pivots_absent':True})
    (out/'sequential-unit-preflight.json').write_text(json.dumps({'first2_linear_matrix':[[str(c) for c in row] for row in lin2.tolist()],
      'first2_determinant':'-16*(1-d^2)*(d^2+5)','origin_minor':'-80','tests':tests,
      'remaining_degree_unit_tests':3,'no_general_algebraic_root':True},indent=2)+'\n')
    roots={};results=[];source_terms=sum(len(p.terms()) for p in HP);peak=source_terms;max_product_bound=0
    def nterms(p):return 0 if p.is_zero else len(p.terms())
    def resident():return source_terms+sum(nterms(v['num']) for v in roots.values())+sum(nterms(v['poly']) for v in results)
    def multiply(A,B,live=0):
        nonlocal peak,max_product_bound
        if A.is_zero or B.is_zero:return s.Poly(0,*R,domain=s.QQ)
        degree=A.total_degree()+B.total_degree();bound=min(nterms(A)*nterms(B),math.comb(degree+4,4))
        total=resident()+live+nterms(A)+nterms(B)+bound
        max_product_bound=max(max_product_bound,bound);peak=max(peak,total)
        assert total<=150000,('STOP_pre_product_expanded_term_bound',total)
        C=A*B;assert resident()+live+nterms(C)<=150000;return C
    def powpoly(A,k,live):
        q=s.Poly(1,*R,domain=s.QQ)
        for j in range(k):q=multiply(q,A,live)
        return q
    def substitute(poly):
        # All non-base variables in this polynomial must already have roots.
        allp=[q for q in vars if q not in R]
        needed={q for q in allp if q in poly.as_expr().free_symbols};assert needed<=roots.keys(),(needed-roots.keys())
        md=mp=0
        for ex,co in poly.terms():
            md=max(md,sum(ex[vars.index(q)]*roots[q]['D'] for q in needed))
            mp=max(mp,sum(ex[vars.index(q)]*roots[q]['P'] for q in needed))
        acc=s.Poly(0,*R,domain=s.QQ)
        for ex,co in poly.terms():
            base_exp=tuple(ex[vars.index(q)] for q in R);term=s.Poly.from_dict({base_exp:co},R,domain=s.QQ)
            ed=ep=0
            for q in needed:
                k=ex[vars.index(q)]
                if k:
                    term=multiply(term,powpoly(roots[q]['num'],k,nterms(acc)),nterms(acc))
                    ed+=k*roots[q]['D'];ep+=k*roots[q]['P']
            unit=powpoly(D,md-ed,nterms(acc)+nterms(term));unit=multiply(unit,powpoly(P,mp-ep,nterms(acc)+nterms(term)),nterms(acc)+nterms(term))
            term=multiply(term,unit,nterms(acc));acc+=term
            assert resident()+nterms(acc)<=150000
        return acc,md,mp
    def cancel_known_units(num,dd,pp):
        while dd:
            q,r=num.div(D)
            if not r.is_zero:break
            num=q;dd-=1
        while pp:
            q,r=num.div(P)
            if not r.is_zero:break
            num=q;pp-=1
        return num,dd,pp
    for k,(q,row) in enumerate(zip(pivots,source_rows)):
        tail,ed,ep=substitute(constant_rows[k]);c0,cd,cp=coefficients[k]
        num=tail.mul_ground(-1/c0);num,ed,ep=cancel_known_units(num,ed+cd,ep+cp)
        assert num.eval(dict.fromkeys(R,0))==0
        roots[q]={'num':num,'D':ed,'P':ep,'row':row,'unit_coefficient':str(s.Poly(H[row],q).nth(1))}
    # The source5 ideal is the exact graph ideal: each row has one target
    # pivot with a checked unit coefficient after all earlier roots.
    identities=[]
    for row in source_rows:
        num,ed,ep=substitute(HP[row]);assert num.is_zero
        identities.append({'row':row,'whole_source_zero_verified':True,'D_clear_power':ed,'P_clear_power':ep})
    rest=[i for i in range(12) if i not in source_rows];residuals=[]
    for row in rest:
        num,ed,ep=substitute(HP[row]);num,ed,ep=cancel_known_units(num,ed,ep)
        if num.is_zero:constant=content=s.Integer(1);primitive=num
        else:
            constant=s.ilcm(*[s.denom(c) for c in num.coeffs()]);integer=num.mul_ground(constant)
            content=s.igcd(*[abs(int(c)) for c in integer.coeffs()]);primitive=integer.mul_ground(1/content)
        results.append({'poly':primitive});residuals.append({'source_H_row':row,'D_denominator_power':ed,'P_denominator_power':ep,
          'QQ_unit_content_numerator':int(content),'QQ_unit_content_denominator':int(constant),
          'actual_core_polynomial':encode(primitive),'whole_identity':'Hrow(roots)=QQ_unit*K/(D^dpow P^ppow)'})
    expanded=resident();assert expanded<=150000
    origin=dict.fromkeys(R,0);K=[r['poly'].as_expr() for r in results]
    assert all(q.subs(origin)==0 for q in K)
    J=s.Matrix([[s.diff(q,z).subs(origin) for z in R] for q in K]);rank=J.rank()
    factored=[str(s.factor(q)) for q in K]
    proof={'source9_variable_order':list(map(str,vars)),'actual_core4_variable_order':list(map(str,R)),
      'selected5_source_H_rows':source_rows,'selected5_pivot_order':list(map(str,pivots)),
      'first2_linear_determinant':'-16*(1-d^2)*(d^2+5)','first2_constant_minor':'-80','unit_tests':tests,
      'five_actual_rational_roots':{str(q):{'numerator':encode(r['num']),'D_power':r['D'],'P_power':r['P'],
        'equation_row':r['row'],'unit_coefficient':r['unit_coefficient']} for q,r in roots.items()},
      'denominator_unit_dictionary':{'D':{'polynomial':'1-d^2','origin':1},'P':{'polynomial':'d^2+5','origin':5}},
      'all5_whole_source_identities':identities,'remaining7_actual_core_polynomials':residuals,
      'actual_core_factored_expressions':factored,'origin_Jacobian7_by4':[[str(c) for c in row] for row in J.tolist()],'origin_rank':rank,
      'whole_E4_presentation':'QQ[[u1..u5,a,d,v2,v6]] / d*(u1..u5,K1..K7)',
      'whole_d_times_selected5_relations_retained':True,'core_quotient_only':'QQ[[a,d,v2,v6]]/(K1..K7)',
      'all5_unit_linear_graph_ideal_proofs_complete':True,
      'statistics':{'stored_expanded_terms':expanded,'max_simultaneous_preproduct_term_guard':peak,'max_single_product_term_bound':max_product_bound},
      'no_GB_localstandardbasis_CM_D0_depth_or_further_elimination':True}
    p=out/'actual4variable7core-polynomials.json';p.write_text(json.dumps(proof,separators=(',',':'))+'\n');assert p.stat().st_size<=32*1024*1024
    peer=ROOT/'reports/AMBIENT_B18_PROPER_RATIONAL_ELIMINATION_20261104_v63_REVIEW/report.json';assert read(peer)['all_passed']
    paths=[Path(__file__),ROOT/'research/notes'/f'{STEM}.md',pp,pr,peer]
    report={'id':STEM,'all_passed':True,'status':'all5_sequential_unit_linear_eliminations_and_actual4variable7core_input',
      'all5_whole_source_identities_verified':True,'remaining_sequential_unit_tests':3,'all_unit_values':[str(s.sympify(t['coefficient']).subs(origins)) for t in tests],
      'denominator_units':{'1-d^2':1,'d^2+5':5},'actual_core_variables':4,'actual_core_polynomials':7,
      'actual_core_factored_expressions':factored,'core_origin_rank':rank,
      'whole_E4_d_times_selected5_and7_residuals_retained':True,
      'GB_CM_D0_or_depth_claimed':False,'output_bytes':p.stat().st_size,**proof['statistics'],
      'input_sha256':{str(p.relative_to(ROOT)):sha(p) for p in paths},
      'evidence_sha256':{p.name:sha(p),'sequential-unit-preflight.json':sha(out/'sequential-unit-preflight.json')},
      'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
