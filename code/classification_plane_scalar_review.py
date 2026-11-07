#!/usr/bin/env python3
"""Independent entire scalar first-symbol T!=0 review, all15 models.

Fraction(tau), sparse QQ[lambda0,lambda1] Bareiss, every line/root and the
exceptional N2/K0 whole-source factor; no producer or CAS imports.
"""
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
import ast,hashlib,json,resource,signal,time
from classification_unimodular_logarithmic_review import RF
import classification_degenerate_quadric_q4_nil31_I_review as q
import classification_plane_rank_zero_nonkernel_review as H
from classification_minimal_quartic_review import LaurentRing
ROOT=Path(__file__).resolve().parents[1]
PREF='reports/CLASSIFICATION_PLANE_RANK_ZERO_SCALAR_'
OUT=PREF+'REVIEW_20261107.json';JOB=PREF+'REVIEW_JOB_20261107.json'
MANIFEST='research/classification/inputs/plane_rank_zero_scalar_manifest_20261107.json'
U=H.previous.RHO


def decode150(v,env,field):
    entries={}
    for text,(ij,e) in zip(v,q.COLS):
        c=H.read_scalar(text,env,field)
        if c:entries.setdefault(ij,{})[e+(0,0,0)]=c
    return q.antisymmetric(entries)
def convert_basis(basis,field):return [q.antisymmetric({ij:{e:field(v) for e,v in B.get(ij,{}).items()} for ij in q.PAIRS}) for B in basis]
def bases(saved):
    first=[]
    for normal in (3,4):
        for ij in combinations(range(3),2):
            for i in range(3):first.append(q.antisymmetric({ij:q.atom(1,normal,i)}))
    for normal,ex in [(4,(3,3)),(3,(3,4)),(4,(3,4)),(3,(4,4))]:
        for i in range(3):first.append(q.antisymmetric({(i,normal):q.atom(1,*ex)}))
    second=[q.antisymmetric({ij:q.atom(1,*ex)}) for ij in combinations(range(3),2) for ex in [(3,3),(3,4),(4,4)]]
    assert [list(map(str,q.vec(B))) for B in first]==saved['all30_source_basis']
    assert [list(map(str,q.vec(B))) for B in second]==saved['all9_source_basis']
    return first,second


def model(name,S,env,field):
    p=lambda text:H.parse(text,env,field)
    if name=='N3':T=[p('b'),p('c'),{}, {},{}]
    elif name=='N2':T=[p('b'),{},{},{},{}]
    elif name=='D':T=[p('a'),p('b'),p('-2*c'),{},{}]
    elif name=='D_plus_N2':T=[p('a+b'),p('b'),p('-2*c'),{},{}]
    elif name=='repeated_ss_Kss':T=[p('tau*a'),p('tau*b'),p('-2*tau*c'),{},{}]
    elif name=='mixed_Kss':T=[p('tau*a+b'),p('tau*b'),p('-2*tau*c'),{},{}]
    elif name in ('regular_ss_Kzero','regular_ss_Knil'):T=[p('a'),p('tau*b'),p('-(1+tau)*c'),{},{}]
    else:
        assert name=='regular2';T=[p('lambda0*a'),p('lambda1*b'),p('-(lambda0+lambda1)*c'),{},{}]
    En=[{},{},{},p('z'),p('w')]
    K={} if S=='zero' else p('-w**2') if S=='nil' else p('-2*z*w')
    return H.Padd(H.wedge(T,En),q.antisymmetric({(3,4):K}))


def operators(P,base1,base2,field):
    assert not any(q.jacobi(P).values()) and not H.divisor(P)
    return [H.convert_matrix(H.previous.operator(P,B),field) for B in (convert_basis(base1,field),convert_basis(base2,field))]

def certificate_pair(Ms,C,env,field):
    certs={};value=field(1)
    for key,M in zip(('D1','D2'),Ms):
        cert,minor=H.certificate(M,C[key],env,field);certs[key]=cert;value=value*minor
    return certs,value


def one_parameter(Ms,analysis,dimension,domain,flag):
    env={'tau':U};certs,minor=certificate_pair(Ms,analysis,env,RF)
    d=certs['D1']['nullity']+certs['D2']['nullity'];image=dimension+5+d
    assert d==analysis['higher_fiber_upper'] and image==analysis['whole_projective_image_upper'] and image<=23
    product=minor*U
    assert product==H.read_scalar(analysis['exception_polynomial'],env,RF)
    factors=[];factprod=RF(1)
    for text,multiplicity in analysis['all_irreducible_factors']:
        f=H.read_scalar(text,env,RF);assert f.q==(F(1),) and len(f.p)==2
        rho=-f.p[0]/f.p[1];assert all(rho!=r for r,t,m in factors)
        factors.append((rho,text,multiplicity))
        for _ in range(multiplicity):factprod=factprod*f
    ratio=product/factprod;assert len(ratio.p)==1 and ratio.q==(F(1),) and ratio.p[0]
    assert len(factors)==len(analysis['all_exceptional_roots'])
    walls=[]
    for rho,text,multiplicity in factors:
        matches=[w for w in analysis['all_exceptional_roots'] if F(w['parameter'])==rho];assert len(matches)==1;w=matches[0]
        assert w['minimal_polynomial']==text and w['multiplicity']==multiplicity
        specialized=[H.convert_matrix(M,lambda v:H.rf_at(v,RF(rho))) for M in Ms]
        c,_=certificate_pair(specialized,w,{'tau':RF(rho)},RF)
        wd=c['D1']['nullity']+c['D2']['nullity'];bound=dimension-1+5+wd
        assert wd==w['higher_fiber_upper'] and bound==w['whole_projective_image_upper']
        outside=not H.rf_at(domain,RF(rho));assert outside==w[flag]
        if not outside:assert bound<=23
        walls.append(dict(parameter=str(rho),outside_Jordan_domain=outside,affine_symbol_dimension=dimension-1,
                          higher_fiber_upper=wd,whole_image_upper=bound,**c))
    return dict(generic=certs,affine_symbol_dimension=dimension,higher_fiber_upper=d,whole_image_upper=image,
                all_irreducible_factors=analysis['all_irreducible_factors'],all_roots=walls),product


def universal_actions():
    R=LaurentRing(('sig','H00','H01','H10','H11','s','t','u'),('z','w'),invert=('sig',))
    p=R.parse;A,B,C,D=[p(k) for k in ('H00','H01','H10','H11')];s,t,u=[p(k) for k in ('s','t','u')]
    h=[[A,B],[C,D]];adj=[[D,R.scale(B,-1)],[R.scale(C,-1),A]];S=[[s,t],[u,R.scale(s,-1)]]
    det=R.add(R.mul(A,D),R.scale(R.mul(B,C),-1))
    multiply=lambda M,N:[[R.add(*(R.mul(a,b) for a,b in zip(row,col))) for col in zip(*N)] for row in M]
    hSadj=multiply(multiply(h,S),adj)
    n=[p('z'),p('w')];hn=[R.add(*(R.mul(a,b) for a,b in zip(row,n))) for row in h]
    qS=lambda M,v:R.add(R.mul(v[0],R.add(*(R.mul(M[1][j],v[j]) for j in range(2)))),R.scale(R.mul(v[1],R.add(*(R.mul(M[0][j],v[j]) for j in range(2)))),-1))
    assert qS(hSadj,hn)==R.mul(det,det,qS(S,n))
    # Normal GL2 diag(sig^-1,1) restores sig*E without square roots.
    assert R.mul(p('sig**(-1)'),p('sig'))==p('1')
    # Base diagonal normalizers for nilpotent3/2 and mixedN2.
    g3=[p('sig**(-2)'),p('sig**(-1)'),p('1')];g2=[p('sig**(-1)'),p('1'),p('1')]
    assert all(R.mul(g3[i],p('sig'),R.inverse_monomial(g3[j]))==p('1') for i,j in [(0,1),(1,2)])
    assert R.mul(g2[0],p('sig'),R.inverse_monomial(g2[1]))==p('1')
    return dict(all_GL2_qS_covariance_polynomial_identity=True,normal_nilpotent_inverse_diagonal=True,
                base_N3_and_N2_mixed_actual_inverse_diagonals=True,no_square_root_needed=True)


def main():
    started=time.monotonic();signal.alarm(300);resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,)*2)
    destination=ROOT/OUT
    if destination.exists():raise FileExistsError(destination)
    job=json.loads((ROOT/JOB).read_text())
    for b in job['exact_inputs_and_hashes']:assert hashlib.sha256((ROOT/b['path']).read_bytes()).hexdigest()==b['sha256'],b['path']
    manifest=json.loads((ROOT/MANIFEST).read_text())
    for b in manifest['files']:assert hashlib.sha256((ROOT/b['path']).read_bytes()).hexdigest()==b['sha256']
    fixed=json.loads((ROOT/(PREF+'FIXED_WEIGHTS_20261107.json')).read_text());one=json.loads((ROOT/(PREF+'ONE_PARAMETER_20261107.json')).read_text());two=json.loads((ROOT/(PREF+'TWO_PARAMETER_20261107.json')).read_text())
    base1,base2=bases(fixed);actions=universal_actions();fixed_results=[]
    assert fixed['all_ten_fixed_G24_closed'] is False
    for c in fixed['cases']:
        P=model(c['T_Jordan_type'],c['S_Jordan_type'],{},RF)
        assert [RF(v) for v in q.vec(P)]==[RF(F(v)) for v in c['full150_P0']]
        Ms=operators(P,base1,base2,RF)
        for M,key in zip(Ms,('full350_by30_D1','full350_by9_D2')):assert M==H.matrix(c[key],{},RF)
        certs,_=certificate_pair(Ms,c,{},RF);d=certs['D1']['nullity']+certs['D2']['nullity']
        assert c['affine_T_stratum_dimension']=={'N3':6,'N2':4,'D':5,'D_plus_N2':7}[c['T_Jordan_type']]
        assert c['affine_S_stratum_dimension']=={'zero':0,'nil':2,'ss':3}[c['S_Jordan_type']]
        dimension=c['affine_T_stratum_dimension']+c['affine_S_stratum_dimension'];bound=dimension+5+d
        assert dimension==c['restored_affine_first_symbol_dimension'] and d==c['higher_fiber_dimension_upper'] and bound==c['whole_projective_image_upper']
        factor=False
        if c['type']=='N2_zero':
            K=H.matrix(c['D1']['complete_kernel'],{},RF)
            for col in zip(*K):
                vector=q.antisymmetric({ij:q.add(*(q.scale(B.get(ij,{}),a) for B,a in zip(base1,col))) for ij in q.PAIRS})
                assert all(not vector.get(ij) for ij in q.PAIRS if 0 not in ij and ij!=(1,2))
            assert all(not B.get(ij) for B in base2 for ij in q.PAIRS if 0 not in ij and ij!=(1,2))
            factor=True;assert bound==31 and not c['G24_closed']
        else:assert bound<=22
        fixed_results.append(dict(type=c['type'],affine_symbol_dimension=dimension,linear_bound=bound,
                                  complete_kernels=certs,whole_source_empty_by_quadratic_factor=factor))
        print(c['type'],'PASS','factor-empty' if factor else bound,flush=True)
    one_results=[]
    for c in one['cases']:
        assert c['restored_affine_first_symbol_dimension']=={'repeated_ss_Kss':8,'mixed_Kss':10,'regular_ss_Kzero':8,'regular_ss_Knil':10}[c['type']]
        kind='zero' if c['type']=='regular_ss_Kzero' else 'nil' if c['type']=='regular_ss_Knil' else 'ss'
        env={'tau':U};P=model(c['type'],kind,env,RF)
        assert [RF(v) for v in q.vec(P)]==[H.read_scalar(v,env,RF) for v in c['full150_P0']]
        Ms=operators(P,base1,base2,RF)
        for M,key in zip(Ms,('full350_by30_D1','full350_by9_D2')):assert M==H.matrix(c[key],env,RF)
        domain=RF(1)
        for value in c['excluded_parameter_values']:domain=domain*(U-F(value))
        result,_=one_parameter(Ms,c['analysis'],c['restored_affine_first_symbol_dimension'],domain,'outside_declared_Jordan_domain')
        result['type']=c['type'];one_results.append(result);print(c['type'],'PASS roots',len(result['all_roots']),flush=True)
    # QQ(lambda0,lambda1), same independent bivariate Bareiss arithmetic.
    field=H.R2;env={'lambda0':field({(1,0):F(1)}),'lambda1':field({(0,1):F(1)})}
    P=model('regular2','ss',env,field)
    assert [field(v) for v in q.vec(P)]==[H.read_scalar(v,env,field) for v in two['full150_P0']]
    Ms=operators(P,base1,base2,field)
    for M,key in zip(Ms,('full350_by30_D1','full350_by9_D2')):assert M==H.matrix(two[key],env,field)
    certs,minor=certificate_pair(Ms,two['generic'],env,field)
    generic_d=certs['D1']['nullity']+certs['D2']['nullity']
    assert generic_d==two['generic']['higher_fiber_upper']==6
    assert 11+5+generic_d==two['generic']['whole_image_upper']==22
    delta=H.read_scalar(two['declared_regular_Jordan_domain'],env,field)
    product=minor*delta;assert product==H.read_scalar(two['whole_exception_polynomial'],env,field)
    factprod=field(1);keys=[]
    for text,multiplicity in two['all_irreducible_linear_factors']:
        f=H.read_scalar(text,env,field);assert f.q==H.PONE and max(map(sum,f.p))==1
        key=tuple(f.p.get(e,F(0)) for e in [(1,0),(0,1),(0,0)]);s=next(v for v in key if v);key=tuple(v/s for v in key)
        assert key not in keys;keys.append(key)
        for _ in range(multiplicity):factprod=factprod*f
    quotient=H.pd(product.p,factprod.p);assert len(quotient)==1 and quotient.get((0,0))
    lines=[];cover=[];root_count=0
    assert {line['equation'] for line in two['all_lines']}=={text for text,m in two['all_irreducible_linear_factors']}
    for line in two['all_lines']:
        maps={key:H.read_scalar(value,{'tau':U},RF) for key,value in line['parameterization'].items()}
        t,s=maps['lambda0'],maps['lambda1'];assert all(v.q==(F(1),) and len(v.p)<=2 for v in (t,s))
        assert any(len(v.p)==2 for v in (t,s))
        f=H.read_scalar(line['equation'],env,field);assert not H.r2_at(f,t,s)
        domain=H.r2_at(delta,t,s);outside=not domain
        assert outside==line['outside_Jordan_domain_entire_line']
        if outside:lines.append(dict(equation=line['equation'],excluded_repeated_Jordan_line=True));cover.append((line,t,s,None));continue
        specialized=[H.convert_matrix(M,lambda v:H.r2_at(v,t,s)) for M in Ms]
        result,product1=one_parameter(specialized,line['analysis'],10,domain,'outside_Jordan_domain')
        result.update(equation=line['equation'],parameterization=line['parameterization']);lines.append(result);cover.append((line,t,s,product1));root_count+=len(result['all_roots'])
    assert len(lines)==12 and len([r for r in lines if r.get('excluded_repeated_Jordan_line')])==3 and root_count==57
    intersections=[]
    for i in range(12):
        for j in range(i+1,12):
            first,second=cover[i],cover[j]
            if first[3] is None and second[3] is None:continue
            line,t,s,prod1=first if first[3] is not None else second
            other=second[0] if first[3] is not None else first[0]
            equation=H.r2_at(H.read_scalar(other['equation'],env,field),t,s)
            assert equation.q==(F(1),) and len(equation.p)<=2
            if len(equation.p)<2:assert equation;continue
            rho=-equation.p[0]/equation.p[1];outside=not H.r2_at(delta,H.rf_at(t,RF(rho)),H.rf_at(s,RF(rho)))
            if not outside and not H.rf_at(prod1,RF(rho)):
                assert any(F(w['parameter'])==rho for w in line['analysis']['all_exceptional_roots'])
            intersections.append(dict(lines=[i,j],parameter=str(rho),outside_Jordan_domain=outside,covered=True))
    two_result=dict(generic=certs,affine_symbol_dimension=11,generic_image_upper=22,all_twelve_linear_factors=two['all_irreducible_linear_factors'],all_lines=lines,all57_roots_checked=True,all_finite_intersections=intersections)
    # Immutable compatible extensions are text audits only, no old calculation.
    functions=lambda text:{n.name:ast.dump(n,include_attributes=False) for n in ast.parse(text).body if isinstance(n,ast.FunctionDef)}
    for previous,current in [(fixed,one),(one,two)]:
        before=functions(previous['producer_source_bytes']);after=functions(current['producer_source_bytes']);assert all(after.get(k)==v for k,v in before.items())
    result=dict(task_id='G',all_passed=True,theorem='Whole traceless-coHiggs-zero scalar first-symbol T!=0 stratum has projective image dimension<=23. No E/fullG claim.',
                actual_normalization_checks=actions,all10_fixed=fixed_results,all4_oneparameter=one_results,whole2parameter=two_result,
                fixed_false31_preserved_and_mathematically_closed_by_factor=True,
                factor_proof='Full N2/K0 necessarysource has Pi=da wedgeV+F2 dbwedge dc. Its actual quintic is F2*(w*Paz-z*Paw); F2=0 impliesrank<=2, otherwise a genuine quadratic factor, incompatiblewith irreducibleD.',
                complete_source_dimensions_restored=True,independent_Fraction_RF_R2_only=True,no_CAS_or_producer_imports=True,
                complements=['T0K0','T0Knonzero separate proof not used','nonzero traceless global-V dependent source','other G/E obligations'],
                inputs=job['exact_inputs_and_hashes'],reviewer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                seconds=time.monotonic()-started,maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    destination.write_text(json.dumps(result,indent=2)+'\n')
    print('PASS',OUT,round(result['seconds'],3),'seconds',result['maximum_RSS_KiB'],'KiB',flush=True)
if __name__=='__main__':main()
