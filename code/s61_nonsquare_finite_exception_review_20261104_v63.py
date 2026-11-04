#!/usr/bin/env python3
"""Independent squarefree/filter certificates and actual28-dimensional L witness."""
import argparse,hashlib,json,resource,signal,subprocess,sys,time
from fractions import Fraction as F
from pathlib import Path
from s61_nonsquare_Q_component_review_20261104_v63 import Program
from s61_nonsquare_raw_H_B_curve_components_review_20261104_v63 import poly
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_RAW_H_B_FINITE_EXCEPTION_REDUCTION_20261104_v63_REVIEW'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def ep(p):return {int(i):F(v) for i,v in p if F(v)}
def ua(p,q):
    r=dict(p)
    for i,c in q.items():
        r[i]=r.get(i,F(0))+c
        if not r[i]:del r[i]
    return r
def us(p,c):return {i:v*c for i,v in p.items() if v*c}
def um(p,q):
    # Clear denominators before convolution: independent integer arithmetic.
    if not p or not q:return {}
    import math
    dp=math.lcm(*(c.denominator for c in p.values()));dq=math.lcm(*(c.denominator for c in q.values()));r={}
    for i,c in p.items():
        for j,v in q.items():r[i+j]=r.get(i+j,0)+int(c*dp)*int(v*dq)
    return {i:F(c,dp*dq) for i,c in r.items() if c}
def main(out,work):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));sys.set_int_max_str_digits(0)
    out.mkdir(parents=True,exist_ok=False);work.mkdir(parents=True,exist_ok=False);b=ROOT/'reports';pre='S61_NONSQUARE_RAW_H_B_';end='_20261104_v63';hashes={}
    def bind(p):hashes[str(p.relative_to(ROOT))]=sha(p)
    def audit(p):
        d=read(p);bind(p)
        for n,v in d.get('input_sha256',{}).items():assert sha(ROOT/n)==v,n;hashes[n]=v
        for n,v in d.get('evidence_sha256',{}).items():assert sha(p.parent/n)==v,n;bind(p.parent/n)
        return d
    for name in ('EXCEPTION_SQUAREFREE','Q_EXCEPTION_NORM_FILTER','L_FINITE_UNIT_FIXED'):audit(b/(pre+name+end)/'report.json')
    for name in ('FINITE_BAD_LOCUS','Q_COMPONENT_FINITE','CURVE_COMPONENTS'):
        assert audit(b/(pre+name+end+'_REVIEW')/'report.json')['all_passed']
    sq=read(b/(pre+'EXCEPTION_SQUAREFREE'+end)/'native-certificate.json');sqsummary=read(b/(pre+'EXCEPTION_SQUAREFREE'+end)/'certificate.json')
    fil=read(b/(pre+'Q_EXCEPTION_NORM_FILTER'+end)/'native-certificate.json')
    primary=read(b/(pre+'FINITE_BAD_LOCUS'+end)/'certificate.json')
    qfield=read(b/(pre+'COMPONENT_EXACT_QUOTIENT_Q'+end)/'certificate.json');qroots=read(b/(pre+'Q_ROOT_ALL_MINORS'+end)/'native-certificate.json')
    P=Program();zero=P.source({});one=P.source({(0,0,0):F(1)})
    up=lambda p:{(0,int(i),0):F(v) for i,v in p if F(v)}
    reg=lambda p:P.source(up(p))
    der=lambda p:[[int(i)-1,str(int(i)*F(v))] for i,v in p if int(i) and F(v)]
    product=lambda regs:_product(P,regs,one)
    degree=lambda p:max(ep(p),default=-1)
    counts=[]
    for index,(source,rec) in enumerate(zip([c['e_polynomial_factors'] for c in primary['all_B_bad_contained_in_three_covers'][1:]],sq['groups'])):
        assert len(source)==len(rec['normalized_input_factors'])==len(rec['original_input_leading_constants'])
        for p,m,lc in zip(source,rec['normalized_input_factors'],rec['original_input_leading_constants']):
            assert F(lc)==ep(p)[degree(p)] and ep(m)[degree(m)]==1
            P.equal(P.op('M',reg(m),reg([[0,lc]])),reg(p))
        prod=product([reg(p) for p in rec['normalized_input_factors']]);P.equal(prod,reg(rec['monic_product_P']))
        assert ep(der(rec['monic_product_P']))==ep(rec['derivative_Pprime'])
        pp,g,s,dp=map(reg,[rec['monic_product_P'],rec['derivative_gcd_G'],rec['squarefree_S_P_over_G'],rec['derivative_Pprime']])
        bez=rec['P_Pprime_Bezout'];P.equal(P.op('A',P.op('M',reg(bez[0]),pp),P.op('M',reg(bez[1]),dp)),g)
        P.equal(P.op('M',g,s),pp);P.equal(P.op('M',g,reg(rec['Pprime_quotient_by_G'])),dp)
        assert ep(der(rec['squarefree_S_P_over_G']))==ep(rec['S_derivative'])
        bez=rec['S_Sprime_Bezout1'];P.equal(P.op('A',P.op('M',reg(bez[0]),s),P.op('M',reg(bez[1]),reg(rec['S_derivative']))),one)
        current=pp
        for step in rec['root_support_removal_steps']:
            factor=reg(step['removed_factor']);after=reg(step['remaining_after']);P.equal(P.op('M',factor,after),current)
            P.equal(P.op('M',factor,reg(step['S_quotient'])),s);current=after
        P.equal(current,one);counts.append(len(rec['root_support_removal_steps']))
        expected=sqsummary['squarefree_L'if index==0 else'squarefree_Q'];assert ep(rec['squarefree_S_P_over_G'])==ep(expected)
    assert [degree(g['squarefree_S_P_over_G']) for g in sq['groups']]==[9,39] and counts==[29,59]
    # Protect all Q denominator fibers before filtering the complement.
    S39=sqsummary['squarefree_Q'];srcE=qfield['finite_e_exception_polynomials']+[[[3,'486'],[2,'243'],[1,'-297'],[0,'-54']]]
    monic=lambda p:[[i,str(c/ep(p)[degree(p)])] for i,c in ep(p).items()]
    E=product([reg(monic(p)) for p in srcE]);P.equal(E,reg(fil['denominator_exception_product']))
    exc,res=reg(fil['protected_exception_gcd']),reg(fil['initial_generic_root_part']);P.equal(P.op('M',exc,res),reg(S39))
    P.equal(P.op('M',exc,reg(fil['denominator_product_quotient_by_exc'])),E)
    bez=fil['protected_gcd_Bezout'];P.equal(P.op('A',P.op('M',reg(bez[0]),reg(S39)),P.op('M',reg(bez[1]),E)),exc)
    ids=qroots['nonzero_remainder_indices'];assert len(ids)==len(fil['all15_generic_norm_gcd_steps'])==15
    for i,step in zip(ids,fil['all15_generic_norm_gcd_steps']):
        norm=[[m[1],v] for m,v in qroots['all22_evaluations'][i]['remainder_quadratic_norm_e']]
        assert ep(step['norm'])==ep(norm);P.equal(res,reg(step['before']));g=reg(step['gcd_after'])
        P.equal(P.op('M',g,reg(step['before_quotient'])),res);P.equal(P.op('M',g,reg(step['norm_quotient'])),reg(norm))
        bez=step['Bezout'];P.equal(P.op('A',P.op('M',reg(bez[0]),res),P.op('M',reg(bez[1]),reg(norm))),g);res=g
    P.equal(res,one);P.equal(res,reg(fil['final_generic_part']));safe=reg(fil['new_safe_squarefree_cover'])
    P.equal(P.op('M',exc,res),safe);P.equal(P.op('M',safe,reg(fil['old_cover_quotient_by_new'])),reg(S39))
    assert degree(fil['protected_exception_gcd'])==degree(fil['new_safe_squarefree_cover'])==37 and degree(fil['initial_generic_root_part'])==2
    # Confirm the residual h1/k1 leading coefficients without collapsing a.
    old=read(b/(pre+'PROJECTION_CURVE_RESIDUAL_INPUT'+end)/'certificate.json')
    for source,record in zip([old['oldh_quotient_h1'],old['newK_quotient_k1']],sqsummary['old_residual_LC_preflight']):
        pp=poly(source);d=max(m[1] for m in pp);lc={m:c for m,c in pp.items() if m[1]==d}
        assert set(lc)=={(0,d,0)} and lc[(0,d,0)]==F(record['LC_is_QQconstant_global_unit'])!=0
        assert {m:c/lc[(0,d,0)] for m,c in pp.items()}==poly(record['monic_normalized'])
    # One fixed ordinary QQ[e] identity program. a never occurs in this phase.
    inp=out/'fixed-univariate-identities.txt';inp.write_text('0\n0\n'+str(len(P.lines))+'\n'+'\n'.join(P.lines)+'\n');assert inp.stat().st_size<32*1024*1024
    cpp=ROOT/'code/s61_nonsquare_Q_component_review_20261104_v63.cpp';exe=work/'verify-univariates'
    subprocess.run(['bash','-c','source /workspace/cloud-setup/activate.sh; c++ -O2 -std=c++17 "$1" -o "$2" $(pkg-config --cflags --libs flint) -lgmp','compile',str(cpp),str(exe)],check=True,timeout=12)
    with(out/'native.log').open('w')as log:
        subprocess.run(['bash','-c','source /workspace/cloud-setup/activate.sh; exec "$1" "$2"','check',str(exe),str(inp)],stdout=log,stderr=subprocess.STDOUT,check=True,timeout=max(1,35-(time.monotonic()-start)))
    native=read(out/'native.log');assert native['all_passed']and native['field_inverse_or_GCD_calls']==0
    univariate_seconds=time.monotonic()-start
    # Independent finite algebra E7[d]/H, entirely over QQ.
    n=read(b/(pre+'L_FINITE_UNIT_FIXED'+end)/'native-certificate.json');assert n['membership_solved']
    aL={2:F(-9,50),1:F(-3,50),0:F(3,25)};current=ep(sqsummary['squarefree_L'])
    for step in n['a_zero_clopen_removal_steps']:
        factor,after=ep(step['removed_factor']),ep(step['remaining_after'])
        assert um(factor,after)==current and um(factor,ep(step['aL_quotient']))==aL;current=after
    modulus=ep(n['remaining_e_modulus']);assert current==modulus and max(modulus)==7 and modulus[7]==1
    bez=list(map(ep,n['modulus_aL_Bezout1']));assert ua(um(modulus,bez[0]),um(aL,bez[1]))=={0:F(1)}
    def reduce_e(p):
        p=dict(p)
        while p and max(p)>=7:
            degree=max(p);c=p[degree];p=ua(p,{degree-7+j:-c*v for j,v in modulus.items()})
        return p
    def em(p,q):return reduce_e(um(p,q))
    onee={0:F(1)};al=reduce_e(aL);ainv=ep(n['aL_inverse']);assert em(al,ainv)==onee
    hm=read(b/(pre+'PAIR_MONIC_COMPLETION'+end)/'certificate.json');H={(m[0]+5,m[1],m[2]):v for m,v in poly(hm['monic_H_laurent']).items()}
    minors=read(b/(pre+'35_MINORS'+end)/'certificate.json');Pi=list(map(poly,minors['localized_ideal_generators_monic_QQ_after_namedunit_and_QQ_only_scaling']));U=poly(minors['FULLU_product'])
    source=[*Pi,U,H];ma=max(m[0] for p in source for m in p);me=max(m[1] for p in source for m in p);md=max(m[2] for p in source for m in p)
    ap=[onee];ee=[onee]
    for j in range(ma):ap.append(em(ap[-1],al))
    for j in range(me):ee.append(em(ee[-1],{1:F(1)}))
    hc=[]
    for j in range(5):
        v={}
        for m,c in H.items():
            if m[2]==j:v=ua(v,us(em(ap[m[0]],ee[m[1]]),c))
        hc.append(v)
    assert hc==list(map(ep,n['Hclear_coefficients']))
    inv5=onee
    for _ in range(5):inv5=em(inv5,ainv)
    h=[em(v,inv5) for v in hc];assert h==list(map(ep,n['H_normalized_monic_coefficients']))and h[4]==onee
    def az():return [{}for _ in range(4)]
    def ao():return [onee,{}, {},{}]
    def aa(p,q):return [ua(x,y)for x,y in zip(p,q)]
    def am(p,q):
        v=[{}for _ in range(7)]
        for i in range(4):
            for j in range(4):v[i+j]=ua(v[i+j],em(p[i],q[j]))
        for j in range(6,3,-1):
            c=v[j]
            for i in range(4):v[j-4+i]=ua(v[j-4+i],us(em(c,h[i]),-1))
        return v[:4]
    dp=[ao()];dd=az();dd[1]=onee
    for _ in range(md):dp.append(am(dp[-1],dd))
    def evaluate(p):
        v=az()
        for (a,e,d),c in p.items():
            t=us(em(ap[a],ee[e]),c);v=aa(v,[em(t,x)for x in dp[d]])
        return v
    images=[evaluate(p)for p in Pi];ui=evaluate(U)
    assert len(images)==22 and images==[[ep(p)for p in v]for v in n['all22_original_source_images']]
    assert ui==list(map(ep,n['U_image'])) and n['dimension_N']==28 and n['base_e_degree']==7
    basis=[]
    for j in range(4):
        for i in range(7):v=az();v[j]={i:F(1)};basis.append(v)
    cols=n['all22N_ideal_matrix_columns'];assert len(cols)==22*28
    nonzero=0
    for source_index,img in enumerate(images):
        for basis_index,v in enumerate(basis):
            actual=am(img,v);assert actual==list(map(ep,cols[28*source_index+basis_index]))
            nonzero+=sum(len(p)for p in actual)
    assert nonzero==n['matrix_nonzero_entries']==17222
    target=ao()
    for _ in range(28):target=am(target,ui)
    assert target==list(map(ep,n['target_U_power_N']))
    witness=[[ep(p)for p in v]for v in n['source22_cofactor_witness']];assert len(witness)==22
    dot=az()
    for c,img in zip(witness,images):dot=aa(dot,am(c,img))
    assert dot==target
    proof=dict(squarefree_support_exact_L9_Q39=True,all_source_products_derivatives_Bezout_and_support_removal_identities=True,
        native_univariate=native,univariate_seconds=univariate_seconds,Q37_denominator_roots_protected=True,
        only_generic_complement_degree2_norm_filtered=True,all15_source_norms_and_exact_GCD_certificates_checked=True,
        old_residual_h1_k1_e_leading_coefficients_QQ_units=True,
        L_a_zero_clopen9_to7_source_bound=True,L_a_inverse_Bezout_checked=True,L_actual_H_monic_and_all22_source_images_checked=True,
        L_all616_ideal_columns_and17222_coefficients_checked=True,L_22_original_cofactors_dot_equals_U28=True,
        L_entire_FULLU_bad_locus_empty=True,remaining_Q_geometric_cover_points_atmost=296,
        old_residual_geometric_bound=1600,combined_geometric_point_bound=1896,
        bound_not_scheme_length_no_globalB3_rawH7_p0_cubic57_quartic239_or_exhaustion=True,
        no_author_solver_new_GCD_point_prime_GB_or_field_inverse=True)
    (out/'independent-proof.json').write_text(json.dumps(proof,indent=2)+'\n')
    for p in (Path(__file__),cpp,ROOT/'code/s61_nonsquare_Q_component_review_20261104_v63.py',ROOT/f'research/notes/{ID}.md'):bind(p)
    report=dict(id=ID,all_passed=True,L_cover_FULLU_empty_proved=True,Q_protected_cover_degree=37,
        combined_geometric_B_bad_point_bound=1896,input_sha256=hashes,evidence_sha256={p.name:sha(p)for p in out.iterdir()},
        seconds=time.monotonic()-start,max_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        scope='Exact squarefree9/39 covers, protectedQ37, whole L finite-cover unit identity. Q37 and old finite residual remain; no global rank or scheme-length claim.')
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(dict(all_passed=True,seconds=report['seconds'],rss=report['max_rss_kib'],report_sha256=sha(out/'report.json'))))
def _product(P,regs,one):
    result=one
    for r in regs:result=P.op('M',result,r)
    return result
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);p.add_argument('--work',required=True,type=Path);a=p.parse_args();main(a.output,a.work)
