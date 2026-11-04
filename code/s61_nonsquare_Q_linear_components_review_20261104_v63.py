#!/usr/bin/env python3
"""Independent whole-closed-fiber QQ algebra certificates for five linear factors."""
import argparse,hashlib,json,resource,signal,sys,time
from fractions import Fraction as F
from pathlib import Path
from s61_nonsquare_finite_exception_review_20261104_v63 import ep,ua,us,um
from s61_nonsquare_raw_H_B_curve_components_review_20261104_v63 import poly
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_RAW_H_B_Q_LINEAR_COMPONENTS_20261104_v63_REVIEW'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main(out):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));sys.set_int_max_str_digits(0)
    out.mkdir(parents=True,exist_ok=False);b=ROOT/'reports';pre='S61_NONSQUARE_RAW_H_B_';end='_20261104_v63';hashes={}
    def bind(p):hashes[str(p.relative_to(ROOT))]=sha(p)
    def audit(p):
        r=read(p);bind(p)
        for n,v in r.get('input_sha256',{}).items():assert sha(ROOT/n)==v,n;hashes[n]=v
        for n,v in r.get('evidence_sha256',{}).items():assert sha(p.parent/n)==v,n;bind(p.parent/n)
        return r
    for name in ['Q_EXCEPTION_FACTOR','Q_FIRST_COMPONENT_UNIT','Q_LINEAR_BATCH']:audit(b/(pre+name+end)/'report.json')
    assert audit(b/(pre+'FINITE_EXCEPTION_REDUCTION'+end+'_REVIEW')/'report.json')['all_passed']
    factor=read(b/(pre+'Q_EXCEPTION_FACTOR'+end)/'certificate.json')
    safe=read(b/(pre+'Q_EXCEPTION_NORM_FILTER'+end)/'native-certificate.json')['new_safe_squarefree_cover']
    assert ep(safe)==ep(factor['source_safe_SF37'])
    factors=[];product={0:F(factor['factor_constant'])};assert product=={0:F(1)}
    for i,r in enumerate(factor['irreducible_factors']):
        p3=poly(r['factor']);assert all(m[0]==m[2]==0 for m in p3) and r['exponent']==1
        p={m[1]:c for m,c in p3.items()};assert p[max(p)]==1
        factors.append(p);product=um(product,p)
        inv=factor['factor_inventory'][i];assert inv['factor_index']==i and inv['degree_e']==max(p) and poly(inv['factor'])==p3
    assert product==ep(safe) and [max(p)for p in factors]==[1,1,1,1,1,6,6,9,11]
    expected=[F(-1),F(2,3),F(-1,6),F(-5,9),F(-8,9)]
    assert [-p.get(0,F(0))/p[1]for p in factors[:5]]==expected
    qsource=read(b/(pre+'COMMON_CURVE_FACTOR'+end)/'certificate.json');Q=poly(qsource['factors'][0]['factor'])
    hsource=read(b/(pre+'PAIR_MONIC_COMPLETION'+end)/'certificate.json');H={(m[0]+5,m[1],m[2]):c for m,c in poly(hsource['monic_H_laurent']).items()}
    msource=read(b/(pre+'35_MINORS'+end)/'certificate.json');Pi=list(map(poly,msource['localized_ideal_generators_monic_QQ_after_namedunit_and_QQ_only_scaling']));U=poly(msource['FULLU_product'])
    folders=[b/(pre+'Q_FIRST_COMPONENT_UNIT'+end)]+[b/(pre+'Q_LINEAR_BATCH'+end)/f'component-{i}'for i in range(1,5)]
    cases=[]
    for index,(evalue,folder)in enumerate(zip(expected,folders)):
        report=audit(folder/'report.json');n=read(folder/'native-certificate.json');assert n['membership_solved']
        # This is Q modulo the selected linear IDEAL, not evaluation of a
        # generic field/root formula at an exceptional numerical sample.
        qa={}
        for (a,e,d),c in Q.items():assert d==0;qa=ua(qa,{a:c*evalue**e})
        assert max(qa)==2 and qa[2]==1
        current=qa
        for step in n['a_zero_clopen_removal_steps']:
            factor_a,remaining=ep(step['removed_factor']),ep(step['remaining_after'])
            assert um(factor_a,remaining)==current and um(factor_a,ep(step['a_quotient']))=={1:F(1)}
            current=remaining
        modulus=ep(n['remaining_a_modulus']);en=n['base_a_degree'];N=n['dimension_N']
        assert current==modulus and max(modulus)==en and modulus[en]==1 and N==4*en
        assert N==[4,4,4,4,8][index]
        if index==0:assert qa=={2:F(1),1:F(26,225)}
        bez=list(map(ep,n['modulus_a_Bezout1']));assert ua(um(modulus,bez[0]),um({1:F(1)},bez[1]))=={0:F(1)}
        def reduce_a(p):
            p=dict(p)
            while p and max(p)>=en:
                degree=max(p);c=p[degree];p=ua(p,{degree-en+i:-c*v for i,v in modulus.items()})
            return p
        def am0(p,q):return reduce_a(um(p,q))
        one={0:F(1)};aval=reduce_a({1:F(1)});ainv=ep(n['a_inverse']);assert am0(aval,ainv)==one
        source=[*Pi,U,H];ma=max(m[0]for p in source for m in p);md=max(m[2]for p in source for m in p)
        ap=[one]
        for _ in range(ma):ap.append(am0(ap[-1],aval))
        hc=[]
        for j in range(5):
            v={}
            for (a,e,d),c in H.items():
                if d==j:v=ua(v,us(ap[a],c*evalue**e))
            hc.append(v)
        assert hc==list(map(ep,n['Hclear_coefficients']))
        inv5=one
        for _ in range(5):inv5=am0(inv5,ainv)
        h=[am0(p,inv5)for p in hc];assert h==list(map(ep,n['H_normalized_monic_coefficients']))and h[4]==one
        def zero():return[{}for _ in range(4)]
        def unit():return[one,{}, {},{}]
        def add(p,q):return[ua(x,y)for x,y in zip(p,q)]
        def multiply(p,q):
            v=[{}for _ in range(7)]
            for i in range(4):
                for j in range(4):v[i+j]=ua(v[i+j],am0(p[i],q[j]))
            for j in range(6,3,-1):
                c=v[j]
                for i in range(4):v[j-4+i]=ua(v[j-4+i],us(am0(c,h[i]),-1))
            return v[:4]
        dp=[unit()];dd=zero();dd[1]=one
        for _ in range(md):dp.append(multiply(dp[-1],dd))
        def image(p):
            v=zero()
            for (a,e,d),c in p.items():
                cc=us(ap[a],c*evalue**e);v=add(v,[am0(cc,x)for x in dp[d]])
            return v
        images=[image(p)for p in Pi];ui=image(U)
        assert images==[[ep(p)for p in v]for v in n['all22_original_source_images']]and ui==list(map(ep,n['U_image']))
        basis=[]
        for j in range(4):
            for i in range(en):v=zero();v[j]={i:F(1)};basis.append(v)
        cols=n['all22N_ideal_matrix_columns'];assert len(cols)==22*N
        nonzero=0
        for source_index,img in enumerate(images):
            for basis_index,v in enumerate(basis):
                actual=multiply(img,v);assert actual==list(map(ep,cols[N*source_index+basis_index]))
                nonzero+=sum(len(p)for p in actual)
        assert nonzero==n['matrix_nonzero_entries']
        target=unit()
        for _ in range(N):target=multiply(target,ui)
        assert target==list(map(ep,n['target_U_power_N']))
        cofs=[[ep(p)for p in v]for v in n['source22_cofactor_witness']];assert len(cofs)==22
        dot=zero()
        for c,img in zip(cofs,images):dot=add(dot,multiply(c,img))
        assert dot==target
        cases.append(dict(factor_index=index,closed_factor_e_value=str(evalue),entire_quadratic_a_fiber_retained_before_known_a_clopen=True,
            remaining_a_degree=en,dimension_N=N,all22_source_images=True,ideal_columns=22*N,nonzero_matrix_entries=nonzero,
            all22_original_cofactors_dot_equals_U_power_N=True,FULLU_empty=True))
    remaining_degrees=[max(p)for p in factors[5:]];assert sum(remaining_degrees)==32
    proof=dict(source_safe_SF37_and_exact_all9_factor_product_checked=True,factor_degrees=[max(p)for p in factors],
        all_five_displayed_linear_factors_whole_closed_fibers_checked=True,cases=cases,
        higher_factor_irreducibility_not_needed_for_safe_remaining_degree32=True,
        remaining_Q_geometric_point_bound=8*32,old_finite_residual_point_bound=1600,total_geometric_bound=1856,
        no_sampling_new_beta_probe_point_prime_GCD_factor_search_GB_author_solver=True,
        no_globalB3_rawH7_p0_cubic57_quartic239_or_scheme_length_claim=True)
    (out/'independent-proof.json').write_text(json.dumps(proof,indent=2)+'\n')
    for p in (Path(__file__),ROOT/'code/s61_nonsquare_finite_exception_review_20261104_v63.py',ROOT/'code/s61_nonsquare_raw_H_B_curve_components_review_20261104_v63.py',ROOT/f'research/notes/{ID}.md'):bind(p)
    report=dict(id=ID,all_passed=True,all_five_Q_linear_components_FULLU_empty_proved=True,remaining_Q_factor_degrees=[6,6,9,11],
        combined_geometric_B_bad_point_bound=1856,input_sha256=hashes,evidence_sha256={p.name:sha(p)for p in out.iterdir()},
        seconds=time.monotonic()-start,max_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        scope='Whole five closed linear Q components eliminated by exact source U^N identities; higher32-degree cover and old finite residual remain. No point sampling or global-rank conclusion.')
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(dict(all_passed=True,seconds=report['seconds'],rss=report['max_rss_kib'],report_sha256=sha(out/'report.json'))))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args();main(a.output)
