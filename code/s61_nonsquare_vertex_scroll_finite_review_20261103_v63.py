#!/usr/bin/env python3
"""Independent vertex/scroll intersection and finite quartic-defect audit."""
import argparse
from fractions import Fraction as Q
import itertools
import json
from pathlib import Path
import resource
import shutil
import signal
import tarfile
import time
import sympy as s
from s61_nonsquare_pair_empty_QQ_review_20261103_v63 import read, sha, decode, add, mul, scale, power, subtract
from s61_nonsquare_ternary_pencil_finite_review_20261103_v63 import determinant

ROOT = Path(__file__).resolve().parents[1]
STEM = 'S61_NONSQUARE_VERTEX_SCROLL_FINITE_20261103_v63_REVIEW'


def main(out):
    started=time.monotonic(); resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000)); signal.alarm(60)
    out.mkdir(parents=True,exist_ok=False)
    parent=ROOT/'reports/S61_NONSQUARE_VERTEX_SCROLL_FINITE_20261103_v63'; author=read(parent/'report.json')
    assert author['frozen'] and author['all_exact_QQ_actual_source_checks_passed']
    archive=ROOT/author['archive']; assert archive.stat().st_size==152942 and sha(archive)==author['archive_sha256']
    snap=out/'immutable-input'; snap.mkdir()
    with tarfile.open(archive) as tar:
        for member in tar:
            assert member.isfile(); dst=snap/member.name; assert dst.resolve().is_relative_to(snap.resolve())
            dst.parent.mkdir(parents=True,exist_ok=True)
            with tar.extractfile(member)as src,dst.open('wb')as target:shutil.copyfileobj(src,target)
    for rel,h in author['input_sha256'].items():assert sha(snap/rel if(snap/rel).exists()else ROOT/rel)==h,rel
    deps=[ROOT/f'reports/{stem}/report.json'for stem in [
        'S61_NONSQUARE_FOUR_QUADRIC_BINARY_REDUCTION_20261103_v63_REVIEW',
        'S61_NONSQUARE_PAIR_EMPTY_QQ_20261103_v63_REVIEW',
        'S61_NONSQUARE_TERNARY_PENCIL_FINITE_20261103_v63_REVIEW',
        'S61_NONSQUARE_QUAD78_CORRECTED_CONDITIONAL_20261103_v63',
        'S61_NONSQUARE_GCD_GATE_PROPORTIONAL_EXCLUSION_20261103_v63_REVIEW',
        'S61_NONSQUARE_FINITE_GCD_EXCEPTION_NAKAYAMA_20261103_v63_REVIEW',
        'S61_NONSQUARE_ZW_PRODUCT_RANK19_20261103_v63_REVIEW']]
    for p in deps:
        rr=read(p);assert rr['all_passed']
        for rel,h in rr['input_sha256'].items():assert sha(ROOT/rel)==h,rel
    assert read(deps[4])['no_proportionality_on_a_P_J_open']
    assert read(deps[5])['finite_V_algebra_without_flatness_reviewed']
    binary=read(snap/'work/nonsquare-v62/four-quadric-binary-certificate-v71/certificate.json')
    pair=read(snap/'work/nonsquare-v62/pair-finite-proof-v76/certificate.json')
    base=snap/'work/nonsquare-v62/vertex-curve-E5-v85'
    profile=read(base/'profile.json'); inter=read(base/'intersection-profile.json'); cert=read(base/'intersection-certificate.json')
    a,e,P,z=s.symbols('a e P z'); vs=(a,e,P); J=5*P-6*e-4; c=s.Integer(6075)
    pd=lambda p:{tuple(ex):Q(str(co))for ex,co in s.Poly(p,*vs,domain=s.QQ).terms()if co}
    cc=[[s.sympify(q,locals={'a':a,'e':e,'P':P})for q in row]for row in binary['C1_C2_four_coordinate_coefficients']]
    alpha,beta,gamma,delta=cc[0];eps,zeta,eta,theta=cc[1]
    assert s.cancel(beta-32*a**3*J**2/243)==0
    R1=s.Poly(s.cancel(P**2*(beta*eps-(zeta-40*z/c)*(alpha-z/c))),z)
    R2=s.Poly(s.cancel(P**2*(beta*(eta-theta*z)-(zeta-40*z/c)*(gamma-delta*z))),z)
    aa=[pd(R1.coeff_monomial(z**k))for k in [2,1,0]];bb=[pd(R2.coeff_monomial(z**k))for k in [2,1,0]]
    assert aa[0]=={(0,0,2):Q(-40,int(c)**2)}
    result=determinant([aa+[{}],[{}]+aa,bb+[{}],[{}]+bb],3)
    Eold=decode(profile['E5_original'],3)
    factor=mul(power(pd(a),9,3),mul(power(pd(J),6,3),power(pd(P),2,3)))
    expected=scale(mul(factor,Eold),Q(-67108864,4749088069944112812890625))
    assert result==expected
    # Exact natural-coordinate substitution, without dividing a new polynomial.
    T={(1,0,0):Q(1)};Pn={(0,1,0):Q(1),(1,0,0):Q(-8)};en={(0,0,1):Q(1,3),(1,0,0):Q(-20,3)}
    images=[mul(T,Pn),en,Pn]
    def compose(p,images,dim):
        result={}
        for ex,co in p.items():
            term={(0,)*dim:co}
            for image,k in zip(images,ex):term=mul(term,power(image,k,dim))
            result=add(result,term)
        return result
    E5=decode(profile['E5_natural'],3);D={(1,0):Q(55),(0,1):Q(-66),(0,0):Q(18)}
    N=scale(mul({(0,1):Q(1),(0,0):Q(2)},power({(0,1):Q(7),(0,0):Q(-11)},2,2)),Q(2))
    lift=lambda p:{(0,ex[0],ex[1]):co for ex,co in p.items()}
    assert E5==subtract(lift(N),mul(T,power(lift(D),2,3)))
    assert compose(Eold,images,3)==mul(power(Pn,3,3),E5)
    assert profile['E5_removed_known_P_power']==3 and Q(profile['E5_natural_QQcontent'])==1
    # Universal homogeneous cross-product: no smoothness assumption is used.
    al,be,ga,de,ep,ze,et,th,ss,tau=s.symbols('al be ga de ep ze et th ss tau')
    row1=[ss*al-tau/c,ss*be,ss*ga-tau*de];row2=[ss*ep,ss*ze-40*tau/c,ss*et-tau*th]
    cross=[s.expand(row1[1]*row2[2]-row1[2]*row2[1]),s.expand(row1[2]*row2[0]-row1[0]*row2[2]),s.expand(row1[0]*row2[1]-row1[1]*row2[0])]
    assert all(s.expand(sum(x*y for x,y in zip(row,cross)))==0 for row in [row1,row2])
    Ac,Xc,Wc=cross; assert s.expand(Wc.subs(ss,0)-40*tau**2/c**2)==0
    phi=[ss*Ac,ss*Xc,ss*Wc,-tau*Wc]
    C1=al*phi[0]+be*phi[1]+ga*phi[2]+de*phi[3];C2=ep*phi[0]+ze*phi[1]+et*phi[2]+th*phi[3]
    assert s.expand(C1-tau*Ac/c)==0 and s.expand(C2-40*tau*Xc/c)==0
    assert s.expand(-phi[3]-tau*Wc)==0
    # Two genuine vertex minors; no B*S or B*U coefficient enters.
    reorder=lambda p:{(ex[0],ex[2],ex[1]):co for ex,co in p.items()}
    rows=[[reorder(decode(p,3))for p in row]for row in pair['natural_rows']]
    equations=[]
    for idx,col in [(0,2),(2,4)]:
        q=reorder(decode(pair['linear_t_equations'][idx],3)); record=pair['linear_equation_normalizations'][idx]
        assert record['pair_coefficient_index']==col
        raw=subtract(mul(rows[0][5],rows[1][col]),mul(rows[1][5],rows[0][col]))
        assert raw==scale(mul(power(T,record['t_unit_power'],3),q),Q(record['QQcontent']))
        assert all(ex[0]<=1 for ex in q); equations.append(q)
    coefft=lambda q,k:{(ex[1],ex[2]):co for ex,co in q.items()if ex[0]==k}
    Fs=[]
    for q,normalization,encoded in zip(equations,inter['normalizations'],inter['generic_projection_polynomials']):
        assert normalization['removed_D_J_N_f']==[0,0,0,0]
        raw=add(mul(power(D,2,2),coefft(q,0)),mul(N,coefft(q,1)))
        f=decode(encoded,2);assert raw==scale(f,Q(normalization['QQcontent']));Fs.append(f)
    assert Fs==[decode(p,2)for p in cert['projection_polynomials']]
    assert list(map(len,Fs))==[26,21]
    content={}
    for (index,degree),weight in zip(cert['coefficient_labels'],cert['coefficient_Bezout']):
        cp={(ex[1],):co for ex,co in Fs[index].items()if ex[0]==degree}
        content=add(content,mul(cp,decode(weight,1)))
    assert content=={(0,):Q(1)}
    def evalv(p,value):
        ans={}
        for (ih,iv),co in p.items():ans=add(ans,{(ih,):co*value**iv})
        return ans
    vv=Q(cert['good_special_v']);special=[evalv(p,vv)for p in Fs]
    assert [max(ex[0]for ex in p)for p in special]==cert['special_H_degrees']==[max(ex[0]for ex in p)for p in Fs]==[3,2]
    u,w=[decode(p,1)for p in cert['special_Bezout']]
    assert add(mul(u,special[0]),mul(w,special[1]))=={(0,):Q(1)}
    branches=[];assert len(cert['D_zero_branches'])==2
    for branch,(vv,hh)in zip(cert['D_zero_branches'],[(Q(-2),Q(-30,11)),(Q(11,7),Q(120,77))]):
        assert Q(branch['v'])==vv and Q(branch['H'])==hh
        assert sum(co*hh**ex[0]*vv**ex[1]for ex,co in D.items())==0
        assert sum(co*hh**ex[0]*vv**ex[1]for ex,co in N.items())==0
        actual=[]
        for q in equations:
            p={}
            for (it,ih,iv),co in q.items():p=add(p,{(it,):co*hh**ih*vv**iv})
            actual.append(p)
        assert actual==[decode(p,1)for p in branch['polynomials']]
        u,w=[decode(p,1)for p in branch['Bezout']]
        assert add(mul(u,actual[0]),mul(w,actual[1]))=={(0,):Q(1)}
        branches.append({'v':str(vv),'H':str(hh),'literal_actual_linear_t_Bezout1':True})
    # Universal low-degree colon lifting, not an invocation of full CI.
    ell=s.symbols('ell');qs=s.symbols('q0:4');ls=s.symbols('l0:4');ks=s.symbols('k01 k02 k03 k12 k13 k23')
    cof=[ell*l for l in ls]
    for k,(i,j)in zip(ks,[(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]):cof[i]+=k*qs[j];cof[j]-=k*qs[i]
    assert s.expand(sum(q*h for q,h in zip(qs,cof))-ell*sum(q*l for q,l in zip(qs,ls)))==0
    evidence=dict(existing_finiteG2_review_frozen_current=True,
        actual_scroll_coefficient_beta_known_unit='32*a^3*J^2/243',actual_fixed_Sylvester_all_coefficients_checked_by_Fraction=True,
        exact_resultant_factor='-67108864*a^9*J^6*P^2*Eold/4749088069944112812890625',
        Eold_equals_P3_E5_under_natural_substitution=True,E5='2*(v+2)*(7*v-11)^2-t*(55H-66v+18)^2',
        universal_P1_cross_product_and_actual_minor_coverage_verified=True,
        cross_product_at_infinity_W='40*tau^2/6075^2',
        integral_coverage='Resultant nonzero => every ratio kernel rank2; cubic basepoint-free P1 map covers C and is birational; degree3 and CM/S1 imply reduced integral.',
        exact_vertex_columns=[2,4,5],chosen_minor_columns=[2,4],no_f_N_D_J_factor_removed_from_projection=True,
        projection_content_Bezout1_and_degree_preserving_v0_Bezout1_verified=True,projection_H_degrees=[3,2],
        D_zero_exhaustive_branches=branches,vertex_proportional_intersection_E5_dimension_at_most0=True,
        G1_failure_finite_on_ALL_aPJ='Proportionality excluded using a/P/J only; rank2 forces DB=N=M=0; constant e4/P28 leading coefficients make the QQ cover finite before FBC-unit/Nakayama stages.',
        raw_quad78_defect_dimension_at_most0_on_ALL_aPJ_via_corrected_criterion=True,
        degree3_colon_from_Koszul_lifts_not_full_CI=True,original_light_rank146_and_weighted152_outside_finite_on_FULLU=True,
        no_old_omitted_Bsquare_syzygy_inference_used=True,no_allpoint_lower_or_fullquartic_or_ambient_claim=True)
    (out/'independent-certificate.json').write_text(json.dumps(evidence,indent=2)+'\n');shutil.rmtree(snap)
    result=dict(id=STEM,all_passed=True,scroll_geometrically_integral_on_aPJ_D_E5_accepted=True,
        vertex_rankdrop_intersection_E5_finite_all_branches_accepted=True,
        corrected_conditional_raw78_defect_dimension_at_most0_on_ALL_aPJ=True,
        G1_finite_cover_needs_no_F_B0_C_unit_stage=True,
        original_light146_weighted152_outside_finite_FULLU_via_degree3_colon=True,
        actual_finite_exception_ranks_and_fullquartic_or_ambient_not_proved=True,
        old_false_omit_Bsquare_syzygy_inference_not_used=True,no_GB_new_field_prime_or_heavy_launch=True,
        seconds=time.monotonic()-started)
    paths=[Path(__file__),ROOT/f'research/notes/{STEM}.md',parent/'report.json',archive]+deps
    result['input_sha256']={str(p.relative_to(ROOT)):sha(p)for p in paths}
    result['evidence_sha256']={'independent-certificate.json':sha(out/'independent-certificate.json')}
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');signal.alarm(0)
    print(json.dumps({k:v for k,v in result.items()if k not in ['input_sha256','evidence_sha256']},indent=2),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
