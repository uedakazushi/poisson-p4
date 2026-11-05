#!/usr/bin/env python3
"""Exact upper619 source clearing; fixed-point checks only, no relation values."""
import argparse, hashlib, json, math, resource, signal, time
from fractions import Fraction as F
from pathlib import Path
import ambient_B18_ell2_T2_mu_20261105_v63 as H
from ambient_B18_ell2_T2_G3_active4_source8_fold_20261105_v63 import count

ROOT=H.ROOT
STEM='AMBIENT_B18_ELL2_T2_UPPER619_SOURCE_CLEARING_20261105_v63'
MIN='reports/AMBIENT_B18_ELL2_T2_MINIMAL6_SOURCE_OPERATOR_20261105_v63/'
SRC='reports/AMBIENT_B18_ELL2_HONEST_T2_SOURCE_OPERATOR_20261105_v63/'
ACT='reports/AMBIENT_B18_ELL2_HONEST_T2_ACTIVE_UNIT_TRANSPORT_20261105_v63/'
MU='reports/AMBIENT_B18_ELL2_T2_MU_20261105_v63/'
UP='reports/AMBIENT_B18_ELL2_T2_UPPER2_CUBIC_OUTPUT_SOURCE_20261105_v63/'

def main(out):
    start=time.monotonic();signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
    resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024))
    out.mkdir(parents=True,exist_ok=False);hs={}
    def field(p,k): hs[p]=H.sha(ROOT/p);return H.selected(p,k)
    def read(p): hs[p]=H.sha(ROOT/p);return json.loads((ROOT/p).read_text())
    def bind(p,k,shape):
        hs[p]=H.sha(ROOT/p)
        return {'path':p,'sha256':hs[p],'field':k,'shape':shape}
    up=read(UP+'fixed-upper2-one-cubic-output-identity-and-source-frontier.json')
    sy=list(map(H.parse,field(MU+'actual741-Y0-fixed-rank-and-T2-minimal-generators.json','actual_selected_Y0_columns102')))
    s=[H.parse(up['s0_on_F2_source102']),H.parse(up['s1_on_F2_source102'])]
    lambdas=up['lambda_two'];assert len(sy)==99 and len(s)==2
    # Positive denominator/absolute-numerator envelopes are computed BEFORE
    # either fixed signed fixture. This is not a remaining-source value query.
    maximum=0
    def mark(q):
        nonlocal maximum
        maximum=max(maximum,q[0].bit_length(),q[1].bit_length());return q
    def exact(x): x=F(x);return mark((x.denominator,abs(x.numerator)))
    def product(a,b):return mark((a[0]*b[0],a[1]*b[1])) if a[1] and b[1] else (1,0)
    def sum_bound(a,b):
        if not a[1]:return b
        if not b[1]:return a
        d=math.lcm(a[0],b[0]);return mark((d,a[1]*(d//a[0])+b[1]*(d//b[0])))
    envelopes=[]
    for row in s:
        for col in sy:
            q=(1,0)
            for i,a in row.items():q=sum_bound(q,product(exact(a),exact(col.get(i,F(0)))))
            envelopes.append(list(q))
    assert maximum<=384
    owned=count([up,sy,s,envelopes,lambdas])+16000
    assert owned<60000
    H.dump(out/'phase.json',{'phase':'source_and_positive_point_fixture_bound_frozen_before_signed_checks','positive_prior_bits':maximum,'owned_numeric_bound':owned,'remaining619_coefficients_queried':False})
    actual=0
    for row in s:
        for col in sy:
            q=F(0)
            for i,a in row.items():
                t=a*col.get(i,F(0));actual=max(actual,abs(t.numerator).bit_length(),t.denominator.bit_length())
                q+=t;actual=max(actual,abs(q.numerator).bit_length(),q.denominator.bit_length())
                assert actual<=maximum
            assert q==0
    # Nodes remain vector-/operator-typed. Existing complete source operators
    # are referenced by their literal field and source index, not re-evaluated.
    refs={
      'original_full_R0_X768':bind(SRC+'actual-honest-T2-459-by2070-factored-source-operator.json','source_X768_nodes',[56,768]),
      'actual_comparison_M':bind(ACT+'actual-honest-T2-132-by1186-factored-source-operator.json','actual_Mnew_node',[26,56]),
      'actual_Y768':bind(ACT+'actual-honest-T2-132-by1186-factored-source-operator.json','all768_cycle_lift_top_nodes',[102,768]),
      'complete_minimal6':bind(MIN+'actual-honest-T2-minimal6-by1060-factored-source-operator.json','complete_minimal_T2_node',[6,1060]),
      'first27_unit':bind(MIN+'actual-honest-T2-minimal6-by1060-factored-source-operator.json','first27_quotient_node',[105,132]),
      'next99_unit':bind(MIN+'actual-honest-T2-minimal6-by1060-factored-source-operator.json','next99_quotient_node',[6,105]),
      'point_chart':bind(MIN+'actual27-plus99-unit-complements-and-source-combinations.json','fixed_bottom30_inverse_columns',[30,30]),
      'original27':bind(MIN+'actual27-plus99-unit-complements-and-source-combinations.json','original27_X_source_columns',[27]),
      'star99':bind(MU+'actual741-Y0-fixed-rank-and-T2-minimal-generators.json','actual_selected_Y0_source_columns',[99]),
      'zero_origin741':bind('reports/AMBIENT_B18_ELL2_T2_MU_INPUT_GATE_20261105_v63/actual-T2-mu-input-source-and-cache-gate.json','all741_zero_origin_source_combinations',[741]),
      'actual619':bind('reports/AMBIENT_B18_ELL2_T2_KNOWN_FIRST_PROFILE_20261105_v63/actual-known-bottom-first-rank-and-corrected-source-kernel.json','all_zero_bottom_and_Y0_source_combinations',[619]),
      'fixed25':bind('reports/AMBIENT_B18_ELL2_T2_KNOWN_FIRST_PROFILE_20261105_v63/phase1-known-G3linear-fixed-source-rank.json','fixed_rank_receipts',[25]),
      'known_linear_NF':bind('reports/AMBIENT_B18_ELL2_T2_KNOWN_LINEAR_INPUT_GATE_20261105_v63/actual-complete-L48-and-degree2-rewrite-source-input.json','all264_NF_root_nodes',[6,44]),
    }
    # Test every bound source field exists without retaining its possibly-large
    # coefficient array. Streaming selected-field readers are used by H.
    for key in ['original27','star99','fixed25']:
        d=field(refs[key]['path'],refs[key]['field']);assert len(d)==refs[key]['shape'][0];del d
    fixture={'lambda_two':lambdas,'s0_s1_two_rows':[list(sorted((i,str(x)) for i,x in row.items())) for row in s],
      'all198_selected99_point_annihilations':True,'positive_envelopes':envelopes,
      'actual_max_bits':actual,'positive_prior_bits':maximum}
    nodes=[]
    def put(kind,**kw):nodes.append({'kind':kind,**kw});return len(nodes)-1
    x=put('external_actual_source_family',binding=refs['original_full_R0_X768'])
    m=put('coherent_comparison_operator',binding=refs['actual_comparison_M'],frame='M^A=M+G2*A;Y^A=Y+A*X;typed_actual_A_retained')
    y=put('coherent_actual_source_lifts',binding=refs['actual_Y768'],source=x,comparison=m)
    piv=put('actual_pivot_source_columns27',source=x,binding=refs['original27'])
    star=put('actual_zero_origin_star_sources99',source=x,kernel741_binding=refs['zero_origin741'],selected_kernel741_indices=refs['star99'])
    c=put('external_actual619_constant_combinations',binding=refs['actual619'],source=x)
    ac=put('actual_moving_pivot27_coordinates',source=c,pivots=piv,unit_binding=refs['first27_unit'],unknown_Yp_retained=True)
    astar=put('actual_moving_pivot27_coordinates',source=star,pivots=piv,unit_binding=refs['first27_unit'],unknown_Yp_retained=True)
    zc=put('subtract_actual_pivots',source=c,pivots=piv,coefficients=ac)
    zs=put('subtract_actual_pivots',source=star,pivots=piv,coefficients=astar)
    ybar=put('actual_lift_minus_Yp_times_moving_coordinates',lifts=y,source=c,coefficients=ac,unknown_Yp_retained=True)
    q=put('fixed25_G3_linear_correction',argument=['PiY0',['degree1',ybar]],binding=refs['fixed25'],full_domain_slots=918)
    b=put('first99_coordinates_of_actual_unit_inverse',unit_binding=refs['next99_unit'],argument=['Ybar_c_minus_G3J_q_and_hc',ybar,q],shape=[99,619],lazy=True)
    e=put('exact_actual_source_residual',formula='E=Zc-Zstar*b',Zc=zc,Zstar=zs,b=b)
    u=put('exact_actual_upper_residual',formula='U=Ybarc-G3J*q-Ybarstar*b',Ybarc=ybar,q=q,b=b,constant_and_linear_zero=True)
    lo=put('exact_actual_lower_residual',formula='B=hc-hstar*b',b=b,constant_and_linear_zero=True)
    a=put('combined_actual_pivot_coefficient',formula='A=ac-astar*b',ac=ac,astar=astar,b=b)
    p2=put('degree2_before_point_pivot_clearing',formula='Xc2-Xp1*a1-Xstar1*b1',source=c,pivots=piv,star=star,a1=['degree1',a],b1=['degree1',b])
    a2=put('point_pivot_projection_of_P2',P2=p2,binding=refs['point_chart'])
    b2=put('point_lower_projection_of_P2',P2=p2,binding=refs['point_chart'])
    h=put('literal_free3_injection',original56_axes=[23,29,30],formula='E2=H0*B2',B2=b2)
    known_cubic=put('retained_degree3_source_expression',formula='Xc3-Xp2*a1-Xstar2*b1-Xp1*A2',source=c,pivots=piv,star=star,a2=a2,b1=['degree1',b])
    outupper=put('two_upper_quadratic_coordinates',formula='lambda_t*(M0*retained_degree3+M1*H0*B2)_[r^2*ell]',lambda_two=lambdas,comparison=m,degree3=known_cubic,E2=h,source_cases=619,ordinary_coefficient=True,values_uncomputed=True)
    ep=out/'exact619-upper-source-clearing-and-order-frontier.json'
    H.dump(ep,{'source_bindings':refs,'operator_nodes':nodes,'upper_two_output_node':outupper,'fixed_point_fixture':fixture,
      'exact_identities':['G2*U=M*E','U0=U1=B0=B1=E0=E1=0','E2=H0*B2','M0*Xp0=0','M0*Xstar1=G2linear*Ystar0','s0*Ystar0=s1*Ystar0=0'],
      'dropped_terms_with_reasons':{'G2quadratic_U1':'U1=0 afteractual27/99/G3clearing','M2_E1':'E1=0','M0_Xp0_A3':'M0Xp0=0','M0_Xstar1_b2':'all198pointannihilationsandlambdaG2linearidentity'},
      'retained_unknowns':['b1 full918-slot source operator','A2=P27P2','Yp0 and allorder unit inverses through exact source definition','M1 in the same coherent normalized frame','honest quotient third r²ell','native L2rr*S1ell+L2rell*S1r'],
      'input_order_frontier':{'Xc3':['r²ell'],'weighted_Xp2_Xstar2':['r²','r*ell'],'Xp1':['r','ell'],'a1_b1':['r','ell'],'A2_B2':['r²','r*ell'],'possible_b1_adjoint_extra_inputs':'retain full918 domain; no plane-support or associated-graded-exactness assumption'},
      'all619_values_uncomputed':True,'no_new_physical_gauge_or_comparison_coefficients':True})
    report={'id':STEM,'all_passed':True,'exact_source_identity_and_fixed_point_fixture_only':True,'all198_star_point_annihilations':True,
      'unknown_Yp_b1_A2_M1_and_native_L2_retained':True,'remaining619_values_computed':False,'new_physical_or_gauge_coefficients':False,
      'positive_point_prior_bits':maximum,'actual_fixed_point_bits':actual,'owned_numeric_bound':owned+count(nodes)+count(refs),
      'input_sha256':hs,'evidence_sha256':{ep.name:H.sha(ep)},'code_sha256':H.sha(Path(__file__)),
      'note_sha256':H.sha(ROOT/'research/notes'/f'{STEM}.md'),'seconds':time.monotonic()-start,'independent_source_review_pending':True}
    assert report['owned_numeric_bound']<60000
    H.dump(out/'report.json',report)
    print(json.dumps({k:report[k] for k in ['all_passed','positive_point_prior_bits','actual_fixed_point_bits','owned_numeric_bound','seconds']}))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output.resolve())
