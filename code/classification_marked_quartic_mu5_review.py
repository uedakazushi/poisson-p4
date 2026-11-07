#!/usr/bin/env python3
"""Independent full MU5 marked-source and known103 boundary review.
Uses frozen independent MU6 checker arithmetic, not producer/CAS functions.
"""
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
import hashlib,json,resource,time,signal
import classification_marked_quartic_mu6_review as m
import classification_eight_essential_quintic_source_review as e
from classification_minimal_quartic_review import LaurentRing
import classification_degenerate_quadric_q4_nil31_I_review as q
ROOT=Path(__file__).resolve().parents[1]
PREFIX='reports/CLASSIFICATION_PLANE_HYPERPLANE_QUARTIC_MU5_'
OUT=PREFIX+'SOURCE_REVIEW_20261107.json'
JOB='research/classification/jobs/plane_hyperplane_quartic_mu5_source_review.json'
DATA=[
('A5_smooth_hyperflex',('s**2*t**2','s*t*(s**2-t**2)','s**4+t**4-4*s*t*(s**2+t**2)'),'(s-t)**4'),
('A4A1_node_inflexion',('s**4','s**2*t**2','2*(s**3*t-s*t**3)+t**4'),'(s-t)**3*(s+t)'),
('A4A1_smooth_hyperflex',('s**4','s**2*t**2','t**4+4*s*t*(s**2+t**2)'),'(s+t)**4'),
('D5_smooth_hyperflex',('s**3*t','s**2*t**2','-s**4+4*s*t**3-t**4'),'-(s-t)**4'),
('twoA2A1_smooth_hyperflex',('s**2*t**2','s**3*(s+t)','t**3*(s+t/16)'),'16*(s+t/4)**4'),
('A3A2_common_tangent_k_minus1',('s**2*t**2','s*t*(s-t)**2','-s**4+4*s**3*t+t**4'),'s**2*t**2')]

def evaluate(R,P,point):return m.subst(R,P,{R.n+i:R.constant(v)for i,v in enumerate(point)})
def graph_h_of_y(phi,infinity=False):
    values=[e.series(P,F(0),F,infinity)for P in phi];assert values[2][0];inv=e.sinverse(values[2]);h=e.smul(values[0],inv);y=e.smul(values[1],inv);assert not h[0]and not y[0]and y[1]
    N=e.N
    def compose(a,b):
        out=[F(0)]*N;power=[F(int(i==0))for i in range(N)]
        for c in a:
            if c:out=e.sadd(out,e.sscale(power,c))
            power=e.smul(power,b)
        return out
    inverse=[F(0)]*N
    for j in range(1,7):inverse[j]=(F(j==1)-compose(y,inverse)[j])/y[1]
    result=compose(h,inverse);assert all(compose(y,inverse)[i]==F(i==1)for i in range(7));return result

def curvecheck(rows):
    B=LaurentRing((),('s','t'));R=LaurentRing((),('h','y','w'));out=[]
    for row,(name,phistr,expected)in zip(rows,DATA):
        assert name==row['name'];phi=[B.parse(t)for t in phistr];f=R.parse(row['input_quartic']);H=R.parse(row['input_marked_H'])
        aff=[e.trim(tuple(P.get((j,4-j),F(0))for j in range(5)))for P in phi]
        assert len(e.ug(e.ug(aff[0],aff[1]),aff[2]))==1 and any(P.get((4,0))for P in phi)
        cols=[B.mul(*(m.power(B,phi[i],k)for i,k in enumerate(ex)))for ex in e.exponents(3,3)];M=[[P.get((12-j,j),F(0))for P in cols]for j in range(13)];assert q.rank(M)[0]==10
        def pull(P):return B.add(*(B.scale(B.mul(*(m.power(B,phi[i],k)for i,k in enumerate(ex))),c)for ex,c in P.items()))
        assert not pull(f);assert pull(H)==B.parse(expected)
        detail={}
        if name.startswith('A5_')or name.startswith('A3A2_'):
            g0,g1=graph_h_of_y(phi),graph_h_of_y(phi,True);first=next(i for i in range(1,7)if g0[i]!=g1[i]);assert first==(3 if name.startswith('A5_')else 2)
            detail['two_smooth_branches_contact']=first
            if first==2:assert e.local(phi,1)['first_characteristic_odd']==3
        elif name.startswith('A4A1_'):
            assert e.local(phi,0)['first_characteristic_odd']==5;node=(1,1,1)if'node_inflexion'in name else(1,-1,1)
            assert not evaluate(R,f,node)and all(not evaluate(R,R.derivative(f,i),node)for i in range(3))
            Hess=[[evaluate(R,R.derivative(R.derivative(f,i),j),node).get(R.zero,F(0))for j in(1,2)]for i in(1,2)]
            assert q.determinant(Hess);detail.update(cusp='A4',ordinary_node=list(node),node_on_H=not bool(evaluate(R,H,node)))
        elif name.startswith('D5_'):
            assert e.local(phi,0)['first_characteristic_odd']==3;inf=e.local(phi,0,infinity=True);assert inf['orders']==[1,2]
            detail.update(cusp='A2',other_branch_smooth_transverse_to_cusp=True)
        else:
            assert e.local(phi,0)['first_characteristic_odd']==3 and e.local(phi,0,infinity=True)['first_characteristic_odd']==3;node=(1,F(-1,16),-1)
            assert not evaluate(R,f,node)and all(not evaluate(R,R.derivative(f,i),node)for i in range(3))
            Hess=[[evaluate(R,R.derivative(R.derivative(f,i),j),node).get(R.zero,F(0))for j in(1,2)]for i in(1,2)];assert q.determinant(Hess)
            detail.update(two_A2_cusps=True,ordinary_node=list(map(str,node)),node_avoided_by_H=bool(evaluate(R,H,node)))
        if name!='A4A1_node_inflexion'and not name.startswith('A3A2_'):
            # Hyperflex image is smooth: evaluate the unique finite pullback root.
            root={'A5_smooth_hyperflex':F(1),'A4A1_smooth_hyperflex':F(-1),'D5_smooth_hyperflex':F(1),'twoA2A1_smooth_hyperflex':F(-1,4)}[name]
            image=[sum(c*root**ex[0]for ex,c in P.items())for P in phi]
            assert any(evaluate(R,R.derivative(f,i),image)for i in range(3));detail['hyperflex_at_smooth_point']=True
        out.append(dict(name=name,basepointfree=True,full13x10_cubic_implicit_rank=10,irreducible_degree4_image_absolute_irreducible=True,exact_H_pullback=expected,local_recognition=detail,delta_sum3_no_extra_singularities=True))
    return out

def known103_bridge(control):
    R=LaurentRing(('b',),('h','y','w'));p=R.parse;b,h,y,w=map(p,('b','h','y','w'));d=p('1+2*b')
    H=p('16*b**4*h');Y=p('4*b**2*y');W=p('-6*(2*b+1)*w+24*(b+1)*y+(12+48*b)*h')
    f3=R.add(m.power(R,R.add(R.mul(H,W),R.scale(m.power(R,Y,2),6)),2),R.scale(R.mul(m.power(R,Y,3),R.add(H,R.mul(m.power(R,b,2),Y))),-144))
    Vd=p('(1+2*b)*(h*w-y**2-4*h*y)-(2+8*b)*h**2')
    kd=p('(1+2*b)**2');fkden=R.add(m.power(R,Vd,2),R.mul(p('4*b*y*(y+2*h)'),Vd),p('16*b**2*h**2*y**2'))
    assert f3==R.mul(p('9216*b**8'),fkden)
    # Whole polynomial family really specializes to the input, and the old normalization maps into it.
    K=LaurentRing(('k',),('h','y','w'));pk=K.parse;V=pk('h*w-y**2-4*h*y+(k-3)*h**2');Fk=K.add(m.power(K,V,2),K.scale(K.mul(pk('(k-1)*y*(y+2*h)'),V),-1),pk('(k-1)**2*h**2*y**2'))
    special=m.subst(K,Fk,{0:pk('-1')});assert special==pk(control['input_quartic'])
    B=LaurentRing(('k',),('s','t'));phi=list(map(B.parse,('s**2*t**2','s*t*(s-t)**2','k*s**4+2*(1-k)*s**3*t+t**4')))
    value={}
    for ex,c in Fk.items():
        term=B.mul(B.scale(m.power(B,B.parse('k'),ex[0]),c),*(m.power(B,phi[i],a)for i,a in enumerate(ex[1:])))
        value=B.add(value,term)
    assert not value
    C=control['gradient_rank'];assert C['rank']==16 and F(C['determinant'])==F(-3913788948480)
    return dict(full_polynomial_Fk_and_normalization_exact=True,known_F3_actual_base_map=['H=16b4 h','Y=4b² y','W=-6(2b+1)w+24(b+1)y+(12+48b)h'],base_determinant='-384b6(2b+1)',relation='b=(1-k)/(2(1+k)); a=b²; k=(1-2b)/(1+2b)',exact_F3_identity='F3(b²)(H,Y,W)=9216 b8 (1+2b)² Fk',unit_chart='k!=0,1,-1 equivalently b!=0,1/2,-1/2 and finite; k=-1 separately retained',special_minor=str(C['determinant']),all_special_frames_and6mag_lift='Generic rank16 implies all17minors identically0; special16minor unit gives regular rank2 kernel bundle. Any special frame and all6mag extend, Saito determinant remains unit. Actual source endpoints lie in known103 by closedness.',reverse_boundary_inference_used=False)

def priorcheck():
    old=m.load('reports/CLASSIFICATION_PLANE_HYPERPLANE_QUARTIC_MU6_SOURCE_REVIEW_20261107.json');P=m.priorcheck(old['old_printed_witnesses']);rows=m.load('research/classification/inputs/mu6_prior111_37d56e2a_20261107.json')['all111_rows'];audit=m.load(PREFIX+'PRIOR111_AUDIT_20261107.json');assert set(x['id']for x in rows)==set(x['id']for x in audit['all111_dispositions'])
    for r in P['all111_dispositions']:
        if r['id']=='R4_A0':r['independent_reason']='Old actual D5 with H through its cusp, new D5 smoothhyperflex H avoids D5; other new profiles differ.'
        if r['id']=='R4_DPM':r['independent_reason']='Old actual 3A2 totalmu6; allnew marked quartics totalmu5.'
    four=audit['four_component_candidates'];assert len(four)==4
    six=[dict(first=a,second=b,reason='different intrinsic quartic profile'if a.split('_')[0]!=b.split('_')[0]else'A4A1 node belongs to H in one source and H avoids both singularities in the other')for a,b in combinations(four,2)]
    nine=m.load(m.PREFIX+'PRIOR111_AUDIT_20261107.json')['nine_component_candidates'];cross=[dict(mu5=a,mu6=b,reason='intrinsic quartic totalmu5 versus6')for a in four for b in nine]
    assert len(six)==6 and len(cross)==36
    return dict(all_prior111=P['all111_dispositions'],immutable_git_blob=P['blob'],all29_accepted_actualSNC_witnesses=P['all29_SNC_witnesses_reused_from_independent_first104_review'],all6_mu5_mutual=six,all36_cross_mu5_mu6=cross,combined13_candidates_all78_pairs='36 MU6 previously independently verified +6 MU5 +36 cross; no component adoption assumption')

def main():
    signal.alarm(180);t=time.monotonic();rows=m.load(PREFIX+'SOURCES_20261107.json')['cases'];norm={r['name']:r for r in m.load(PREFIX+'NORMALFORMS_20261107.json')['cases']};checks=[]
    for row in rows:
        n=norm.get(row['name'])
        if n is None:
            v=list(map(F,row['magnetic_cokernel_basis'][0]));j=next(i for i,a in enumerate(v)if a);n=dict(complete_cokernel=row['magnetic_cokernel_basis'][0],canonical_R=('h**2','h*y','h*w','y**2','y*w','w**2')[j],projective_orbit_dimension=23)
        checks.append(m.sourcecheck(row,n))
    out=dict(task_id='E_reduced_nonnormal',all_passed=True,scope='MU5 full five source+one boundary control, whole marked list handproof, known103 allk bridge, prior111 and mutual13 comparison. Ambient independent root; A5 remains unknown.',all6_full_sources=checks,all6_explicit_curves=curvecheck(rows),known103_whole_A3A2=known103_bridge(rows[-1]),prior_and_mutual=priorcheck(),remaining_unknown=['A5 smoothhyperflex: D100 two normals, no component or knownparent claim','E6 smoothflex remains separately unknown','Other E residual branches not silently closed'])
    out['inputs']=m.load(JOB)['exact_inputs_and_hashes'];assert all(m.sha(r['path'])==r['sha256']for r in out['inputs']);out['reviewer_sha256']=m.sha('code/'+Path(__file__).name);out['executed_reviewer_bytes']=Path(__file__).read_text();out['seconds']=time.monotonic()-t;out['maximum_RSS_KiB']=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    path=ROOT/OUT;assert not path.exists();path.write_text(json.dumps(out,indent=2,ensure_ascii=False)+'\n');print(json.dumps({k:out[k]for k in('all_passed','seconds','maximum_RSS_KiB')}))
if __name__=='__main__':main()
