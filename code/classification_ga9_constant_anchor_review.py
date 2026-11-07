#!/usr/bin/env python3
"""Independent complete Ga9 c1 anchor support/component/distinctness review.

No producer, CAS, Singular, NumPy, SymPy, legacy or producer-helper imports.
Only separate-author standard-library Fraction polynomial/linear arithmetic.
"""
import argparse
from fractions import Fraction as Q
import hashlib
import itertools as it
import json
from pathlib import Path
import resource
import subprocess
import time

from classification_minimal_quartic_review import LaurentRing
from classification_projected_veronese_zero_review import Ring,plus,times,scale,derivative,monomials,sparse_rref,det

ROOT=Path(__file__).resolve().parents[1]
PAIRS=list(it.combinations(range(5),2));TRIPLES=list(it.combinations(range(5),3));QUADS=list(it.combinations(range(5),4))
COLS=[(ij,e)for ij in PAIRS for e in monomials(5,2)]
ROWS=[(quad,e)for quad in QUADS for e in monomials(5,4)]


def anti(P):
    out={ij:P.get(ij,{})for ij in PAIRS}
    for i,j in PAIRS:out[j,i]=scale(out[i,j],-1)
    return out


def jacobi(P,R):
    strict={tri:plus(*(times(P.get((a,l),{}),derivative(P.get((b,c),{}),l))for a,b,c in [tri,(tri[1],tri[2],tri[0]),(tri[2],tri[0],tri[1])]for l in range(5)))for tri in TRIPLES}
    return {quad:plus(*(scale(times(R.var(R.names[quad[k]]),strict[quad[:k]+quad[k+1:]]),(-1)**k)for k in range(4)))for quad in QUADS}


def vec(P):return [P.get(ij,{}).get(e,Q(0))for ij,e in COLS]


def jvec(J):return [J[quad].get(e,Q(0))for quad,e in ROWS]


def transpose(M):return [list(row)for row in zip(*M)]


def certify(M):
    basis,_,rows=sparse_rref(M);columns=sorted(basis);rank=len(columns);rowset=rows[:rank]
    value=det([[M[i][j]for j in columns]for i in rowset]);assert value
    return dict(shape=[len(M),len(M[0])],rank=rank,rows=rowset,columns=columns,determinant=str(value),
        full_matrix_sha256=hashlib.sha256(json.dumps(M,default=str,separators=(',',':')).encode()).hexdigest())


def source_reconstruction(data):
    names=('r','b','alpha','p','q','A','D','G','H','k');R=LaurentRing(names,('h','x','y','z','w'),invert=('r',));p=R.parse
    assert list(names)==data['parameters']
    phi=p('r*y**4/4+A*y**3+(b-1/4)*y**2*w+D*y**2+(alpha-((1-3*b)*A+q*(b-1/4))/r)*y*w-(1-b)*(b-1/4)*w**2/(2*r)+G*y+H*w')
    Vy=p('-alpha-b*y');Vw=p('p+q*y+b*w+r*y**2');u=p('x-y*z')
    A2=R.add(phi,R.mul(Vy,R.derivative(phi,2)),R.mul(Vw,R.derivative(phi,4)),R.scale(p('k'),-1))
    assert all(sum(e[-5:])<=2 for e in A2)
    six=[R.add(R.derivative(phi,4),R.mul(p('y'),Vy)),R.add(u,A2,R.scale(R.mul(p('z'),Vy),-1)),
         R.add(R.scale(R.derivative(phi,2),-1),R.mul(p('y'),Vw),p('z')),R.scale(Vy,-1),p('1'),Vw]
    assert six==[p(t)for t in data['universal_brackets']]
    raw=anti({ij:value for ij,value in zip(it.combinations(range(1,5),2),six)})
    div=[R.add(*(R.derivative(raw.get((i,j),{}),i)for i in range(1,5)))for j in range(1,5)]
    assert div==[p('y'),{},p('1'),{}]
    pf=R.add(R.mul(six[0],six[5]),R.scale(R.mul(six[1],six[4]),-1),R.mul(six[2],six[3]))
    assert pf==R.add(p('k'),R.scale(u,-1),R.scale(phi,-1))==p(data['actual_pfaffian'].replace(' u ',' (x-y*z) '))
    strict=[R.add(*(R.mul(raw.get((a,l),{}),R.derivative(raw.get((b,c),{}),l))for a,b,c in [tri,(tri[1],tri[2],tri[0]),(tri[2],tri[0],tri[1])]for l in range(1,5)))for tri in it.combinations(range(1,5),3)]
    assert all(not value for value in strict)
    # The whole actual quintic belongs to the cube of the plane ideal(h,y).
    quintic={}
    for e,c in pf.items():
        f=list(e);assert f[10]==0;f[10]=5-sum(f[11:]);assert f[10]>=1
        quintic[tuple(f)]=c;assert f[10]+f[12]>=3
    anchor={name:Q(1 if name=='r'else 2 if name=='b'else 0)for name in names}
    def evaluate(poly,point=anchor):
        out={}
        for e,c in poly.items():
            for name,n in zip(names,e[:10]):c*=point[name]**n
            spatial=e[10:];hpower=2-sum(spatial[1:]);assert spatial[0]==0 and hpower>=0
            exponent=(hpower,)+spatial[1:];out[exponent]=out.get(exponent,Q(0))+8*c
        return {e:c for e,c in out.items()if c}
    point=anti({ij:evaluate(v)for ij,v in zip(it.combinations(range(1,5),2),six)})
    directions=[anti({ij:evaluate(derivative(v,k))for ij,v in zip(it.combinations(range(1,5),2),six)})for k in range(10)]
    assert len(sparse_rref([vec(v)for v in directions])[0])==10
    expected=['14*h*w-2*y**2','8*h*x+8*y*z+35*w**2','8*h*z-12*w*y','16*h*y','8*h**2','16*h*w+8*y**2']
    X=Ring(('h','x','y','z','w'))
    assert all(point[ij]==X.parse(t)for ij,t in zip(it.combinations(range(1,5),2),expected))
    assert all(not v for v in jacobi(point,X).values())
    previous=anti({ij:evaluate(v,anchor|{'k':Q(1)})for ij,v in zip(it.combinations(range(1,5),2),six)})
    inverse_translation=[X.var('h'),plus(X.var('x'),X.var('h')),X.var('y'),X.var('z'),X.var('w')]
    translated={}
    for ij in PAIRS:
        result={}
        for exponent,c in previous[ij].items():
            term=X.const(c)
            for value,power in zip(inverse_translation,exponent):term=times(term,X.power(value,power))
            result=plus(result,term)
        translated[ij]=result
    assert all(translated[ij]==point[ij]for ij in PAIRS)
    return point,directions,dict(full_ten_parameter_source_reconstructed_from_potential=True,all_six_saved_brackets_exactly_match=True,
        whole_source_parameter_domain='Gm(r) x A9(b,alpha,p,q,A,D,G,H,k)',whole_source_irreducible=True,
        all_four_strict_Jacobi_zero=True,whole_cofactor=['y','0','1','0'],whole_source_Pfaffian='k-(x-y*z)-phi',
        whole_actual_quintic_in_plane_ideal_cube=True,triple_plane='h=y=0',previous_k1_anchor_exact_x_minus_h_translation_to_k0=True,all_ten_actual_source_derivatives_retained=True,
        source_completeness_dependency='Previously independently accepted complete c=1 source reduction; the present checker reconstructs its actual tensor and all parameter derivatives, rather than reclassifying all Ga9 c types.')


def frame_review(point,directions,saved):
    X=Ring(('h','x','y','z','w'));Dcols=[]
    for ij,e in COLS:
        variation=dict(point);variation[ij]=plus(point[ij],{e:Q(1)});variation[ij[::-1]]=scale(variation[ij],-1)
        Dcols.append(jvec(jacobi(variation,X)))
    D=transpose(Dcols);E=[];G=[]
    for a,b in it.product(range(5),repeat=2):
        v={};o={}
        for i,j in PAIRS:
            v[i,j]=plus(times(X.var(X.names[i]),X.var(X.names[b]))if j==a else {},scale(times(X.var(X.names[j]),X.var(X.names[b])),-1)if i==a else {})
            o[i,j]=plus(point.get((b,j),{})if i==a else {},point.get((i,b),{})if j==a else {},scale(times(X.var(X.names[b]),derivative(point[i,j],a)),-1))
        E.append(vec(v));G.append(vec(o))
    actual=E+G+[vec(point)];Fcols=actual+[vec(v)for v in directions];F=transpose(Fcols)
    assert all(sum(a*b for a,b in zip(row,col))==0 for row in D for col in Fcols)
    certificates={name:certify(M)for name,M in [('D',D),('Euler',transpose(E)),('GL5_Euler_scalar',transpose(actual)),('whole_source',F)]}
    assert [certificates[name]['rank']for name in certificates]==[98,24,47,48]
    Ncols=[[Q(v.get(str(i),'0'))for i in range(150)]for v in saved['normal_vectors']]
    assert len(Ncols)==4 and all(sum(a*b for a,b in zip(row,col))==0 for row in D for col in Ncols)
    assert len(sparse_rref(transpose(Fcols+Ncols))[0])==52
    U=saved['implicit_columns'];rs=saved['implicit_rows'];assert len(U)==len(rs)==len(set(U))==len(set(rs))==98
    implicit_minor=det([[D[i][j]for j in U]for i in rs]);assert implicit_minor
    Lcols=[[Q(i==j)for i in range(150)]for j in U]+Ncols
    Lcert=certify(transpose(Lcols));directsum=certify(transpose(Fcols+Lcols))
    assert Lcert['rank']==102 and directsum['rank']==150
    weights=(0,4,1,3,2);character=-3
    relative=[sum(a*b for a,b in zip(ex,weights))-weights[i]-weights[j]-character for(i,j),ex in COLS]
    rowweights=[sum(a*b for a,b in zip(ex,weights))-sum(weights[i]for i in quad)-2*character for quad,ex in ROWS]
    assert (min(relative),max(relative))==(-4,10)and(min(rowweights),max(rowweights))==(-4,16)
    assert all(relative[j]==0 for j,a in enumerate(vec(point))if a)
    assert all(not a or rowweights[i]==relative[j]for i,row in enumerate(D)for j,a in enumerate(row))
    assert all(len({relative[i]for i,a in enumerate(col)if a})<=1 for col in Fcols)
    normalweights=[]
    for col in Ncols:
        ws={relative[i]for i,a in enumerate(col)if a};assert len(ws)==1;normalweights.append(next(iter(ws)))
    assert normalweights==[1,2,2,3]==saved['all_four_normals_weights']
    return Ncols,U,rs,relative,rowweights,dict(certificates=certificates,whole_D_times_full_source_zero=True,
        all_four_full_normal_vectors_in_kernel=True,F_plus_normals_rank=52,whole_kernel_dimension=52,
        implicit_full98_minor=str(implicit_minor),true_cut=Lcert,actual_F_plus_true_L=directsum,
        actual_independent_linear_cut_count=48,ambient_coordinate_weights=list(weights),fixed_tensor_character=character,
        normal_positive_weights=normalweights,all_full_350_by150_weights_compatible=True,whole_F_and_L_weight_stable=True,
        ambient_implicit_weight_interval=[min(relative),max(relative)],all_equation_weight_interval=[min(rowweights),max(rowweights)],
        complete_implicit_bound=10,complete_residual_bound=16)


def full_map_review(point,Ncols,U,rs,relative,rowweights,saved):
    R=Ring(('h','x','y','z','w','t0','t1','t2','t3'));T=Ring(('t0','t1','t2','t3'));nw=[1,2,2,3]
    coefficients=[{}for _ in range(150)]
    for k,col in enumerate(Ncols):
        e=tuple(int(i==k)for i in range(4))
        for j,c in enumerate(col):
            if c:coefficients[j]=plus(coefficients[j],{e:c})
    seen=set();assert Q(saved['common_delta'])==1
    for row in saved['complete_weight10_implicit_polynomials']:
        j,d=row['column'],row['weight'];assert j in U and j not in seen;seen.add(j)
        f=T.parse(row['numerator'].replace('^','**'));assert f and relative[j]==d<=10
        assert all(sum(e)>=2 and sum(a*b for a,b in zip(e,nw))==d for e in f)
        coefficients[j]=plus(coefficients[j],f)
    assert len(seen)==78
    P={ij:{e+(0,0,0,0):c for e,c in point[ij].items()}for ij in PAIRS}
    for poly,(ij,e)in zip(coefficients,COLS):P[ij]=plus(P[ij],{e+t:c for t,c in poly.items()})
    P=anti(P);J=jacobi(P,R)
    residual=[]
    for quad,e in ROWS:residual.append({powers[5:]:c for powers,c in J[quad].items()if powers[:5]==e})
    assert len(residual)==350 and all(not residual[i]for i in rs)
    expected={row['row']:T.parse(row['numerator'].replace('^','**'))for row in saved['complete_all350_remaining_equations_through_weight16']}
    assert len(expected)==179 and {i:p for i,p in enumerate(residual)if p}==expected
    assert all(sum(a*b for a,b in zip(e,nw))==rowweights[i]for i,p in enumerate(residual)for e in p)
    assert all(not p or 0<rowweights[i]<=16 for i,p in enumerate(residual))
    return residual,dict(all78_nonzero_implicit_coordinate_polynomials_verified=True,all98_implicit_coordinates_in_true_cut=True,
        implicit_correction_starts_at_ordinary_degree_two=True,all350_original_Jacobi_exactly_reconstructed=True,
        all98_selected_equations_identically_zero=True,all179_literal_nonzero_residuals_exactly_match=True,
        all_zero_and_negative_weight_equations_retained=True,all_residuals_exactly_weighted_homogeneous=True,
        formal_IFT_unique_complete_polynomial_map=True,
        proof='The exact98 minor and zero constant/linear implicit correction prove this polynomial map is the unique formal IFT solution of the selected98 original equations. Its equivariance and positive four-normal weights give finite all-orders bounds10/16. All350 original equations, including every zero row, were checked as polynomial identities; this is a true complete cut, not a finite jet.',
        full_map_term_count=sum(map(len,coefficients)),complete_residual_term_count=sum(map(len,residual)))


def ideal_review(residual,support,certificate):
    T=Ring(('t0','t1','t2','t3'));equations=[residual[row['row']]for row in support['complete_all350_remaining_equations_through_weight16']]
    assert len(equations)==certificate['original_complete_equation_count']==179
    outputs=[];targets=[]
    for packet in certificate['four_lifted_polynomials']:
        f=T.parse(packet['polynomial'].replace('^','**'));terms={int(k):T.parse(v.replace('^','**'))for k,v in packet['original_equation_coefficients'].items()}
        assert plus(*(times(g,equations[i])for i,g in terms.items()))==f
        degree=packet['weighted_degree'];assert all(sum(a*b for a,b in zip(e,[1,2,2,3]))==degree for e in f)
        lead=max(f,key=lambda e:(sum(e),tuple(-i for i in reversed(e))))
        assert {lead:Q(1)}==T.parse(packet['leading_monomial'].replace('^','**'))
        outputs.append(dict(normal_variable=packet['normal_variable'],weighted_degree=degree,
            original_179_equation_indices=sorted(terms),all_original_ideal_membership_coefficients_checked=True,
            leading_monomial=packet['leading_monomial']));targets.append(f)
    assert targets[1:]==[T.parse(t)for t in ['t1**4','t2**4','t3**3']]
    restriction={e:c for e,c in targets[0].items()if not any(e[1:])}
    coefficient=Q(514114689789605014659184048028741872800)
    assert restriction=={(5,0,0,0):coefficient}
    assert all(sum(e)>0 for f in residual for e in f)
    return dict(all_four_lifted_polynomials=outputs,all_four_explicit_ideal_identities_verified=True,
        pure_powers=['t1^4','t2^4','t3^3'],last_polynomial_restriction_coefficient=str(coefficient),last_polynomial_restriction='nonzero*t0^5',
        exact_reduced_true_cut_support='origin',complete_Groebner_basis_not_assumed=True,
        support_proof='Any common zero has t1=t2=t3=0 from three literal powers; the fourth polynomial then forces t0=0. All original equations vanish at the origin. Each polynomial is an explicitly verified polynomial combination of the complete original179 residuals.',
        local_uniqueness_or_smoothness_claimed=False)


def old_quartic_flags():
    R=Ring(('rho','h','x','y','z','w'));p=R.parse;h,x,y,z,w=[R.var(a)for a in ('h','x','y','z','w')]
    examples={
      'R4_H_QUARTIC_ESSENTIAL5':['h**2-x**2+x*z-x*w','-x**2','0','3*x*y-4*x*z+z**2-3*z*w','x**2+x*y-2*x*z+2*x*w+z**2/2-2*z*w+w**2','h**2+x**2+x*w'],
      'R4_MINIMAL_M':['h**2+3*x**2-6*x*z','0','-2*x**2','3*x*z-18*z**2','6*x**2-2*x*y-3*x*w+30*z*w','h**2+6*x*z'],
      'R4_A0':['h**2+3*x*z/4','0','z**2+3*x**2/4','-z**2/4','x**2+2*x*z+3*z**2-3*x*y/4+z*w/4','h**2+x*z/4']}
    # D+ from the printed common normal form at rho=sqrt(3).  Divergence
    # cancels rho identically, so no algebraic-field assumption is needed.
    c=p('rho/18');d=p('1/2+rho/18');a=p('-rho/6');b=p('-rho/216');u=p('1/6');v=p('-1/36');one=R.const(1)
    generic=[plus(times(plus(one,scale(d,-1)),x,z),scale(times(v,x,x),-1),scale(times(u,z,z),-1)),{},
        plus(times(a,z,z),scale(times(u,x,z),-4),scale(times(plus(c,scale(d,-3),R.const(2)),x,x),Q(1,2))),
        plus(times(b,x,x),scale(times(v,x,z),-4),scale(times(plus(c,d),z,z),-Q(1,2))),
        plus(p('x**2+2*x*z+3*z**2'),times(plus(scale(d,2),scale(c,-1),R.const(-1)),x,y),times(c,z,w),scale(times(u,y,z),2),scale(times(v,x,w),2)),
        plus(times(d,x,z),times(v,x,x),times(u,z,z))]
    generic[0]=plus(generic[0],times(h,h));generic[5]=plus(generic[5],times(h,h));examples['R4_DPM']=generic
    answer={}
    for label,texts in examples.items():
        six=[p(t)if isinstance(t,str)else t for t in texts];P=anti({ij:t for ij,t in zip(it.combinations(range(1,5),2),six)})
        div=[plus(*(derivative(P.get((i,j),{}),1+i)for i in range(5)))for j in range(5)]
        expected=[{}, {},p('6*x')if label=='R4_MINIMAL_M'else z,{},{}if label=='R4_MINIMAL_M'else x]
        assert div==expected
        N=[[field.get(tuple(int(k==j+1)for k in range(6)),Q(0))for j in range(5)]for field in div]
        assert all(not row[0]for row in N)
        rank=len(sparse_rref(N)[0]);assert rank==(1 if label=='R4_MINIMAL_M'else 2)
        columns_of_dual=N # Rows of N are columns of N^t.
        assert len(sparse_rref(columns_of_dual+[[Q(i==0)for i in range(5)]])[0])==rank+1
        answer[label]=dict(actual_directly_recomputed_divergence=['0','0','6*x'if label=='R4_MINIMAL_M'else'z','0','0'if label=='R4_MINIMAL_M'else'x'],
            full_homogeneous_N_matrix=[[str(a)for a in row]for row in N],rank=rank,h_not_in_dual_image=True,
            A0_or_DPM_type_not_assumed=True,arithmetic='DPM rho coefficient cancels identically; no CAS or sqrt3 approximation')
    return answer


def distinctness(prior,flags):
    rows=prior['components'];assert len(rows)==prior['component_count']==97 and len({r['id']for r in rows})==97
    # Bind the immutable accepted97 Git blob, never the mutable live catalogue
    # which will legitimately advance after a new component is accepted.
    frozen=subprocess.check_output(['git','cat-file','blob',prior['catalogue_git_blob']],cwd=ROOT)
    assert hashlib.sha256(frozen).hexdigest()==prior['catalogue_sha256']
    assert hashlib.sha1(('blob '+str(len(frozen))+'\0').encode()+frozen).hexdigest()==prior['catalogue_git_blob']
    catalogue=json.loads(frozen);original={r['id']:r for r in catalogue['components']}
    assert len(original)==97 and set(original)=={r['id']for r in rows}
    assert all(all(original[r['id']].get(k)==v for k,v in r.items())for r in rows)
    records=[];same=[]
    for row in rows:
        label=row['id'];dim=row['projective_dimension'];rank=row['generic_rank']
        if rank==2:reason='generic rank2 differs from the new generic rank4'
        elif dim!=23:reason='different full component dimension'
        elif row['meets_SNC_log_symplectic_locus']:reason='new whole closure has a plane of multiplicity>=3; SNC divisors do not'
        elif label in flags:
            assert row.get('generic_pgl5_orbit_dimension')==23
            same.append(label);reason='old dense23 PGL orbit has h not in im(N*); whole new c1 source has h in im(N*)'
        else:
            f=row['generic_actual_degeneracy_divisor'];assert f
            signature=(f.get('irreducible_factor_degrees'),f.get('corresponding_multiplicities'))
            assert signature!=([1,4],[1,1]),label
            reason='different generic factor degrees or multiplicities; proper factor multiplication images separate from irreducible h*F4'
        records.append(dict(id=label,reason=reason))
    assert set(same)==set(flags)
    return dict(all_prior97_rows_checked=True,all_projection_fields_exactly_match_immutable_committed_git_blob=True,prior_catalogue_sha256=prior['catalogue_sha256'],prior_catalogue_blob=prior['catalogue_git_blob'],
        frozen_snapshot_metadata_head=prior['catalogue_git_head'],all97_individual_reasons=records,same23_Hquartic_cases_directly_checked=flags,
        new_whole_source_cofactor='N=y*d_x+h*d_z',new_N_rank=2,new_N_squared_zero=True,new_h_in_dual_image=True,
        intrinsic_flag_proof='The unique linear divisor factor H defines h. Its h-Casimir homogeneous lift is unique: E wedge L also h-Casimir forces L=cE. Its divergence cofactor is therefore intrinsic up to conjugacy/nonzero scale. Every new source point has N*=sending dz to dh. The old dense23 orbits computed from actual tensors have no h in the dual image. Equality of component closures would make their dense constructible images intersect, contradicting this intrinsic flag.',
        triple_plane_closed_condition_proof='Quintics D in I_L^3 for some plane L form the proper image of the projective bundle over Gr(3,5). Every whole-source D lies in (h,y)^3. A reduced SNC divisor has at most two branches at a generic codimension2 point, so cannot have multiplicity>=3 along a plane. Include the zero quintic in the affine cone for rank-drop boundaries.',
        no_local_uniqueness_or_smoothness_claim=True,distinct_from_all97=True)


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',default='reports/CLASSIFICATION_GA9_CONSTANT_ANCHOR_REVIEW_20261107.json');args=ap.parse_args()
    resource.setrlimit(resource.RLIMIT_CPU,(239,240));resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,)*2);started=time.monotonic()
    paths=['reports/CLASSIFICATION_HYPERPLANE_QUARTIC_NIL211_SHORT1_CONSTANT_C_20261107.json','reports/CLASSIFICATION_GA9_CONSTANT_ANCHOR_SUPPORT_20261107.json','reports/CLASSIFICATION_GA9_CONSTANT_ANCHOR_CERTIFICATE_20261107.json','reports/CLASSIFICATION_GA9_CONSTANT_PRIOR97_20261107.json']
    source,support,certificate,prior=[json.loads((ROOT/p).read_text())for p in paths]
    point,directions,source_record=source_reconstruction(source['nil_short_translation_constant_c'])
    N,U,rs,relative,rw,frame=frame_review(point,directions,support)
    residual,map_record=full_map_review(point,N,U,rs,relative,rw,support)
    ideal=ideal_review(residual,support,certificate)
    differences=distinctness(prior,old_quartic_flags())
    imports=['code/classification_minimal_quartic_review.py','code/classification_quadratic_cofactor_review.py','code/classification_projected_veronese_zero_review.py',str(Path(__file__).relative_to(ROOT)),
        'legacy/v06/nilpotent22.tex','legacy/v07/minimal_thickness.tex','legacy/v07/nil22_v7.tex','legacy/v07/dplus_fullp4.tex']
    result=dict(task_id='E_reduced_nonnormal',all_passed=True,no_producer_or_CAS_import=True,whole_source=source_record,
        complete_frame=frame,complete_true_IFT_map=map_record,true_cut_support=ideal,all97_distinctness=differences,
        component_implication=dict(raw_local_dimension_upper=48,actual_whole_source_rank_lower=48,whole_irreducible_source_closure_raw_dimension=48,
            projective_component_dimension=23,whole_c1_source_and_all_closure_boundaries_covered=True,local_uniqueness_or_smoothness_claimed=False,
            proof='The actual whole ten-source/PGL/Euler/scale image has rank48 at the anchor and irreducible domain Gm(r)xA9 times irreducible groups. The genuine102 slice is defined by48 independent linear cuts and has support0 by explicit ideal identities. Krull bounds raw local dimension by48, matching the real image lower bound. Hence the whole irreducible c1 source image closure is an ambient component of dimension48, projective23. No assumption of reduced or unique local component, or of uniform b-field IFT, is needed.'),
        scope='The whole c1 Ga9 source produces one new23 ambient component relative to the audited prior97. Other c types and all Ga9 support are not exhausted.',
        input_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in paths+imports},seconds=time.monotonic()-started,max_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    dest=ROOT/args.output
    if dest.exists():raise FileExistsError(str(dest))
    dest.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:result[k]for k in ['all_passed','seconds','max_rss_kib']}))

if __name__=='__main__':main()
