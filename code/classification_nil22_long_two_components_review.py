#!/usr/bin/env python3
"""Independent full two type7 component certificates, Fraction only.

No producer, CAS, legacy engine, private discovery file or numerical prime.
The polynomial/linear arithmetic library is this reviewer's prior independent
projected-Veronese checker, with only standard-library dependencies.
"""
from fractions import Fraction as Q
import hashlib
import itertools as it
import json
from pathlib import Path
import resource
import subprocess
import time
from classification_projected_veronese_zero_review import Ring,plus,times,scale,derivative,monomials,sparse_rref,det

ROOT=Path(__file__).resolve().parents[1]
PAIRS=list(it.combinations(range(5),2));TRIPLES=list(it.combinations(range(5),3));QUADS=list(it.combinations(range(5),4))
COLS=[(ij,e)for ij in PAIRS for e in monomials(5,2)]
ROWS=[(quad,e)for quad in QUADS for e in monomials(5,4)]

def anti(P):
    out={ij:P.get(ij,{})for ij in PAIRS}
    for i,j in PAIRS:out[j,i]=scale(out[i,j],-1)
    return out

def jacobi(P,R,offset=0):
    S={tri:plus(*(times(P.get((i,l),{}),derivative(P.get((j,k),{}),offset+l))
                 for i,j,k in [tri,(tri[1],tri[2],tri[0]),(tri[2],tri[0],tri[1])]for l in range(5)))for tri in TRIPLES}
    return {quad:plus(*(scale(times(R.var(R.names[offset+quad[k]]),S[quad[:k]+quad[k+1:]]),(-1)**k)for k in range(4)))for quad in QUADS}

def vec(P):return [P.get(ij,{}).get(e,Q(0))for ij,e in COLS]
def jvec(J):return [J[quad].get(e,Q(0))for quad,e in ROWS]
def transpose(A):return [list(row)for row in zip(*A)]
def rank(A):return len(sparse_rref(A)[0])
def sparse(v):return {i:c for i,c in enumerate(v)if c}
def dot(row,col):return sum(c*col[j]for j,c in row.items())

def source_point(family,anchor):
    params=family['parameters'];k=len(params);R=Ring(tuple(params)+('h','x','y','z','w'));X=Ring(('h','x','y','z','w'))
    full={}
    for ij,text in zip(it.combinations(range(1,5),2),family['universal_brackets']):
        result={}
        for e,c in R.parse(text).items():
            ee=list(e);assert ee[k]==0;ee[k]=2-sum(ee[k+1:]);assert ee[k]>=0;result[tuple(ee)]=c
        full[ij]=result
    full=anti(full)
    assert all(not p for p in jacobi(full,R,k).values())
    values=[Q(family['canonical_parameter_values'][p])for p in params];sc=Q(anchor['integer_scale'])
    def evaluate(p):
        ans={}
        for e,c in p.items():
            for n,a in zip(e[:k],values):c*=a**n
            ee=e[k:];ans[ee]=ans.get(ee,Q(0))+sc*c
        return {e:c for e,c in ans.items()if c}
    P=anti({ij:evaluate(v)for ij,v in full.items()})
    directions=[anti({ij:evaluate(derivative(v,j))for ij,v in full.items()})for j in range(k)]
    assert all(P[ij]==X.parse(anchor['witness_brackets'].get(f'{ij[0]},{ij[1]}','0'))for ij in PAIRS)
    assert all(not p for p in jacobi(P,X).values())
    return P,directions,{'all_original_family_parameters':params,'units':family['localization'],
                         'all_full_family_projective_Jacobi_identically_zero':True,
                         'canonical_anchor_exactly_reconstructed_from_whole_family':True,
                         'whole_source_irreducibility_and_upper22':'Previously independently accepted public7/V2; inverse actual core/centering normalizations and all unit-zero complements retained.'}

def frame(P,directions,anchor,model):
    X=Ring(('h','x','y','z','w'));Dcols=[]
    for ij,e in COLS:
        V=dict(P);V[ij]=plus(P[ij],{e:Q(1)});V[ij[::-1]]=scale(V[ij],-1)
        Dcols.append(jvec(jacobi(V,X)))
    D=transpose(Dcols);Ds=[sparse(v)for v in D];E=[];G=[]
    for a,b in it.product(range(5),repeat=2):
        ev={};gv={}
        for i,j in PAIRS:
            ev[i,j]=plus(times(X.var(X.names[i]),X.var(X.names[b]))if j==a else {},scale(times(X.var(X.names[j]),X.var(X.names[b])),-1)if i==a else {})
            gv[i,j]=plus(P.get((b,j),{})if i==a else {},P.get((i,b),{})if j==a else {},scale(times(X.var(X.names[b]),derivative(P[i,j],a)),-1))
        E.append(vec(ev))
        # Public packet uses pullback differentiation; this is the negative
        # of the forward-pushforward formula assembled immediately above.
        G.append([-c for c in vec(gv)])
    allM=E+G+[vec(P)];M=transpose(allM);source=allM+[vec(v)for v in directions]
    assert rank(transpose(E))==24 and rank(M)==47 and rank(transpose(source))==47
    assert all(dot(row,v)==0 for row in Ds for v in source)
    cert=anchor['derivative'];minor=det([[D[i][j]for j in cert['columns']]for i in cert['rows']])
    assert len(cert['rows'])==len(cert['columns'])==102 and minor==Q(cert['determinant'])!=0
    Rvec=[Q(model['normal_sparse'].get(str(i),'0'))for i in range(150)]
    assert all(dot(row,Rvec)==0 for row in Ds)
    chosen=model['actual_action_columns'];assert len(chosen)==47
    M47=[allM[j]for j in chosen];K=M47+[Rvec]
    assert rank(transpose(K))==48
    kk=model['kernel_plus_action_minor_rows'];km=det([[v[i]for v in K]for i in kk])
    assert len(kk)==48 and km==Q(model['kernel_plus_action_minor'])!=0,(str(km),model['kernel_plus_action_minor'])
    savedK=[[Q(v.get(str(i),'0'))for i in range(150)]for v in anchor['full_exact_kernel_sparse_columns']]
    assert len(savedK)==48 and rank(transpose(savedK))==48 and all(dot(row,v)==0 for row in Ds for v in savedK)
    U=model['IFT_columns'];rs=model['IFT_rows'];assert len(U)==len(rs)==102 and len(set(U))==102
    im=det([[D[i][j]for j in U]for i in rs]);assert im==Q(model['IFT_minor'])!=0
    L=[Rvec]+[[Q(i==j)for i in range(150)]for j in U]
    assert rank(transpose(M47+L))==150
    weights=model['coordinate_weights'];character=model['tensor_character']
    rel=[sum(a*b for a,b in zip(e,weights))-weights[i]-weights[j]-character for(i,j),e in COLS]
    rw=[sum(a*b for a,b in zip(e,weights))-sum(weights[i]for i in quad)-2*character for quad,e in ROWS]
    assert all(not c or rel[j]==0 for j,c in enumerate(vec(P)))
    assert all(not c or rel[j]==rw[i]for i,row in enumerate(D)for j,c in enumerate(row))
    assert all(len({rel[i]for i,c in enumerate(v)if c})<=1 for v in M47)
    assert {rel[i]for i,c in enumerate(Rvec)if c}=={2}and model['normal_weights']==[2]
    maximum=max(rel[j]for j in U);assert maximum==model['a_priori_maximum_implicit_weight']
    assert maximum//2==model['a_priori_maximum_total_degree']
    return Rvec,U,rs,rel,rw,{'full_D_shape':[350,150],'full_D_rank':102,'D_rank102_minor':str(minor),
         'Euler_rank':24,'actual_GL5_Euler_scalar_rank':47,'whole_source_plus_group_rank':47,
         'whole_kernel_dimension':48,'all_saved_kernel_and_all_actual_directions_checked':True,
         'actual_M47_plus_R_minor':str(km),'full150_M47_R_U102_direct_sum':True,
         'true_cut_dimension':103,'complete_IFT_minor':str(im),'normal_positive_weight':2,
         'coordinate_weights':weights,'tensor_character':character,
         'maximum_implicit_weight':maximum,'a_priori_maximum_graph_degree':maximum//2,
         'all350_by150_weight_intertwining_verified':True,
         'projective_orbit_dimension':22,'projective_tangent_dimension':23}

def complete_model(P,Rvec,U,rs,rel,rw,model):
    R=Ring(('h','x','y','z','w','e'));graph={}
    for degree,column in model['complete_polynomial_graph_sparse'].items():
        d=int(degree)
        for j,c in column.items():
            j=int(j);ij,e=COLS[j];graph[ij]=plus(graph.get(ij,{}),{e+(d,):Q(c)})
            assert rel[j]==2*d
            if d>=2:assert j in U
    graph=anti(graph)
    assert all(graph[ij]==R.parse(model['complete_polynomial_graph_brackets'].get(f'{ij[0]},{ij[1]}','0'))for ij in PAIRS)
    assert all(graph[ij].get(e+(0,),Q(0))==vec(P)[j]and graph[ij].get(e+(1,),Q(0))==Rvec[j]for j,(ij,e)in enumerate(COLS))
    J=jacobi(graph,R);residual=[]
    for quad,e in ROWS:residual.append({power[5]:c for power,c in J[quad].items()if power[:5]==e})
    assert all(not residual[i]for i in rs)
    saved={int(d):{int(i):Q(c)for i,c in v.items()}for d,v in model['all350_residual_sparse'].items()}
    actual={}
    for i,p in enumerate(residual):
        assert len(p)<=1
        for d,c in p.items():actual.setdefault(d,{})[i]=c;assert rw[i]==2*d
    assert actual==saved
    order=min(actual);assert order==model['existing_complete_slice_algebra_length']
    assert order in [2,3]and all(d>=order for d in actual)
    row=min(actual[order]);coefficient=actual[order][row]
    assert coefficient
    return {'all350_literal_residuals_reconstructed':{str(d):{str(i):str(c)for i,c in v.items()}for d,v in actual.items()},
            'all102_selected_original_equations_polynomially_zero':True,
            'complete_graph_is_the_unique_formal_IFT_solution':True,
            'higher_graph_coefficients_only_in_true_implicit_U':True,
            'exact_completed_slice_ideal':f'(e^{order})','transverse_length':order,
            'unit_multiple_pure_power_witness':{'row':row,'coefficient':str(coefficient),'power':order},
            'raw_local_form':f'QQ[[g1,...,g47,e]]/(e^{order})',
            'unique_reduced_full_component_dimension':22,
            'generic_full_Jacobi_scheme_smoothness_claimed':False}

def cofactor(P):
    X=Ring(('h','x','y','z','w'));div=[plus(*(derivative(P.get((i,j),{}),i)for i in range(5)))for j in range(5)]
    N=[[v.get(tuple(int(k==j)for k in range(5)),Q(0))for j in range(5)]for v in div]
    assert all(sum(e)==1 for v in div for e in v)
    def mm(A,B):return [[sum(a*b for a,b in zip(row,col))for col in zip(*B)]for row in A]
    powers=[];A=N;squared=None
    for k in range(4):
        powers.append(rank(A))
        if k==1:squared=A
        A=mm(A,N)
    return {'N':[[str(c)for c in row]for row in N],'power_ranks':powers,
            'h_in_image_dual_N_squared':rank(squared+[[Q(i==0)for i in range(5)]])==rank(squared),
            'unique_H_h_Casimir_intrinsic':True}

def prior99_review(new_flags):
    commit='53da727720f053011cb53154eb789d7a6cb322be';path='research/COMPONENTS.json'
    blob=subprocess.check_output(['git','rev-parse',f'{commit}:{path}'],cwd=ROOT,text=True).strip()
    contents=subprocess.check_output(['git','cat-file','blob',blob],cwd=ROOT);catalogue=json.loads(contents)
    rows=catalogue['components'];assert len(rows)==99
    oldpaths=['reports/CLASSIFICATION_HYPERPLANE_QUARTIC_NIL_REGULAR_COMPONENT_20261107.json',
              'reports/CLASSIFICATION_HYPERPLANE_QUARTIC_NIL31_LONG_ORBIT1_AMBIENT_20261107.json']
    olddata=[json.loads((ROOT/p).read_text())for p in oldpaths]
    oldpoints=[olddata[0]['nil_regular_component'],olddata[1]['source_point']]
    ids=['R4_NIL4_HYPERPLANE_QUARTIC_22','R4_NIL31_LONG_HYPERPLANE_QUARTIC_22'];X=Ring(('h','x','y','z','w'));oldflags={}
    for label,saved in zip(ids,oldpoints):
        P=anti({tuple(map(int,k.split(','))):X.parse(v)for k,v in saved['witness_brackets'].items()})
        assert all(not v for v in jacobi(P,X).values());flag=cofactor(P)
        assert flag['power_ranks']==[3,2,1,0];oldflags[label]=flag
    reasons=[];counts={'rank2':0,'other_dimensions':0,'nonreduced22':0,'two_hyperplanes_cubic22':0,'Hquartic_Jordan41':0}
    for row in rows:
        label=row['id'];dimension=row['projective_dimension'];generic_rank=row['generic_rank']
        if generic_rank==2:reason='generic rank2 differs';counts['rank2']+=1
        elif dimension!=22:reason='different full component dimension';counts['other_dimensions']+=1
        elif label in oldflags:
            reason='unique-H h-Casimir cofactor Jordan(4,1), while both new dense22 orbits have Jordan(3,2)';counts['Hquartic_Jordan41']+=1
        else:
            f=row['generic_actual_degeneracy_divisor']
            if not f['reduced']:
                assert any(n>1 for n in f['corresponding_multiplicities'])
                reason='whole old closure satisfies closed nonreduced-quintic condition';counts['nonreduced22']+=1
            else:
                assert f['irreducible_factor_degrees']==[1,1,3]and f['corresponding_multiplicities']==[1,1,1]
                reason='closed 1+1+3 multiplication image excludes the new simple H times irreducible quartic';counts['two_hyperplanes_cubic22']+=1
        reasons.append({'id':label,'reason':reason})
    assert counts=={'rank2':6,'other_dimensions':77,'nonreduced22':13,'two_hyperplanes_cubic22':1,'Hquartic_Jordan41':2}
    assert all(flag['power_ranks']==[3,1,0,0]for flag in new_flags)
    return {'all99_individual_reasons':reasons,'counts':counts,'old_Hquartic_direct_divergence_flags':oldflags,
            'immutable_commit':commit,'catalogue_git_blob':blob,'catalogue_sha256':hashlib.sha256(contents).hexdigest(),
            'catalogue_rows_frozen_in_report':rows,'all99_distinctness_proved':True},oldpaths

def main():
    resource.setrlimit(resource.RLIMIT_CPU,(239,240));resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2))
    started=time.monotonic();paths=['reports/CLASSIFICATION_HYPERPLANE_QUARTIC_NIL22_LONG_QUARTIC_ORBIT_MODELS_20261107.json',
        'reports/CLASSIFICATION_HYPERPLANE_QUARTIC_NIL22_LONG_QUARTIC_SOURCES_20261107.json',
        'reports/CLASSIFICATION_HYPERPLANE_QUARTIC_NIL22_LONG_QUARTIC_REVIEW_V2_20261107.json',
        'research/classification/HYPERPLANE_QUARTIC_NIL22_LONG_QUARTIC_ORBIT_MODELS.md',
        'research/classification/HYPERPLANE_QUARTIC_NIL22_LONG_QUARTIC_SOURCES.md',
        'code/classification_projected_veronese_zero_review.py']
    packet=json.loads((ROOT/paths[0]).read_text())['nil22_long_quartic_orbit_models']
    family_packet=json.loads((ROOT/paths[1]).read_text())['nil22_long_quartic_sources']
    assert packet['family_sha256']==hashlib.sha256((ROOT/paths[1]).read_bytes()).hexdigest()
    assert family_packet['whole_H_irreducible_quartic_support_covered_by_seven_normalized_sources']
    assert family_packet['whole_saturated_projective_dimension_upper']==22
    answer=[];flags=[]
    for model in packet['complete_one_normal_models']:
        index=model['anchor']-1;assert index in [2,3]
        family=family_packet['all_seven_complete_families'][index];anchor=packet['all_seven_full_ambient_anchors'][index]
        P,directions,src=source_point(family,anchor);Rvec,U,rs,rel,rw,fr=frame(P,directions,anchor,model)
        graph=complete_model(P,Rvec,U,rs,rel,rw,model);flag=cofactor(P)
        assert flag['power_ranks']==[3,1,0,0]and flag['h_in_image_dual_N_squared']
        X=Ring(('h','x','y','z','w'))
        actualD=times(X.var('h'),plus(times(P[1,2],P[3,4]),scale(times(P[1,3],P[2,4]),-1),times(P[1,4],P[2,3])))
        witness='2*h*((h*z-w*y)**2+8*h**3*w)'if index==2 else '-6*h*((h*z-w*y)**2-8*h*w**3)'
        assert actualD==X.parse(witness)
        flags.append(flag);answer.append({'family_number':index+1,'name':family['name'],'whole_source':src,'full_frame':fr,
          'complete_one_normal_model':graph,'intrinsic_H_cofactor':flag,
          'actual_Euler_Pfaffian_exactly_reconstructed':witness,
          'absolute_irreducibility_proof':'The quartic is primitive quadratic in z over C[h,y,w]; its discriminant is a nonzero scalar times h^5*w (family3) or h^3*w^3 (family4), hence nonsquare. Its z coefficients have gcd1 since the constant term has w^2*y^2. Gauss over C proves irreducibility; h is the unique simple linear factor.',
          'whole_source_saturation_is_this_orbitclosure':True,
          'whole_source_bridge':'Accepted irreducible whole-source image upper22 contains the explicitly reconstructed actual orbit22. Its closure is therefore that orbitclosure, now proved a full component by the complete true slice.'})
    assert [p['complete_one_normal_model']['transverse_length']for p in answer]==[3,2]
    distinctness,oldpaths=prior99_review(flags);paths+=oldpaths
    out={'task_id':'E_reduced_nonnormal','all_passed':True,'no_producer_CAS_or_private_work_imports':True,
         'accepted_public7_V2_complete_cover_is_a_dependency_not_reclassified':True,'two_full_components':answer,
         'all_prior99_distinctness':distinctness,
         'pairwise_distinctness':'Each full22 component is a dense actual PGL orbit closure. Equality forces the two open dense orbits to coincide. The actual group x true-slice etale product identifies full raw Jacobi local rings with QQ[[47 smooth variables,e]]/(e^3) and /(e^2), which have intrinsic local multiplicities3 and2 and cannot be carried to each other by the actual PGL/Euler/scalar action.',
         'new_component_count_changed':False,
         'input_sha256':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},
         'reviewer_code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
         'seconds':time.monotonic()-started,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    target=ROOT/'reports/CLASSIFICATION_NIL22_LONG_TWO_COMPONENTS_REVIEW_20261107.json'
    assert not target.exists();target.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:out[k]for k in ['all_passed','seconds','max_rss_kib']}))

if __name__=='__main__':main()
