#!/usr/bin/env python3
"""Independent complete21-source G review, using only prior Fraction libraries.

Mutable producer/jobs are read from archived bytes only. No producer/CAS/GB
imports. Full350/30/9, all walls, cokernels, actual shear slices and all700
residual equations are reconstructed, including a monic gcd certificate.
"""
import ast
from fractions import Fraction as F
import hashlib,json,resource,signal,time
from pathlib import Path
import classification_plane_rank_zero_rank_two_review as r
import classification_plane_global_nilpotent_nonzero_b_review as b1review
import classification_plane_high_m_rank_at_least_two_review as h
from classification_unimodular_logarithmic_review import ud

ROOT=Path(__file__).resolve().parents[1]
PFX='reports/CLASSIFICATION_PLANE_RANK_ZERO_GLOBAL_'
JOB='research/classification/jobs/plane_global_lowparameter_closure_review.json'
OUT='reports/CLASSIFICATION_PLANE_GLOBAL_LOWPARAMETER_CLOSURE_REVIEW_20261107.json'
def read(suffix):return json.loads((ROOT/(PFX+suffix+'_20261107.json')).read_text())
def poly(s,names):return {e:h.constant(c)for e,c in r.parse(s,tuple(names),len(names)).items()}
def scalar(s):
    def visit(node):
        if isinstance(node,ast.Constant)and type(node.value)is int:return r.RF(node.value)
        if isinstance(node,ast.Name):assert node.id=='tau';return r.RHO
        if isinstance(node,ast.UnaryOp):return -visit(node.operand)if isinstance(node.op,ast.USub)else visit(node.operand)
        if isinstance(node,ast.BinOp):
            a,b=visit(node.left),visit(node.right)
            if isinstance(node.op,ast.Add):return a+b
            if isinstance(node.op,ast.Sub):return a-b
            if isinstance(node.op,ast.Mult):return a*b
            if isinstance(node.op,ast.Div):return a/b
            if isinstance(node.op,ast.Pow):
                k=h.constant(b);assert k.denominator==1;k=int(k)
                if k<0:a=r.RF(1)/a;k=-k
                out=r.RF(1)
                for _ in range(k):out=out*a
                return out
        raise ValueError(ast.dump(node))
    return visit(ast.parse(str(s),mode='eval').body)
def matrix(M):return [[scalar(v)for v in row]for row in M]
def certificate(M,C):
    rank=C['rank'];n=len(M[0]);k=n-rank;rr,cc=C['minor_rows'],C['minor_columns']
    assert C['nullity']==k and len(rr)==len(cc)==rank
    minor=r.det([[M[i][j]for j in cc]for i in rr]);assert minor and minor==scalar(C['minor_determinant'])
    K=matrix(C['complete_kernel']);assert len(K)==n and all(len(row)==k for row in K)
    assert all(not sum((a*b for a,b in zip(row,col)if a and b),r.RF())for row in M for col in zip(*K))
    assert not k or r.rank(K)==k
    return dict(shape=[len(M),n],rank=rank,nullity=k,full_kernel_checked=True,nonzero_minor=r.encode(minor))
def mat(M):return [[h.constant(r.scalar(v))for v in row]for row in M]
def mulmat(A,B):return [[sum(a*b for a,b in zip(row,col))for col in zip(*B)]for row in A]
def eye(n):return [[F(i==j)for j in range(n)]for i in range(n)]
def value(poly_,n):return {e+(0,)*n:h.constant(c)for e,c in poly_.items()}
def pad(poly_,n):return {e+(0,)*n:c for e,c in poly_.items()}


def model(name):
    if name.startswith('NIL_B1'):return b1review.models(name)[0]
    a,b,c,z,w=[r.atom(1,i)for i in range(5)];t=r.RHO
    H=[a,b,r.q.scale(c,-2),{},{}];N=[b,{},{},{},{}];J=[b,c,{},{},{}];J2=[c,{},{},{},{}]
    commonT=[c,{},{},{},{}];commonV=[{},c,{},{},{}]
    def plus(A,B):return [r.q.add(a,b)for a,b in zip(A,B)]
    if name.startswith('SS_'):
        if name=='SS_NILSQUARE':T,V,K=J,J2,r.q.scale(r.q.mul(z,w),t)
        else:T,V,K=commonT,commonV,r.q.scale(r.q.mul(z,w),int(name.endswith('K1')))
        normal=[{},{},{},z,r.q.scale(w,-1)]
    else:
        normal=[{},{},{},{},z]
        if name=='NIL_B0_NILREG':T,V,k=plus(J,J2),J,t
        elif 'TSQUARE'in name:T,V,k=J2,J,int(name.endswith('K1'))
        elif 'VSQUARE'in name:T,V,k=J,J2,int(name.endswith('K1'))
        elif 'COMMON'in name:T,V,k=commonT,commonV,int(name.endswith('K1'))
        elif name=='NIL_B0_MIX_TN':T,V,k=N,plus(H,N),t
        elif name=='NIL_B0_MIX_VH':T,V,k=plus(H,N),H,t
        elif name=='NIL_B0_MIX_TN_VH':T,V,k=N,H,t
        elif name=='NIL_B0_MIX_VN':T,V,k=plus(H,N),N,t
        elif name=='NIL_B0_MIX_TH':T,V,k=H,plus(H,N),t
        elif name.startswith('NIL_B0_MIX_TH_VN'):T,V,k=H,N,int(name.endswith('K1'))
        else:raise ValueError(name)
        K=r.q.scale(r.q.mul(z,z),k)
    return r.padd(r.wedge(T,[{},{},{},z,w]),r.wedge(V,normal),r.q.antisymmetric({(3,4):K}))


def operators(saved,basis1,basis2):
    records=[];cache={}
    for c in saved['cases']:
        name=c['type'];P=model(name)
        assert P==r.q.antisymmetric({ast.literal_eval(ij):r.parse(v)for ij,v in c['normal_form_P0'].items()})
        assert r.divisor(P)==r.parse(c['actual_D'])
        A,B,S=r.operator(P,basis1),r.operator(P,basis2),r.shear(P,basis1)
        assert A==r.matrix(c['complete_D1_350_by30'])and B==r.matrix(c['complete_D2_350_by9'])and S==r.matrix(c['complete_actual_shear_30_by6'])
        a=c['analysis'];checks={k:certificate(M,a[k])for k,M in [('D1',A),('D2',B),('actual_shear',S)]};walls=[]
        if c['varying_parameter']:
            product=r.RF(1)
            for k in checks:product=product*scalar(a[k]['minor_determinant'])
            assert product==scalar(a['exception_polynomial'])
            factored=r.RF(1);roots=set()
            for f,m in a['all_irreducible_factors']:
                f=scalar(f);assert len(f.p)==2 and len(f.q)==1;roots.add(-f.p[0]/f.p[1])
                for _ in range(m):factored=factored*f
            ratio=product/factored;assert ratio and len(ratio.p)==len(ratio.q)==1
            assert roots=={F(w['parameter'])for w in a['all_exceptional_roots']}
            for wall in a['all_exceptional_roots']:
                t=F(wall['parameter']);certs={k:certificate(r.specialization(M,t),wall[k])for k,M in [('D1',A),('D2',B),('actual_shear',S)]}
                walls.append(dict(parameter=str(t),certificates=certs,m2_excluded=wall.get('excluded_multiplicity_higher',False)))
        records.append(dict(type=name,generic=checks,walls=walls,affine_symbol_dimension=c['affine_first_symbol_dimension_upper']))
        cache[name]=(P,A,B,S,c)
    return records,cache


def tensor(P0,basis1,c1,basis2,c2,n):
    def spatial(P):return {ij:{e[:5]+(0,)*n:h.constant(c)for e,c in P.get(ij,{}).items()}for ij in r.PAIRS}
    out=spatial(P0)
    for B,coefficient in list(zip(basis1,c1))+list(zip(basis2,c2)):
        B=spatial(B);coeff={(0,)*5+e:c for e,c in coefficient.items()}
        for ij in r.PAIRS:out[ij]=r.q.add(out[ij],r.q.mul(B[ij],coeff))
    return r.q.antisymmetric(out)
def coeffs(M,n):return [r.q.add(*(r.q.scale(h.variable(n,j),c)for j,c in enumerate(row)))for row in M]
def polar(basis1,C,n):
    P=tensor({},basis1,coeffs(C,n),[],[],n)
    return h.jacobi_and_divisor(P,n)[0]
def cokernel(phi,B):return all(not sum(F(c)*h.constant(B[int(i)][j])for i,c in phi)for j in range(9))
def combine(phi,Q):return r.q.add(*(r.q.scale(Q[int(i)],F(c))for i,c in phi))


def onecuts(saved,cache,basis1):
    records=[]
    for c in saved['cases']:
        name=c['type'];P,A,B,S,original=cache[name];wall=next(w for w in original['analysis']['all_exceptional_roots']if w['parameter']=='0')
        assert r.matrix(c['full_weight1_kernel'])==r.matrix(wall['D1']['complete_kernel'])
        assert r.matrix(c['full_weight2_kernel'])==r.matrix(wall['D2']['complete_kernel'])
        Q=polar(basis1,mat(c['full_weight1_kernel']),8);names=[f'u{i}'for i in range(8)]
        assert Q==[poly(s,names)for s in c['full_Q2_350']]
        phi=c['cokernel_functional'];B0=r.specialization(B,F(0))
        assert cokernel(phi,B0)and combine(phi,Q)==poly(c['obstruction_polynomial'],names)and combine(phi,Q)
        rows,cols=c['weight2_pivot_rows'],c['weight2_pivot_columns'];minor=[[h.constant(B0[i][j])for j in cols]for i in rows]
        assert mulmat(minor,mat(c['weight2_pivot_inverse']))==eye(6)
        assert c['PGL_image_upper']==9+7+3+5==24
        records.append(dict(type=name,full350_Q2_and_fullD2_cokernel=True,nonzero_obstruction=True,whole_image_upper=24))
    return records


def full_comp(c,cache,basis1,basis2):
    P,A,B,S,original=cache[c['type']];R=mat(c['complete_transverse_weight1_30_by_n']);n=len(R[0]);K=mat(original['analysis']['D2']['complete_kernel'])
    assert r.matrix(c['full_first_weight_kernel'])==r.matrix(original['analysis']['D1']['complete_kernel'])
    assert len(R)==30 and n in (6,8)
    columns=c['actual_shear_selected_columns'];L=[[h.constant(S[i][j])for j in columns]+R[i]for i in range(30)]
    assert len(columns)==c['actual_shear_rank']
    minor=[L[i]for i in c['shear_plus_slice_basis_rows']]
    assert minor==mat(c['shear_plus_slice_basis_minor'])and mulmat(minor,mat(c['shear_plus_slice_basis_inverse']))==eye(len(minor))
    assert all(not sum(a*b for a,b in zip(row,col))for row in A for col in zip(*R))
    assert len(minor)==original['analysis']['D1']['nullity']
    Q=polar(basis1,R,n);names=[f'u{i}'for i in range(n)]
    assert Q==[poly(v,names)for v in c['full_Q2_350']]
    rows,cols=c['D2_pivot_rows'],c['D2_pivot_columns'];I=mat(c['D2_pivot_inverse']);M=[[h.constant(B[i][j])for j in cols]for i in rows]
    assert mulmat(M,I)==eye(6)
    particular=[{}for _ in range(9)]
    for j,row in zip(cols,I):particular[j]=r.q.add(*(r.q.scale(Q[i],-v)for i,v in zip(rows,row)))
    residual=[r.q.add(q,r.q.add(*(r.q.scale(p,h.constant(v))for p,v in zip(particular,row))))for q,row in zip(Q,B)]
    gens=[poly(s,names)for s in c['necessary_compatibility_generators']]
    assert len(gens)==len(c['obstruction_cokernel_functionals'])
    for phi,g in zip(c['obstruction_cokernel_functionals'],gens):assert cokernel(phi,B)and combine(phi,Q)==g
    mons=sorted(set().union(*(set(p)for p in residual+gens)))
    rm=[[p.get(e,F())for e in mons]for p in residual];gm=[[p.get(e,F())for e in mons]for p in gens]
    assert r.rank(rm)==r.rank(gm)==r.rank(rm+gm)
    return dict(P=P,R=R,n=n,names=names,particular=particular,K=K,gens=gens)


def gcd_monic(Fpoly,Gpoly,names):
    """A constant leading coefficient makes this exact specialization conclusive.

    Every nonconstant common divisor is monic in u0 up to a scalar; thus its
    positive u0 degree survives every specialization of the other variables.
    A coprime specialized pair excludes all such factors over C as well as Q.
    """
    x=names.index('u0');degree=max(e[x]for e in Fpoly)
    leading={e[:x]+e[x+1:]:c for e,c in Fpoly.items()if e[x]==degree}
    assert set(leading)=={(0,)*(len(names)-1)}and leading[(0,)*(len(names)-1)]
    values={name:F(name in ('u1','u2','u6'))for name in names if name!='u0'}
    def specialize(P):
        ans=[F()]* (max((e[x]for e in P),default=0)+1)
        for e,c in P.items():
            for i,k in enumerate(e):
                if i!=x:c*=values[names[i]]**k
            ans[e[x]]+=c
        while ans and not ans[-1]:ans.pop()
        return tuple(ans)
    a,b=specialize(Fpoly),specialize(Gpoly);original=(a,b);remainders=[]
    while b:
        quotient,remainder=ud(a,b);remainders.append(dict(quotient=list(map(str,quotient)),remainder=list(map(str,remainder))));a,b=b,remainder
    assert len(a)==1 and a[0]
    return dict(variable='u0',leading_coefficient=str(next(iter(leading.values()))),
                exact_specialization={k:str(v)for k,v in values.items()},
                univariate_pair=[list(map(str,p))for p in original],Euclidean_chain=remainders,gcd_one_over_C=True)


def branch(case,data,basis1,basis2):
    n0=data['n'];names=data['names']+['v0','v1','v2'];n=len(names)
    images=[h.variable(n,j)for j in range(n)]
    for key,value_ in case['substitution'].items():images[names.index(key)]=poly(value_,names)
    c1=[h.compose(p,images)for p in coeffs(data['R'],n)]
    c2=[h.compose(r.q.add(pad(p,3),*(r.q.scale(h.variable(n,n0+j),k)for j,k in enumerate(row))),images)
        for p,row in zip(data['particular'],data['K'])]
    for g in data['gens']:assert not h.compose(pad(g,3),images)
    P1=tensor({},basis1,c1,[],[],n);P2=tensor({},[],[],basis2,c2,n)
    for saved,P in [('complete_weight1_coefficients',P1),('complete_weight2_coefficients',P2)]:
        assert P==r.q.antisymmetric({ast.literal_eval(ij):poly(v,list(r.NAMES)+names)for ij,v in case[saved].items()})
    P=tensor(data['P'],basis1,c1,basis2,c2,n);J,D=h.jacobi_and_divisor(P,n)
    J3=[];J4=[]
    for row,((quad,ex),j)in enumerate(zip(r.q.ROWS,J)):
        weight=ex[3]+ex[4]-sum(i>=3 for i in quad)
        if weight not in (3,4):assert not j
        J3.append(j if weight==3 else{});J4.append(j if weight==4 else{})
    assert J3==[poly(v,names)for v in case['all350_weight3_Jacobi']]
    assert J4==[poly(v,names)for v in case['all350_weight4_Jacobi']]
    assert J3[65]
    if case['type']=='SS_COMMON_K0':certificate=gcd_monic(J3[65],J3[66],names)
    else:certificate=dict(nonzero_original_weight3_row=65,polynomial=case['all350_weight3_Jacobi'][65])
    return dict(type=case['type'],linear_branch=case.get('linear_branch'),all700_rows_and_fullP2_kernel=True,certificate=certificate)


def main():
    resource.setrlimit(resource.RLIMIT_CPU,(235,236));resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,)*2);signal.alarm(240)
    start=time.monotonic();job=json.loads((ROOT/JOB).read_text())
    for row in job['exact_inputs_and_hashes']:assert hashlib.sha256((ROOT/row['path']).read_bytes()).hexdigest()==row['sha256'],row['path']
    manifest=json.loads((ROOT/'research/classification/inputs/plane_rank_zero_global_vector_lowparameter_packet_manifest_20261107.json').read_text())
    archive=json.loads((ROOT/'research/classification/inputs/plane_rank_zero_global_lowparameter_frozen_inputs_20261107.json').read_text())
    archived={p['path']:p for p in archive['snapshots']};assert len(archived)==9
    for row in manifest['files']:
        content=archived[row['path']]['bytes_utf8'].encode()if row['path']in archived else(ROOT/row['path']).read_bytes()
        assert hashlib.sha256(content).hexdigest()==row['sha256']
        if row['path']in archived:assert archived[row['path']]['sha256']==row['sha256']
    low=read('LOWPARAMETER_WEIGHTS');assert low['all_passed']is False
    basis1,basis2,audit=r.source_bases(low);ops,cache=operators(low,basis1,basis2)
    cuts=onecuts(read('ONECUT'),cache,basis1)
    compatibility=read('COMPATIBILITY');datas={c['type']:full_comp(c,cache,basis1,basis2)for c in compatibility['cases']}
    records=[]
    for c in read('RESIDUAL')['cases']:
        if 'linear_branch'in c:records.append(branch(c,datas[c['type']],basis1,basis2))
    assert len(records)==7
    ts=read('TSQUARE');data=datas['NIL_B0_TSQUARE_K0'];names=data['names'];p=lambda s:poly(s,names)
    g=data['gens'];L=p('-u2+3*u3+u4');d=p('u2-u4')
    assert g[1]==p('-12*u5**2')and g[3]==r.q.scale(r.q.mul(p('u3'),L),F(1,3))
    assert g[5]==r.q.add(p('14*u0*u5/3'),r.q.scale(r.q.mul(L,r.q.add(r.q.scale(d,2),p('9*u3'))),F(1,6)))
    assert r.q.mul(L,L)==r.q.add(r.q.scale(r.q.mul(p('u3'),L),F(15,2)),r.q.scale(r.q.mul(L,r.q.add(r.q.scale(d,2),p('9*u3'))),-F(1,2)))
    ts_case=dict(ts,substitution=ts['linear_substitution'],type='NIL_B0_TSQUARE_K0')
    records.append(branch(ts_case,data,basis1,basis2))
    assert poly(ts['nonzero_witness']['polynomial'],names+['v0','v1','v2'])==poly('u3*(2*u0**2-v1)',names+['v0','v1','v2'])
    nonlinear_fixed={'SS_COMMON_K0','NIL_B0_TSQUARE_K0','NIL_B0_COMMON_K0'}
    for row in ops:
        if row['type']in nonlinear_fixed:continue
        B=row['affine_symbol_dimension'];C=row['generic']
        assert B+C['D1']['nullity']+C['D2']['nullity']+5<=24
        for wall in row['walls']:
            if wall['m2_excluded']:continue
            C=wall['certificates'];bound=B-1+C['D1']['nullity']+C['D2']['nullity']+5
            if bound>24:assert row['type']in {'SS_NILSQUARE','NIL_B0_NILREG'}and wall['parameter']=='0'and bound==25
    out=dict(task_id='G',all_passed=True,producer_or_CAS_imports=False,
             exact_input_bindings=job['exact_inputs_and_hashes'],archived_mutable_files_checked=9,
             complete_source_quotient=audit,whole21_operators_and_all_walls=ops,
             both_onecuts=cuts,all_common_and_Tsquare_branches=records,
             reduced_Tsquare_support_equal_u5_L_zero=True,all21_whole_images_at_most24=True,
             old_partial_false_flags_preserved=True,GB_or_lift_dependency=False,
             complements=['generic two/three parameter sources','high-m restored symbols not separately accepted','fullG/E exhaustion'],
             elapsed_seconds=time.monotonic()-start,peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    (ROOT/OUT).write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:out[k]for k in ['all_passed','elapsed_seconds','peak_rss_kib']}))


if __name__=='__main__':main()
