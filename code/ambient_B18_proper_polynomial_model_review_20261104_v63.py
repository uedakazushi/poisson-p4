#!/usr/bin/env python3
"""Own physical CA/CB formula and exact proper-source14/17 model audit.
No producer imports, source jet, Groebner basis or local elimination.
"""
import gzip,hashlib,itertools,json,resource,signal,time
from pathlib import Path
import sympy as S
ROOT=Path(__file__).resolve().parents[1]
ID='AMBIENT_B18_PROPER_POLYNOMIAL_MODEL_20261104_v63_REVIEW'
def load(p):return json.loads(gzip.decompress(p.read_bytes())if p.suffix=='.gz'else p.read_bytes())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def mat(a):return S.Matrix([[S.Rational(c)for c in row]for row in a])
def em(a):return [[str(c)for c in row]for row in a.tolist()]
def check_parent(stem,paths):
    p=ROOT/'reports'/stem;r=load(p/'report.json');assert r['all_passed'];paths.add(p/'report.json')
    for rel,h in r.get('input_sha256',{}).items():assert sha(ROOT/rel)==h,rel;paths.add(ROOT/rel)
    for rel,h in r.get('evidence_sha256',{}).items():assert sha(p/rel)==h,rel;paths.add(p/rel)
    return p,r
def main(out):
    started=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));out.mkdir(exist_ok=False,parents=True);paths=set()
    pb,pr=check_parent('AMBIENT_B18_PROPER_POLYNOMIAL_MODEL_PREFLIGHT_20261104_v63',paths)
    mb,mr=check_parent('AMBIENT_B18_PROPER_POLYNOMIAL_MODEL_20261104_v63',paths)
    p=load(pb/'literal-proper-source14-preflight.json');q=load(mb/'actual-proper-source14-equations17.json')
    h,x,y,z,w=S.symbols('h x y z w');physical=(h,x,y,z,w);four=(0,1,3,4)
    pairs=list(itertools.combinations(range(5),2));mons=list(itertools.combinations_with_replacement(range(5),2));pairindex={v:i for i,v in enumerate(pairs)}
    def vecY(Y):
        v=S.zeros(150,1)
        for j in four:
            pol=S.Poly(S.expand(Y[j]),*physical);ij=tuple(sorted((2,j)));sgn=1 if j>2 else-1
            for k,(u,vv)in enumerate(mons):v[15*pairindex[ij]+k]=sgn*pol.coeff_monomial(physical[u]*physical[vv])
        return v
    # Independently derive canonical ey-wedge vector fields AFTER the physical h/x swap.
    # CA divergence before canonicalization is -f_h; CB is -2*f_h.
    # Euler subtraction therefore adds respectively f_h*X/5 and 2*f_h*X/5.
    def sourceY(f,H,cb):
        fh=S.diff(f,h)
        if cb:return S.Matrix([-f+2*h*fh/5,2*x*fh/5,0,-S.diff(H,w)-z*fh/10,S.diff(H,z)-w*fh/10])
        return S.Matrix([-f+h*fh/5,x*fh/5,0,-S.diff(H,w)+z*fh/5,S.diff(H,z)+w*fh/5])
    cubic=[S.prod((x,h,z,w)[j]for j in mon)for mon in itertools.combinations_with_replacement(range(4),3)if any(j>=2 for j in mon)]
    binary=[x*x,x*h,h*h];assert len(cubic)==16
    AY=[sourceY(f,0,False)for f in binary]+[sourceY(0,H,False)for H in cubic]
    BY=[sourceY(f,0,True)for f in binary]+[sourceY(0,H,True)for H in cubic]
    AC=S.Matrix.hstack(*map(vecY,AY));BC=S.Matrix.hstack(*map(vecY,BY))
    assert AC==mat(p['all150_CA_fixed_flag19']) and BC==mat(p['all150_CB_fixed_flag19'])
    Y0=S.Matrix([x*x,0,0,h*x,h*h+z*z]);P0=S.zeros(5)
    for j in four:P0[2,j]=Y0[j];P0[j,2]=-Y0[j]
    point=vecY(Y0);assert point==S.Matrix([S.Rational(c)for c in p['actual_P0_physical150']])
    op=ROOT/'reports/AMBIENT_B18_ACTUAL_INTERSECTION_ORBIT_PRODUCT_20261103_v63/certificate.json.gz';O=load(op);paths.add(op)
    assert point==S.Matrix([S.Rational(c)for c in O['physical_point150']])
    rows=p['zero_coordinate_rows36'];ar=p['CA_constant_left_inverse_rows19'];br=p['CB_constant_left_inverse_rows19']
    A36=AC[rows,:];B36=BC[rows,:];ai=mat(p['CA_left_inverse19']);bi=mat(p['CB_left_inverse19'])
    assert ai*A36[ar,:]==bi*B36[br,:]==S.eye(19)
    ap=S.Matrix([S.Rational(c)for c in p['actual_CA_source19_center']]);bp=S.Matrix([S.Rational(c)for c in p['actual_CB_source19_center']])
    assert AC*ap==BC*bp==point
    rest=[i for i in range(36)if i not in br];N=S.zeros(17,36);inj=S.zeros(36,17)
    for k,r in enumerate(rest):
        N[k,r]=1;inj[r,k]=1;response=B36[r,:]*bi
        for j,c in enumerate(br):N[k,c]-=response[j]
    assert N==mat(p['CB_annihilator17_by36']) and N*B36==S.zeros(17,19) and N[:,rest]==S.eye(17)
    # Standard physical GL4 action, independently applied to literal P0.
    def vector(P):
        ans=[]
        for i,j in pairs:
            pol=S.Poly(S.expand(P[i,j]),*physical)
            ans.extend(pol.coeff_monomial(physical[u]*physical[v])for u,v in mons)
        return S.Matrix(ans)
    def action(i,j):
        M=S.zeros(5);M[i,j]=1;v=M*S.Matrix(physical)
        return vector(M*P0+P0*M.T-sum((v[k]*P0.diff(physical[k])for k in range(5)),S.zeros(5)))
    gl4=[(i,j)for i in four for j in four];orbits={ij:action(*ij)for ij in gl4}
    orbit=S.Matrix.hstack(*[orbits[ij]for ij in gl4]);assert orbit==mat(p['all150_GL4_orbit16']) and orbit.rank()==15
    assert sum((orbits[i,i]for i in four),S.zeros(150,1))==-point
    par=[tuple(ij)for ij in p['parabolic_positions11']];bad=[(1,0),(1,3),(1,4),(0,3),(0,4)]
    assert len(par)==11 and sorted(par+bad)==sorted(gl4)
    for i,j in par:assert not(i==1 and j!=1)and not(i==0 and j not in(0,1))
    para=S.Matrix.hstack(*[ai*orbits[ij][[rows[k]for k in ar],0]for ij in par])
    assert para==mat(p['CA_parabolic_coefficient_matrix19_by11']) and AC*para==S.Matrix.hstack(*[orbits[ij]for ij in par]) and para.rank()==10
    chosen=p['chosen_parabolic10_column_indices'];fixed=p['CA_fiber10_fixed_coefficient_indices'];free=p['CA_free_fiber9_coefficient_indices']
    assert sorted(fixed+free)==list(range(19)) and len(fixed)==10 and len(free)==9
    pg=para[:,chosen];assert pg[fixed,:].det()==S.Rational(p['parabolic_fiber10_minor'])
    # Under I+tM, flag covectors have derivative -bad(M); the restoring frame
    # has the opposite bad action. Therefore the actual source gauge is diag(-I5,M10).
    tt=S.Symbol('group_t')
    for k,(i,j)in enumerate(bad):
        elementary=S.zeros(5);elementary[i,j]=1;gg=S.eye(5)+tt*elementary;ggi=S.eye(5)-tt*elementary
        assert gg*ggi==S.eye(5)
        ell=S.eye(5)[1,:]*ggi;ell=ell/ell[1]
        kk=S.eye(5)[0,:]*ggi;kk=kk-kk[1]*ell;kk=kk/kk[0]
        actual_flags=[ell[0],ell[3],ell[4],kk[3],kk[4]]
        flagder=S.Matrix([S.diff(f,tt).subs(tt,0)for f in actual_flags]);assert flagder==-S.eye(5)[:,k]
        correcting=S.eye(5);correcting[1,0],correcting[1,3],correcting[1,4],correcting[0,3],correcting[0,4]=actual_flags
        assert ((correcting*gg).diff(tt)).subs(tt,0)==S.zeros(5)
    gauge=S.diag(-S.eye(5),pg[fixed,:]);assert gauge==mat(p['source_gauge15_flag5_then_fiber10_Jacobian']) and gauge.det()==-S.Rational(2,3)
    g15=bad+[par[i]for i in chosen];assert [list(ij)for ij in g15]==p['actual_G15_generator_positions']
    cutrows=p['actual_common_cut15_rows'];orbit15=S.Matrix.hstack(*[orbits[ij]for ij in g15]);cut=orbit15[cutrows,:]
    assert cut==mat(p['actual_common_cut15_response']) and cut.det()==32
    ep=ROOT/'reports/AMBIENT_B18_ZERO_NORMALIZER_TENSOR_MODEL_20261104_v63/actual-normalizer-tensor-free4-DAG.json';E=load(ep);paths.add(ep)
    assert rows==E['coordinate_rows36'] and cutrows==E['common_cut_zero_rows15']
    # Own exact coordinate frame and its inverse. No source-jet construction.
    a,b,c,d,e=S.symbols('a b c d e');flags=(a,b,c,d,e);G=S.eye(5)
    G[0,3]=d;G[0,4]=e;G[1,0]=a;G[1,3]=b;G[1,4]=c
    Gi=S.eye(5);Gi[0,3]=-d;Gi[0,4]=-e;Gi[1,0]=-a;Gi[1,3]=a*d-b;Gi[1,4]=a*e-c
    assert G*Gi==S.eye(5) and G.det()==1
    # The saved frame entries are polynomials, not rational scalar matrices.
    assert G==S.Matrix(p['relative_CB_frame5_by5']) and Gi==S.Matrix(p['relative_CB_inverse_frame5_by5'])
    Ys=[Y0]+[AY[i]for i in free];T=S.zeros(150,10)
    for j,Y in enumerate(Ys):
        moved=(G*Y.subs(dict(zip(physical,Gi*S.Matrix(physical))),simultaneous=True)).applyfunc(S.expand)
        assert S.expand(sum(S.diff(moved[i],physical[i])for i in four))==0
        T[:,j]=vecY(moved)
    F17=(N*T[rows,:]).applyfunc(S.expand);R19=(bi*T[[rows[i]for i in br],:]).applyfunc(S.expand)
    assert (T[rows,:]-B36*R19-inj*F17).applyfunc(S.expand)==S.zeros(36,10)
    rec=S.zeros(150,36)
    for i,r in enumerate(rows):rec[r,i]=1
    for i,r in enumerate(E['eliminated_zero_rows4']):rec[r,:]=S.Matrix([E['eliminated_four_coordinate_reconstruction'][i]])
    assert (T-rec*T[rows,:]).applyfunc(S.expand)==S.zeros(150,10)
    assert (T-BC*R19-rec*inj*F17).applyfunc(S.expand)==S.zeros(150,10)
    def decode(p,variables):return S.expand(sum(S.Rational(n,den)*S.prod(v**k for v,k in zip(variables,m))for m,n,den in p))
    for saved,actual in [(q['all17_by10_polynomial_coefficients'],F17),(q['all36_by10_transformed_CA_coefficients'],T[rows,:]),(q['all19_by10_CB_coefficient_readbacks'],R19)]:
        assert len(saved)==actual.rows
        for i in range(actual.rows):
            for j in range(10):assert S.expand(decode(saved[i][j],flags)-actual[i,j])==0
    v=S.symbols('v0:9');variables=tuple(v)+flags
    eq=[S.expand(F17[i,0]+sum(F17[i,j+1]*v[j]for j in range(9)))for i in range(17)]
    for i in range(17):assert decode(q['equations14variable_sparse'][i],variables)==eq[i]
    origin=dict.fromkeys(variables,0);J=S.Matrix([[S.diff(f,variable).subs(origin)for variable in variables]for f in eq])
    assert J==mat(q['Jacobian17_by14'])==mat(p['polynomial_Jacobian17_by14_at_origin']) and J.rank()==5
    jr=q['five_unit_pivot_rows'];jc=q['five_unit_pivot_columns'];assert J[jr,jc].det()==S.Rational(8,3)
    assert [str(variables[i])for i in jc]==q['pivot_variable_names']
    # Full coordinate Npair ideal and its exact17 source substitutions.
    ca=ap.copy()
    for i,j in enumerate(free):ca[j]+=v[i]
    A,B,C=ca[:3,0];I=[C,B-A*a,b,c,d,e]
    assert list(map(str,I))==q['Npair_source_six_generators']==['v1','a + v0','b','c','d','e']
    subst={v[1]:0,v[0]:-a,b:0,c:0,d:0,e:0}
    assert all(S.expand(f.subs(subst,simultaneous=True))==0 for f in eq)
    il=S.Matrix([[S.diff(f,variable).subs(origin)for variable in variables]for f in I]);assert il==mat(q['Npair_conormal6_by14']) and il.rank()==6
    assert il.col_join(J[jr,:]).rank()==6
    assert q['explicit_free_flag_generator_after_5_unit_IFT']=='d'
    con=J[jr,:].col_join(S.eye(14)[variables.index(d),:]);cols=q['Npair_Nakayama6_minor_columns'];minor=con[:,cols].det()
    assert minor==S.Rational(q['Npair_Nakayama6_minor'])!=0
    for gate in [ROOT/'reports/AMBIENT_B18_QUADRATIC_IMAGE_REDUCTION_20261103_v63_REVIEW/report.json',
      ROOT/'reports/AMBIENT_B18_ACTUAL_INTERSECTION_ORBIT_PRODUCT_20261103_v63_REVIEW.json',
      ROOT/'reports/AMBIENT_B18_ZERO_NORMALIZER_TENSOR_MODEL_20261104_v63_REVIEW/report.json']:
        assert load(gate)['all_passed'];paths.add(gate)
    proof={'all_passed':True,'own_source_formulas':'Physical ey-wedge canonical CA/CB vectors, using divergence -f_h/-2f_h and Euler correction; no producer imports.',
      'own_all150_CA_CB_19_source_coefficients':True,'own_P0_and_both_source_centers':True,'own_CB19_leftinverse_complete_annihilator17':True,
      'own_full_GL4_orbit_rank':15,'own_parabolic_rank':10,'own_source_gauge15_minor':str(gauge.det()),'own_same_physical_cut15_minor':str(cut.det()),
      'GL4center_affine_scalar_included_no_projectivization':True,'own_all360_transformed36_and_all150_image_readback_checks':True,
      'own_all170_parameter_polynomials_and17_full_equations':True,'own_Jacobian_rank':5,'own_Jacobian_minor':str(J[jr,jc].det()),
      'literal_Npair_ideal':list(map(str,I)),'own_all17_Npair_substitutions_zero':True,'own_Npair_conormal6_minor':str(minor),'explicit_remaining_flag':'d',
      'allorder_argument':'The two smooth proper zero sources have dimension24 each at their accepted P0 normalization points. The15-parameter GL4 section and CA source gauge have unit derivative diag(-I5,M10), giving a formal product of the full48-source germ with its33-dimensional gauge slice. The equality incidence ideal is GL4-invariant, so pulls back from that slice. The same15 group parameters have physical-cut derivative determinant32, hence the actual physical zero cut is formally isomorphic to the source-gauge incidence. The four positive cuts vanish on this fixed locus. CB19 coefficients are eliminated by a constant left inverse with exact kernel annihilator17, yielding the full14-variable17-equation ideal. This is E4, not D0.',
      'flag_factor_argument':'Every Fi lies in the coordinate ideal I=(v1,a+v0,b,c,d,e), as proved by exact substitution. The five selected Fi together with d span its six independent conormal classes, with the saved nonzero6 minor. Nakayama gives I=(selected5 Fi,d) in the completed local ring. After the5-unit IFT, all12 remaining Fi are multiples of the actual residual flag d. No rational elimination or CM conclusion is inferred.',
      'E4_not_D0':True,'new_GB_sourcejets_local_elimination_CM_depth_claims':False}
    (out/'independent-proof.json').write_text(json.dumps(proof,indent=2)+'\n')
    paths.update([Path(__file__),ROOT/'research/notes'/(ID+'.md')])
    report={'id':ID,'all_passed':True,'source_prefight_and_exact14_17_model_independently_checked':True,'E4_not_D0':True,
      'input_sha256':{str(z.relative_to(ROOT)):sha(z)for z in sorted(paths)},'evidence_sha256':{'independent-proof.json':sha(out/'independent-proof.json')},
      'seconds':time.monotonic()-started,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'frozen':True}
    assert sum(f.stat().st_size for f in out.iterdir())<32*1024*1024;signal.alarm(0)
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2),flush=True)
if __name__=='__main__':
    import argparse
    pa=argparse.ArgumentParser(description=__doc__);pa.add_argument('--output',type=Path,default=ROOT/'reports'/ID)
    main(pa.parse_args().output.resolve())
