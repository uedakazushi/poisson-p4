#!/usr/bin/env python3
"""Independent physical tensor-factor Lie jets, actual57 constraints and genericV proof audit."""
import argparse
from fractions import Fraction as Q
import gzip
import hashlib
import itertools
import json
import math
from pathlib import Path
import resource
import signal
import subprocess
import sys
import tempfile
import time
import sympy as s
from sympy.polys.matrices import DomainMatrix

ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_CA_V1_ODD_ORDER_20261103_v63_REVIEW'
PAIRS=list(itertools.combinations(range(5),2));MONS=list(itertools.combinations_with_replacement(range(5),2))


def read(path):
    with (gzip.open(path,'rt')if path.suffix=='.gz'else path.open())as f:return json.load(f)


def sha(path):
    h=hashlib.sha256()
    with path.open('rb')as f:
        for block in iter(lambda:f.read(1<<20),b''):h.update(block)
    return h.hexdigest()


def matrix(rows):return s.Matrix([[s.Rational(v)for v in row]for row in rows])
def mm(A,B):return(DomainMatrix.from_Matrix(A).convert_to(s.QQ)*DomainMatrix.from_Matrix(B).convert_to(s.QQ)).to_Matrix()
def rank(M):return len(DomainMatrix.from_Matrix(M).convert_to(s.QQ).rref()[1])
def rational(v):return Q(int(s.numer(v)),int(s.denom(v)))


def tensor_factor_operator(a,b):
    # Two-index exterior action and the separate quadratic-monomial contragredient action.
    W=[];V=[]
    for out,(i,j)in enumerate(PAIRS):
        for old,(k,l)in enumerate(PAIRS):
            c=int(i==a)*(int(b==k and j==l)-int(b==l and j==k))
            c+=int(j==a)*(int(i==k and b==l)-int(i==l and b==k))
            if c:W.append((out,old,c))
    for old,mon in enumerate(MONS):
        terms={}
        for occurrence in range(2):
            if mon[occurrence]==a:
                word=list(mon);word[occurrence]=b;new=MONS.index(tuple(sorted(word)));terms[new]=terms.get(new,0)-1
        V.extend((old,new,c)for new,c in terms.items()if c)
    return W,V


def act(op,point):
    result=[Q(0)]*150;W,V=op
    for out,old,c in W:
        for mon in range(15):result[15*out+mon]+=c*point[15*old+mon]
    for old,new,c in V:
        for pair in range(10):result[15*pair+new]+=c*point[15*pair+old]
    return result


def add_scaled(target,vector,c):
    if c:
        for i,v in enumerate(vector):
            if v:target[i]+=c*v


def jet(point,fiber,operators,q,N):
    # Ordered representation exponentials on150 coefficients, independent of the5x5 coordinate-pushforward producer.
    current=[[Q(0)]*150 for _ in range(N+1)];current[0]=point[:]
    for j,F in enumerate(fiber):
        for degree in range(1,N+1):add_scaled(current[degree],F,q[20+j][degree])
    for j in range(19,-1,-1):
        if not any(q[j][:N+1]):continue
        result=[v[:]for v in current];term=[v[:]for v in current]
        for power in range(1,N+1):
            differentiated=[act(operators[j],v)for v in term]
            following=[[Q(0)]*150 for _ in range(N+1)]
            for degree in range(N+1):
                for k in range(1,degree+1):add_scaled(following[degree],differentiated[degree-k],q[j][k]/power)
            term=following
            for degree in range(N+1):add_scaled(result[degree],term[degree],Q(1))
        current=result
    return [[current[d][i]for d in range(N+1)]for i in range(150)]


def main(out):
    start=time.monotonic();out.mkdir(parents=True,exist_ok=False)
    resource.setrlimit(resource.RLIMIT_AS,(2200000000,2200000000));signal.alarm(60)
    v2=ROOT/'reports/AMBIENT_CA_V1_ODD_ORDER_V2_20261103_v63'
    supplemental=ROOT/'reports/AMBIENT_CA_V1_ACTUAL_ODD_RANK9_20261103_v63'
    checkpoint=ROOT/'reports/AMBIENT_CA_V1_ODD_ORDER_CHECKPOINT_20261103_v63/checkpoint.json.gz'
    hessian=ROOT/'reports/AMBIENT_CA_V1_COEFFICIENT_HESSIAN_20261103_v63/hessian.json.gz'
    critpath=ROOT/'reports/AMBIENT_B18_CA_CRITICAL_SCHEME_20261103_v63/certificate.json.gz'
    deps=[ROOT/'reports/AMBIENT_CA_GLOBAL_FLAG_MODEL_20261103_v63_REVIEW.json',
          ROOT/'reports/AMBIENT_B18_QUADRATIC_IMAGE_REDUCTION_20261103_v63_REVIEW/report.json',
          ROOT/'reports/AMBIENT_B18_CA_N_CRITICAL_SCHEMES_20261103_v63_REVIEW.json',supplemental/'report.json']
    checked={}
    for path in deps:
        d=read(path);assert d['all_passed']
        for name,want in d.get('input_sha256',{}).items():assert sha(ROOT/name)==want;checked[name]=want
        for name,want in d.get('evidence_sha256',{}).items():assert sha(path.parent/name)==want;checked[str((path.parent/name).relative_to(ROOT))]=want
    report=read(v2/'report.json');assert report['all_passed']and report['exact_normal_conductor_order']==1
    assert sha(checkpoint)==report['checkpoint_sha256']and sha(hessian)==report['hessian_sha256']
    assert sha(ROOT/'code/ambient_CA_V1_odd_order_v2_20261103_v63.py')==report['source_sha256']
    for name,want in report['evidence_sha256'].items():assert sha(v2/name)==want
    ck,helper,critical,setup,jets=[read(p)for p in [checkpoint,hessian,critpath,v2/'setup.json.gz',v2/'jets.json.gz']]
    C=matrix(critical['fixed_fiber_linear_basis150_by20']);assert C==matrix(ck['fixedCA150_by20'])
    X=s.symbols('h x y z w');h,x,y,z,w=X
    P0=[Q(0)]*150
    fields=[x*x,0,0,h*x,h*h+z*z]
    for k,(i,j)in enumerate(PAIRS):
        expr=(fields[j]if i==2 else 0)-(fields[i]if j==2 else 0)
        polynomial=s.Poly(expr,*X)
        for l,(a,b)in enumerate(MONS):P0[15*k+l]=rational(polynomial.coeff_monomial(X[a]*X[b]))
    assert C*s.Matrix(critical['center_fixed20'])==s.Matrix(P0)
    # The ENTIRE V curve, not a direction jet, lies in the actual critical source scheme.
    v=s.Symbol('v');sub={s.Symbol('u'):0,s.Symbol('w'):0,s.Symbol('v'):v}
    curve=[s.expand(s.sympify(critical['inverse_substitution'][f's{i}']).subs(sub))for i in range(10)]
    ss=s.symbols('s0:10')
    for expression in critical['full150_critical_Schur_column']:
        assert s.expand(s.sympify(expression).subs(dict(zip(ss,curve))))==0
    delta=s.zeros(20,1)
    for i,value in zip(critical['affine_slice_free_indices'],curve):delta[i]=value
    curvepoint=s.Matrix(P0)+C*delta;point=list(map(Q,ck['point150']))
    assert list(curvepoint.subs(v,1))==list(map(s.Rational,point))and curvepoint.subs(v,0)==s.Matrix(P0)
    allops=[tensor_factor_operator(*divmod(i,5))for i in range(25)]
    O=s.Matrix.hstack(*[s.Matrix(act(op,point))for op in allops]);assert rank(O)==20
    gi=ck['group20_GL_indices'];extra=ck['extra9_fiber_indices'];operators=[allops[i]for i in gi]
    F=[list(map(rational,C[:,i]))for i in extra]
    J=O[:,gi].row_join(C[:,extra]);assert J==matrix(ck['full150_by29_target_Jacobian'])and rank(J)==28
    assert rank(J[:,:20])==20
    assert J[:,10]==s.Matrix(point)and gi[10]==12 #noncentral E22, never centralGL scaling
    flags=critical['flag_positions'];whole=C.row_join(O[:,[5*i+j for i,j in flags]])
    piv=list(DomainMatrix.from_Matrix(C.T).rref()[1]);left=C[piv,:].inv();source=s.zeros(29,25)
    for k,(i,j)in enumerate(itertools.product(range(5),repeat=2)):
        if [i,j]in flags:source[20+flags.index([i,j]),k]=1
        else:source[:20,k]=mm(left,O[piv,k])
    assert mm(whole,source)==O
    frame=source[:,gi].row_join(s.eye(29)[:,extra]);assert frame==matrix(ck['normalizer_source29_frame'])and frame.det()!=0
    kappa=list(map(Q,ck['kernel29']));kap=s.Matrix(list(map(s.Rational,kappa)))
    assert J*kap==s.zeros(150,1)
    LP=[list(map(rational,O[:,i]))for i in gi]
    LF=[[act(op,F[j])for j in range(9)]for op in operators]
    LL={(i,j):act(operators[i],LP[j])for i in range(20)for j in range(i,20)}
    columns=[]
    for j in range(20):
        column=[Q(0)]*150
        for i in range(20):add_scaled(column,LL[min(i,j),max(i,j)],kappa[i])
        for f in range(9):add_scaled(column,LF[j][f],kappa[20+f])
        columns.append(column)
    for f in range(9):
        column=[Q(0)]*150
        for i in range(20):add_scaled(column,LF[i][f],kappa[i])
        columns.append(column)
    H=s.Matrix.hstack(*[s.Matrix(list(map(s.Rational,c)))for c in columns])
    assert H==matrix(helper['H_kappa_original150_by29'])
    jc=setup['source_r28_columns'];S=s.eye(29)[:,jc];T=kap.row_join(S);assert T.det()!=0
    Hn=mm(H,T);assert Hn==matrix(setup['full150_by29_mixed_Hessian'])
    Jr=J[:,jc];Y=setup['actual_target_Y28'];rg=setup['orbit20_raw_rows'];branch=setup['actual_branch29th_row'];selected=setup['actual_odd_selected8_rows']
    assert len(Y)==28 and Y[:20]==rg and len(selected)==8
    assert rank(J[rg,:20])==20
    assert selected==[27,28,29,35,60,65,67,68]
    inv=Jr[Y,:].inv();W=Hn[:,0]/2
    def rem(M):return M-mm(Jr,mm(inv,M[Y,:]))
    wc=rem(W);weight=wc[branch];assert weight!=0 and list(map(str,wc))==setup['corrected_quadratic150']
    Jn=s.zeros(150,1).row_join(Jr);odd=Hn.row_join(Jr);base=Y+[branch]
    block=Jn[rg,:].row_join(s.zeros(20,28)).col_join(odd[base+selected,:])
    assert block.shape==(57,57)and block==matrix(setup['full57_linear_block'])
    inverse=matrix(setup['inverse57']);assert mm(block,inverse)==mm(inverse,block)==s.eye(57)
    pure=[[Q(0)]*4 for _ in range(29)]
    for i,c in enumerate(kappa):pure[i][1]=c
    baseline=jet(point,F,operators,pure,3)
    assert all(row[0]==point[i]and row[1]==0 and row[2]==rational(W[i])for i,row in enumerate(baseline))
    Z=s.Matrix([s.Rational(row[3])for row in baseline])
    rhs=s.Matrix([-s.Rational(baseline[i][2])for i in rg]+[-s.Rational(baseline[i][3])for i in base+selected])
    solution=mm(inverse,rhs);de=mm(T,solution[:29,:]);do=mm(S,solution[29:,:])
    stage=jets['stages'][0];assert len(jets['stages'])==1 and stage['source_degree']==3
    assert list(map(str,de))==stage['source_even_coeff']and list(map(str,do))==stage['source_odd_coeff']
    finalsource=[[Q(v)for v in row[:4]]for row in jets['final_source29_series']]
    final=jet(point,F,operators,finalsource,3)
    assert final==[[Q(v)for v in row]for row in stage['all150_output_series']]
    assert all(final[i][2]==0 for i in rg)and all(final[i][3]==0 for i in base+selected)
    remaining={i:str(row[3])for i,row in enumerate(final)if row[3]};assert remaining=={63:'800/27',64:'-40/27'}
    assert mm(T.inv(),s.Matrix([s.Rational(row[1])for row in finalsource]))==s.eye(29)[:,0]
    assert mm(T.inv(),do)[0]==0
    # Actual trace-quadratic odd linearization, in raw target coordinates Y plusbranch.
    mixed=Hn[:,1:];ZZ=rem(Z-mm(mixed,mm(inv,W[Y,:])));EE=rem(mm(mixed,inv))
    def normal(M):return M-mm(wc/weight,M[branch,:])
    BY=normal(EE);Bnu=normal(ZZ)/weight;L=mm(Jr[branch,:],inv)
    A=(BY-mm(Bnu,L)).row_join(Bnu)
    ac=read(supplemental/'certificate.json.gz');assert A==matrix(ac['all150_actual_odd_linear'])
    ar,aq=ac['nine_rows'],ac['nine_columns'];minor=A[ar,aq]
    determinant=s.Rational(DomainMatrix.from_Matrix(minor).convert_to(s.QQ).det())
    assert determinant==s.Rational(ac['nine_minor'])==s.Rational('-70000000/59049')
    assert mm(A[:,aq],matrix(ac['all_column_factorization']))==A
    assert rank(rem(Hn))==9 and rem(Hn)==matrix(ac['whole_physical_critical_cotangent'])
    regular=s.eye(28)[:20,:].col_join(BY[selected,:]);regdet=s.Rational(DomainMatrix.from_Matrix(regular).convert_to(s.QQ).det())
    assert regdet==s.Rational(100,6561)==s.Rational(setup['regular_odd8_plus_orbit20_minor'])
    assert Z==matrix(ac['source_pure_kernel3'])
    # Fresh author57 run is only a replay; the independent physical computations are above.
    with tempfile.TemporaryDirectory(prefix='CA-V1-review-')as tmp:
        dest=Path(tmp)/'V2'
        subprocess.run([sys.executable,str(ROOT/'code/ambient_CA_V1_odd_order_v2_20261103_v63.py'),'--output',str(dest),'--checkpoint',str(checkpoint),'--hessian',str(hessian)],capture_output=True,text=True,check=True,timeout=30)
        for name in ['setup.json.gz','jets.json.gz']:assert(dest/name).read_bytes()==(v2/name).read_bytes()
    proper=read(ROOT/'reports/AMBIENT_B18_QUADRATIC_IMAGE_REDUCTION_20261103_v63/report.json')
    assert proper['whole_flag_fibers']['CA']['scheme_length']==2
    # Independent polynomial check of all81 proposed model relations; wholekernel uses the Koszul proof in the note.
    nu=s.Symbol('nu');b=[nu,*s.symbols('s1:9')];yy=s.symbols('Y0:9');tt=s.Symbol('t')
    relations=[b[j]*yy[i]-b[i]*yy[j]for i,j in itertools.combinations(range(9),2)]
    relations +=[yy[i]*yy[j]-nu*b[i]*b[j]for i in range(9)for j in range(i,9)]
    assert len(relations)==81
    for relation in relations:
        image=s.expand(relation.subs({y:bi*tt for y,bi in zip(yy,b)}))
        assert s.rem(image,tt*tt-nu,tt)==0
    (out/'own-jet.json').write_text(json.dumps(dict(all150_final_output_series=[[str(v)for v in row]for row in final],remaining_odd_cubic=remaining,critical_V_curve_s10=list(map(str,curve)),trace_quadratic_unit=str(weight),regular8_minor=str(regdet),actual9_minor=str(determinant)),indent=2)+'\n')
    companion=ROOT/'code/ambient_CA_generic_V_completion_20261103_v63.py';companionnote=ROOT/'research/notes/AMBIENT_CA_GENERIC_V_FULL_COMPLETION_20261103_v63.md'
    paths=[Path(__file__),ROOT/f'research/notes/{STEM}.md',*deps,checkpoint,hessian,critpath,v2/'report.json',v2/'setup.json.gz',v2/'jets.json.gz',supplemental/'certificate.json.gz',companion,companionnote,ROOT/'research/notes/AMBIENT_CA_V1_ODD_ORDER_RESULT_20261103_v63.md']
    result=dict(id=STEM,all_passed=True,current_dependency_hashes_checked=len(checked),
        independent_exterior10_and_quadratic15_tensor_factor_physical_Lie_action=True,
        helper_Hessian_treated_as_author_contribution_not_prior_independent_review=True,
        actual_V_curve_all150_critical_polynomials_zero_and_P0_V1_endpoints=True,
        full_source_frame29_invertible_GL20_target_rank20_and_scalar_column10_P=True,
        all150_by29_J_and_ordered_exponential_mixed_Hessian_rebuilt=True,
        actual57_block_and_both_inverse_products_verified=True,actual_regular8_with_orbit20_minor=str(regdet),
        all150_pure_kernel_and_corrected_finite_source_jets_independently_replayed=True,
        actual_odd_constraint_rows=selected,remaining_odd_cubic={'63':'800/27','64':'-40/27'},
        full_actual_odd_linear_rank=9,exact9minor=str(determinant),full_critical_cotangent_rank=9,
        all_author_setup_and_jets_fresh_byte_equal=True,raw_Y_slice_not_substituted_for_actual8odd=True,
        checkpoint_I_plus_E_firstjet_recipe_not_used_for_higher_diagonal_jets=True,
        fullproper_P0_fiber2_and_dense_V1_orbit_transport_to_fullfiber2_accepted=True,
        normalization_conclusion_after_nonzero_odd_I_not_assumed_from_proper_model=True,
        actual_trace_quadratic_image_R_plus_tI_and_critical_R_I_nu_identification_accepted=True,
        smooth_codim9_J_and_actual_I_eq_J_by_Nakayama_accepted=True,
        actual_CA_generic_V_fullcompletion='C[[g20,s8,t^2,s1*t,...,s8*t,t^3]]',
        complete81_kernel_Koszul_products_proof_and_conductor_depth_multiplicity_accepted=True,
        cone_dimension=29,cone_embedding_dimension=38,cone_depth=21,multiplicity=2,
        projective_dimension=28,projective_embedding_dimension=37,projective_depth=20,
        fulltype_scope='CA V1 and its dense V orbit; notCB,notB18origin,notactualCA_CBintersection',
        no_GB_or_symbolic_multivariate_Hessian_or_new_highorder_probe=True,
        input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},evidence_sha256={'own-jet.json':sha(out/'own-jet.json')},seconds=time.monotonic()-start)
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');signal.alarm(0)
    print('INDEPENDENT CA V1 ACTUAL57/ODD9/FULL GENERIC V COMPLETION PASS',result['seconds'],flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
