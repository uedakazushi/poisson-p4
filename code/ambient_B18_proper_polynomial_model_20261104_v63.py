#!/usr/bin/env python3
"""ONE actual proper-source14/17 polynomial input; no local GB/standard basis."""
import argparse,hashlib,itertools,json,math,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
import sympy as s
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_PROPER_POLYNOMIAL_MODEL_20261104_v63'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def encode_poly(p,vars):
    if p==0:return []
    return [[list(m),int(s.numer(c)),int(s.denom(c))] for m,c in s.Poly(p,*vars).terms() if c]
def matrix_strings(M):return [[str(c) for c in row] for row in M.tolist()]
def main(out):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False)
    pp=ROOT/'reports/AMBIENT_B18_PROPER_POLYNOMIAL_MODEL_PREFLIGHT_20261104_v63/literal-proper-source14-preflight.json';pr=pp.parent/'report.json'
    P=read(pp);assert read(pr)['all_passed']
    AC=s.Matrix(P['all150_CA_fixed_flag19']);BC=s.Matrix(P['all150_CB_fixed_flag19']);point=s.Matrix(P['actual_P0_physical150'])
    N=s.Matrix(P['CB_annihilator17_by36']);Bi=s.Matrix(P['CB_left_inverse19']);rows=P['zero_coordinate_rows36'];br=P['CB_constant_left_inverse_rows19']
    free=P['CA_free_fiber9_coefficient_indices'];CAs=s.Matrix.hstack(point,AC[:,free]);assert CAs.shape==(150,10)
    Gs=s.Matrix(P['source_gauge15_flag5_then_fiber10_Jacobian']);cut=s.Matrix(P['actual_common_cut15_response'])
    assert Gs.det()==-s.Rational(2,3) and cut.det()==32
    assert Bi*BC[[rows[j] for j in br],:]==s.eye(19) and N*BC[rows,:]==s.zeros(17,19)
    h,x,y,z,w=s.symbols('h x y z w');phys=(h,x,y,z,w);a,b,c,d,e=s.symbols('a b c d e');flags=(a,b,c,d,e)
    G=s.Matrix(P['relative_CB_frame5_by5']);Gi=s.Matrix(P['relative_CB_inverse_frame5_by5']);assert G.det()==1 and G*Gi==s.eye(5)
    pairs=list(itertools.combinations(range(5),2));mons=list(itertools.combinations_with_replacement(range(5),2));pi={p:i for i,p in enumerate(pairs)}
    T=s.zeros(150,10)
    for col in range(10):
        Y=s.zeros(5,1)
        for j in [0,1,3,4]:
            pair=tuple(sorted((2,j)));sgn=1 if j>2 else-1
            Y[j]=sum(sgn*CAs[15*pi[pair]+k,col]*phys[u]*phys[v] for k,(u,v) in enumerate(mons))
        moved=(G*Y.subs(dict(zip(phys,Gi*s.Matrix(phys))),simultaneous=True)).applyfunc(s.expand)
        assert moved[2]==0 and s.expand(sum(s.diff(moved[j],phys[j]) for j in [0,1,3,4]))==0
        for j in [0,1,3,4]:
            poly=s.Poly(moved[j],h,x,z,w);assert poly.is_zero or all(sum(m)==2 for m in poly.monoms())
            pair=tuple(sorted((2,j)));sgn=1 if j>2 else-1
            for k,(u,v) in enumerate(mons):
                if 2 in (u,v):continue
                T[15*pi[pair]+k,col]=s.expand(sgn*poly.coeff_monomial(phys[u]*phys[v]))
    T36=T[rows,:];R19=(Bi*T36[br,:]).applyfunc(s.expand);F17=(N*T36).applyfunc(s.expand)
    rest=[i for i in range(36) if i not in br];inj=s.zeros(36,17)
    for k,r in enumerate(rest):inj[r,k]=1
    assert (T36-BC[rows,:]*R19-inj*F17).applyfunc(s.expand)==s.zeros(36,10)
    ep=ROOT/'reports/AMBIENT_B18_ZERO_NORMALIZER_TENSOR_MODEL_20261104_v63/actual-normalizer-tensor-free4-DAG.json';E=read(ep)
    rec=s.zeros(150,36)
    for i,r in enumerate(rows):rec[r,i]=1
    for j,r in enumerate(E['eliminated_zero_rows4']):rec[r,:]=s.Matrix([E['eliminated_four_coordinate_reconstruction'][j]])
    assert (T-rec*T36).applyfunc(s.expand)==s.zeros(150,10)
    assert (T-BC*R19-rec*inj*F17).applyfunc(s.expand)==s.zeros(150,10)
    v=s.symbols('v0:9');eq=[s.expand(F17[i,0]+sum(F17[i,j+1]*v[j] for j in range(9))) for i in range(17)]
    zero=dict.fromkeys(flags,0);J=s.zeros(17,14)
    for i in range(17):
        for j in range(9):J[i,j]=F17[i,j+1].subs(zero)
        for j,flag in enumerate(flags):J[i,9+j]=s.diff(F17[i,0],flag).subs(zero)
    assert J==s.Matrix(P['polynomial_Jacobian17_by14_at_origin']) and J.rank()==5
    jr=P['Jacobian_pivot_rows'];jc=P['Jacobian_pivot_columns'];assert J[jr,jc].det()==s.Rational(8,3)
    variables=tuple(v)+flags;pivvars=[variables[j] for j in jc]
    degrees=[s.Poly(eq[i],*pivvars).total_degree() for i in jr]
    # Explicit six-generator Npair FLAG ideal, directly from the raw source
    # identity. B,C are independent free fiber coefficients, A remains a unit.
    ca=s.Matrix(P['actual_CA_source19_center'])
    src=ca.copy()
    for j,k in enumerate(free):src[k]+=v[j]
    AA,BB,CC=src[:3,0];assert 1 in free and 2 in free and AA.subs(dict.fromkeys(v,0))==-1
    I=[CC,BB-AA*a,b,c,d,e]
    substitute={v[free.index(2)]:-ca[2],v[free.index(1)]:AA*a-ca[1],b:0,c:0,d:0,e:0}
    assert all(s.expand(q.subs(substitute,simultaneous=True))==0 for q in eq)
    Ilin=s.Matrix([[s.diff(q,t).subs(dict.fromkeys(variables,0)) for t in variables] for q in I]);assert Ilin.rank()==6
    assert Ilin.col_join(J[jr,:]).rank()==6
    choices=[]
    for j in [10,11,12,13]:
        if j in jc:continue
        line=s.zeros(1,14);line[0,j]=1
        if J[jr,:].col_join(line).rank()==6:choices.append(j)
    assert choices
    flagcolumn=choices[0];explicit_flag=str(variables[flagcolumn]);assert explicit_flag in ['b','c','d','e']
    conormal=J[jr,:].col_join(s.eye(14)[flagcolumn,:]);rc=list(conormal.rref()[1]);iminor=s.factor(conormal[:,rc].det());assert iminor
    def em(M):return [[encode_poly(q,flags) for q in row] for row in M.tolist()]
    fe,te,re=em(F17),em(T36),em(R19)
    allpolys=[p for mat in [fe,te,re] for row in mat for p in row]
    terms=sum(len(p) for p in allpolys);assert terms<=150000
    fterms=sum(len(p) for row in fe for p in row);assert fterms<=42840
    maxdeg=max((sum(m) for p in allpolys for m,n,dn in p),default=0);assert maxdeg<=5
    allco=[s.Rational(n,dn) for p in allpolys for m,n,dn in p]
    fco=[s.Rational(n,dn) for row in fe for p in row for m,n,dn in p]
    den=s.ilcm(*[s.denom(c) for c in fco]);assert 15%den==0
    allDen=s.ilcm(*[s.denom(c) for c in allco])
    maxnum=max(abs(int(c*15)) for c in fco);assert maxnum<=int(P['future_budget']['conservative_integer_numerator_abs_bound'])
    source={'flag_variables':list(map(str,flags)),'CA_fiber_variables':list(map(str,v)),'coefficient_monomial_format':'[five flag exponents, numerator, denominator]',
      'all17_by10_polynomial_coefficients':fe,'all36_by10_transformed_CA_coefficients':te,
      'all19_by10_CB_coefficient_readbacks':re,'source_fixed_P0_file':str(pp.relative_to(ROOT)),
      'all36_image_plus_annihilator_reconstruction_identities':True,'all150_zero_divergence_reconstruction_and_CB_image_identities':True,
      'complete_equation_rule':'Fi=coefficient[i][0]+sum_j coefficient[i][j+1]*v_j',
      'equations14variable_sparse':[encode_poly(q,variables) for q in eq],
      'source_gauge15_minor':'-2/3','same_actual_common_cut15_minor':'32','Jacobian17_by14':matrix_strings(J),
      'five_unit_pivot_rows':jr,'five_unit_pivot_columns':jc,'five_unit_minor':'8/3',
      'pivot_variable_names':list(map(str,pivvars)),
      'pivot_distribution':{'CA_fiber':sum(j<9 for j in jc),'relative_flag':sum(j>=9 for j in jc)},
      'selected_equations_degrees_in_all5_pivot_variables':degrees,
      'selected5_jointly_linear_in_pivots':all(n<=1 for n in degrees),
      'Npair_source_six_generators':list(map(str,I)),'Npair_substitution_all17_polynomials_zero':True,
      'Npair_conormal6_by14':matrix_strings(Ilin),'explicit_free_flag_generator_after_5_unit_IFT':explicit_flag,
      'Npair_Nakayama6_minor_columns':rc,'Npair_Nakayama6_minor':str(iminor),
      'exact_formal_Npair_identity':'I_Npair=(five selected polynomial equations, explicit_free_flag)',
      'all12_remaining_after5IFT_have_actual_common_flag_factor':True,
      'remaining12_rational_polynomials_or_unit_solve_constructed':False,
      'model_is_normalizer_tensor_E4_not_actual_D0':True,'no_GB_standardbasis_newjets_or_CM_depth_assertion':True,
      'statistics':{'F_coefficient_terms':fterms,'all_stored_coefficient_matrix_terms':terms,'max_flag_degree':maxdeg,
        'max_total_equation_degree':max(s.Poly(q,*variables).total_degree() for q in eq),'coefficient_common_denominator':str(den),'all_readback_common_denominator':str(allDen),'max_numerator_when_scaled_to_den15':maxnum,'coefficient_numerator_bit_bound':int(maxnum).bit_length()}}
    p=out/'actual-proper-source14-equations17.json';p.write_text(json.dumps(source,separators=(',',':'))+'\n');assert p.stat().st_size<=32*1024*1024
    paths=[Path(__file__),ROOT/'research/notes'/f'{STEM}.md',pp,pr,ep]
    rep={'id':STEM,'all_passed':True,'status':'actual_complete_14_variable17_polynomial_E4_source_model_built',
      'actual_polynomial_equations':17,'variables':14,'linear_CA_fiber_variables':9,'parameter_coefficient_matrix':[17,10],
      'all36_and150_source_readback_image_identities_checked':True,'source_gauge15_minor':'-2/3','physical_cut15_minor':'32','linear5_minor':'8/3',
      'pivot_distribution':source['pivot_distribution'],'pivot_variables':source['pivot_variable_names'],
      'pivot_equations_jointly_linear':source['selected5_jointly_linear_in_pivots'],
      'explicit_Npair_flag_factor_after5IFT':explicit_flag,'Npair_ideal_and_flag_factor_proved_by_source_substitution_and_Nakayama':True,
      'rational12_equations_GB_localstandardbasis_higherjets_D0_or_CM_claimed':False,
      **source['statistics'],'output_bytes':p.stat().st_size,
      'input_sha256':{str(p.relative_to(ROOT)):sha(p) for p in paths},'evidence_sha256':{p.name:sha(p)},
      'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    (out/'report.json').write_text(json.dumps(rep,indent=2)+'\n');print(json.dumps(rep))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
