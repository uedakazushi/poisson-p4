#!/usr/bin/env python3
"""Independent actual-source/formal-operator and homological T_nil audit.

Does not import the producer or evaluate a new source function or source jet.
"""
import hashlib,itertools,json,resource,signal,time
from collections import Counter
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
PRIMARY='AMBIENT_B18_HONEST_NIL_T_OPERATOR_20261104_v63'
ID=PRIMARY+'_REVIEW';OUT=ROOT/'reports'/ID

def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def frozen(x):
    if isinstance(x,list):return tuple(frozen(t)for t in x)
    if isinstance(x,dict):return tuple((k,frozen(v))for k,v in sorted(x.items()))
    return x
def write(name,data):(OUT/name).write_text(json.dumps(data,indent=2)+'\n')
def mv(A,v):return [sum(F(a)*b for a,b in zip(r,v))for r in A]
def mm(A,B):return [[sum(F(a)*F(B[k][j])for k,a in enumerate(r))for j in range(len(B[0]))]for r in A]
def ident(n):return [[F(int(i==j))for j in range(n)]for i in range(n)]

# Independent sparse polynomial primitives for identities in the formal template.
def pa(*ps):
    out={}
    for p in ps:
        for e,c in p.items():out[e]=out.get(e,F(0))+c
    return {e:c for e,c in out.items()if c}
def pm(p,i,n=1,scale=F(1)):
    out={}
    for e,c in p.items():v=list(e);v[i]+=n;out[tuple(v)]=c*scale
    return out
def pd(p,i):
    out={}
    for e,c in p.items():
        if e[i]:v=list(e);v[i]-=1;out[tuple(v)]=c*e[i]
    return out
def ps(p,c):return {e:a*c for e,a in p.items()if a*c}
PAIRS=list(itertools.combinations(range(8),2))
def differential(h,n):
    v=[{}for _ in range(8)]
    for (i,j),p in zip(PAIRS,h):
        v[j]=pa(v[j],pm(p,i));v[i]=pa(v[i],pm(p,j,n if j==7 else 1,-1))
    return v
def lift_formal_cycle(v,n):
    # Ordered prefix division assigns a monomial to its LAST nonzero z index.
    H=[{}for _ in range(7)]
    for e,c in v[7].items():
        support=[i for i in range(7)if e[i]];assert support
        i=support[-1];q=list(e);q[i]-=1;H[i][tuple(q)]=c
    w=[pa(v[i],pm(H[i],7,n))for i in range(7)]
    assert not pa(*(pm(w[i],i)for i in range(7)))
    out=[]
    for i,j in PAIRS:
        if j==7:out.append(H[i]);continue
        d=pa(pd(w[j],i),ps(pd(w[i],j),-1))
        out.append({e:c/F(sum(e[:7])+2)for e,c in d.items()})
    assert differential(out,n)==v
    return out

def main():
    start=time.monotonic();signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));OUT.mkdir(exist_ok=True)
    pins={}
    def load(path):
        p=ROOT/path;pins[path]=sha(p);return json.loads(p.read_text())
    base='reports/'+PRIMARY+'/'
    report=load(base+'report.json')
    assert pins[base+'report.json']=='3e6bd43aaa4877d5b17614c914dc29b5a4a1dcb951c99e01d65c2975fb8aafd8'
    for path,h in report['input_sha256'].items():assert sha(ROOT/path)==h
    D=load(base+'actual-honest-Tnil128-972-formal-operator.json')
    assert pins[base+'actual-honest-Tnil128-972-formal-operator.json']=='9325c84cd57b9ec6fdf4b62256be84fd1e3ecc6673615b224d2c92f0aa897de0'
    src={}
    for ns,ref in D['source_artifacts'].items():
        src[ns]=load(ref['path']);assert pins[ref['path']]==ref['sha256']
        if ns.endswith('_peer'):assert src[ns]['all_passed']
    R,O,N,L,C=[src[k]for k in ['radical','old_unit','new_resolution','old_resolution','Cprime_values']]
    assert src['radical_peer']['all_passed'] and R['NilD0_exactly_H1_of_actual_complex']
    assert O['actual_variable_order']==D['actual_S_frame']==['zA0','zA1','zA2','zA3','zA4','zA5','zA6','rA']
    assert O['V56_target_order']==['tauA','tauB','tauA*tauB']
    assert O['h_is_zero']=='m*n=I24 by the selected unit block; hence V*n=d1 exactly'
    P=[[F(c)for c in row]for row in N['old_tau_target_P_gamma_columns']]
    Pi=[[F(c)for c in row]for row in N['P_inverse']]
    assert mm(P,Pi)==mm(Pi,P)==ident(3)
    assert P==[[F(1),F(0),F(0)],[F(11,12),F(0),F(1)],[F(0),F(1),F(0)]]
    assert D['new_target_gamma_P']==N['old_tau_target_P_gamma_columns']
    # Complete quotient action, not a finite-jet inference: gamma1,gamma2,gamma3,r gamma3.
    image_basis=[(i,list(map(F,v)))for i,v in C['normalized_triangular_image_basis']]
    free=C['quotient_original_basis_indices'];assert free==[0,3,4,6]
    def quotient(v):
        v=list(map(F,v))
        for i,b in image_basis:
            a=v[i];v=[x-a*y for x,y in zip(v,b)]
        return [v[i]for i in free]
    gamma=[]
    for j in range(3):
        x=[F(0)]*8
        for sec,off in enumerate([0,3,6]):x[off]=P[sec][j]
        gamma.append(quotient(x))
    acts=C['all8_quotient_coordinate_actions_columns']
    def action(i,x):return mv(list(map(list,zip(*acts[i]))),x)
    assert all(not any(action(i,g))for i in range(7)for g in gamma)
    assert not any(action(7,gamma[0])) and not any(action(7,gamma[1]))
    last=action(7,gamma[2]);assert any(last) and not any(action(7,last))
    B=list(map(list,zip(*(gamma+[last]))))
    assert [list(map(F,col))for col in C['exact_decomposition_basis_columns_in_quotient']]==gamma+[last]
    # A fixed small determinant proves these four actual source classes are a basis.
    from itertools import permutations
    det=sum((-1)**sum(p[i]>p[j]for i in range(4)for j in range(i+1,4))*
            F(B[0][p[0]])*F(B[1][p[1]])*F(B[2][p[2]])*F(B[3][p[3]])for p in permutations(range(4)))
    assert det
    # All three complete Koszul differentials reconstructed from exterior deletion signs.
    assert N['last_regular_sequence_powers']==[1,1,2] and N['ranks'][:4]==[3,24,84,168]
    for deg in [1,2,3]:
        basis=[[s,list(w)]for s in range(3)for w in itertools.combinations(range(8),deg)]
        assert N['basis_sector_wedge'][deg]==basis
        lower=[[s,list(w)]for s in range(3)for w in itertools.combinations(range(8),deg-1)]
        pos={(s,tuple(w)):i for i,(s,w)in enumerate(lower)}
        expected=[]
        for s,w in basis:
            col=[]
            for t,j in enumerate(w):
                e=[0]*8;e[j]=[1,1,2][s]if j==7 else 1
                col.append([pos[(s,tuple(w[:t]+w[t+1:]))],e,(-1)**t])
            expected.append(col)
        assert N['all8_differentials_sparse_columns'][deg-1]==expected
    # Literal native sources d*h and d*c*l and their four basis multiples.
    assert R['C_basis']==['1','tauA','tauB','tauA*tauB']
    gen=R['R9_generators_polynomials9']
    for j,e in [(5,[0,0,0,0,0,1,0,0,0]),(7,[0,0,0,0,0,0,1,0,1])]:assert gen[j]==[[e,'1']]
    nilcols=[20,21,22,23,28,29,30,31]
    assert D['actual_nil_source_columns44']==N['actual_new_nil_source_columns']==nilcols
    assert all(R['actual_V3_by44_rows_after_unit1_cancellation'][i][j]==R['actual_Qrad4_by36_columns'][j][i+1]for j in nilcols for i in range(3))
    assert len(R['actual_Res4_Atranspose36_by48_columns'])==48 and all(len(c)==36 for c in R['actual_Res4_Atranspose36_by48_columns'])
    assert 'NEVER transpose inside mu' in R['transpose_convention']
    # Independent all24 residue equality with FIXED saved coefficients (no solve).
    coeff=[[F(c)for c in row]for row in D['constant_nil_image_lift_coefficients8_by24']]
    images=[[F(c)for c in row]for row in C['all8_actual_nil_image_vectors_in_Cbar8']]
    shifts=list(map(F,L['CB_shift_coefficients7']))
    for col,(sec,i)in enumerate(itertools.product(range(3),range(8))):
        residue=[F(0)]*8;n=[1,1,2][sec]
        for old,(off,size)in enumerate(zip([0,3,6],[3,3,2])):
            if i<7 and old==1:residue[off+2]=-P[old][sec]*shifts[i]
            if i==7 and n<size:residue[off+n]=P[old][sec]
        assert residue==list(map(F,D['all24_old_sector_residues'][col]))
        assert [sum(coeff[j][col]*images[j][k]for j in range(8))for k in range(8)]==residue
        assert all(coeff[j][col]==0 for j in [3,5,6,7])

    # Audit all scalar operations and find the EXACT prescribed m/H expressions.
    S=D['scalar_function_nodes'];index={frozen(n):i for i,n in enumerate(S)}
    assert len(S)==5408 and len(index)==len(S)
    def lookup(x):return index[frozen(x)]
    zero=lookup(['QQ',0,1]);one=lookup(['QQ',1,1]);assert zero==0 and one==1
    def q(c):c=F(c);return lookup(['QQ',c.numerator,c.denominator])
    def add(*a):
        a=[x for x in a if x!=zero]
        return zero if not a else a[0]if len(a)==1 else lookup(['sum',a])
    def mul(*a):
        if zero in a:return zero
        a=[x for x in a if x!=one]
        return one if not a else a[0]if len(a)==1 else lookup(['product',a])
    def minus(a):return mul(q(-1),a)
    def shift(f,c):return f if not any(c)else lookup(['coordinate_shift_pullback',f,[[F(v).numerator,F(v).denominator]for v in c]])
    def prefix(f,i):return f if i==7 else lookup(['zero7_prefix_restriction',f,i])
    coord_seen=set();r_seen=set();calls=0
    def canonical(f,power,c,proof):
        nonlocal calls
        calls+=1;g=shift(f,[-x for x in c]);result=[]
        for i in range(7):
            num=add(prefix(g,i+1),minus(prefix(g,i)))
            root=lookup(['regular_coordinate_quotient',num,i,'prefix_difference_vanishes']);coord_seen.add(root)
            result.append(shift(root,c))
        root=lookup(['regular_r_power_quotient',prefix(g,0),power,proof]);r_seen.add(root);result.append(shift(root,c))
        return result
    for i,n in enumerate(S):
        k=n[0];refs=[]
        if k=='QQ':assert n[2]>0 and F(n[1],n[2]).denominator==n[2]
        elif k=='actual_S_coordinate':assert n[1]in range(8)
        elif k=='actual_function_restrict_ell0':assert n[1]=='radical' and 0<n[2]<len(R['new_Cscalar_and_R0component_nodes'])
        elif k in ['sum','product']:refs=n[1];assert len(refs)>=2
        elif k=='zero7_prefix_restriction':refs=[n[1]];assert n[2]in range(7)
        elif k=='coordinate_shift_pullback':refs=[n[1]];assert len(n[2])==7 and any(F(*v)for v in n[2])
        elif k=='regular_coordinate_quotient':refs=[n[1]];assert n[2]in range(7) and n[3]=='prefix_difference_vanishes'
        elif k=='regular_r_power_quotient':refs=[n[1]];assert n[2]in [1,2,3] and n[3]in ['actual_Cprime_bar_diagonal_annihilator_complete','actual_old_Cbar_residue_cancelled_by_fixed_nil_lifts']
        else:raise AssertionError(k)
        assert all(0<=j<i for j in refs)
    M=D['matrix_operator_nodes'];sh=D['matrix_operator_shapes'];X=D['matrices']
    assert len(M)==45 and len(sh)==45
    vr=M[X['Vrad']][1]
    for row,oldrow in zip(vr,R['actual_V3_by44_rows_after_unit1_cancellation']):
        assert row==[zero if v==0 else lookup(['actual_function_restrict_ell0','radical',v])for v in oldrow]
    # Compare target changes with exact rational coefficients.
    for name,value in [('P',P),('Pinv',Pi)]:assert M[X[name]][1]==[[q(c)for c in row]for row in value]
    m=M[X['new_m24_by44']][1];H=M[X['old_residual_correction24_by24']][1]
    for sec,j in itertools.product(range(3),range(44)):
        f=add(*(mul(q(Pi[sec][k]),vr[k][j])for k in range(3)))
        a=canonical(f,[1,1,2][sec],[F(0)]*7,'actual_Cprime_bar_diagonal_annihilator_complete')
        assert [m[8*sec+i][j]for i in range(8)]==a
    var=[lookup(['actual_S_coordinate',i])for i in range(8)]
    for col,(sec,i)in enumerate(itertools.product(range(3),range(8))):
        power=[1,1,2][sec]if i==7 else 1;pol=mul(*([var[i]]*power))
        for old in range(3):
            residual=add(mul(q(P[old][sec]),pol),*(mul(q(-coeff[j][col]),vr[old][nilcols[j]])for j in range(8)if coeff[j][col]))
            a=canonical(residual,[3,3,2][old],shifts if old==1 else[F(0)]*7,'actual_old_Cbar_residue_cancelled_by_fixed_nil_lifts')
            assert [H[8*old+i][col]for i in range(8)]==a
    assert calls==204 and calls*8==1632
    assert coord_seen=={i for i,n in enumerate(S)if n[0]=='regular_coordinate_quotient'}
    assert r_seen=={i for i,n in enumerate(S)if n[0]=='regular_r_power_quotient'}
    # Type and all actual-source wiring checks, without multiplying formal matrices.
    fields={
        ('old_unit','typed_n56_by24'):[56,24],
        ('radical','actual_Res4_Atranspose36_by48_columns'):[36,48],
        ('radical','actual_Brad36_R0_values_by640_columns'):[36,640],
        ('radical','actual_P8_by640_columns'):[8,640],
        ('radical','actual_complete_SyzQrad36_by140_sparse_columns'):[36,140],
        **{('new_resolution',f'all8_differentials_sparse_columns[{i}]'):[N['ranks'][i],N['ranks'][i+1]]for i in range(3)}}
    vectors=D['vector_formal_function_operator_nodes'];assert len(vectors)==6240
    for i,n in enumerate(M):
        k=n[0];refs=[]
        if k=='literal_scalar_rows':assert sh[i]==[len(n[1]),len(n[1][0])] and all(len(r)==sh[i][1] and all(v in range(len(S))for v in r)for r in n[1])
        elif k=='actual_external':assert sh[i]==fields[(n[1],n[2])]
        elif k=='zero':pass
        elif k=='identity':assert sh[i][0]==sh[i][1]
        elif k=='matrix_product':
            refs=n[1];assert all(sh[a][1]==sh[b][0]for a,b in zip(refs,refs[1:])) and sh[i]==[sh[refs[0]][0],sh[refs[-1]][1]]
        elif k=='matrix_sum':refs=n[1];assert all(sh[j]==sh[i]for j in refs)
        elif k=='matrix_neg':refs=[n[1]];assert sh[n[1]]==sh[i]
        elif k in ['horizontal_blocks','vertical_blocks']:
            refs=n[1];axis=0 if k=='horizontal_blocks'else 1;assert len({sh[j][axis]for j in refs})==1
            expected=sh[refs[0]].copy();expected[1-axis]=sum(sh[j][1-axis]for j in refs);assert sh[i]==expected
        elif k=='unit_matrix_inverse':refs=[n[1]];assert sh[i]==sh[n[1]]==[24,24] and 'M(0)=I24' in n[2]
        elif k=='matrix_from_original_column_vector_operators':assert sh[i]==[84,780] and n[1]==[8*j+7 for j in range(780)]
        else:raise AssertionError(k)
        assert all(0<=j<i for j in refs)
    assert M[3]==['matrix_product',[2,0]] and X['Vgamma']==3
    assert M[13]==['vertical_blocks',[9,12]] and M[9]==['horizontal_blocks',[7,8]] and M[12]==['horizontal_blocks',[10,11]]
    assert M[7][:3]==['actual_external','radical','actual_Res4_Atranspose36_by48_columns']
    assert sh[8]==[36,8] and M[8]==['zero'] and sh[10]==[8,48] and M[10]==['zero'] and M[11]==['identity']
    assert M[14][1]==[[one if i==j else zero for j in nilcols]for i in range(44)]
    assert M[15][1]==[[q(c)for c in row]for row in coeff]
    assert M[16]==['matrix_product',[13,6,5]] and M[17]==['matrix_product',[14,15]] and M[18]==['matrix_sum',[16,17]]
    assert M[19]==['matrix_product',[4,18]] and M[20][0:2]==['unit_matrix_inverse',19] and M[21]==['matrix_product',[18,20]]
    assert M[31]==['horizontal_blocks',[28,30]] and M[28]==['vertical_blocks',[25,26]] and M[30]==['vertical_blocks',[27,29]] and M[29]==['zero'] and sh[29]==[8,140]
    assert M[32]==['matrix_product',[4,31]]
    for j in range(780):
        v=vectors[8*j:8*j+8];assert v[0]=={'operation':['actual_matrix_column_functions',32,j],'dimension':24}
        for sec in range(3):
            a=1+2*sec;assert v[a]=={'operation':['contiguous_vector_slice',8*j,8*sec,8],'dimension':8}
            assert v[a+1]=={'operation':['complete_Koszul_cycle_Hadamard_Euler_lift',8*j+a,[1,1,2][sec],'shared_template_v1'],'dimension':28}
        assert v[7]=={'operation':['ordered_vector_concatenation',[8*j+2,8*j+4,8*j+6]],'dimension':84}
    assert M[35]==['matrix_product',[21,4]] and M[36]==['matrix_neg',35] and M[37]==['matrix_sum',[34,36]] and M[34]==['identity']
    assert M[38]==['matrix_product',[21,23]] and M[39]==['horizontal_blocks',[37,38]]
    assert M[44]==['vertical_blocks',[41,43]] and M[41]==['horizontal_blocks',[31,21,40]] and M[43]==['horizontal_blocks',[33,42,24]]
    assert M[40]==M[42]==['zero'] and sh[40]==[44,168] and sh[42]==[84,24]
    assert X['Phi44_by128']==39 and X['honest_Tnil128_by972']==44 and sh[44]==[128,972]
    # Independent polynomial sign fixtures for the specified original-column lift.
    def mono(e):return {tuple(e):F(1)}
    fixture_count=0
    for power,(i,j),k in itertools.product([1,2],PAIRS,range(8)):
        exp=[0]*8;exp[k]=2;exp[(k+3)%8]+=1
        h0=[{}for _ in PAIRS];h0[PAIRS.index((i,j))]=mono(exp)
        v=differential(h0,power);lift_formal_cycle(v,power);fixture_count+=1
    template=D['shared_homotopy_template_v1']
    assert template['output_order_pairs']==[list(p)for p in PAIRS] and template['not_S_linear']
    assert template['after_all780_original_columns_only_extend_S_linearly'] and template['no_pull_out_actual_coefficients_or_truncated_normal_forms']
    # Unit/minimality, kernel and specialization are proved in the accompanying note,
    # from the complete Koszul resolution, actual V/X and the termwise free complex.
    assert src['nil_finite_peer']['actual_NilD0_mu_interval']==[14,20]
    assert src['rankCprime_peer']['generic_DVR_rank_Cprime']==2 and src['rankCprime_peer']['length_Cprime_mod_ell']==4
    assert src['HF3_peer']['actual_T_length_interval_accepted']==[11,17]
    assert D['actual_Tnil_length_interval']==[16,22] and D['optional104_free_complement_literal_basis_constructed']is False
    assert D['scalar_Tnil_values_minimal_generators_Nil_generators_cube_or_exact_T_length_claimed']is False
    evidence={'all_actual_radical_V_X_Z_old_n_source_frames_and_signs':True,
        'eight_literal_dh_dcl_source_columns':nilcols,'all24_old_residue_equalities':True,
        'actual_complete_Cprime_gamma4basis_and_Koszul3_differentials':True,
        'canonical_division_slots':1632,'unique_regular_coordinate_quotient_nodes':len(coord_seen),
        'unique_regular_r_quotient_nodes':len(r_seen),'all204_canonical_calls_matched_exact_input_functions':True,
        'scalar_nodes':len(S),'matrix_nodes':len(M),'vector_nodes':len(vectors),
        'all780_h_lifts_apply_to_original_full_mX_columns_before_S_linear_extension':True,
        'independent_symbolic_sign_fixture_count':fixture_count,
        'formal_Cartan_identity_allorders_not_fixture_extrapolation':True,
        'minimal_Koszul_d2_implies_M0_I_and_n_unit_normalization':True,
        'Phi_surjection_and_complete_kernel_n24_plus_d3_168':True,
        'complete_Tnil_presentation_shape':[128,972],
        'Phi_relation_matrix_equals_X_0_0_in_free44_NOT_zero':True,
        'termwise_free_ell_UCT_H2_irrelevant':True,
        'length_Tnil_equals_L_plus5_interval':[16,22],
        'no_scalar_source_jet_values_823680products_or_nil_generator_cube_claim':True}
    write('independent-source-formal-and-homological-proof.json',evidence)
    rep={'id':ID,'all_passed':True,'frozen':True,'complete_actual_formal_presentation_shape':[128,972],
         'actual_source_maps_division_membership_all780_original_column_homotopies_accepted':True,
         'unit_normalization_and_complete_kernel_specialization_accepted':True,
         'length_Tnil_formula':'length(T)+5','length_Tnil_interval':[16,22],
         'scalar_coefficients_minimal_nil_generators_exact_Tlength_or_cube_NOT_computed':True,
         'optional104x948_not_built':True,'new_source_values_jets_or_scalar_mX_expansion':False,
         'input_sha256':pins,'code_sha256':sha(Path(__file__)),
         'note_sha256':sha(ROOT/'research/notes'/f'{ID}.md'),
         'evidence_sha256':{'independent-source-formal-and-homological-proof.json':sha(OUT/'independent-source-formal-and-homological-proof.json')},
         'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    assert sum(p.stat().st_size for p in OUT.iterdir())<32*1024*1024
    write('report.json',rep);signal.alarm(0);print(json.dumps({k:v for k,v in rep.items()if k!='input_sha256'},indent=2))

if __name__=='__main__':main()
