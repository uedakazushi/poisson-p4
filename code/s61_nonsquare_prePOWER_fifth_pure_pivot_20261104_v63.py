#!/usr/bin/env python3
"""Typed fifth pure pivot, global heavy coverage, and bordered-six support only."""
import argparse, hashlib, importlib.util, itertools, json, math, resource, signal, time
from fractions import Fraction as Q
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
ID='S61_NONSQUARE_PREPOWER_FIFTH_PURE_PIVOT_20261104_v63'
def read(p): return json.loads(p.read_text())
def sha(p):
    with p.open('rb') as f: return hashlib.file_digest(f,'sha256').hexdigest()
def dump(p,d): p.write_text(json.dumps(d,indent=2)+'\n')

def formal_check():
    # A separate universal Laurent ring: v and all entries c,r,D are independent.
    def atom(x,k=1): return {((x,k),):1}
    one={():1}; zero={}
    def add(a,b,sign=1):
        z=a.copy()
        for m,c in b.items():
            z[m]=z.get(m,0)+sign*c
            if not z[m]: del z[m]
        return z
    def mul(a,b):
        z={}
        for m,c in a.items():
            for n,d in b.items():
                q=dict(m)
                for x,k in n:q[x]=q.get(x,0)+k
                e=tuple(sorted((x,k)for x,k in q.items()if k))
                z[e]=z.get(e,0)+c*d
                if not z[e]:del z[e]
        return z
    def mm(a,b):
        out=[]
        for row in a:
            rr=[]
            for j in range(len(b[0])):
                s={}
                for k,c in enumerate(row):s=add(s,mul(c,b[k][j]))
                rr.append(s)
            out.append(rr)
        return out
    def eye(n):return [[one if i==j else zero for j in range(n)]for i in range(n)]
    v=atom('v');iv=atom('v',-1);r=[atom('r'+str(j))for j in range(12)];c=[atom('c'+str(i))for i in range(5)]
    D=[[atom('d%d_%d'%(i,j))for j in range(12)]for i in range(5)]
    A=[[v]+r]+[[c[i]]+D[i]for i in range(5)]
    L=eye(6);Li=eye(6);L[0][0]=iv;Li[0][0]=v
    U=eye(13);Ui=eye(13)
    for i in range(5):L[i+1][0]=mul({():-1},mul(c[i],iv));Li[i+1][0]=c[i]
    for j in range(12):U[0][j+1]=mul({():-1},mul(r[j],iv));Ui[0][j+1]=mul(r[j],iv)
    for a,b,n in [(L,Li,6),(Li,L,6),(U,Ui,13),(Ui,U,13)]:assert mm(a,b)==eye(n)
    T=[[add(D[i][j],mul(mul(c[i],r[j]),iv),-1)for j in range(12)]for i in range(5)]
    target=[[one]+[zero]*12]+[[zero]+T[i]for i in range(5)]
    assert mm(mm(L,A),U)==target
    return dict(all_four_inverse_products=True,all78_diagonalization_entries=True,
                formal_variables=78,universal_T_coefficients=120,new_parameter_coefficients=0)

def main(out):
    start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60)
    out.mkdir(parents=True,exist_ok=False);pins={}
    def load(st,name='certificate.json'):
        p=ROOT/'reports'/st/name;pins[str(p.relative_to(ROOT))]=sha(p);return read(p)
    def gate(st,key='all_passed'):
        d=load(st,'report.json');assert d[key];return d
    direct='S61_NONSQUARE_PREPOWER_DIRECT_P4_H5_SCHUR_20261104_v63'
    dop=load(direct,'operator.json');gate(direct+'_REVIEW')
    ent='S61_NONSQUARE_PREPOWER_U4_T15_SCHUR_ENTRY_20261104_v63'
    e=load(ent,'native.json');assert not e['numerator_is_zero']and e['numerator_terms']==7022
    # Keep only metadata; the accepted 7022 coefficients are referenced, not duplicated.
    del e
    binding=load(ent,'source-binding.json')
    binding={k:binding[k]for k in ['Mark_entry','scalar_relation']}
    gate('S61_NONSQUARE_PREPOWER_U4_T15_MODP_POLYNOMIAL_20261104_v63_REVIEW')
    assert binding['Mark_entry']==[0,9]and dop['context']['Q9_eight_other_pure_tails']
    for st,key in [
      ('S61_NONSQUARE_RAW_H_P_BRANCH_CHAIN_20261104_v63_REVIEW','all_passed'),
      ('S61_NONSQUARE_RAW_H_GENERIC_RANK_20261104_v63_REVIEW','all_passed'),
      ('S61_NONSQUARE_THETA_RANK_DROP_CURVES_20261104_v63_REVIEW','all_passed'),
      ('S61_NONSQUARE_THETA_P80_REPLACEMENT_20261104_v63_REVIEW','ALL_P80_full_cubic57_proved'),
      ('S61_NONSQUARE_THETA_P100_REPLACEMENT_20261104_v63_REVIEW','ALL_P100_full_cubic57_proved'),
      ('S61_NONSQUARE_RAW_H_B_K13_CUBIC_REGULAR_REP_20261104_v63_REVIEW','whole_actual_Ka13_Bbad_full_cubic57')]:gate(st,key)
    coverage=dict(scope='All geometric points of FULLU; not equality of heavy Fitting schemes.',
      cover=['D(p)','V(g1)','V(Theta)'],p_factor='750000p=g1*Theta',comaximal='1=Theta/B^2+C*g2*g1/B^2',
      Dp_and_g1_rank='rawHrank=4+rank(B7x3)',Theta_rank='rawHrank6 on P80*P100=0, rank7 elsewhere',
      Theta_B_selected_minor_is_global_unit=True,
      rawH_rank_below7_exact_geometric_support=['actual Bbad13','Theta P80','Theta P100'],
      all_three_exceptional_sets_already_full_cubic57=True,
      unresolved_actual_fullcubic_defect_contained_in_H7=True,
      H7_pure_source='local free kernel10 -> pure17; quotient by unimodular Q9 is local free9 -> pure16',
      projected_full_rank_only_sufficient_for_actual_full57=True)
    nodes=[]
    def n(op,sh,args=(),**kw):
        assert all(0<=i<len(nodes)for i in args);i=len(nodes);nodes.append(dict(id=i,op=op,shape=sh,args=list(args),**kw));return i
    def ref(exp,sh):return n('accepted_external_export',sh,reference=dict(path=f'reports/{direct}/operator.json',sha256=pins[f'reports/{direct}/operator.json'],export=exp,node=dop['exports'][exp]))
    def sel(a,rr,cc):return n('submatrix',[len(rr),len(cc)],[a],rows=rr,columns=cc)
    def prod(a,b):assert nodes[a]['shape'][1]==nodes[b]['shape'][0];return n('matrix_product',[nodes[a]['shape'][0],nodes[b]['shape'][1]],[a,b])
    def neg(a):return n('negative',nodes[a]['shape'],[a])
    def scale(c,a):return n('central_scalar_multiple',nodes[a]['shape'],[c,a])
    def eye(k):return n('identity',[k,k])
    def zero(m,k):return n('zero',[m,k])
    def block(rows):
        h=[nodes[r[0]]['shape'][0]for r in rows];w=[nodes[x]['shape'][1]for x in rows[0]]
        for i,row in enumerate(rows):
            for j,x in enumerate(row):assert nodes[x]['shape']==[h[i],w[j]]
        return n('block_matrix',[sum(h),sum(w)],[x for r in rows for x in r],blocks=rows)
    E=ref('residual6x13',[6,13]);v=sel(E,[0],[9]);iv=n('certified_chart_scalar_inverse',[1,1],[v],proof=binding['scalar_relation'],new_chart_factor='Saux7022')
    jj=[j for j in range(13)if j!=9];rr=list(range(1,6));r=sel(E,[0],jj);c=sel(E,rr,[9]);D=sel(E,rr,jj)
    T=n('difference',[5,12],[D,scale(iv,prod(c,r))]);Tclear=n('difference',[5,12],[scale(v,D),prod(c,r)])
    L=block([[iv,zero(1,5)],[neg(scale(iv,c)),eye(5)]]);Li=block([[v,zero(1,5)],[c,eye(5)]])
    U=block([[eye(1),neg(scale(iv,r))],[zero(12,1),eye(12)]]);Ui=block([[eye(1),scale(iv,r)],[zero(12,1),eye(12)]])
    Pc=n('column_permutation',[13,13],order=[9]+jj);Pci=n('transpose',[13,13],[Pc])
    S15=block([[eye(9),zero(9,6)],[zero(6,9),L]]);S15i=block([[eye(9),zero(9,6)],[zero(6,9),Li]])
    V13=prod(Pc,U);V13i=prod(Ui,Pci)
    V22=block([[eye(9),zero(9,13)],[zero(13,9),V13]]);V22i=block([[eye(9),zero(9,13)],[zero(13,9),V13i]])
    fs=prod(S15,ref('full_source15',[15,15]));fsi=prod(ref('full_source15_inverse',[15,15]),S15i)
    ft=prod(ref('full_target22',[22,22]),V22);fti=prod(V22i,ref('full_target22_inverse',[22,22]))
    old25=prod(sel(L,rr,list(range(6))),ref('surviving6_source_original25',[6,25]))
    remaining_labels=[dop['context']['remaining13_original32_seed_labels'][j]for j in jj]
    assert remaining_labels==[1,2,3,4,6,7,8,9,14,17,18,19]
    operator=dict(nodes=nodes,exports=dict(residual5x12=T,cleared_residual5x12=Tclear,source6=L,source6_inverse=Li,
      target13=V13,target13_inverse=V13i,full_source15=fs,full_source15_inverse=fsi,full_target22=ft,full_target22_inverse=fti,
      surviving5_source_original25=old25),chart='FULLU D(R434*Delta5*Saux)',
      remaining_current15=[8,9,11,12,14],remaining_original25=[18,19,21,22,24],remaining_original32=remaining_labels,
      extra8_Q9_tails_retained=remaining_labels[:8],new_parameter_coefficients=0,
      identities=['L6*E6x13*Pcol*U13=diag(1,T5x12)','all four inverse products are identities',
       'full_source15*M15x22*full_target22=diag(I10,T5x12)',
       'I6(E)=I5(T) on stated chart; this is a sufficient-source criterion, not converse actual57',
       'surviving5_source_original25=L6[1:6,:]*previous_surviving6_original25'],
      unresolved_complements=['p=0 (covered only for heavy rank by g1/Theta analysis)','Delta5=0','R434=0','Saux=0'])
    dump(out/'operator.json',operator);dump(out/'coverage.json',coverage);dump(out/'formal-check.json',formal_check())
    # Only source supports. No parameter product, determinant, or scalar cancellation is evaluated.
    sp=importlib.util.spec_from_file_location('bounds',ROOT/'code/s61_nonsquare_prePOWER_h22_kernel9_support_preflight_20261104_v63.py');b=importlib.util.module_from_spec(sp);sp.loader.exec_module(b)
    old=load('S61_NONSQUARE_PREPOWER_H22_KERNEL9_SUPPORT_PREFLIGHT_20261104_v63')
    q=load('S61_NONSQUARE_CUBIC_Q2_L_UNIT_COMPLETION_20261104_v63');cp=load('S61_NONSQUARE_CUBIC_ZW25_UNITS_20261104_v63','checkpoint.json');pd=load('S61_NONSQUARE_CUBIC_ZW25_PAIR_DUAL_PREFLIGHT_20261104_v63')
    h=load('S61_NONSQUARE_PURE_KERNEL12_SOURCE_PREFLIGHT_20261104_v63','prePOWER-h9-source.json');fixed=load('S61_NONSQUARE_PREPOWER_H22_PURE3_FIXED_MINOR_20261104_v63','native.json')
    def literal_terms(x):
        if isinstance(x,list):
            if x and all(isinstance(t,list)and len(t)==2 and isinstance(t[0],list)and len(t[0])==3 and all(isinstance(a,int)for a in t[0])and isinstance(t[1],str)for t in x):return len(x)
            return sum(map(literal_terms,x))
        return sum(map(literal_terms,x.values()))if isinstance(x,dict)else 0
    # Entire source-object storage, including unused checkpoint/cofactor duplicates.
    source_terms=sum(literal_terms(x)for x in [q,cp,pd,h,fixed])
    assert source_terms<150000
    cached={x['old_DAG_root']:b.atom(x['Laurent_polynomial'])for x in fixed['node_values']};values=[]
    del fixed
    def rf(nd):
        num,den=nd
        if len(den)==1:assert den[0][0][1:]==[0,0]and Q(den[0][1]);return b.shift(b.atom(num),-den[0][0][0])
        dd={tuple(m):Q(c)for m,c in den};bb={tuple(m):Q(c)for m,c in b.BP};ratio=dd[next(iter(bb))]/bb[next(iter(bb))]
        assert dd.keys()==bb.keys()and ratio and all(dd[m]==ratio*bb[m]for m in bb)
        x=b.atom(num);return None if x is None else(x[0],x[1],1)
    for i,a in enumerate(q['scalar_DAG']['nodes']):
        if i in cached:x=cached[i]
        elif a[0]=='z':x=None
        elif a[0]=='one':x=b.ONE
        elif a[0]=='checkpoint':x=rf(cp['matrix25x32'][a[1]][a[2]])
        elif a[0]=='accepted_pair_second':x=rf(pd['Lframe_dual6x2'][a[1]][1])
        elif a[0]=='new_unit_atom':x=rf(q['scalar_DAG']['new_unit_atoms_num_den'][a[1]])
        elif a[0]=='m':x=b.mm(values[a[1]],values[a[2]])
        else:x=None if a[0]=='s'and a[1]==a[2]else b.aa(values[a[1]],values[a[2]])
        values.append(x)
    cols=q['new22_target_original32_columns'][:17];M=[[values[i]for i in row[:17]]for row in q['remaining15x22_scalar_roots']]
    for lab,p in zip(h['original_target_labels'],h['source_h9']):M[0][cols.index(lab)]=b.atom(p)
    piv=[22,11,13,25];P=[[M[i][cols.index(j)]for j in piv]for i in [0,1,3,6]]
    def detbound(A):
        states={0:b.ONE}
        for row in A:
            nxt={}
            for mask,x in states.items():
                for j,y in enumerate(row):
                    if not(mask>>j&1):nxt[mask|1<<j]=b.aa(nxt.get(mask|1<<j),b.mm(x,y))
            states=nxt
        return states[(1<<len(A))-1]
    def fromenc(x):return None if x['zero']else(x['support_upper'],tuple(tuple(z)for z in x['Laurent_exponent_box']),x['known_B_denominator_power'])
    Bhat=[]
    for profile in old['six_staged_coupled_row_bounds']:
        src={int(i):fromenc(x)for i,x in profile['source_cofactor_support_bounds'].items()}
        row=[]
        for j in range(17):
            z=None
            for i,x in src.items():z=b.aa(z,b.mm(x,M[i+1][j]))
            row.append(z)
        Bhat.append(row)
    candidates=[]
    for ii in range(1,6):
        for lab in remaining_labels:
            labels=piv+[15,lab]
            A=[[M[i][cols.index(j)]for j in labels]for i in [0,1,3,6]]+[[Bhat[i][cols.index(j)]for j in labels]for i in [0,ii]]
            z=detbound(A)
            input_sum=sum(b.size(x)for row in A for x in row)
            candidates.append(dict(remaining_row=ii-1,compact_source=old['six_staged_coupled_row_bounds'][ii]['compact_source_row'],original32_target=lab,
              heavy_kind=old['six_staged_coupled_row_bounds'][ii]['kind'],direct_bordered6_support=b.enc(z),
              input_entry_support_sum=input_sum,output_plus_input_support_upper=input_sum+b.size(z),
              scalar_evaluation_authorized=False,resource_gate='Support-only; does not certify lifetime, coefficient height, or serialization.'))
    candidates.sort(key=lambda x:(x['output_plus_input_support_upper'],x['direct_bordered6_support']['support_upper'],x['remaining_row'],x['original32_target']))
    costs=dict(source_only=True,all60_candidates=candidates,smallest=candidates[0],
      all60_output_support_upper=sum(x['direct_bordered6_support']['support_upper']for x in candidates),
      scalar_source_explicit_terms_loaded_upper=source_terms,new_actual_parameter_terms=0,
      direct_identity='S0t*Sij-Sit*S0j=D4*detborder6(P4,Bhat0,Bhati; piv4,t,j)',
      identity_is_universal_polynomial_after_D4_localization_then_domain_descent=True,
      next_value_job_NOT_authorized=True,
      caveat='No use of the 7022 polynomial in products; all60 values remain unexpanded. All B poles are tracked.')
    dump(out/'support-preflight.json',costs)
    note=ROOT/f'research/notes/{ID}.md';pins[str(note.relative_to(ROOT))]=sha(note);pins[str(Path(__file__).relative_to(ROOT))]=sha(Path(__file__))
    artifact=sum(p.stat().st_size for p in out.iterdir());assert artifact<32*1024*1024
    report=dict(id=ID,all_passed=True,frozen=True,source_and_type_only=True,typed_nodes=len(nodes),
      actual_new_scalar_coefficients=0,actual_same_cokernel_shape=[5,12],all8_Q9_tails_retained=True,
      global_heavy_bad_support_already_full57=True,global_cubic57_NOT_proved=True,
      no_new_point_prime_minor_rank_GB=True,scalar_source_explicit_terms_loaded_upper=source_terms,
      smallest_candidate=candidates[0],all60_output_support_upper=costs['all60_output_support_upper'],
      input_sha256=pins,evidence_sha256={p.name:sha(p)for p in out.iterdir()},raw_artifact_bytes_before_report=artifact,
      seconds=time.monotonic()-start)
    dump(out/'report.json',report);signal.alarm(0)
    print(json.dumps({k:v for k,v in report.items()if k not in ['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True);main(ap.parse_args().output.resolve())
