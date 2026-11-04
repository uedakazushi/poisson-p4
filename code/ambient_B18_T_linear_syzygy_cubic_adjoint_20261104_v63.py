#!/usr/bin/env python3
"""Two cubic quotient covectors and X2-only source contractions, no X2 query."""
import argparse,itertools,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
from ambient_B18_comparison_T_constant_rank_20261104_v63 import read,sha,rec,encode,addscaled
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_T_LINEAR_SYZYGY_CUBIC_ADJOINT_20261104_v63'
def dump(out,name,x):(out/name).write_text(json.dumps(x,separators=(',',':'))+'\n')
def put(v,k,a):
    z=v.get(k,F(0))+a
    if z:v[k]=z
    else:v.pop(k,None)
def enc(v):return [[list(k) if isinstance(k,tuple) else k,str(a)] for k,a in sorted(v.items()) if a]
def main(out):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
    out=out.resolve();out.mkdir(parents=True,exist_ok=False)
    paths={
      'rewrite':'reports/AMBIENT_B18_T_KNOWN_QUADRATIC_CUBIC_CLOSURE_20261104_v63/linear-module-rewrite-certificate.json',
      'cubic':'reports/AMBIENT_B18_T_KNOWN_QUADRATIC_CUBIC_CLOSURE_20261104_v63/actual-known-quadratic-cubic-closure-certificate.json',
      'linear':'reports/AMBIENT_B18_COMPARISON_T_LINEAR_QUADRATIC_20261104_v63/phase1-actual-five-target-linear-relations.json',
      'mu':'reports/AMBIENT_B18_COMPARISON_T_CONSTANT_RANK_20261104_v63/actual-linear-cycle-muT-certificate.json',
      'm0':'reports/AMBIENT_B18_FULLV56_COMPARISON_UNIT_20261104_v63/actual-fullV56-m0-and-typed-comparison.json',
      'pi':'reports/AMBIENT_B18_T_QUADRATIC_M1_CANCELLATION_PREFLIGHT_20261104_v63/actual-projected-three-class-X2-operator-preflight.json',
      'syzygies':'reports/AMBIENT_B18_T_LINEAR_SYZYGY_PURE_R3_20261104_v63/actual-linear-syzygies-and-pure-r3-source-values.json',
      'a2':'reports/AMBIENT_B18_T_CUBIC_X2_ONLY_PREFLIGHT_20261104_v63/actual119-cubic-X2-only-source-operator-preflight.json',
      'chain':'reports/AMBIENT_B18_NORMALIZER_CHAIN_LIFT_640_NEGATIVE_BOUNDARY_20261104_v63/actual-UP-QB-640-typed-operator.json',
      'C0':'reports/AMBIENT_B18_T0_HONEST_C0_20261104_v63/honest-C0-D0-D4-DAG.json'}
    hashes={p:sha(ROOT/p) for p in paths.values()};dump(out,'phase.json',{'phase':'fixed_input_hashes_before_adjoints','input_sha256':hashes})
    D={k:read(ROOT/p) for k,p in paths.items()};NF,C,L,M,V,P,Y,A,B,H=[D[k] for k in paths]
    assert sha(ROOT/paths['a2'])=='b611397190163530969212be236a544cf95b3227816d2ec095c332f6649fdeb7'
    mons2=list(itertools.combinations_with_replacement(range(8),2));mi2={m:i for i,m in enumerate(mons2)}
    mons3=list(itertools.combinations_with_replacement(range(8),3));mi3={m:i for i,m in enumerate(mons3)}
    assert NF['monomials3_120']==[list(m) for m in mons3]
    rank=C['complete_reduced_rank_certificate'];rb=list(map(rec,rank['basis_normalized_columns']));lead={min(v):i for i,v in enumerate(rb)}
    standard=NF['standard_global600_positions'];missing=[i for i in range(len(standard)) if i not in lead];assert missing==[3,4]
    quotient=[]
    for j in range(5):
        v={j:F(1)}
        for i in sorted(lead):
            if i in v:addscaled(v,rb[lead[i]],-v[i])
        assert set(v)<=set(missing);quotient.append({k:v.get(i,F(0)) for k,i in enumerate(missing) if v.get(i)})
    Pi=[{} for _ in range(2)]
    for index,step,v in NF['all600_rewrite_steps']:
        for j,c in rec(v).items():
            for k,a in quotient[j].items():put(Pi[k],index,c*a)
    assert Pi[1]=={599:F(1)}
    # Exact annihilation of every known cubic relation, before lifting its
    # functional through the constant comparisons.
    for col in C['all1424_reduced_cubic_columns']:
        v=rec(col);assert all(sum((v.get(j,F(0))*a for j,a in d.items()),F(0))==0 for d in [{j:quotient[j].get(k,F(0)) for j in range(5)} for k in range(2)])
    qnorm=list(map(rec,L['normalized84_unitcycle_Q1_columns_sparse']))
    beta=[]
    for d in Pi:
        b={}
        for p,q in enumerate(qnorm):
            for m,pair in enumerate(mons2):
                z=-sum((a*d.get(120*(j//8)+mi3[tuple(sorted(pair+(j%8,)))],F(0)) for j,a in q.items()),F(0))
                if z:b[36*p+m]=z
        beta.append(b)
    wb=list(map(rec,M['complete_rank_certificate']['basis_normalized_columns']));wl={min(v):j for j,v in enumerate(wb)};wcoords=[]
    for j in range(84):
        v={j:F(1)};c={}
        while v:
            i=min(v);p=wl[i];a=v[i];c[p]=c.get(p,F(0))+a;addscaled(v,wb[p],-a)
        wcoords.append(c)
    wedge=M['literal_wedge_basis84'];m0=[list(map(F,row)) for row in V['all24_by56_m0']]
    pirows=[]
    for row in P['pi_free32_rows']:
        pirows.append({P['free_kernel_rows32'][j]:c for j,c in rec(row).items()})
    kappa=[];aadj=[]
    for t,d in enumerate(Pi):
        kadj={}
        # Direct five-component projection.
        for index,c in d.items():
            target,m=divmod(index,120)
            for src,a in pirows[target].items():put(kadj,120*src+m,c*a)
        a={}
        for w,(sector,pair) in enumerate(wedge):
            i,j=pair
            for m,pair2 in enumerate(mons2):
                c=sum((b*beta[t].get(36*p+m,F(0)) for p,b in wcoords[w].items()),F(0))
                if not c:continue
                if j==7:terms=[(8*sector+7,i,1)]
                else:terms=[(8*sector+j,i,1),(8*sector+i,j,-1)]
                for outrow,var,sign in terms:
                    mon=tuple(sorted(pair2+(var,)));k=sum(v<7 for v in mon);den=k if j==7 else k+1
                    z=c*F(sign*mon.count(var),den);put(a,120*outrow+mi3[mon],z)
        for index,c in a.items():
            row,m=divmod(index,120)
            for src,b in enumerate(m0[row]):
                if b:put(kadj,120*src+m,c*b)
        kappa.append(kadj);aadj.append(a)
    assert kappa[1]=={6719:F(1)}
    dump(out,'phase1-cubic-Pi2-and-source-adjoint.json',{'monomials2_36':mons2,'monomials3_120':mons3,
      'two_quotient_dual_rows600':list(map(encode,Pi)),'two_normalized_wedge_quadratic_adjoints':list(map(encode,beta)),
      'two_m0Z3_quadratic_homotopy_adjoints':list(map(encode,aadj)),
      'two_Xsource3_adjoint_rows56_times120':list(map(encode,kappa)),
      'homotopy_Euler_by_zdegree':'wij=(di aj-dj ai)/(k+1),wi7=di a7/k; degree3 input',
      'Euler_formula_valid_on_actual_cycles_mod_known_quadratic_times_variables':True,
      'only_X2_after_linear_coefficient_contraction':True,
      'Pi_supports':[[list(mons3[i%120]) for i in sorted(v)] for v in Pi],
      'source_adjoint_supports':[[list(mons3[i%120]) for i in sorted(v)] for v in kappa]})
    dump(out,'phase.json',{'phase':'exact_cubic_duals_and_source_homotopy_adjoints_complete','Pi_nonzeros':list(map(len,Pi)),'source_adjoint_nonzeros':list(map(len,kappa))})
    source37=[rec(r['actual768_constant_source_combination']) for r in Y['actual37_e5_r2_source_coefficients']]
    assert [encode(v) for v in source37]==A['actual_selected37_original768_constant_combinations']
    syz=[rec(r['coefficients296_sparse']) for r in Y['all119_linear_coefficient_syzygies']]
    c1=[]
    for c in syz:
        rows=[{} for _ in range(8)]
        for index,a in c.items():j,v=divmod(index,8);addscaled(rows[v],source37[j],a)
        c1.append(rows)
    bcols=[dict(c) for c in B['B48_by640_sparse_columns']];bops=B['operator_scalar_nodes']
    def key(j,row,pair):
        p=mi2[pair]
        if j>=640:return (2,j-640,row,p) if row<48 else None
        if row<48:
            root=bcols[j].get(row)
            if root is None:return None
            n=bops[root]
            if n[0]=='c':return None
            return (1,n[1],p) if n[0]=='actual_C0_scalar_function' else (0,root,p)
        root=H['honest_D0_relation_columns'][j][row-48+1]
        return None if H['DAG_scalar_nodes'][root][0]=='c' else (1,root,p)
    # Structural labels before allocating Fraction functionals: no source
    # values are assigned. It gives a tight output-support upper bound.
    symbolic=[];upper=0
    for rows in c1:
        per=[]
        for kap in kappa:
            labels=set()
            for index in kap:
                src,m=divmod(index,120);mon=mons3[m]
                for v in set(mon):
                    pair=list(mon);pair.remove(v);pair=tuple(pair)
                    for j in rows[v]:
                        q=key(j,src,pair)
                        if q is not None:labels.add(q)
            upper+=len(labels);per.append(labels)
        symbolic.append(per)
    forecast=upper+sum(len(v) for rows in c1 for v in rows)+sum(map(len,kappa))+6000
    dump(out,'preflight.json',{'two_cubic_dual_nonzeros':list(map(len,Pi)),'two_source_adjoint_nonzeros':list(map(len,kappa)),
      'all119_source_linear_combinations_QQ_support':sum(len(v) for rows in c1 for v in rows),
      'folded_scalar_root_functional_structural_support_bound':upper,'TOTAL_bound':forecast,'cap':150000,
      'no_new_source_X2_value_or_jet':True})
    assert forecast<150000,('SOURCE_ADJOINT_SUPPORT_GUARD',forecast)
    xp1=list(map(rec,A['actual_Xp1_sparse27_columns']));a2=list(map(rec,A['all119_actual_a2_sparse27_times36']))
    folded=[];neededpairs=set();rootkinds={0:set(),1:set(),2:set()};support=0
    for number,rows in enumerate(c1):
        fs=[];cor=[]
        for t,kap in enumerate(kappa):
            f={};known=F(0)
            for index,c in kap.items():
                src,m=divmod(index,120);mon=mons3[m]
                for v in set(mon):
                    pair=list(mon);pair.remove(v);pair=tuple(pair)
                    for j,a in rows[v].items():
                        q=key(j,src,pair)
                        if q is not None:put(f,q,c*a)
                    # Xp1*a2 contributes exactly once per distinct first
                    # variable; ordinary coefficient convention, no Hessians.
                    for p,first in enumerate(xp1):known-=c*first.get(8*src+v,F(0))*a2[number].get(36*p+mi2[pair],F(0))
            assert set(f)<=symbolic[number][t]
            fs.append(f);cor.append(known)
            for q in f:
                neededpairs.add(mons2[q[-1]]);rootkinds[q[0]].add(q[1])
            support+=len(f)+bool(known)
        folded.append({'syzygy119_index':number,'two_scalar_root_X2_functionals':list(map(enc,fs)),
          'known_Xp1_a2_correction':list(map(str,cor))})
    assert support+sum(len(v) for rows in c1 for v in rows)<150000
    b11={(0,7),(1,7),(4,7),(5,7),(7,7)}
    extra=sorted(neededpairs-b11)
    dump(out,'phase2-linear-syzygy-source-root-contractions.json',{'pairs36_order':mons2,
      'all119_linear_coefficient_original768_combinations':[[encode(v) for v in rows] for rows in c1],
      'all119_two_cubic_source_functionals':folded,'root_kind0_typedB_kind1_C0_kind2_native':True,
      'actual_requested_X2_pairs':sorted(neededpairs),'requested_pairs_outside_B11':extra,
      'known_Xp1_a2_correction_included':True,'native_individual_X2_not_queried':True,
      'new_source_coefficients_or_jets_evaluated':False})
    dump(out,'phase.json',{'phase':'FINAL_SOURCE_PREFLIGHT','needed_pairs':sorted(neededpairs),'extra_pairs_outside_B11':extra})
    note=ROOT/'research/notes'/f'{STEM}.md'
    report={'id':STEM,'all_passed':True,'two_cubic_dual_QQ_support':list(map(len,Pi)),
      'two_source_homotopy_adjoint_QQ_support':list(map(len,kappa)),
      'requested_X2_pairs':sorted(neededpairs),'extra_pairs_outside_B11':extra,
      'typed_B_root_count':len(rootkinds[0]),'original_C0_root_count':len(rootkinds[1]),'native_column_count':len(rootkinds[2]),
      'actual_folded_functional_QQ_support':support,'no_new_source_function_jet_or_rank':True,
      'actual_cubic_relation_values_or_total_T_length_computed':False,
      'input_sha256':hashes,'code_sha256':sha(Path(__file__)),'note_sha256':sha(note),
      'evidence_sha256':{p.name:sha(p) for p in out.iterdir()},'seconds':time.monotonic()-start,
      'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
    dump(out,'report.json',report);print(json.dumps({k:report[k] for k in ['id','all_passed','two_cubic_dual_QQ_support','two_source_homotopy_adjoint_QQ_support','extra_pairs_outside_B11','actual_folded_functional_QQ_support','seconds']}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);main(p.parse_args().output)
