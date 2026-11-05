#!/usr/bin/env python3
"""Exact pure-r³ source functional for619 constant quadratic-kernel relations."""
import argparse,itertools,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
from ambient_B18_comparison_T_constant_rank_20261104_v63 import read,sha,rec,encode,addscaled
from ambient_B18_comparison_T_next3_source_preflight_20261104_v63 import receipt_combos
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_T_CONSTANT_QUADRATIC_KERNEL_R3_ADJOINT_20261104_v63'
def dump(out,name,x):(out/name).write_text(json.dumps(x,separators=(',',':'))+'\n')
def put(v,k,a):
    c=v.get(k,F(0))+a
    if c:v[k]=c
    else:v.pop(k,None)
def enc(v):return [[list(k) if isinstance(k,tuple) else k,str(a)] for k,a in sorted(v.items()) if a]
def main(out):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out=out.resolve();out.mkdir(parents=True,exist_ok=False)
    paths={
      'boundary':'reports/AMBIENT_B18_T_CUBIC_BOUNDARY_SOURCE3_PREFLIGHT_20261104_v63/phase1-corrected-d3-pure-r3-channel.json',
      'gates':'reports/AMBIENT_B18_T_CUBIC_BOUNDARY_SOURCE3_PREFLIGHT_20261104_v63/phase2-known-Schur-second-r-gates.json',
      'low':'reports/AMBIENT_B18_COMPARISON_T_LINEAR_QUADRATIC_20261104_v63/actual-linear-quadratic-T-certificate.json',
      'linear':'reports/AMBIENT_B18_COMPARISON_T_LINEAR_QUADRATIC_20261104_v63/phase1-actual-five-target-linear-relations.json',
      'mu':'reports/AMBIENT_B18_COMPARISON_T_CONSTANT_RANK_20261104_v63/actual-linear-cycle-muT-certificate.json',
      'pi':'reports/AMBIENT_B18_T_QUADRATIC_M1_CANCELLATION_PREFLIGHT_20261104_v63/actual-projected-three-class-X2-operator-preflight.json',
      'm0':'reports/AMBIENT_B18_FULLV56_COMPARISON_UNIT_20261104_v63/actual-fullV56-m0-and-typed-comparison.json',
      'source620':'reports/AMBIENT_B18_COMPARISON_T_NEXT3_SOURCE_PREFLIGHT_20261104_v63/zero-linear620-source-combinations.json',
      'rank1':'reports/AMBIENT_B18_T_ACTUAL_PROJECTED_QUADRATIC_RANK_20261104_v63/actual-three-by620-quadratic-image-certificate.json',
      'chain':'reports/AMBIENT_B18_NORMALIZER_CHAIN_LIFT_640_NEGATIVE_BOUNDARY_20261104_v63/actual-UP-QB-640-typed-operator.json',
      'C0':'reports/AMBIENT_B18_T0_HONEST_C0_20261104_v63/honest-C0-D0-D4-DAG.json'}
    hashes={p:sha(ROOT/p) for p in paths.values()};dump(out,'phase.json',{'phase':'sources_bound_before619_adjoints','input_sha256':hashes})
    D={k:read(ROOT/p) for k,p in paths.items()};B,G,L,L1,M,P,V,S,R,C,H=[D[k] for k in paths]
    assert B['all119_linear_syzygies_annihilated_so_theta_is_well_defined']
    assert all(F(r['corrected_boundary_pure_r3'])==0 for r in B['all168_corrected_boundaries'])
    basis=list(map(rec,L['quadratic_certificate']['basis_normalized_columns']));lead={min(v):j for j,v in enumerate(basis)}
    theta_basis=list(map(F,B['theta_on_normalized177_quadratic_basis']));theta={}
    for i in range(180):
        v={i:F(1)};a=F(0)
        for k in sorted(lead):
            if k in v:
                c=v[k];j=lead[k];a+=c*theta_basis[j];addscaled(v,basis[j],-c)
        if a:theta[i]=a
    assert all(sum((theta.get(i,F(0))*c for i,c in b.items()),F(0))==a for b,a in zip(basis,theta_basis))
    mons=list(itertools.combinations_with_replacement(range(8),2));mi={m:i for i,m in enumerate(mons)}
    star_rr=[F(r['known_rr']) for r in G['normalized_Ystar84_source55_rr']];assert len(star_rr)==84
    qnorm=list(map(rec,L1['normalized84_unitcycle_Q1_columns_sparse']))
    nbeta={}
    for p,q in enumerate(qnorm):
        for v in range(8):
            a=-sum((c*theta.get(36*(i//8)+mi[tuple(sorted((i%8,v)))],F(0)) for i,c in q.items()),F(0))
            if v==7:a+=star_rr[p]
            if a:nbeta[8*p+v]=a
    wb=list(map(rec,M['complete_rank_certificate']['basis_normalized_columns']));wl={min(v):j for j,v in enumerate(wb)};wc=[]
    for i in range(84):
        v={i:F(1)};b={}
        while v:
            k=min(v);j=wl[k];a=v[k];b[j]=b.get(j,F(0))+a;addscaled(v,wb[j],-a)
        wc.append(b)
    pirows=[{P['free_kernel_rows32'][j]:c for j,c in rec(p).items()} for p in P['pi_free32_rows']]
    adj={}
    for i,c in theta.items():
        t,m=divmod(i,36)
        for src,a in pirows[t].items():put(adj,36*src+m,c*a)
    m0=[list(map(F,row)) for row in V['all24_by56_m0']];aa={}
    for w,(sector,pair) in enumerate(M['literal_wedge_basis84']):
        i,j=pair
        for v in range(8):
            c=sum((a*nbeta.get(8*p+v,F(0)) for p,a in wc[w].items()),F(0))
            if not c:continue
            terms=[(8*sector+7,i,1)] if j==7 else [(8*sector+j,i,1),(8*sector+i,j,-1)]
            for row,k,sign in terms:
                mon=tuple(sorted((v,k)));zd=sum(t<7 for t in mon);den=zd if j==7 else zd+1
                a=c*F(sign*mon.count(k),den);put(aa,36*row+mi[mon],a)
    for i,c in aa.items():
        row,m=divmod(i,36)
        for src,a in enumerate(m0[row]):
            if a:put(adj,36*src+m,c*a)
    dump(out,'phase1-theta-and-pure-r3-second-adjoint.json',{'theta_extension180':encode(theta),
      'theta_plus_Ystar2_beta_quadratic_source_adjoint56_times36':encode(adj),
      'normalized_wedge_linear_adjoint':encode(nbeta),
      'formula':'rawP55,3 - Xp55,2rr*a1r - [theta(Qproxy2)+Ystar55,2rr*beta1r]',
      'Theta_source_argument':'Xc2-Xp1*a1','all_d3_quadratic_ambiguities_killed_by_theta':True})
    norm,kernel,labels=receipt_combos(R['complete_fixed_rank_certificate'],620);assert len(norm)==1 and len(kernel)==619
    bcols=[dict(v) for v in C['B48_by640_sparse_columns']];ops=C['operator_scalar_nodes'];xp1=list(map(rec,P['actual_pivot_X1_by27_sparse56_times8']))
    source620=S['all620'];all620=[]
    for row in source620:
        c=rec(row['combination_in_actual768']);a1=rec(row['linear_Xp27_correction_row8j_plus_variable']);p3={};xp2={};x2={};known=F(0)
        for j,a in c.items():
            if j<640:put(p3,H['honest_D0_relation_columns'][j][8],a)
        for p in range(27):
            a=a1.get(8*p+7,F(0))
            if a:
                for root,b in G['normalized_Xp27_source55_rr'][p]['rawP8_root_coefficients']:put(xp2,root,-a*F(b))
        for index,a in adj.items():
            src,m=divmod(index,36);u,v=mons[m]
            for p,first in enumerate(xp1):
                known+=a*a1.get(8*p+u,F(0))*first.get(8*src+v,F(0))
                if u!=v:known+=a*a1.get(8*p+v,F(0))*first.get(8*src+u,F(0))
            for j,b in c.items():
                if j>=640:
                    if src<48:put(x2,(2,j-640,src,m),-a*b)
                elif src<48:
                    root=bcols[j].get(src)
                    if root is None:continue
                    n=ops[root]
                    if n[0]=='c':continue
                    key=(1,n[1],m) if n[0]=='actual_C0_scalar_function' else (0,root,m);put(x2,key,-a*b)
                else:
                    root=H['honest_D0_relation_columns'][j][src-48+1]
                    if H['DAG_scalar_nodes'][root][0]!='c':put(x2,(1,root,m),-a*b)
        # Literal source constants have no third or second coefficient.
        p3={i:c for i,c in p3.items() if H['DAG_scalar_nodes'][i][0]!='c'}
        xp2={i:c for i,c in xp2.items() if H['DAG_scalar_nodes'][i][0]!='c'}
        all620.append((p3,xp2,x2,known))
    output=[];support=0;pairs=set();p3roots=set();xp2roots=set()
    for i,c in enumerate(kernel):
        p3={};xp2={};x2={};known=F(0)
        for j,a in c.items():
            q,r,s,t=all620[j];addscaled(p3,q,a);addscaled(xp2,r,a);addscaled(x2,s,a);known+=a*t
        support+=len(c)+len(p3)+len(xp2)+len(x2)+bool(known);assert support<150000
        pairs.update(mons[k[-1]] for k in x2);p3roots.update(p3);xp2roots.update(xp2)
        output.append({'constant_quadratic_kernel_index619':i,'fixed_zero_receipt_source620':labels[i],
          'combination_in620':encode(c),'rawP8_last_cubic_roots':encode(p3),
          'Xp_second_r2_roots':encode(xp2),'quadratic_source_root_correction':enc(x2),
          'known_Xp1_a1_contraction':str(known)})
    dump(out,'phase2-all619-pure-cubic-source-functionals.json',{'pairs36_order':mons,'all619':output,
      'actual619_is_kernel_of_new_quadratic_rank1_image':True,'actual_required_X2_pairs':sorted(pairs),
      'raw_P3_root_count':len(p3roots),'Xp_unknown_second_root_ids':sorted(xp2roots),
      'exact_four_Xp_roots_are_old_minus_new_g':[[i,H['DAG_scalar_nodes'][i]] for i in sorted(xp2roots)],
      'all_ordinary_division_second_and_third_values_still_symbolic':True,
      'no_new_source3_or5_jet_or_value':True,'expanded_saved_QQ_support':support})
    note=ROOT/'research/notes'/f'{STEM}.md';dump(out,'phase.json',{'phase':'FINAL_SOURCE_ADJOINT','pairs':sorted(pairs),'source_support':support})
    report={'id':STEM,'all_passed':True,'actual_constant_quadratic_kernel_dimension':619,
      'theta180_nonzeros':len(theta),'second_source_adjoint_nonzeros':len(adj),'required_X2_pairs':sorted(pairs),
      'raw_P3_root_count':len(p3roots),'symbolic_Xp2_root_count':len(xp2roots),'expanded_QQ_support':support,
      'no_source3_or5_jet_or_actual619_values_computed':True,'input_sha256':hashes,
      'code_sha256':sha(Path(__file__)),'note_sha256':sha(note),
      'evidence_sha256':{p.name:sha(p) for p in out.iterdir()},'seconds':time.monotonic()-start,
      'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
    dump(out,'report.json',report);print(json.dumps({k:report[k] for k in ['id','all_passed','theta180_nonzeros','second_source_adjoint_nonzeros','required_X2_pairs','raw_P3_root_count','expanded_QQ_support','seconds']}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);main(p.parse_args().output)
