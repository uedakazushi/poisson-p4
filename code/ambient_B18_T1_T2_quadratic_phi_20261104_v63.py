#!/usr/bin/env python3
"""ONE cached-source2 contraction and exact quadratic functional on Q1/Q2.

No new source jet/IFT. Unknown old regular second derivatives and g first
derivatives are formal until their coefficients or directional gate vanish.
"""
import argparse
from fractions import Fraction as F
import hashlib,json
from pathlib import Path
import resource,signal,time
from ambient_B18_top_graph_telescoping_20261104_v63 import independent_rows,linear_solve
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_T1_T2_QUADRATIC_PHI_20261104_v63'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def enc(M):return [[str(v) for v in r] for r in M]
def main(out):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));out.mkdir(parents=True,exist_ok=False)
    sp=ROOT/'reports/AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63/actual-source-coefficients.json'
    fp=ROOT/'reports/AMBIENT_B18_FULL_RELATIVE_INPUT_20261104_v63/full-relative186-DAG.json'
    hp=ROOT/'reports/AMBIENT_B18_T3_HONEST10_PRESENTATION_20261104_v63/honest10-presentation-DAG.json'
    up=ROOT/'reports/AMBIENT_B18_T3_HONEST10_UNIT_REDUCTION_20261104_v63/unit-reduced-presentation-DAG.json'
    jp=ROOT/'reports/AMBIENT_B18_T3_FULL9_FIRSTJET_20261104_v63/full9-firstjet-certificate.json'
    qp=ROOT/'reports/AMBIENT_B18_T1_T2_FIRSTJET_20261104_v63/Q1-Q2-input-and-firstjets.json'
    bp=ROOT/'reports/AMBIENT_B18_HONEST9_TO_N8_20261104_v63/honest9-to-N8-DAG.json'
    zp=ROOT/'reports/AMBIENT_B18_T3_ZERO_STRUCTURAL_20261104_v63/report.json'
    old,full,honest,unit,jet,qdata,BN=map(read,[sp,fp,hp,up,jp,qp,bp]);assert read(zp)['all_passed']
    oldnodes=old['DAG_scalar_nodes'];nodes=unit['DAG_scalar_nodes'];N=len(oldnodes);mon=list(map(tuple,jet['source2_monomials']))
    var={n[1]:i for i,n in enumerate(oldnodes) if n[0]=='v'};xids=[var[f'x{i}'] for i in old['physical_zero_base_rows9']];xi=old['physical_zero_base_rows9'].index(61)
    bindings={};line={};YA1=None
    for name in ['CA','CB']:
        src=jet['both57_source2'][name];sol={i:list(map(F,v)) for i,v in src['all57_source2jets']}
        L=[[sol[i][1+j] for j in range(9)] for i in xids]
        y,_=linear_solve(L,[[F(int(i==xi))] for i in range(9)]);y1=[r[0] for r in y]
        def quadratic(v):return sum(v[k]*y1[a]*y1[b] for k,(a,b) in enumerate(mon[10:],10))
        y,_=linear_solve(L,[[-quadratic(sol[i])] for i in xids]);y2=[r[0] for r in y]
        b={i:[F(0),sum(v[1+j]*y1[j] for j in range(9)),sum(v[1+j]*y2[j] for j in range(9))+quadratic(v)] for i,v in sol.items()}
        assert all(b[i]==[F(0),F(int(j==xi)),F(0)] for j,i in enumerate(xids))
        for i,v in b.items():
            if i in bindings:assert bindings[i]==v
            bindings[i]=v
        line[name]={'physical_x61_line_Y_first':list(map(str,y1)),'physical_x61_line_Y_second':list(map(str,y2)),
                    'all57_source_coefficient_line_jets':[[i,list(map(str,v))] for i,v in b.items()]}
        if name=='CA':YA1=y1
    original=[]
    for n in oldnodes:
        o=n[0]
        if o=='c':v=[F(n[1],n[2]),F(0),F(0)]
        elif o=='v':v=bindings[var[n[1]]]
        elif o=='a':v=[a+b for a,b in zip(original[n[1]],original[n[2]])]
        elif o=='n':v=[-a for a in original[n[1]]]
        elif o=='m':
            a,b=original[n[1]],original[n[2]];v=[a[0]*b[0],a[0]*b[1]+a[1]*b[0],a[0]*b[2]+a[1]*b[1]+a[2]*b[0]]
        elif o=='i':
            a=original[n[1]];assert a[0];v=[1/a[0],-a[1]/a[0]**2,a[1]**2/a[0]**3-a[2]/a[0]**2]
        else:raise ValueError(o)
        original.append(v)
    # Weighted opaque monomials: t=0, gamma1..4 weight1, theta coefficients weight2.
    gnodes=[q['node'] for q in full['all4_known_regular_g_divisions']];rnodes=[q['node'] for q in full['all_coordinate_divisions']]
    theta={i:5+j for j,i in enumerate(rnodes+gnodes)};weight={0:1,**{i:1 for i in range(1,5)},**{j:2 for j in theta.values()}}
    def degree(m):return sum(weight[i] for i in m)
    def add(a,b):
        c=dict(a)
        for m,v in b.items():
            c[m]=c.get(m,F(0))+v
            if not c[m]:c.pop(m)
        return c
    def scale(a,c):return {m:v*c for m,v in a.items() if v*c}
    def mul(a,b,limit=2):
        c={}
        for m,v in a.items():
            for n,w in b.items():
                p=tuple(sorted(m+n))
                if degree(p)<=limit:
                    c[p]=c.get(p,F(0))+v*w
                    if not c[p]:c.pop(p)
        return c
    def inverse(a,limit):
        c=a.get((),F(0));assert c;h=add(a,{():-c});v=add({():1/c},scale(h,-1/c**2))
        return add(v,scale(mul(h,h),1/c**3)) if limit==2 else v
    def base(v,limit):return {m:a for m,a in zip([(),(0,),(0,0)],v) if a and len(m)<=limit}
    def hline(v):return [F(v[0]),sum(F(v[1+j])*YA1[j] for j in range(9))]
    known2={};known1={}
    for i in rnodes:
        h=hline(jet['all784_CB_h_firstjets_in_YA'][str(i)]);p=base(h,1);known1[i]=p;known2[i]=add(p,{(theta[i],):F(1)})
    for i in gnodes:
        j=gnodes.index(i);known1[i]={(1+j,):F(1)};known2[i]={(1+j,):F(1),(theta[i],):F(1)}
    for key,v in jet['all792_CA_h_firstjets'].items():known1[int(key)]=base(hline(v),1)
    gtail=full['CB_g_restriction_node'];known1[gtail]=known2[gtail]={}
    def evaluate(roots,known,limit):
        needed=set();stack=roots[:]
        while stack:
            i=stack.pop()
            if i in needed:continue
            needed.add(i)
            if i<N or i in known:continue
            n=nodes[i]
            if n[0] in ['a','m']:stack.extend(n[1:])
            elif n[0] in ['n','i']:stack.append(n[1])
            elif n[0]!='c':raise AssertionError(('uncontrolled',i,n))
        V={}
        for i in sorted(needed):
            if i<N:v=base(original[i],limit)
            elif i in known:v=known[i]
            else:
                n=nodes[i];o=n[0]
                if o=='c':v={():F(n[1],n[2])} if n[1] else{}
                elif o=='a':v=add(V[n[1]],V[n[2]])
                elif o=='n':v=scale(V[n[1]],-1)
                elif o=='m':v=mul(V[n[1]],V[n[2]],limit)
                elif o=='i':v=inverse(V[n[1]],limit)
                else:raise ValueError(o)
            V[i]=v
        return V
    active=unit['remaining_D3_relation_original_indices'];positions=sorted({i for family in qdata['families'].values() for i in family['honest9_column_source_positions_in_Q3_reduced']})
    oldidx={i:active[i] for i in positions}
    raw={i:(honest['all36_Koszul_kernel_columns_R0_10'][k][0] if k<36 else full['actual_Qn_R0_squared_generator_matrices']['3'][k-36][0]) for i,k in oldidx.items()}
    c5={i:(honest['all36_Koszul_kernel_columns_R0_10'][k][5] if k<36 else honest['all2034_Q3_columns_lifted_R0_10'][k-36][5]) for i,k in oldidx.items()}
    step=unit['all_constant_unit_pivot_steps'][0];pe=step['relation_before_projection'][0];u=step['pivot_node']
    E2=evaluate(list(raw.values())+[pe],known2,2);E1=evaluate(list(c5.values())+[u],known1,1)
    assert E2[pe].get((),F(0))==0;u0=E1[u].get((),F(0));assert u0==F(10,3)
    p1={m:v for m,v in E2[pe].items() if degree(m)==1};p2={m:v for m,v in E2[pe].items() if degree(m)==2};u1={m:v for m,v in E1[u].items() if degree(m)==1}
    cfirst=[F(0)]+[F(t['all8_coordinate_divisions'][0]['origin']) for t in BN['all64_regular_coordinate_divisions']]
    assert len(cfirst)==9 and [F(c[0]) for c in BN['B_N_origin_matrix_columns']]==[F(1)]+[F(0)]*8
    phis={};records=[]
    for i in positions:
        c0=E1[c5[i]].get((),F(0));c1={m:v for m,v in E1[c5[i]].items() if degree(m)==1}
        a2={m:v for m,v in E2[raw[i]].items() if degree(m)==2}
        reduced=add(add(a2,scale(p2,-c0/u0)),add(scale(mul(p1,c1),-1/u0),scale(mul(p1,u1),c0/u0**2)))
        first=list(map(F,jet['all2069_known_firstjet90_columns'][i]));vel=[sum(first[10*r+1+j]*YA1[j] for j in range(9)) for r in range(9)]
        assert vel[0]==0
        correction=sum(cfirst[j]*vel[j] for j in range(1,9));phi=add(reduced,{(0,0):correction} if correction else{})
        unknown={m:v for m,v in phi.items() if m!=(0,0)}
        exact=not unknown;phis[i]=phi.get((0,0),F(0)) if exact else None
        records.append({'source_position_in_Q3_reduced':i,'old_original_relation_index':oldidx[i],
                        'raw_E0_node':raw[i],'old_component5_node':c5[i],
                        'reduced_old_E0_quadratic_formal':[[list(m),str(v)] for m,v in sorted(reduced.items())],
                        'all9_reduced_first_line_coefficients':list(map(str,vel)),
                        'BN_firstrow_correction':str(correction),'phi_formal':[[list(m),str(v)] for m,v in sorted(phi.items())],
                        'all_unknown_first_g_and_second_regular_sensitivities_zero':exact,'exact_phi':str(phis[i]) if exact else None})
    results={}
    for n,fam in qdata['families'].items():
        ids=fam['honest9_column_source_positions_in_Q3_reduced'];exact=all(phis[i] is not None for i in ids)
        if exact:
            A=[[F(v) for v in row]+[phis[i]] for i,row in zip(ids,fam['all_exact_firstjet90_columns'])]
            piv=independent_rows(A);rank=len(piv);assert rank in [43,44]
            basis=[A[i] for i in piv];coords=independent_rows([list(r) for r in zip(*basis)])
            _,det=linear_solve([[row[j] for j in coords] for row in basis],[[F(int(i==j)) for j in range(rank)] for i in range(rank)])
        else:A=[];piv=[];rank=None;coords=[];det=None
        results[n]={'all_phi_exact':exact,'all_phi':[str(phis[i]) if phis[i] is not None else None for i in ids],
                    'all_augmented91_columns':enc(A),'augmented_rank':rank,'independent_column_indices':piv,
                    'minor_coordinate_indices':coords,'minor_determinant':str(det) if det else None,
                    'module_result':'ZERO' if rank==44 else 'CYCLIC_NONZERO' if rank==43 else 'UNDETERMINED'}
    data={'physical_line':'all8actualCohenphysicalcoords0,oldrho61=v EXACTLY','both_family_exact_second_coordinate_contractions':line,
          'old_regular_second_derivative_formal_ids':[[i,j] for i,j in theta.items()],'four_g_first_formal_ids':[[i,1+j] for j,i in enumerate(gnodes)],
          'BN_firstrow_linear_v_coefficients_from_FULL_composed_H_image':list(map(str,cfirst)),
          'Schur_unit_origin':str(u0),'pivot_even_first_formal':[[list(m),str(v)] for m,v in p1.items()],
          'all_relation_phi_records':records,'families':results,
          'no_unknown_derivative_assigned_zero':True,'no_nil_kappa_or_inverse_coordinate_jet_needed':True}
    encoded=json.dumps(data,separators=(',',':'))+'\n';assert len(encoded.encode())<32*1024*1024
    (out/'quadratic-phi-certificate.json').write_text(encoded)
    paths=[Path(__file__),sp,fp,hp,up,jp,qp,bp,zp]
    rep={'id':STEM,'all_passed':True,'status':'cached_source2_exact_quadratic_functional',
         'Q1_all_phi_exact':results['1']['all_phi_exact'],'Q2_all_phi_exact':results['2']['all_phi_exact'],
         'Q1_augmented_rank':results['1']['augmented_rank'],'Q2_augmented_rank':results['2']['augmented_rank'],
         'T1_result_under_accepted_minimal44_theorem':results['1']['module_result'],
         'T2_result_under_accepted_minimal44_theorem':results['2']['module_result'],
         'full_B18_depth_or_T0_not_solved':True,'no_new_source_IFT_jets_Honest2_matrix_GB_or_points':True,
         'input_sha256':{str(p.relative_to(ROOT)):sha(p) for p in paths},
         'evidence_sha256':{'quadratic-phi-certificate.json':sha(out/'quadratic-phi-certificate.json')},
         'certificate_bytes':len(encoded.encode()),'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    (out/'report.json').write_text(json.dumps(rep,indent=2)+'\n');print(json.dumps({k:rep[k] for k in ['all_passed','Q1_all_phi_exact','Q2_all_phi_exact','Q1_augmented_rank','Q2_augmented_rank','T1_result_under_accepted_minimal44_theorem','T2_result_under_accepted_minimal44_theorem','seconds','max_rss_kib']}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
