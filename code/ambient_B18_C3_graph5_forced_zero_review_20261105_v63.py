#!/usr/bin/env python3
"""Independent finite homogeneous/source proof; no graph-root readback."""
import hashlib,itertools,json,pathlib,resource,signal,time
from fractions import Fraction as Q
ROOT=pathlib.Path(__file__).resolve().parents[1]
ID='AMBIENT_B18_C3_GRAPH5_FORCED_ZERO_20261105_v63_REVIEW'
PRI='AMBIENT_B18_C3_GRAPH5_FORCED_ZERO_20261105_v63'
PRIMARY='a35b929755230acc1621df59d25d21a55d2617892925a1b2a027d7b5b1451d8b'
def sha(p):
    h=hashlib.sha256()
    with pathlib.Path(p).open('rb') as f:
        for b in iter(lambda:f.read(65536),b''):h.update(b)
    return h.hexdigest()
def read(p):return json.loads(pathlib.Path(p).read_text())
def save(p,x):pathlib.Path(p).write_text(json.dumps(x,indent=2)+'\n')
def field(t,k):
    p=t.index(':',t.index(json.dumps(k)))+1
    while t[p].isspace():p+=1
    return json.JSONDecoder().raw_decode(t,p)[0]
def plus(a,b,c=Q(1)):
    d=dict(a)
    for m,v in b.items():
        d[m]=d.get(m,Q(0))+c*v
        if not d[m]:del d[m]
    return d
def product(a,b):
    d={}
    for m,v in a.items():
        for n,w in b.items():
            k=tuple(i+j for i,j in zip(m,n));d[k]=d.get(k,Q(0))+v*w
    return {m:v for m,v in d.items() if v}
def composition(rows,x,y):
    d={}
    for (a,b),c in rows:
        v={(0,0):Q(c)}
        for _ in range(a):v=product(v,x)
        for _ in range(b):v=product(v,y)
        d=plus(d,v)
    return d
def normal(p,co):
    # Direct monomial rewriting, not the author's closed NF formula.
    p={m:c for m,c in p.items() if m[1]<3}
    while any(m[0]>=3 for m in p):
        m=max(m for m in p if m[0]>=3);c=p.pop(m);r,l=m
        for j in [3,2,1]:
            n=(r-j,l+j)
            if n[1]<3:p=plus(p,{n:-c*co[j]/co[0]})
    return p
def vector(p,basis):return [p.get(m,Q(0)) for m in basis]
def main():
    t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30))
    out=ROOT/'reports'/ID;out.mkdir(exist_ok=False)
    try:
        rp=ROOT/'reports'/PRI/'report.json';assert sha(rp)==PRIMARY;pr=read(rp);assert pr['all_passed']
        pins={str(rp.relative_to(ROOT)):PRIMARY}
        for p,h in pr['input_sha256'].items():assert sha(ROOT/p)==h;pins[p]=h
        for p,h in pr['evidence_sha256'].items():
            q=ROOT/'reports'/PRI/p;assert sha(q)==h;pins[str(q.relative_to(ROOT))]=h
        E=read(ROOT/'reports'/PRI/'actual-cover-NF-and-five-triangular-cubic-quartic-contradictions.json')
        def source(suffix):return next(ROOT/p for p in pr['input_sha256'] if p.endswith(suffix))
        for p in pr['input_sha256']:
            if p.endswith('_REVIEW/report.json'):assert read(ROOT/p)['all_passed']
        rank=read(ROOT/'reports/AMBIENT_B18_COMPARISON_C_GENERIC_DVR_RANK_20261104_v63_REVIEW/report.json')
        assert rank['actual_comparison_C_generic_DVR_rank']==4 and rank['actual_mu_D_C']==8
        C2=read(source('actual-physical12-and-C2-right-exact-quotient.json'))
        assert C2['actual_C2_dimension']==14 and C2['ell_action_rank']==6
        rankC=rank['actual_comparison_C_generic_DVR_rank'];mu=14-6;tors=mu-rankC
        exp1=2*tors-(14-2*rankC);ge2=tors-exp1
        lower=3*rankC+exp1+2*ge2;assert (mu,tors,exp1,ge2,lower)==(8,4,2,2,18)
        pt=source('actual-special-fiber20-coefficients.json').read_text()
        arows=field(pt,'gA_plane3');brows=field(pt,'gB_own_plane3');yb=field(pt,'YB_in_YA_linear');del pt
        r={(1,0):Q(1)};s={(1,0):Q(1),(0,1):Q(1,5)}
        rb=plus({m:Q(yb[7][7])*c for m,c in r.items()},{m:Q(yb[7][8])*c for m,c in s.items()})
        sb=plus({m:Q(yb[8][7])*c for m,c in r.items()},{m:Q(yb[8][8])*c for m,c in s.items()})
        polys=[composition(arows,r,s),composition(brows,rb,sb)]
        curves=[[p.get((3-i,i),Q(0)) for i in range(4)] for p in polys]
        assert curves[0][0]==Q(125,24) and curves[1][0]==Q(-15125,144)
        nf3=[];nf4=[];block=[]
        for k,co in enumerate(curves):
            assert co==[Q(c) for c in E['two_block_source_curve_NF_receipts'][k]['curve3_r3_r2ell_rell2_ell3']]
            n3=[vector(normal({(3-i,i):Q(1)},co),[(2,1),(1,2)]) for i in range(4)]
            n4=[vector(normal({(4-i,i):Q(1)},co),[(2,2)])[0] for i in range(5)]
            q3=[[n3[j][i] for j in range(4)] for i in range(2)]
            assert [[str(c) for c in v] for v in q3]==E['two_block_source_curve_NF_receipts'][k]['NF3_to_r2ell_rell2']
            assert list(map(str,n4))==E['two_block_source_curve_NF_receipts'][k]['NF4_to_r2ell2']
            assert not normal(polys[k],co)
            assert not normal(product(polys[k],r),co) and not normal(product(polys[k],{(0,1):Q(1)}),co)
            assert normal({(1,2):Q(1)},co)=={(1,2):Q(1)}
            assert normal({(2,2):Q(1)},co)=={(2,2):Q(1)}
            # Completeness is witnessed by unit triangular relation blocks, not only annihilation.
            block.append({'degree3_relation_minor':str(co[0]),'degree4_relation_minor':str(co[0]**2),
              'degree3_free2_and_degree4_free1':True})
            nf3.append(q3);nf4.append(n4)
        N3=[[Q(0)]*12 for _ in range(4)];N4=[[Q(0)]*15 for _ in range(2)]
        for i in range(2):
            for j in range(2):N3[2*i+j][4*i:4*i+4]=nf3[i][j]
            N4[i][5*i:5*i+5]=nf4[i]
        assert [[str(c) for c in v] for v in N3]==E['whole_NF3_rows4_columns12']
        assert [[str(c) for c in v] for v in N4]==E['whole_NF4_rows2_columns15']
        raw=read(source('actual-ell2-product33-and-raw45-actions.json'));basis=raw['actual33_basis']
        frame=read(source('actual-first21-quadratic-frame-and-twelve-cubic-zero-sensitivity-proof.json'))
        fm={x['physical_row']:x for x in frame['all12_verified_source_records']}
        tauA=read(source('nine-complete-tauA-source-kernel-and-unpaired-cubic-cache-adjoints.json'))
        K=[[Q(c) for c in row] for row in tauA['complete_constant_kernel12_by9']]
        rows=[63,64,68,73,74];cr={x['physical_row']:x for x in C2['physical12_coefficients']}
        mixed_idx=[basis.index([i+1,a,b]) for a,b in [(1,0),(0,1)] for i in range(7)]
        premises=[]
        for row in rows:
            j=tauA['physical_row_order'].index(row);kc=tauA['free_columns9'].index(j)
            assert all(K[i][kc]==int(i==j) for i in range(12))
            assert all(not Q(c) for c in fm[row]['bB_second_plane']+fm[row]['bA_first_z7'])
            de=dict(cr[row]['delta']);ba=dict(cr[row]['bA'])
            assert all(not Q(de.get(i,'0')) for i in mixed_idx)
            assert all(not Q(de.get(basis.index(m),'0')) for m in [[0,3,0],[0,2,1]])
            mx=[[basis[i],ba[i]] for i in mixed_idx if i in ba and Q(ba[i])]
            assert mx and mx==next(x for x in E['same_five_actual_source_premises'] if x['physical_row']==row)['actual_bA2_H_mixed_terms_retained']
            premises.append({'row':row,'same_free_tauA_source_ej':True,'unknown_zz_terms_not_dropped':True,
              'actual_mixed_bA2_H_terms_retained':mx,'tauB_H_mixed_carry':'r/ell*RH has eB plane2=delta_i2_plane=0',
              'tauA_eA_carry':'Only exact RA clears A, so eA cubic=delta3; all B/H carry remains opaque.'})
        # Fresh universal minor calculation over QQ[d,u,v], no value assigned to d.
        d={(1,0,0):Q(1)};u={(0,1,0):Q(1)};v={(0,0,1):Q(1)};z={}
        p=[z,d,u,v];q=[z,z,z,d]
        cubic=plus(product(p[1],q[3]),product(q[1],p[3]),-1)
        aB=curves[1][1]/curves[1][0]
        rp0=d;rp1=plus(v,u,-aB);rq0=z;rq1=d
        quartic=plus(product(rp0,rq1),product(rq0,rp1),-1)
        assert cubic==quartic=={(2,0,0):Q(1)}
        coverHF=[3,6,6,4,2];assert sum(coverHF)==21 and 21-2-2==17<lower
        proof={'all_exact_source_prior_peer_and_primary_hashes_match':True,
          'actual_module':'C=coker(D0->E4) over SAME R0, C3=C/ell3C; not an algebra/direct splitting claim',
          'DVR_lower_derivation':{'rank':rankC,'C2_dimension':14,'C2_ell_rank':6,'C_mod_ell_dimension':mu,
            'torsion_summands':tors,'exactly_exponent1':exp1,'atleast_exponent2':ge2,'C3_lower':lower},
          'both_actual_curve_compositions_and_all_NF3_NF4_entries_independently_rewritten':True,
          'relation_block_completeness':block,'whole_NF3_4x12_and_NF4_2x15_match':True,
          'known_cover_HF':coverHF,'known_cover_length':21,'cover_only_not_actual_module_split':True,
          'all5_actual_same_source_premises':premises,
          'same_delta3_coefficient_on_diagonal':'Existing raw33 literal Taylor duals give delta_r3=delta_r2ell=0. Modell3, both A and B source components are d*r*ell2.',
          'all_opaque_B_carry_retained_formal_u_v':True,
          'actual_cubic_classes_mod_known_cover':['(0,d,u,v)','(0,0,0,d)'],
          'actual_r_multiple_quartics':['(d,-(gB_r2ell/gB_r3)*u+v)','(0,d)'],
          'fresh_universal_cubic_and_quartic_fixed_minors':'d^2 in QQ[d,u,v]',
          'source_and_filtered_dimension_argument':'If d!=0, actual initial submodule contains2 independent degree3 classes and2 independent r-multiples in degree4 of the length21 cover. Hence actual length<=17, contradicting actualDVR lower18.',
          'all5_actual_delta_r_ell2_zero_forced':True,
          'all5_graph_tauB_cubic_classes_zero':True,
          'ell3_coefficient_remains_unknown_not_zero_filled':True,
          'no_graph5_scalar_readbacks_or_IFT_gauge_jet_or_rank_search':True,
          'no_graph5_V1_or_V2_execution_needed':True,
          'no_actual_C3_length_or_all_cubic_or_torsion_exponent_determination':True,
          'static_numeric_envelope':50000,'hard_numeric':150000,'seconds':60,'address_bytes':2<<30,'hard_RAW_bytes':32<<20}
        save(out/'independent-triangular-source-and-DVR-dimension-proof.json',proof)
        report={'id':ID,'all_passed':True,'frozen':True,'source_only':True,'primary_report_sha256':PRIMARY,
          'input_sha256':pins,'code_sha256':sha(__file__),
          'evidence_sha256':{'independent-triangular-source-and-DVR-dimension-proof.json':sha(out/'independent-triangular-source-and-DVR-dimension-proof.json')},
          'all5_graph_coefficients_and_cubic_classes_zero_without_readback_accepted':True,
          'full_mixed_H_carry_and_same_tauA_tauB_source_diagonal_accepted':True,
          'no_actual_C3_length_or_torsion_exponent_claim':True,'no_author_producer_imports':True,'elapsed_seconds':time.monotonic()-t}
        save(out/'report.json',report);signal.alarm(0)
        print(json.dumps({'all_passed':True,'report':str(out/'report.json'),'sha256':sha(out/'report.json'),'seconds':report['elapsed_seconds']}))
    except Exception as e:
        save(out/'failure.json',{'error':type(e).__name__+': '+str(e),'no_positive_claim':True});raise
if __name__=='__main__':main()
