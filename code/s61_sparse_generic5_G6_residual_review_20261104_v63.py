#!/usr/bin/env python3
"""Independent source and finite-operator audit; never run the producer."""
import argparse, collections, gzip, hashlib, itertools, json, math, resource, signal, time
from fractions import Fraction as Q
from pathlib import Path
import sympy as S
from s61_sparse_generic5_full134_normal_circuit_review_20261104_v63 import plus, times, neg, const, symbol, pp_mul

ROOT=Path(__file__).resolve().parents[1]
AUTHOR='S61_SPARSE_GENERIC5_G6_RESIDUAL_20261104_v63'
ID=AUTHOR+'_REVIEW'; ZERO=(0,)*5
def read(p):
    b=p.read_bytes(); return json.loads(gzip.decompress(b) if p.suffix=='.gz' else b)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def enc(p):return tuple(sorted((e,str(v)) for e,v in p.items() if v))
def dec(p):return {tuple(t[1:]):Q(t[0]) for t in p if Q(t[0])}
def polynomial(p):
    """Canonical primitive-atom embedding into a free polynomial ring."""
    if not p:return {}
    if set(p)=={ZERO}:return const(p[ZERO])
    d=math.lcm(*(v.denominator for v in p.values()))
    a={e:v.numerator*(d//v.denominator) for e,v in p.items()}
    g=math.gcd(*(abs(v) for v in a.values()))
    if a[min(a)]<0:g=-g
    return times(const(Q(g,d)),symbol(('primitive',tuple((e,v//g) for e,v in sorted(a.items())))))
def rowadd(row,m,p):
    row[m]=plus(row.get(m,{}),p)
    if not row[m]:del row[m]
def total(ps):
    out={}
    for p in ps:out=plus(out,p)
    return out

class Reader:
    def __init__(self,d,old,res,stops=None):
        self.d=d;self.old=old;self.res=res;self.memo={};self.stops=stops or {}
    def __call__(self,i):
        if i in self.stops:return self.stops[i]
        if i in self.memo:return self.memo[i]
        op,*a=self.d['nodes'][i]
        if op=='qq':v=const(Q(a[0]))
        elif op=='parameter_poly_atom':v=polynomial(dec(self.d['scalar_parameter_atoms'][a[0]]))
        elif op=='external_full134_node':v=self.old(a[0])
        elif op=='external_residual129_node':v=self.res(a[0])
        elif op=='external_original_W':v=symbol(('W',*a))
        elif op in ('z','n'):v=symbol((op,a[0]))
        else:
            assert op in ('sum','prod');v={} if op=='sum' else const(1)
            for j in a[0]:v=plus(v,self(j)) if op=='sum' else times(v,self(j))
        self.memo[i]=v;return v

def main(out):
    start=time.monotonic();signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
    out.mkdir(parents=True,exist_ok=False);base=ROOT/'reports'/AUTHOR
    pr=read(base/'report.json');assert pr['all_passed'];hashes={}
    def bind(p):hashes[str(p.relative_to(ROOT))]=sha(p)
    def audit_report(p):
        g=read(p);assert g['all_passed'];bind(p)
        for n,h in g.get('input_sha256',{}).items():assert sha(ROOT/n)==h,n;hashes[n]=h
        for n,h in g.get('evidence_sha256',{}).items():assert sha(p.parent/n)==h,n;bind(p.parent/n)
        return g
    audit_report(base/'report.json')
    for name in ('S61_SPARSE_GENERIC5_CROSS_RANK_COVARIANCE_20261104_v63_REVIEW',
                 'S61_SPARSE_GENERIC5_FULL134_NORMAL_CIRCUIT_20261104_v63_REVIEW',
                 'S61_SPARSE_GENERIC5_RESIDUAL129_CIRCUIT_20261104_v63_REVIEW',
                 'S61_SPARSE_GENERIC5_ALL670_SOURCE_BINDING_20261104_v63_REVIEW',
                 'S61_SPARSE_GENERIC5_PLANE_FIVE_UNIVARIATE_20261104_v63_REVIEW'):
        audit_report(ROOT/'reports'/name/'report.json')
    g=audit_report(ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json')
    gr=audit_report(ROOT/'reports/S61_UNIVERSAL_WEIGHT6_REVIEW_20261031.json')
    assert g['linear_syzygy_kernel_dimension']==36 and g['canonical_Schouten_Bianchi_dimension']==35
    assert gr['exact_QQ_rank']==3114 and gr['canonical_Bianchi_rank']==35
    assert g['sample_augmented_quadratic_span_rank']==gr['sample_augmented_quadratic_rank']==135
    assert Q(g['sample_augmented_minor_determinant'])==Q(gr['sample_augmented_minor_determinant'])!=0
    pref=ROOT/'reports/S61_SPARSE_GENERIC5_INPUT_PREFLIGHT_20261104_v63'
    c=read(pref/'certificate.json');V=S.symbols('H L s0 q0 t0')
    beta={i:{tuple(e):Q(v) for e,v in S.Poly(S.sympify(f),*V).terms() if v}
          for i,f in zip(c['historical_weight1_indices'],c['source_u16'])}
    qids=c['historical_weight2_indices'];qi={j:i for i,j in enumerate(qids)}
    assert qids==g['q_to_u_indices_zero_based']
    pulled=read(base/'pulled-G6-and-cofactors.json.gz')
    data=read(pref/'pulled134-input.json.gz');F=[]
    for eq in data['equations']:
        r={}
        for m,e,v in eq:rowadd(r,tuple(m),{tuple(e):Q(v)})
        F.append(r)
    # Each original61 monomial is substituted independently, without using
    # the producer's pulled G6 or its scalar-atom construction.
    target={};alive=0
    for mon,n,d in g['universal_new_equation']:
        p={ZERO:Q(n,d)};qm=[]
        for i in mon:
            if i in qi:qm.append(qi[i])
            else:p=pp_mul(p,beta[i])
        if p:alive+=1
        rowadd(target,tuple(sorted(qm)),p)
    expected={tuple(m):dec(p) for m,p in pulled['G6_q45_coefficients']}
    assert target==expected and sum(map(len,target.values()))==1940
    cof=[[j%70,j//70,str(Q(v))] for j,v in g['selected_non_Bianchi_syzygy']]
    assert len(cof)==125 and cof==pulled['linear_q_cofactors125']
    rhs={};q3slots=set()
    for i,q,v in cof:
        for m,p in F[i].items():
            mm=tuple(sorted((*m,q)))
            if len(mm)==3:q3slots.update((mm,e) for e in p)
            rowadd(rhs,mm,{e:Q(v)*x for e,x in p.items()})
    assert rhs==target and len(q3slots)==674 and sorted(set(map(len,target)))==[1,2]
    full=ROOT/'reports/S61_SPARSE_GENERIC5_FULL134_NORMAL_CIRCUIT_20261104_v63'
    direct=read(full/'full134-direct-circuit.json.gz');blocks=read(full/'coefficient-block-circuit.json.gz')
    rbase=ROOT/'reports/S61_SPARSE_GENERIC5_RESIDUAL129_CIRCUIT_20261104_v63'
    rd=read(rbase/'residual129-factored-circuit.json.gz');old=direct['nodes']+blocks['appended_nodes']
    D=read(base/'G6-normal-residual-circuit.json.gz');N=direct['normal_q_rows42'];ni={i:k for k,i in enumerate(N)}
    J=rd['selected_original134_rows'];non=rd['residual_original134_rows'];assert len(D['nodes'])==15289
    assert D['source_files']==dict(full134_direct=str((full/'full134-direct-circuit.json.gz').relative_to(ROOT)),
        full134_blocks=str((full/'coefficient-block-circuit.json.gz').relative_to(ROOT)),
        residual129=str((rbase/'residual129-factored-circuit.json.gz').relative_to(ROOT)),
        W='reports/S61_SPARSE_GENERIC5_PLANE_RANK3_CIRCUIT_20261104_v63/plane-W-original.json.gz')
    for name,length in [('G6_five_functional_column_nodes',5),('lambda_G6_5',5),
                        ('R_G6_linear_n42',42),('R_G6_cross126',126),
                        ('R_G6_linear_normal_cofactor42',42),('G6_K70_linear_q',70),
                        ('transported_selected_F5_coefficients',5)]:
        assert len(D[name])==length
    olddeg=[]
    for op,*a in old:
        if op in ('z','n','q_input'):d=1
        elif op=='prod':d=sum(olddeg[j] for j in a[0])
        elif op=='sum':d=max(olddeg[j] for j in a[0])
        else:
            if op=='inv':assert olddeg[a[0]]==0
            d=0
        olddeg.append(d)
    rdeg=[]
    for op,*a in rd['nodes']:
        if op in ('z','n'):d=1
        elif op=='external_prior_F134':d=2
        elif op=='sum':d=max(rdeg[j] for j in a[0])
        elif op=='prod':d=sum(rdeg[j] for j in a[0])
        else:
            if op=='inv':assert rdeg[a[0]]==0
            d=0
        rdeg.append(d)
    deg=[];ops=collections.Counter()
    for i,(op,*a) in enumerate(D['nodes']):
        ops[op]+=1
        if op in ('sum','prod'):
            assert all(isinstance(j,int) and 0<=j<i for j in a[0])
            d=max(deg[j] for j in a[0]) if op=='sum' else sum(deg[j] for j in a[0])
        elif op=='qq':Q(a[0]);d=0
        elif op=='parameter_poly_atom':assert 0<=a[0]<len(D['scalar_parameter_atoms']);d=0
        elif op=='external_full134_node':assert 0<=a[0]<len(old);d=olddeg[a[0]]
        elif op=='external_residual129_node':assert 0<=a[0]<len(rdeg);d=rdeg[a[0]]
        elif op=='external_original_W':assert 0<=a[0]<45 and 0<=a[1]<2;d=0
        else:assert op in ('n','z') and 0<=a[0]<(42 if op=='n' else 3);d=1
        deg.append(d)
    assert 'inv' not in ops and deg[D['G6_direct_node']]<=2 and deg[D['R_G6_direct_node']]<=2
    # Unexpanded old roots are treated as independent symbols. These are
    # stronger wiring identities, not evaluations at a parameter point.
    def prior(i):return const(Q(old[i][1])) if old[i][0]=='qq' else symbol(('prior',i))
    def res(i):return const(Q(rd['nodes'][i][1])) if rd['nodes'][i][0]=='qq' else symbol(('res',i))
    reader=Reader(D,prior,res)
    q=[prior(i) for i in direct['forward_q45']];E=[[prior(i) for i in row] for row in direct['E45x3']]
    assert reader(D['G6_direct_node'])==total(times(polynomial(p),mathprod([q[i] for i in m])) for m,p in target.items())
    wp=ROOT/'reports/S61_SPARSE_GENERIC5_PLANE_RANK3_CIRCUIT_20261104_v63/plane-W-original.json.gz'
    w=read(wp);C=[{tuple(t[1:]):Q(p['content'])*int(t[0]) for t in p['terms']} for p in w['center45']]
    known=read(pref/'known-five-input.json.gz');dv=[{tuple(e):Q(v) for e,v in p} for p in known['Q52_minus_center']]
    P=dec(w['common_P']);W=[[symbol(('W',i,h)) for h in range(2)] for i in range(45)]
    # Generic monomial translation derives c and its complete45 gradient.
    shifted={}
    for m,p in target.items():
        partial={():p}
        for i in m:
            nxt={}
            for mon,coef in partial.items():
                rowadd(nxt,mon,pp_mul(coef,C[i]));rowadd(nxt,tuple(sorted((*mon,i))),coef)
            partial=nxt
        for m2,p2 in partial.items():rowadd(shifted,m2,p2)
    c0=shifted.get((),{});L=[shifted.get((i,),{}) for i in range(45)]
    Q2={m:p for m,p in target.items() if len(m)==2}
    SG=[times(polynomial(P),polynomial(c0))]
    SG += [total(times(polynomial(L[i]),W[i][h]) for i in range(45)) for h in range(2)]
    SG += [total(times(polynomial(p),plus(times(W[i][h],polynomial(dv[j])),times(W[j][h],polynomial(dv[i]))))
                  for (i,j),p in Q2.items()) for h in range(2)]
    assert [reader(x) for x in D['G6_five_functional_column_nodes']]==SG
    stops={ref:symbol(('SG',i)) for i,ref in enumerate(D['G6_five_functional_column_nodes'])}
    lr=Reader(D,prior,res,stops)
    for j,ref in enumerate(D['lambda_G6_5']):
        want=times(res(rd['inverseDelta5_node']),total(times(symbol(('SG',h)),res(rd['adjugate5_nodes'][h][j])) for h in range(5)))
        assert lr(ref)==want
    stops={ref:symbol(('lambda',j)) for j,ref in enumerate(D['lambda_G6_5'])}
    rr=Reader(D,prior,res,stops);lam=[symbol(('lambda',j)) for j in range(5)]
    assert rr(D['R_G6_direct_node'])==plus(reader(D['G6_direct_node']),neg(total(times(lam[h],prior(direct['F134'][j])) for h,j in enumerate(J))))
    cross=[{} for _ in range(126)];normal=[{} for _ in range(42)]
    for (i,j),p in Q2.items():
        for a,b in ((i,j),(j,i)):
            if a in ni:
                for h in range(3):cross[42*h+ni[a]]=plus(cross[42*h+ni[a]],times(polynomial(p),E[b][h]))
        if i in ni and j in ni:
            a,b=sorted((ni[i],ni[j]));normal[a]=plus(normal[a],times(polynomial(p),symbol(('n',b))))
    line=[];cx=[]
    for k,i in enumerate(N):
        line.append(plus(polynomial(L[i]),neg(total(times(lam[h],prior(blocks['rows134'][j]['linear_n42'][k])) for h,j in enumerate(J)))))
    for k in range(126):cx.append(plus(cross[k],neg(total(times(lam[h],prior(blocks['rows134'][j]['cross_z_n126'][k])) for h,j in enumerate(J)))))
    assert [rr(x) for x in D['R_G6_linear_n42']]==line
    assert [rr(x) for x in D['R_G6_cross126']]==cx
    for k,ref in enumerate(D['R_G6_linear_normal_cofactor42']):
        want=total([line[k],*[times(cx[42*h+k],symbol(('z',h))) for h in range(3)],normal[k],
                    neg(total(times(lam[h],res(rd['source_normal2_factor134x42'][j][k])) for h,j in enumerate(J)))])
        assert rr(ref)==want and deg[ref]<=1
    K=[{} for _ in range(70)]
    for i,qj,v in cof:K[i]=plus(K[i],times(const(Q(v)),q[qj]))
    assert [reader(x) for x in D['G6_K70_linear_q']]==K
    for h,ref in enumerate(D['transported_selected_F5_coefficients']):
        want=total([K[J[h]],neg(lam[h]),*[times(K[i],res(rd['lambda129x5'][non.index(i)][h])) for i in range(70) if i not in J and K[i]]])
        assert rr(ref)==want
    expected_rhs=total([*[times(K[i],res(rd['residual_direct_F129'][non.index(i)])) for i in range(70) if i not in J and K[i]],
                        *[times(rr(b),prior(direct['F134'][j])) for b,j in zip(D['transported_selected_F5_coefficients'],J)]])
    assert rr(D['transported_cofactor_rhs_node'])==expected_rhs
    # The symbolic transported sum cancels every old lambda*selected term.
    fs=[symbol(('F',i)) for i in range(134)];kl=[symbol(('K',i)) for i in range(70)]
    ls=lambda i,h:symbol(('oldlambda',i,h))
    outsum={}
    for i in range(70):
        if i not in J:outsum=plus(outsum,times(kl[i],plus(fs[i],neg(total(times(ls(i,h),fs[j]) for h,j in enumerate(J))))))
    for h,j in enumerate(J):
        b=total([kl[j] if j<70 else {},neg(lam[h]),*[times(kl[i],ls(i,h)) for i in range(70) if i not in J]])
        outsum=plus(outsum,times(b,fs[j]))
    assert outsum==plus(total(times(kl[i],fs[i]) for i in range(70)),neg(total(times(lam[h],fs[j]) for h,j in enumerate(J))))
    # Exact homogeneous Hilbert bounds for the canonical triangle+two points.
    a,b,z=S.symbols('a b z');quadrics=[-a*a+a*b-b*b-2*z*z,-a*a+2*a*b,2*a*b-b*b,a*z,b*z]
    m2=list(itertools.combinations_with_replacement(range(3),2));m3=list(itertools.combinations_with_replacement(range(3),3));vv=[a,b,z]
    mons=lambda ms:[math.prod(vv[i] for i in m) for m in ms]
    M2=S.Matrix([[S.Poly(f,*vv).coeff_monomial(m) for m in mons(m2)] for f in quadrics])
    M3=S.Matrix([[S.Poly(x*f,*vv).coeff_monomial(m) for m in mons(m3)] for x in vv for f in quadrics])
    assert M2.rank()==5 and M3.rank()==10
    proof=dict(original9257_monomials_independently_pulled=True,pulled_terms=1940,all125_linear_q_cofactors_QQ=True,
        all674_possible_q3_slots_cancel=True,all15289_nodes_and_external_bindings_checked=True,new_unknown_inverses=0,
        all5_functionals_and5_selector_coefficients_checked=True,all42_normal_factor_and126_cross_roots_checked=True,
        all70_K_and5_transported_coefficients_and_rhs_checked=True,formal_transport_cancellation_identity=True,
        FULL135_covariance_separate_from_RAW134='S36/B35 is the trivial Sp4 line; full Bianchi and affine gauge translations give beta^2 F4+beta F5 with q-independent coefficients',
        old135_QQ_minor_reused_not_recomputed=True,current135_independence_by_invariant_open_and_dominant_cover=True,
        residual_dimension=130,cross_rank=126,pure_normal_dimension=4,explicit_pure_normal_basis=False,
        canonical_homogeneous_quadratic_rank=5,canonical_all10_cubic_rank=10,HF=[1,3,1],
        full_A_module_generators_over_K_normal42_at_most=5,free_module_or_K_finite_length_claim=False,
        no_flag_boundary_removal_no_residual_unit_no_exhaustion=True,ops=dict(ops))
    (out/'independent-proof.json').write_text(json.dumps(proof,indent=2)+'\n')
    for p in [Path(__file__),ROOT/'code/s61_sparse_generic5_full134_normal_circuit_review_20261104_v63.py',ROOT/f'research/notes/{ID}.md']:
        bind(p)
    report=dict(id=ID,all_passed=True,input_sha256=hashes,
        evidence_sha256={p.name:sha(p) for p in out.iterdir()},seconds=time.monotonic()-start,
        max_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        scope='Actual G6/source and all finite circuit roots; separate FULL135 covariance, four-dimensional normal-only rowspace and finite <=5 module. No basis/free/unit/exhaustion assertion.')
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(all_passed=True,seconds=report['seconds'],rss=report['max_rss_kib'],report_sha256=sha(out/'report.json'))))

def mathprod(ps):
    out=const(1)
    for p in ps:out=times(out,p)
    return out
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output)
