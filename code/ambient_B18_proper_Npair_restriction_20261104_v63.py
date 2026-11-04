#!/usr/bin/env python3
"""ONE exact d=0 core-ideal restriction and five constant-unit eliminations."""
import argparse,hashlib,json,math,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='AMBIENT_B18_PROPER_NPAIR_RESTRICTION_20261104_v63'
CAP=150000;METRICS={'largest_polynomial_support':0,'multiplications':0}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def norm(p):
    p={m:v for m,v in p.items()if v};METRICS['largest_polynomial_support']=max(METRICS['largest_polynomial_support'],len(p));assert len(p)<=CAP;return p
def add(p,q):
    r=dict(p)
    for m,c in q.items():r[m]=r.get(m,F(0))+c
    return norm(r)
def scale(p,c):return norm({m:v*c for m,v in p.items()})
def mul(p,q,flag=None):
    METRICS['multiplications']+=1;r={}
    for m,a in p.items():
        for n,b in q.items():
            e=tuple(i+j for i,j in zip(m,n))
            if flag is None or e[flag]<=1:r[e]=r.get(e,F(0))+a*b
    return norm(r)
def power(p,n,dim,flag=None):
    r={(0,)*dim:F(1)}
    for _ in range(n):r=mul(r,p,flag)
    return r
def variable(k,n):return{tuple(int(i==k)for i in range(n)):F(1)}
def decode(p):return norm({tuple(m):F(a,b)for m,a,b in p})
def encode(p):return[[list(m),c.numerator,c.denominator]for m,c in sorted(p.items(),reverse=True)]
def substitute(p,images,flag=None):
    dim=len(next(iter(images[0])))if images[0]else len(next(iter(next(v for v in images if v))))
    caches=[{0:{(0,)*dim:F(1)},1:v}for v in images];out={}
    for m,c in p.items():
        v={(0,)*dim:c}
        for j,k in enumerate(m):
            if k:
                if k not in caches[j]:caches[j][k]=power(images[j],k,dim,flag)
                v=mul(v,caches[j][k],flag)
        out=add(out,v)
    return out
def pretty(p,names):
    if not p:return'0'
    terms=[]
    for m,c in sorted(p.items(),reverse=True):
        mon='*'.join(name+(f'^{k}'if k>1 else'')for name,k in zip(names,m)if k)
        body=(str(abs(c))+'*'if abs(c)!=1 and mon else str(abs(c))if not mon else'')+mon
        terms.append(('-'if c<0 else'+'if terms else'')+body)
    return''.join(terms)
def main(out):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
    out=out.resolve();out.mkdir(parents=True,exist_ok=False);hashes={}
    def read(p):hashes[str(p.relative_to(ROOT))]=sha(p);return json.loads(p.read_text())
    def audit(p):
        r=read(p);assert r['all_passed']
        for n,h in r.get('input_sha256',{}).items():assert sha(ROOT/n)==h,n;hashes[n]=h
        for n,h in r.get('evidence_sha256',{}).items():assert sha(p.parent/n)==h,n;hashes[str((p.parent/n).relative_to(ROOT))]=h
        return r
    b=ROOT/'reports/AMBIENT_B18_PROPER_RATIONAL_ELIMINATION_20261104_v63';audit(b/'report.json')
    P=read(b/'actual12-flag-divided-polynomials.json');assert P['remaining9_variable_order']==['a','d','v2','v3','v4','v5','v6','v7','v8']
    sourcefile=ROOT/P['original_source_file'];source=read(sourcefile);audit(sourcefile.parent/'report.json')
    assert source['CA_fiber_variables']+source['flag_variables']==P['original14_variable_order']
    h9=[decode(r['divided_primitive_polynomial'])for r in P['all12_actual_remaining_and_divided_records']]
    h8=[{m[:1]+m[2:]:c for m,c in h.items()if m[1]==0}for h in h9]
    assert len(h8)==12 and all(p for p in h8)
    names=['a','v2','v3','v4','v5','v6','v7','v8'];z8=(0,)*8;retained=[0,1,5];target=['a','v2','v6']
    planned=[(0,3,F(10)),(4,2,F(8)),(1,4,F(80)),(2,6,F(16)),(3,7,F(192))]
    pre=dict(source_H12_sha256=sha(b/'actual12-flag-divided-polynomials.json'),closed_flag='d=0',D_at_flag=1,
        ring_variables8=names,input_terms=sum(map(len,h8)),input_max_support=max(map(len,h8)),
        fixed_pivot_rows=[r for r,j,c in planned],fixed_pivot_names=[names[j]for r,j,c in planned],
        retained_variables=target,term_cap=CAP,hard_seconds=60,memory_bytes=2_000_000_000,
        no_algebraic_roots_GB_decomposition_or_multiplicity_job=True)
    (out/'preflight.json').write_text(json.dumps(pre,indent=2)+'\n')
    # Independently bind H_i(d=0) to the exact coefficient of d in the original
    # source substitution. D^-m=1 mod d², so this is an exact coefficient identity.
    a,d,v2,v3,v4,v5,v6,v7,v8=[variable(j,9)for j in range(9)];zero={}
    root1={'v0':add(scale(a,-1),scale(d,2)),
        'b':scale(mul(d,add(a,scale(v2,-1)),1),F(1,2)),
        'e':mul(d,v6,1),
        'v1':add(scale(mul(d,add(a,scale(v2,-1)),1),F(1,4)),scale(mul(d,v6,1),2)),
        'c':mul(d,add(scale(mul(v6,add(a,v2),1),F(1,2)),scale(v3,-1)),1)}
    for name,entry in P['five_explicit_rational_roots'].items():
        num={m:c for m,c in decode(entry['numerator']).items()if m[1]<=1};assert num==root1[name]
    assert P['original14_variable_order']==['v0','v1','v2','v3','v4','v5','v6','v7','v8','a','b','c','d','e']
    images=[root1['v0'],root1['v1'],v2,v3,v4,v5,v6,v7,v8,a,root1['b'],root1['c'],d,root1['e']]
    original=[decode(q)for q in source['equations14variable_sparse']];assert len(original)==17
    pulled=[substitute(q,images,1)for q in original];assert all(not q for q in pulled[:5])
    for i,(record,h)in enumerate(zip(P['all12_actual_remaining_and_divided_records'],h9)):
        assert record['original_equation_row']==5+i
        expected=scale({tuple(e+(j==1)for j,e in enumerate(m)):c for m,c in h.items()if m[1]==0},F(record['primitive_content_unit'],record['constant_denominator']))
        assert pulled[i+5]==expected
    # Every selected equation is exactly a CONSTANT unit times its graph
    # equation after the preceding substitutions. No nonlinear root choice.
    substitutions=[variable(j,8)for j in range(8)];steps=[];roots={}
    for row,pivot,unit in planned:
        reduced=substitute(h8[row],substitutions)
        assert all(m[pivot]<=1 for m in reduced)
        coefficient={tuple(e-(j==pivot)for j,e in enumerate(m)):c for m,c in reduced.items()if m[pivot]==1}
        assert coefficient=={z8:unit}
        rest={m:c for m,c in reduced.items()if m[pivot]==0};root=scale(rest,-1/unit)
        assert all(all(not m[j]for j in range(8)if j not in retained)for m in root)
        graph=add(variable(pivot,8),scale(root,-1));assert reduced==scale(graph,unit)
        steps.append(dict(source_H_index=row,original_source_equation_row=5+row,pivot=names[pivot],leading_constant_unit=str(unit),
            reduced_equation=encode(reduced),polynomial_root=encode(root),triangular_graph_identity_verified=True))
        substitutions[pivot]=root;roots[names[pivot]]=root
    values=[substitute(p,substitutions)for p in h8]
    chosen=[r for r,j,u in planned];assert all(not values[r]for r in chosen)
    residualrows=[i for i in range(12)if i not in chosen];assert residualrows==list(range(5,12))
    residuals=[]
    for i in residualrows:
        p=values[i];assert all(all(not m[j]for j in range(8)if j not in retained)for m in p)
        q={tuple(m[j]for j in retained):c for m,c in p.items()}
        residuals.append(dict(source_H_index=i,original_source_equation_row=5+i,
            exact_polynomial3=encode(q),expression=pretty(q,target),terms=len(q),degree=max(map(sum,q),default=-1)))
    stored=sum(len(p)for p in h8)+sum(len(r)for r in roots.values())+sum(r['terms']for r in residuals)
    assert stored<=CAP and math.prod(u for r,j,u in planned)==19660800
    result=dict(source_H_file=str((b/'actual12-flag-divided-polynomials.json').relative_to(ROOT)),
        specialized8_variable_order=names,all12_exact_H_d0_polynomials=[encode(p)for p in h8],
        original17_source_linear_d_coefficient_bindings_verified=True,
        source_coefficient_meaning='Each H_i(d=0) equals the exact d coefficient of Fi after the five original unit roots, with its recorded rational primitive scaling; D=1-d² has no linear coefficient.',
        triangular5_elimination=steps,all5_final_selected_equations_zero=True,all_polynomial_roots={name:encode(p)for name,p in roots.items()},
        final3_variable_order=target,seven_remaining_actual_polynomials=residuals,
        quotient_identity='QQ[a,v2,...,v8]/(H0|d0,...,H11|d0) is isomorphic to QQ[a,v2,v6]/(the seven saved residuals); all eliminated pivots have nonzero QQ constant leading coefficients.',
        completed_quotient_identity_same=True,
        correct_geometric_scope='This is (S/(H0,...,H11))/(d), the residual core restricted to the Npair divisor. It is NOT E4/(d); E4 is S/(d*H), whose quotient by d is regular and imposes none of the H equations.',
        original_proper_source_to_E4_transfer_independent_review_still_required_for_geometric_promotion=True,
        no_GB_decomposition_multiplicity_depth_socle_or_D0_claim=True,total_saved_polynomial_terms=stored,metrics=METRICS)
    p=out/'actual-Npair-three-variable-polynomials.json';p.write_text(json.dumps(result,indent=2)+'\n')
    for p in[Path(__file__),ROOT/f'research/notes/{ID}.md']:hashes[str(p.relative_to(ROOT))]=sha(p)
    report=dict(id=ID,all_passed=True,all12_source_d0_coefficients_bound_to_original17=True,
        five_constant_unit_eliminations=True,root_choice_or_unknown_denominator_used=False,
        remaining_variables3=target,remaining_equations7=True,polynomial_degrees=[r['degree']for r in residuals],
        nonzero_residual_polynomials=sum(bool(r['terms'])for r in residuals),constant_pivot_product=19660800,
        source_transfer_peer_pending=True,no_GB_decomposition_multiplicity_or_depth_claim=True,
        input_sha256=hashes,evidence_sha256={p.name:sha(p)for p in out.iterdir()},seconds=time.monotonic()-start,
        max_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    p=out/'report.json';p.write_text(json.dumps(report,indent=2)+'\n');assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
    print(json.dumps(dict(all_passed=True,seconds=report['seconds'],source_terms=stored,remaining_expressions=[r['expression']for r in residuals],report_sha256=sha(p))))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args();main(a.output)
