#!/usr/bin/env python3
"""Independent fixed Q-field certificate checks plus all22 literal source evaluations."""
import argparse, hashlib, itertools, json, math, resource, signal, subprocess, time
from fractions import Fraction as F
from pathlib import Path
from s61_nonsquare_raw_H_B_curve_components_review_20261104_v63 import poly,add,mul,scale,power,Z
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_RAW_H_B_Q_COMPONENT_FINITE_20261104_v63_REVIEW'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
class Program:
    def __init__(self):self.lines=[];self.registers=0;self.exceptions=set();self.cache={}
    def put(self,s,register=True):
        self.lines.append(s)
        if register:i=self.registers;self.registers+=1;return i
    def e_poly(self,p):return ' '.join([str(len(p))]+[f'{i} {F(c)}' for i,c in p])
    def coeff(self,c):
        parts=[str(len(c))]
        for r in c:
            a=r['a_power'];q=r['fraction'];nu=q['numerator_e'];de=q['denominator_e'];assert de
            deg=max(int(i) for i,v in de)
            if deg:
                lc=next(F(v) for i,v in de if int(i)==deg)
                self.exceptions.add(tuple(sorted((int(i),str(F(v)/lc)) for i,v in de)))
            parts += [str(a),self.e_poly(nu),self.e_poly(de)]
        return ' '.join(parts)
    def dp(self,p):
        key=json.dumps(p,sort_keys=True,separators=(',',':'))
        if key not in self.cache:self.cache[key]=self.put('P '+str(len(p))+' '+' '.join(str(r['d_power'])+' '+self.coeff(r['coefficient']) for r in p))
        return self.cache[key]
    def source(self,p):
        by={}
        for (a,e,d),v in p.items():by.setdefault(d,{}).setdefault(a,[]).append([e,str(v)])
        return self.dp([dict(d_power=d,coefficient=[dict(a_power=a,fraction=dict(numerator_e=es,denominator_e=[[0,'1']])) for a,es in aa.items()]) for d,aa in by.items()])
    def op(self,o,*x):return self.put(o+' '+' '.join(map(str,x)))
    def equal(self,a,b):self.put(f'E {a} {b}',False)
    def degree(self,a,d):self.put(f'DEG {a} {d}',False)
    def nz(self,a):self.put(f'NZ {a}',False)
    def summ(self,items,zero):
        result=zero
        for i in items:result=self.op('A',result,i)
        return result
def main(out,work):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
    out.mkdir(parents=True,exist_ok=False);work.mkdir(parents=True,exist_ok=False);hashes={};base=ROOT/'reports';pre='S61_NONSQUARE_RAW_H_B_';end='_20261104_v63'
    def bind(p):hashes[str(p.relative_to(ROOT))]=sha(p)
    def audit(p):
        d=read(p);bind(p)
        for n,h in d.get('input_sha256',{}).items():assert sha(ROOT/n)==h,n;hashes[n]=h
        for n,h in d.get('evidence_sha256',{}).items():assert sha(p.parent/n)==h,n;bind(p.parent/n)
        return d
    oldgate=base/(pre+'CURVE_COMPONENTS'+end+'_REVIEW')/'report.json';assert audit(oldgate)['all_passed']
    normgate=audit(base/(pre+'ADDITIONAL_NORM_PENCIL'+end+'_REVIEW')/'report.json')
    assert normgate['all_passed'] and normgate['source_monic_norm_pencil_and_residual_cover_accepted']
    assert normgate['finite_residual_geometric_point_bound']==1600
    qbase=base/(pre+'COMPONENT_EXACT_QUOTIENT_Q'+end);qb=audit(qbase/'report.json');assert qb['all_whole_field_identities_checked']
    rbase=base/(pre+'Q_ROOT_ALL_MINORS'+end);rb=audit(rbase/'report.json');assert rb['all_independent_integer_evaluations_memberships_norms_checked']
    cert=read(qbase/'certificate.json');rootcert=read(rbase/'native-certificate.json');iv=read(rbase/'independent-verification.json')
    inputs=read(base/(pre+'CURVE_COMPONENTS_INPUT'+end)/'certificate.json');rec=next(r for r in inputs['component_records'] if r['factor_name']=='Q')
    assert cert['factor_source']==rec['factor'] and cert['source_reduced_polynomials']==[p['remainder'] for p in rec['polynomials']]
    Q=poly(rec['factor']);q0={m:c for m,c in Q.items() if m[0]==0};q1={(0,m[1],0):c for m,c in Q.items() if m[0]==1}
    P=Program();zero=P.source({});one=P.source({Z:F(1)});sources=[P.source(poly(p)) for p in cert['source_reduced_polynomials']]
    final=cert['final'];G=P.dp(final['common_gcd']);P.degree(G,1)
    for a,b in zip(sources,final['all3_source_quotients_by_gcd']):P.equal(a,P.op('M',P.dp(b),G))
    P.equal(P.summ([P.op('M',P.dp(b),a) for a,b in zip(sources,final['Bezout3'])],zero),G)
    steps=0;firstg=None
    for index,t in enumerate(cert['exact_Euclid_traces']):
        r0,r1=(sources[0],sources[1]) if index==0 else (firstg,sources[2]);initial=(r0,r1)
        for s in t['steps']:
            lc=P.dp([dict(d_power=0,coefficient=s['divisor_LC'])]);inv=P.dp([dict(d_power=0,coefficient=s['divisor_LC_inverse'])])
            P.equal(P.op('M',lc,inv),one);P.equal(P.op('LC',r1),lc);P.degree(r1,s['divisor_degree'])
            quotient,remainder=P.dp(s['quotient']),P.dp(s['remainder']);P.equal(r0,P.op('A',P.op('M',quotient,r1),remainder))
            assert not s['remainder'] or max(x['d_power'] for x in s['remainder'])<s['divisor_degree']
            r0,r1=r1,remainder;steps+=1
        P.equal(r1,zero);gg=P.dp(t['gcd_monic']);P.equal(P.op('LC',gg),one)
        P.equal(r0,P.op('M',P.op('LC',r0),gg));P.equal(P.summ([P.op('M',P.dp(c),a) for a,c in zip(initial,t['source2_Bezout'])],zero),gg)
        if index==0:firstg=gg
        else:P.equal(gg,G)
    assert steps==6 and not cert['U_support_removal_steps'] and not cert['generic_common_roots_all_U_boundary']
    # Coefficient c in G=d+c equals the literal N/D, with no extension inverse.
    n={(1,0,0):F(18000),(0,2,0):F(6480),(0,1,0):F(3420),Z:F(-980)}
    den={(0,3,0):F(486),(0,2,0):F(243),(0,1,0):F(-297),Z:F(-54)}
    assert n==poly(iv['root_numerator_N']) and den==poly(iv['root_denominator_D'])
    constant=next(r for r in final['common_gcd'] if r['d_power']==0)
    cr=P.dp([constant]);P.equal(P.op('M',cr,P.source(den)),P.source(n))
    negative=P.op('N',cr);U=P.dp(final['U_reduced_on_component']);P.nz(P.op('EV',U,negative))
    Uraw=poly(read(base/(pre+'35_MINORS'+end)/'certificate.json')['FULLU_product'])
    P.equal(U,P.source(Uraw)) # original polynomial interpreted in the exact quadratic field
    listed={tuple(sorted((int(i),str(F(v))) for i,v in p)) for p in cert['finite_e_exception_polynomials']}
    assert P.exceptions==listed and len(listed)==19 and sum(max(i for i,c in x) for x in listed)==320
    header=lambda p:P.e_poly([[m[1],str(v)] for m,v in sorted(p.items())])
    inp=out/'fixed-field-identities.txt';inp.write_text(header(q0)+'\n'+header(q1)+'\n'+str(len(P.lines))+'\n'+'\n'.join(P.lines)+'\n')
    assert inp.stat().st_size<32*1024*1024
    (out/'preflight.json').write_text(json.dumps(dict(serialized_Q_field_terms=2935,fixed_field_instructions=len(P.lines),input_bytes=inp.stat().st_size,all22_source_terms=2963,field_inversions_GCD_or_Euclid_to_compute=0),indent=2)+'\n')
    cpp=ROOT/'code/s61_nonsquare_Q_component_review_20261104_v63.cpp';exe=work/'verify-Q'
    subprocess.run(['bash','-c','source /workspace/cloud-setup/activate.sh; c++ -O2 -std=c++17 "$1" -o "$2" $(pkg-config --cflags --libs flint) -lgmp','compile',str(cpp),str(exe)],check=True,timeout=12)
    with (out/'native.log').open('w') as log:
        subprocess.run(['bash','-c','source /workspace/cloud-setup/activate.sh; exec "$1" "$2"','check',str(exe),str(inp)],stdout=log,stderr=subprocess.STDOUT,check=True,timeout=max(1,45-(time.monotonic()-start)))
    native=read(out/'native.log');assert native['all_passed'] and native['field_inverse_or_GCD_calls']==0
    # Original three-variable, denominator-cleared polynomial checks happen
    # BEFORE reducing modulo Q; these also verify the saved division quotients.
    minors=read(base/(pre+'35_MINORS'+end)/'certificate.json');A=list(map(poly,minors['localized_ideal_generators_monic_QQ_after_namedunit_and_QQ_only_scaling']))
    hc=read(base/(pre+'PAIR_MONIC_COMPLETION'+end)/'certificate.json');H={(m[0]+5,m[1],m[2]):c for m,c in poly(hc['monic_H_laurent']).items()}
    onepoly={Z:F(1)};np=[onepoly];ds=[onepoly]
    for j in range(7):np.append(mul(np[-1],scale(n,-1)));ds.append(mul(ds[-1],den))
    nonzero=[]
    def check_source(source,record):
        degree=max(m[2] for m in source);value={}
        for j in range(degree+1):
            coeff={(m[0],m[1],0):v for m,v in source.items() if m[2]==j}
            value=add(value,mul(coeff,mul(np[j],ds[degree-j])))
        assert value==poly(record['cleared_evaluation'])
        quotient,remainder=poly(record['Q_quotient']),poly(record['remainder'])
        assert add(mul(quotient,Q),remainder)==value and all(m[0]<=1 and m[2]==0 for m in remainder)
        r0={m:v for m,v in remainder.items() if m[0]==0};r1={(0,m[1],0):v for m,v in remainder.items() if m[0]==1}
        conjugate=add(add(r0,scale(mul(q1,r1),-1)),{(1,m[1],0):-v for m,v in r1.items()})
        assert conjugate==poly(record['remainder_conjugate'])
        norm=add(add(mul(r0,r0),scale(mul(q1,mul(r0,r1)),-1)),mul(q0,mul(r1,r1)))
        assert norm==poly(record['remainder_quadratic_norm_e']) and all(m[0]==m[2]==0 for m in norm)
        assert add(mul(remainder,conjugate),scale(mul(poly(record['norm_membership_Q_quotient']),Q),-1))==norm
        return remainder,norm
    assert len(A)==len(rootcert['all22_evaluations'])==22
    for i,(source,record) in enumerate(zip(A,rootcert['all22_evaluations'])):
        remainder,norm=check_source(source,record)
        if remainder:assert norm;nonzero.append(i)
    ur,un=check_source(Uraw,rootcert['localizer_U_evaluation']);hr,hn=check_source(H,rootcert['Hclear_evaluation'])
    assert not hr and not hn and ur and un
    expected=[0,2,4,6,7,8,9,11,12,13,14,16,18,19,21]
    assert nonzero==rootcert['nonzero_remainder_indices']==rb['nonzero_indices']==expected
    proof=dict(independent_Q_field_exact_arithmetic_without_inverse_GCD_or_Euclid=True,native=native,
        all6_saved_Euclid_divisions_LC_inverses_and_Bezout_identities=True,all3_source_divisibilities=True,
        generic_common_gcd='d+20*(900*a+324*e^2+171*e-49)/(27*(3*e-2)*(6*e+1)*(e+1))',
        U_at_generic_root_nonzero=True,all19_finite_field_denominator_exceptions_retained=True,
        all22_original_source_cleared_evaluations_quotients_remainders_conjugates_norms_checked=True,
        nonzero_actual_minor_remainder_indices=nonzero,Hclear_at_root_zero=True,Uclear_and_norm_nonzero=True,
        Q_generic_full_B_bad_curve_excluded=True,
        FULLU_B_bad_finite_or_empty_using_prior_norm_pencil_curve_plus_finite_residual_decomposition=True,
        all_denominator_D_norm_and_old_h_over_G_k_over_G_exceptions_kept=True,
        empty_or_allFULLU_rank3_H7_cubic57_quartic239_or_exhaustion_not_claimed=True)
    (out/'independent-proof.json').write_text(json.dumps(proof,indent=2)+'\n')
    for p in (Path(__file__),cpp,ROOT/'code/s61_nonsquare_raw_H_B_curve_components_review_20261104_v63.py',ROOT/f'research/notes/{ID}.md'):bind(p)
    report=dict(id=ID,all_passed=True,input_sha256=hashes,evidence_sha256={p.name:sha(p) for p in out.iterdir()},seconds=time.monotonic()-start,
        max_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        scope='Exact Q field/source certificates and all22 actual source minors; with accepted L and prior curve+finite decomposition, whole FULLU B-rankbad is finite or empty, never asserted empty.')
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(dict(all_passed=True,seconds=report['seconds'],rss=report['max_rss_kib'],report_sha256=sha(out/'report.json'))))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--work',type=Path,required=True);a=p.parse_args();main(a.output,a.work)
