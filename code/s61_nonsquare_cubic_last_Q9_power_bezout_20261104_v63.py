#!/usr/bin/env python3
"""Fixed conic/power proof and finite a^29 Bezout DAG, without a degree search."""
import argparse
import hashlib
import itertools
import json
import math
import resource
import signal
import time
from pathlib import Path

import sympy as sy
from sympy.polys.rings import ring

ROOT=Path(__file__).resolve().parents[1]
ID='S61_NONSQUARE_CUBIC_LAST_Q9_POWER_BEZOUT_20261104_v63'


def sha(p):
    h=hashlib.sha256()
    with p.open('rb')as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()


def compositions(total,n):
    if n==1:return[(total,)]
    return[(i,)+v for i in range(total+1)for v in compositions(total-i,n-1)]


class DAG:
    def __init__(self,R):
        self.R=R;self.atoms=[];self.ai={};self.nodes=[['z'],['one']];self.ni={('z',):0,('one',):1};self.degrees=[0,0]
    def node(self,n,deg):
        key=tuple(n)
        if key not in self.ni:
            assert len(self.nodes)<200000
            self.ni[key]=len(self.nodes);self.nodes.append(n);self.degrees.append(deg)
        return self.ni[key]
    def c(self,p):
        p=self.R(p)
        if not p:return 0
        if p==self.R.one:return 1
        key=tuple(sorted((m,str(c))for m,c in p.items()))
        if key not in self.ai:self.ai[key]=len(self.atoms);self.atoms.append(p)
        return self.node(['c',self.ai[key]],max(sum(m)for m in p))
    def add(self,a,b):
        if not a:return b
        if not b:return a
        return self.node(['a',a,b],max(self.degrees[a],self.degrees[b]))
    def mul(self,a,b):
        if not a or not b:return 0
        if a==1:return b
        if b==1:return a
        return self.node(['m',a,b],self.degrees[a]+self.degrees[b])
    def prod(self,vals):
        z=1
        for v in vals:z=self.mul(z,v)
        return z
    def power(self,x,n):return self.prod([x]*n)


def main(out):
    began=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(57)
    out.mkdir(parents=True,exist_ok=False);read=lambda p:json.loads(p.read_text())
    b=ROOT/'reports/S61_NONSQUARE_CUBIC_LAST_Q9_PREFLIGHT_20261104_v63'
    u=ROOT/'reports/S61_NONSQUARE_CUBIC_ZW25_UNITS_20261104_v63'
    pd=ROOT/'reports/S61_NONSQUARE_CUBIC_ZW25_PAIR_DUAL_PREFLIGHT_20261104_v63'
    inp=ROOT/'reports/S61_NONSQUARE_CUBIC_ZW25_AUGMENTED_INPUT_20261104_v63/input.json'
    paths=[Path(__file__),b/'report.json',b/'certificate.json',u/'checkpoint.json',u/'report.json',
           pd/'certificate.json',pd/'report.json',inp,
           ROOT/'reports/S61_NONSQUARE_CUBIC_Q2_L_UNIT_COMPLETION_20261104_v63_REVIEW/report.json']
    br,bc,cp,ur,pc,pr,data=[read(p)for p in paths[1:8]]
    assert sha(paths[2])==br['certificate_sha256']and sha(paths[3])==ur['evidence_sha256']['checkpoint.json']
    assert sha(paths[5])==pr['certificate_sha256']and sha(inp)==cp['source_input_sha256']
    assert read(paths[8])['all_passed']
    R,a,e,d=ring('a,e,d',sy.QQ);K=R.to_field()
    decp=lambda ts:R.from_dict({tuple(m):sy.QQ(c)for m,c in ts})
    dec=lambda q:K(decp(q[0]))/K(decp(q[1]))
    encp=lambda p:[[list(m),str(c)]for m,c in p.items()]
    enc=lambda q:[encp(q.numer),encp(q.denom)]
    poly=lambda q:q.numer.exquo(q.denom)if q.denom.is_ground else (_ for _ in()).throw(AssertionError('nonpolynomial atom'))
    Y,y1,y2,y3=ring('y1,y2,y3',K);ys=[y1,y2,y3]
    pairmons=[(2,0,0),(1,1,0),(1,0,1),(0,2,0),(0,1,1),(0,0,2)]
    ym=lambda m: y1**m[0]*y2**m[1]*y3**m[2]
    initial=[[K(decp(q))for q in row]for row in data['matrix25x32'][:10]]
    pair=[[initial[i][j]for j in range(20,26)]for i in [8,9]]
    assert pair==[[dec(q)for q in row]for row in pc['Lframe_pair2x6']]
    dual=[dec(row[1])for row in pc['Lframe_dual6x2']]
    assert sum((pair[0][i]*dual[i]for i in range(6)),K.zero)==0
    assert sum((pair[1][i]*dual[i]for i in range(6)),K.zero)==1
    assert all((K(a*a)*v).denom.is_ground for v in dual)
    P=sum((Y(c)*ym(m)for c,m in zip(pair[0],pairmons)),Y.zero)
    RR=sum((Y(c)*ym(m)for c,m in zip(pair[1],pairmons)),Y.zero)
    ell=Y(initial[5][20])*y1+Y(initial[5][21])*y2+Y(initial[5][22])*y3
    for row,mult in zip([5,6,7],ys):
        assert [initial[row][j]for j in range(20,26)]==[(ell*mult).get(m,K.zero)for m in pairmons]
    assert P.get((2,0,0))==K(-5)/4 and P.get((1,0,1))==K(-5)/2 and not P.get((0,0,2),K.zero)
    assert K(-5)/4*K.zero-(K(-5)/4)**2==K(-25)/16
    # Fixed source transcript: only the ten quadratic rows, never the author scan.
    V=[[K.one if i==j else K.zero for j in range(10)]for i in range(10)]
    for st in cp['all_completed_unit_steps']:
        ri,cj=st['pivot_row'],st['pivot_column']
        if ri>=10:
            assert not any(i<10 for i,_ in st['eliminations']);continue
        pivot=dec(st['pivot']);V[ri]=[x/pivot for x in V[ri]]
        for i,co in st['eliminations']:
            if i<10:
                co=dec(co);V[i]=[x-co*y for x,y in zip(V[i],V[ri])]
    assert all(not V[4][i]for i in range(10)if i not in[0,4])
    cc=dec(bc['c_s_dot_v']);vf=[x-cc*y for x,y in zip(V[9],V[5])]
    assert vf[9]==1
    m=sum((Y(vf[i])*ys[i-5]for i in [5,6,7]),Y.zero)
    n=sum((Y(vf[i])*ys[i-1]for i in [1,2,3]),Y.zero);lam=vf[8]
    hvals=[poly(dec(v['actual_coefficient']))for v in bc['actual_nine_coefficients']]
    puremons=[(2,1,0),(1,2,0),(1,1,1),(1,0,2),(0,2,1),(0,1,2),(0,0,3)]
    AA=sum((Y(K(hvals[i]))*ym(mon)for i,mon in enumerate(puremons)),Y.zero)
    v0=K(6*(18*e+7))/K(5*a)
    BB=Y(K(hvals[7]))*y1*y2+Y(K(hvals[8]))*y1*y3+Y(v0*K(hvals[7])/9-2*K(hvals[8]))*y3*y3
    assert AA==RR*m+P*n
    assert BB==RR+ell*m+Y(lam)*P
    # Coordinate change with only nonzero rational determinant.
    S=Y(K(-5)/2)*y1+Y(K(9*e+4)/5)*y2;TT=y2
    W=y3-S/5-Y(K(3*e+8)/25)*TT
    q=K(625*a+108*e**2+156*e+48)/125
    assert P==S*W+Y(q)*TT*TT
    U,s,t,w=ring('S,T,W',K)
    yi=[U(K(-2)/5)*(s-U(K(9*e+4)/5)*t),t,w+s/5+U(K(3*e+8)/25)*t]
    def substitute(f):
        z=U.zero
        for mon,co in f.items():
            v=U(co)
            for im,ex in zip(yi,mon):v*=im**ex
            z+=v
        return z
    def nf(f):
        z=U.zero
        for (i,j,k),co in f.items():
            h=min(i,k);z+=U(co*(-q)**h)*s**(i-h)*t**(j+2*h)*w**(k-h)
        return z
    assert substitute(P)==s*w+U(q)*t*t
    Rbar=nf(substitute(RR))
    basis2=[s*s,s*t,t*t,t*w,w*w];mons2=[(2,0,0),(1,1,0),(0,2,0),(0,1,1),(0,0,2)]
    rv=[poly(Rbar.get(mon,K.zero))for mon in mons2]
    assert Rbar==sum((U(K(x))*mon for x,mon in zip(rv,basis2)),U.zero)
    # Transport full accepted pair functional; kill P before normalizing R.
    Tbasis=[S*S,S*TT,TT*TT,TT*W,W*W]
    tv=[sum((f.get(mon,K.zero)*dual[i]for i,mon in enumerate(pairmons)),K.zero)for f in Tbasis]
    Tpoly=[poly(K(a*a)*x)for x in tv]
    assert sum((x*y for x,y in zip(Tpoly,rv)),R.zero)==a*a
    # The nine coefficient generators of Rbar^2 mod(SW+qT²).
    Z,qq,r0,r1,r2,r3,r4=ring('q,r0,r1,r2,r3,r4',sy.ZZ);rs=[r0,r1,r2,r3,r4]
    Js=[r0*r0,2*r0*r1,r1*r1+2*r0*r2,2*r1*r2-2*qq*r0*r3,
        r2*r2-2*qq*r1*r3+2*qq*qq*r0*r4,2*r2*r3-2*qq*r1*r4,
        r3*r3+2*r2*r4,2*r3*r4,r4*r4]
    def evaluni(f):
        z=R.zero
        for ex,co in f.items():
            v=R(int(co))*poly(q)**ex[0]
            for r,k in zip(rv,ex[1:]):v*=r**k
            z+=v
        return z
    basis4=[(4,0,0),(3,1,0),(2,2,0),(1,3,0),(0,4,0),(0,3,1),(0,2,2),(0,1,3),(0,0,4)]
    Jactual=[evaluni(f)for f in Js]
    square=nf(Rbar*Rbar)
    assert [poly(square.get(mon,K.zero))for mon in basis4]==Jactual
    # Polynomial matrix a*J=A9*h; no source tail is thrown away.
    columns=[]
    for i in range(9):
        if i<7:expr=-Y(K(a))*ell*ym(puremons[i])
        elif i==7:expr=RR*(Y(K(a))*y1*y2+Y(K(2*(18*e+7))/15)*y3*y3)
        else:expr=RR*(Y(K(a))*y1*y3-Y(K(2*a))*y3*y3)
        ff=nf(substitute(expr));columns.append([poly(ff.get(mon,K.zero))for mon in basis4])
    A9=[[columns[j][i]for j in range(9)]for i in range(9)]
    assert all(sum((A9[i][j]*hvals[j]for j in range(9)),R.zero)==a*Jactual[i]for i in range(9))
    assert nf(substitute(RR*BB-ell*AA))==square
    # Fixed monic rewrite identities, not a Groebner basis computation.
    ruleindices=[0,8,2,6,4]
    head=[r0*r0,r4*r4,r1*r1,r3*r3,r2*r2]
    heads=[next(iter(f))for f in head]
    def witness(power):
        rem=power.copy();cofs=[Z.zero for _ in range(9)];steps=0
        while rem:
            chosen=None
            for mon,co in sorted(rem.items(),reverse=True):
                for gi,hmon in zip(ruleindices,heads):
                    if all(x>=y for x,y in zip(mon,hmon)):
                        chosen=(mon,co,gi,hmon);break
                if chosen:break
            assert chosen is not None
            mon,co,gi,hmon=chosen
            fac=Z.from_dict({tuple(x-y for x,y in zip(mon,hmon)):co})
            rem-=fac*Js[gi];cofs[gi]+=fac;steps+=1
            assert steps<5000
        assert sum((c*g for c,g in zip(cofs,Js)),Z.zero)==power
        return cofs,steps
    nilbounds=[2,4,6,4,2]
    powercofs=[];rewritecounts=[]
    for r,bound in zip(rs,nilbounds):
        cofs,count=witness(r**bound);powercofs.append(cofs);rewritecounts.append(count)
    assert sum(x-1 for x in nilbounds)==13
    # Build only a coefficient circuit; do not expand the degree~hundreds cofactors.
    dag=DAG(R);qr=dag.c(poly(q));rroots=[dag.c(x)for x in rv];Troots=[dag.c(x)for x in Tpoly]
    def uninode(f):
        z=0
        for ex,co in sorted(f.items()):
            z=dag.add(z,dag.prod([dag.c(int(co)),dag.power(qr,ex[0])]+[dag.power(r,k)for r,k in zip(rroots,ex[1:])]))
        return z
    witnessroots=[[uninode(f)for f in row]for row in powercofs]
    powersT=[[dag.power(x,k)for k in range(15)]for x in Troots]
    powersR=[[dag.power(x,k)for k in range(15)]for x in rroots]
    acc=[0]*9;terms=[];coefsum=0
    for ex in compositions(14,5):
        chosen=next(i for i,bound in enumerate(nilbounds)if ex[i]>=bound)
        coefficient=math.factorial(14)//math.prod(math.factorial(k)for k in ex);coefsum+=coefficient
        leftover=list(ex);leftover[chosen]-=nilbounds[chosen]
        root=dag.prod([dag.c(coefficient)]+[powersT[i][ex[i]]for i in range(5)]+[powersR[i][leftover[i]]for i in range(5)])
        for k in range(9):acc[k]=dag.add(acc[k],dag.mul(root,witnessroots[chosen][k]))
        terms.append([list(ex),chosen,coefficient,root])
    assert len(terms)==3060 and coefsum==5**14
    Hroots=[]
    for j in range(9):
        root=0
        for k in range(9):root=dag.add(root,dag.mul(acc[k],dag.c(A9[k][j])))
        Hroots.append(root)
    f_roots=[dag.c(x)for x in hvals];rhs=0
    for h,f in zip(Hroots,f_roots):rhs=dag.add(rhs,dag.mul(h,f))
    lhs=dag.c(a**29)
    sourceproof={'raw_q5L_cofactor_after_old_steps':enc(vf[9]),
                 'q5_m_coefficients':[enc(vf[i])for i in [5,6,7]],
                 'q0_n_coefficients':[enc(vf[i])for i in [1,2,3]],'q0_lambda':enc(lam),
                 'A_equals_Rm_plus_Pn_all_cubic_coefficients':True,
                 'B_equals_R_plus_ellm_plus_lambdaP_all_quadratic_coefficients':True,
                 'third_L_coefficient_recovery':'B2=v0*h22/9-2*h25',
                 'first_GL6_preserves_x0_subspace':True,
                 'second_GL3_NOT_assumed_to_preserve_x0':'Use explicit F=row9-c*row5+(c*d5-d9)*row4, its seven pure-y coefficients, and v.F_L=0 with v2=-9.'}
    certificate={'coefficient_domain':'FULLU; polynomial core identity over QQ[a,e,d], no physical C0 extension asserted',
                 'source_binding':sourceproof,
                 'P_and_R_and_ell_y':[[[list(mon),enc(co)]for mon,co in f.items()]for f in [P,RR,ell]],
                 'P_constant_symmetric2minor':'-25/16',
                 'canonical_coordinates':['S=-5*y1/2+(9e+4)*y2/5','T=y2','W=y3-S/5-(3e+8)*T/25'],
                 'q':encp(poly(q)),'P_equals_SW_plus_qT2':True,
                 'Rbar_normal_basis':['S2','ST','T2','TW','W2'],'r5':[encp(x)for x in rv],
                 'transported_pair_functional_t5':[enc(x)for x in tv],
                 'polynomial_T5_a2_times_t':[encp(x)for x in Tpoly],
                 'literal_sum_Ti_ri_equals_a2':True,
                 'degree4_normal_basis':basis4,
                 'universal_q_r_generator9':[encp(x)for x in Js],
                 'universal_power_bounds':nilbounds,
                 'universal_power_witness5x9':[ [encp(f)for f in row]for row in powercofs],
                 'fixed_rewrite_step_counts':rewritecounts,
                 'all5_universal_power_identities_literal_verified':True,
                 'a_J_equals_matrix_times_actual_h9':[[encp(x)for x in row]for row in A9],
                 'all9_source_membership_polynomial_identities_verified':True,
                 'actual_h9_original_target_labels':bc['query_only_original32_target_columns'],
                 'actual_h9':[encp(x)for x in hvals],
                 'scalar_DAG':{'atoms_polynomial_QQ_a_e_d':[encp(x)for x in dag.atoms],'nodes':dag.nodes},
                 'universal_power_witness_scalar_roots':witnessroots,
                 'multinomial3060_terms_exp_choice_integer_multiplier_root':terms,
                 'universal_J_cofactor_roots':acc,
                 'actual_h9_cofactor_roots':Hroots,
                 'identity_left_root_a29':lhs,'identity_right_root_sum_Hi_hi':rhs,
                 'identity_proof':'Each3060term uses a literal saved r_i^bound witness. Multinomial coefficient enumeration gives (sum Ti ri)^14=a28. Substitute literal aJ=A9*h and multiply once by a. Therefore sum Hi*hi=a29 without expanded Hi.',
                 'actual_9row_dual':'cofactor_root_i / a^29',
                 'polynomial_cofactor_degree_upper_bounds':[dag.degrees[x]for x in Hroots],
                 'no_target_source_transport_or_free_completion_executed':True}
    (out/'certificate.json').write_text(json.dumps(certificate,separators=(',',':'))+'\n')
    report={'id':ID,'all_actual_source_coordinate_pair_and_membership_checks_passed':True,
            'all5_explicit_universal_power_witnesses_passed':True,
            'fixed_pigeonhole_exponent':14,'exact_multinomial_term_count':3060,
            'literal_finite_DAG_Bezout_target':'a^29','actual_row_dual_denominator':'a^29',
            'all_core_atoms_polynomial_no_unknown_factor_inverted':True,
            'cofactor_polynomial_degree_upper_bound':max(dag.degrees[x]for x in Hroots),
            'scalar_nodes':len(dag.nodes),'core_polynomial_atoms':len(dag.atoms),
            'no_GB_point_degree_guess_Gaussian_or_automatic_transport':True,
            'no_full57_full239_or_physical_C0_extension_claim':True,
            'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},
            'certificate_sha256':sha(out/'certificate.json'),
            'seconds':time.monotonic()-began,'frozen':True,'independent_review_pending':True}
    assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');signal.alarm(0)
    print(json.dumps(report,indent=2),flush=True)


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True)
    main(ap.parse_args().output.resolve())
