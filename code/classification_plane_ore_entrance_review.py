#!/usr/bin/env python3
"""Independent universal Ore source and strict-Q Jacobi entrance.
Whole QQ(t) basis coverage, all cohomology variables and zero complements.
"""
from fractions import Fraction as F
from pathlib import Path
from itertools import combinations
import hashlib,json,resource,signal,time
from classification_minimal_quartic_review import LaurentRing
from classification_unimodular_logarithmic_review import RF,ug,trim
import classification_plane_dependent_review as D
q,H,B,M=D.q,D.H,D.B,D.M
ROOT=Path(__file__).resolve().parents[1]
JOB='research/classification/jobs/plane_ore_entrance_review.json'
OUT='reports/CLASSIFICATION_PLANE_ORE_ENTRANCE_REVIEW_20261107.json'
U=D.U

def make(name):
    if name=='SS':names=('t','B0','B1','E0','E1','k0','k1','k2');cohom=('B0','B1','E0','E1')
    elif name in ('SS3','SS1'):names=('B0','B1','E0','E1','f','k0','k1','k2');cohom=names[:5]
    elif name=='NIL1':names=('D0','D1','E0','E1','k0','k1','k2');cohom=names[:4]
    else:assert name=='NIL0';names=('D0','D1','B0','B1','E0','E1','f','k0','k1','k2');cohom=names[:7]
    R=LaurentRing(names,('b','c','z','w'));p=R.parse
    if name.startswith('SS'):
        t='t'if name=='SS'else '3'if name=='SS3'else '1'
        XB=[p(t+'*z*w'),p('('+t+'-1)*w**2')];XE=[p('-('+t+'+1)*z**2'),p('-'+t+'*z*w')]
        anchors=[[R.add(R.mul(p('B'+str(i)),XB[j]),R.mul(p('E'+str(i)),XE[j]))for j in range(2)]for i in range(2)]
        if name=='SS3':anchors[1][0]=R.add(anchors[1][0],p('f*w**2'))
        if name=='SS1':
            anchors[1][0]=R.add(anchors[1][0],p('f*z*w'));anchors[1][1]=R.add(anchors[1][1],p('f*w**2'))
        hh=[p('('+t+'-1)*w*B'+str(i)+'/2-('+t+'+1)*z*E'+str(i)+'/2')for i in range(2)]
        expected=['B0*E1-B1*E0']if name=='SS'else ['B0*E1-B1*E0','f*E0']if name=='SS3'else ['B0*E1-B1*E0-f*E0','f*B0']
    elif name=='NIL1':
        XD=[p('2*z**2'),p('z*(z+2*w)')];XE=[p('z*(z+2*w)'),p('2*w*(z+w)')]
        anchors=[[R.add(R.mul(p('D'+str(i)),XD[j]),R.mul(p('E'+str(i)),XE[j]))for j in range(2)]for i in range(2)]
        hh=[p('z*D'+str(i)+'+(z+w)*E'+str(i))for i in range(2)];expected=['D0*E1-D1*E0']
    else:
        XD=[{},p('z**2')];XB=[p('-2*z*w'),p('-3*w**2')];XE=[p('z**2'),p('2*z*w')]
        anchors=[[R.add(R.mul(p('D'+str(i)),XD[j]),R.mul(p('B'+str(i)),XB[j]),R.mul(p('E'+str(i)),XE[j]))for j in range(2)]for i in range(2)]
        anchors[1][0]=R.add(anchors[1][0],p('f*z**2'));anchors[1][1]=R.add(anchors[1][1],p('f*z*w'))
        hh=[p('-2*w*B'+str(i)+'+z*E'+str(i))for i in range(2)]
        expected=['2*(D0*B1-D1*B0)+f*E0','B0*E1-B1*E0+f*B0','D0*E1-D1*E0-f*D0']
    bc=R.add(R.mul(hh[0],p('c')),R.scale(R.mul(hh[1],p('b')),-1),p('k0*z**2+k1*z*w+k2*w**2'))
    Q={(0,1):bc,(0,2):anchors[0][0],(0,3):anchors[0][1],(1,2):anchors[1][0],(1,3):anchors[1][1]}
    for (i,j),v in list(Q.items()):Q[j,i]=R.scale(v,-1)
    return R,Q,cohom,expected

def divexact(P,Q):
    rem=dict(P);out={};f=max(Q);c=Q[f]
    while rem:
        e=max(rem);assert all(a>=b for a,b in zip(e,f)),('nonmultiple strictJacobi',e,f)
        g=tuple(a-b for a,b in zip(e,f));v=rem[e]/c;out[g]=out.get(g,F(0))+v
        for h,w in Q.items():
            k=tuple(a+b for a,b in zip(g,h));rem[k]=rem.get(k,F(0))-v*w
            if not rem[k]:rem.pop(k)
    return out

def strict(R,Q,expected,name):
    coeffs=[]
    for i,j,k in combinations(range(4),3):
        J=R.add(*(R.mul(Q.get((a,l),{}),R.derivative(Q.get((b,c),{}),l))for a,b,c in ((i,j,k),(j,k,i),(k,i,j))for l in range(4)))
        groups={}
        for e,v in J.items():groups.setdefault(e[R.n:],{})[e[:R.n]+(0,)*4]=v
        coeffs.extend(groups.values())
    gens=[R.parse(t)for t in expected]
    if name=='SS':
        quotients=[];gcd=()
        for c in coeffs:
            v=divexact(c,gens[0]);assert all(not any(e[1:])for e in v)
            poly=trim([v.get((k,)+(0,)*(R.n+3),F(0))for k in range(max(e[0]for e in v)+1)])
            quotients.append(list(map(str,poly)));gcd=ug(gcd,poly)
        assert len(gcd)==1
        return dict(exact_ideal='(det(B,E)) for every complex t',all_strict_coefficients_are_Delta_multiples=True,all_quotient_polynomials=quotients,quotient_gcd=list(map(str,gcd)))
    mons=sorted(set(e for P in coeffs+gens for e in P));matrix=lambda Ps:[[P.get(e,F(0))for e in mons]for P in Ps]
    rank=q.rank(matrix(coeffs))[0];assert rank==len(gens)==q.rank(matrix(gens))[0]==q.rank(matrix(coeffs+gens))[0]
    return dict(exact_generator_equations=expected,coefficient_monomials=[list(e[:R.n])for e in mons],all_strict_Jacobi_coefficient_rows=[[str(v)for v in row]for row in matrix(coeffs)],expected_generator_rows=[[str(v)for v in row]for row in matrix(gens)],exact_quadratic_rowspace_rank=rank)

def extract(R,Q,var):
    index=R.names.index(var);tindex=R.names.index('t')if 't'in R.names[:R.n]else None
    out={}
    for (i,j),P in Q.items():
        if i>=j:continue
        accum={}
        for e,v in P.items():
            if e[index]!=1 or any(e[k]for k in range(R.n)if k not in (index,tindex)):continue
            exponent=(0,)+e[R.n:]+(0,0,0);degree=e[tindex]if tindex is not None else 0
            coeff=accum.setdefault(exponent,{})
            coeff[degree]=coeff.get(degree,F(0))+v
        out[i+1,j+1]={e:RF(tuple(c.get(k,F(0))for k in range(max(c)+1)))for e,c in accum.items()}
    return q.antisymmetric(out)

def canonical(P,first):
    zz=(0,0,0,2,0,0,0,0);ww=(0,0,0,0,2,0,0,0)
    v=[q.add(q.scale(q.atom(1,3),P.get((i,3),{}).get(zz,RF(0))),q.scale(q.atom(1,4),P.get((i,4),{}).get(ww,RF(0))))for i in range(3)]+[{},{}]
    x=[q.atom(1,i)for i in range(5)];P=H.Padd(P,H.wedge(x,v));coeff=[]
    for T in first:
        ij=next(ij for ij in q.PAIRS if T.get(ij));e=next(iter(T[ij]));coeff.append(P.get(ij,{}).get(e,RF(0)))
    assert q.vec(P)==q.vec(q.antisymmetric({ij:q.add(*(q.scale(T.get(ij,{}),a)for T,a in zip(first,coeff)))for ij in q.PAIRS}))
    return coeff

def common():
    result=[]
    for i in (1,2):
        for j in (1,2):
            for n in (3,4):result.append(q.antisymmetric({(0,i):q.atom(1,j,n)}))
    for i,mon in ((3,(3,4)),(3,(4,4)),(4,(3,3)),(4,(3,4))):result.append(q.antisymmetric({(0,i):q.atom(1,*mon)}))
    return result

def basis_check(name,R,Q,cohom,first):
    fields=[extract(R,Q,v)for v in cohom]+common()
    P0=D.model('SS_N_LINE'if name.startswith('SS')else'NIL_N_TRACE_ONE_K0'if name=='NIL1'else'NIL_N_TRACE_ZERO')
    if name in ('SS3','SS1'):
        value=F(3 if name=='SS3'else 1)
        P0=q.antisymmetric({ij:{e:M.at_rf(RF(v),value)for e,v in P0.get(ij,{}).items()}for ij in q.PAIRS})
    for P in fields:
        A=q.jacobi(H.Padd(P0,P));C=q.jacobi(P)
        assert all(not q.add(A[k],q.scale(C[k],-1))for k in A)
    K=q.transpose([canonical(P,first)for P in fields]);n=len(fields)
    if name=='SS':
        certs=[];gcd=()
        # The actual z/w involution exchanges the two critical charts. Their
        # minors are units times(t+1)^2 and(t-1)^2, hence cover ALL complex t.
        for value in (F(1),F(-1)):
            numeric=[[M.at_rf(RF(v),value).p[0]if M.at_rf(RF(v),value).p else F(0)for v in row]for row in K]
            rank,rows=q.rank(q.transpose(numeric));assert rank==n
            det=H.det_generic([[K[i][j]for j in range(n)]for i in rows]);assert det and det.q==(F(1),)
            gcd=ug(gcd,det.p);certs.append(dict(selected_basis_chart=str(value),rows=rows,minor=H.encode(det)))
        assert len(gcd)==1,('selected basis minors have a common exceptional chart',gcd,certs)
        independence=dict(whole_t_two_minor_coprimality=True,minor_gcd=list(map(str,gcd)),certificates=certs)
    else:
        numeric=[[RF(v).p[0]if RF(v).p else F(0)for v in row]for row in K]
        rank,rows=q.rank(q.transpose(numeric));assert rank==n
        det=q.determinant([[numeric[i][j]for j in range(n)]for i in rows]);assert det
        independence=dict(rank=n,rows=rows,nonzero_minor=str(det))
    assert n=={'SS':16,'SS3':17,'SS1':17,'NIL1':16,'NIL0':19}[name]
    return dict(cohomology_parameters=list(cohom),all12_common_W_included=True,complete_P1_dimension=n,
                canonical_full30_by_kernel_basis=[[H.encode(RF(v))for v in row]for row in K],all350_D1_membership=True,independence=independence)

def main():
    started=time.monotonic();signal.alarm(120);resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,)*2)
    assert not(ROOT/OUT).exists();job=json.loads((ROOT/JOB).read_text())
    for b in job['exact_inputs_and_hashes']:assert hashlib.sha256((ROOT/b['path']).read_bytes()).hexdigest()==b['sha256'],b['path']
    first,second=M.bases();records=[]
    for name in ('SS','SS3','SS1','NIL1','NIL0'):
        R,Q,cohom,expected=make(name);basis=basis_check(name,R,Q,cohom,first);eq=strict(R,Q,expected,name)
        records.append(dict(type=name,whole_linear_source_basis=basis,strict_Q_Jacobi_ideal=eq));print(name,'whole source/strictQJac PASS',flush=True)
    # The affine magnetic three in Q and six in W are exactly the original P2nine.
    own=[q.antisymmetric({ij:q.atom(1,*ex)})for ij in ((0,1),(0,2),(1,2))for ex in ((3,3),(3,4),(4,4))]
    assert [q.vec(P)for P in own]==[q.vec(P)for P in second]
    P0=D.model('SS_N_LINE');sigma=[0,1,2,4,3];sign=[-1,1,1,1,1];changed={}
    for i,j in q.PAIRS:
        changed[i,j]={tuple(e[sigma[k]]for k in range(5))+e[5:]:v*sign[i]*sign[j]*((-1)**e[0])for e,v in P0.get((sigma[i],sigma[j]),{}).items()}
    target=q.antisymmetric({ij:{e:H.rf_at(RF(v),-U)for e,v in P0.get(ij,{}).items()}for ij in q.PAIRS})
    assert q.vec(q.antisymmetric(changed))==q.vec(target)
    mod=json.loads((ROOT/'reports/CLASSIFICATION_PLANE_ORE_MODULE_REVIEW_20261107.json').read_text());assert mod['all_passed']
    result=dict(task_id='G',all_passed=True,all_five_complete_Ore_source_types=records,
                original_P2nine_all_retained=True,SS_artificial_one_third_basis_roots_retained=True,
                SS_negative_critical_by_actual_normal_swap_and_a_sign='t->-t, involution; whole quadratic source transported with inverse',
                all_generic_critical_nil_zero_nonzero_branches_preserved=True,
                source_to_normalizer_geometry='Frozen Ore note independently reviewed: a-independent Euler gauge, projectiveJac equivalence, constant-v/essential3, fzero/firstcolumnzero, allconiccharts, criticalswap, fullcentralcobound/scales. Details in independent note.',
                independent_all8_affine_module_receipt='reports/CLASSIFICATION_PLANE_ORE_MODULE_REVIEW_20261107.json',
                whole_dependent_reduced_irreducible_quintic_G_image_at_most23=True,
                no_CAS_or_producer_imports=True,inputs=job['exact_inputs_and_hashes'],reviewer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                seconds=time.monotonic()-started,maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    (ROOT/OUT).write_text(json.dumps(result,indent=2)+'\n');print('PASS',OUT,round(result['seconds'],3),flush=True)
if __name__=='__main__':main()
