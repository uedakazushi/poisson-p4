#!/usr/bin/env python3
"""One exact factored134 substitution, coefficient blocks, no large field expansion."""
import argparse,gzip,hashlib,itertools,json,math,resource,signal,time
from collections import Counter,defaultdict
from fractions import Fraction as Q
from functools import lru_cache
from pathlib import Path
from q_new_family_v_zero_degree26_support_preflight_20261104_v63 import read,sha
from s61_sparse_generic5_input_preflight_20261104_v63 import add,mul,scale
ROOT=Path(__file__).resolve().parents[1]
ID='S61_SPARSE_GENERIC5_FULL134_NORMAL_CIRCUIT_20261104_v63'
ZERO=(0,0,0,0,0)
class Stop(Exception):pass
class Pool:
    def __init__(self):self.nodes=[];self.cache={};self.qq={};self.deg=[];self.truth=[];self.atoms=[];self.atomcache={};self.z=self.const(0);self.o=self.const(1)
    def node(self,op,args,degree=0,truth=None):
        key=(op,tuple(args))
        if key in self.cache:return self.cache[key]
        i=len(self.nodes);self.cache[key]=i;self.nodes.append([op,*args]);self.deg.append(degree);self.truth.append(truth);return i
    def const(self,x):
        x=Q(x);i=self.node('qq',[str(x)],0,bool(x));self.qq[i]=x;return i
    def prod(self,*args):
        c=Q(1);ids=[]
        for i in args:
            if i in self.qq:c*=self.qq[i]
            else:ids.append(i)
        if not c:return self.z
        if not ids:return self.const(c)
        if c!=1:ids.append(self.const(c))
        ids.sort()
        if len(ids)==1:return ids[0]
        tr=False if any(self.truth[i]is False for i in ids)else(True if all(self.truth[i]is True for i in ids)else None)
        return self.node('prod',[tuple(ids)],sum(self.deg[i]for i in ids),tr)
    def summ(self,args):
        terms=defaultdict(Q);constant=Q(0)
        for i in args:
            if i in self.qq:constant+=self.qq[i];continue
            c=Q(1);base=i
            if self.nodes[i][0]=='prod':
                rest=[]
                for j in self.nodes[i][1]:
                    if j in self.qq:c*=self.qq[j]
                    else:rest.append(j)
                base=self.prod(*rest)
            terms[base]+=c
        ids=[self.prod(self.const(c),i)for i,c in terms.items()if c]
        if constant:ids.append(self.const(constant))
        ids.sort()
        if not ids:return self.z
        if len(ids)==1:return ids[0]
        return self.node('sum',[tuple(ids)],max(self.deg[i]for i in ids),None)
    def atom(self,p):
        if not p:return self.z
        if set(p)=={ZERO}:return self.const(p[ZERO])
        den=math.lcm(*(c.denominator for c in p.values()));integers={m:c.numerator*(den//c.denominator)for m,c in p.items()};g=math.gcd(*(abs(v)for v in integers.values()));sg=1 if integers[min(integers)]>0 else -1;content=Q(g*sg,den);key=tuple((m,v//(g*sg))for m,v in sorted(integers.items()))
        if key not in self.atomcache:self.atomcache[key]=len(self.atoms);self.atoms.append([[str(v),*m]for m,v in key])
        i=self.node('parameter_poly_atom',[self.atomcache[key]],0,True);return self.prod(self.const(content),i)
    def inv(self,i):
        assert self.deg[i]==0 and self.truth[i]is not False
        if i in self.qq:return self.const(1/self.qq[i])
        return self.node('inv',[i],0,True)
def nativepoly(p):return {tuple(t[1:]):Q(p['content'])*int(t[0])for t in p['terms']}
def enc(p):return [[str(c),*e]for e,c in sorted(p.items())]
def stats(pool,refs):
    refs=list(refs);t=Counter(pool.truth[i]for i in refs);return dict(entries=len(refs),proven_zero=t[False],proven_nonzero=t[True],undetermined=t[None],nonzero_lower=t[True],nonzero_upper=len(refs)-t[False])
def universal_checks(piv,normal):
    # Scalar formal polynomials: tuple of E labels -> integer. No actual E expansion.
    def ap(a,b):
        c=dict(a)
        for m,v in b.items():
            c[m]=c.get(m,0)+v
            if not c[m]:del c[m]
        return c
    def mp(a,b):
        c={}
        for m,x in a.items():
            for n,y in b.items():k=tuple(sorted(m+n));c[k]=c.get(k,0)+x*y
        return {m:x for m,x in c.items()if x}
    one={():1};EF=lambda i,j:one if i in piv and piv.index(i)==j else({}if i in piv else{(3*i+j,):1})
    T=[]
    for i in range(45):
        row={j:EF(i,j)for j in range(3)if EF(i,j)}
        if i in normal:row[3+normal.index(i)]=one
        T.append(row)
    B=[]
    for i in range(3):B.append({piv[i]:one})
    for i in normal:B.append({i:one}|{piv[j]:{m:-v for m,v in EF(i,j).items()}for j in range(3)})
    for X,Y in[(T,B),(B,T)]:
        for i in range(45):
            row={}
            for k,a in X[i].items():
                for j,b in Y[k].items():row[j]=ap(row.get(j,{}),mp(a,b))
            assert {j:a for j,a in row.items()if a}=={i:one}
    pairs=list(itertools.combinations_with_replacement(range(45),2));checks=46
    for i,j in pairs:
        direct={}
        for a,x in T[i].items():
            for b,y in T[j].items():k=tuple(sorted((a,b)));direct[k]=ap(direct.get(k,{}),mp(x,y))
        built={}
        def put(k,v):built[k]=ap(built.get(k,{}),v)
        for a in range(3):
            for b in range(a,3):put((a,b),ap(mp(EF(i,a),EF(j,b)),mp(EF(i,b),EF(j,a)))if a!=b else mp(EF(i,a),EF(j,b)))
        if i in normal:
            for a in range(3):put((a,3+normal.index(i)),EF(j,a))
        if j in normal:
            for a in range(3):put((a,3+normal.index(j)),EF(i,a))
        if i in normal and j in normal:put(tuple(sorted((3+normal.index(i),3+normal.index(j)))),one)
        assert {k:v for k,v in direct.items()if v}=={k:v for k,v in built.items()if v};checks+=1
    return dict(all1081_qdegree_le2_basis_substitution_identities=True,basis_identities=checks,both45x45_inverse_products_identity=True,inverse_entries_checked=4050,E_pivot_I3_used=True,no_parameter_E_expansion_or_sample=True)
def main(out):
    start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));signal.signal(signal.SIGALRM,lambda *_:(_ for _ in()).throw(Stop('60s wall cap')));signal.alarm(60);out.mkdir(parents=True,exist_ok=False)
    source=ROOT/'reports/S61_SPARSE_GENERIC5_PLANE_RANK3_CIRCUIT_20261104_v63';pref=ROOT/'reports/S61_SPARSE_GENERIC5_INPUT_PREFLIGHT_20261104_v63';tail=ROOT/'reports/S61_SPARSE_GENERIC5_SMALL_TAIL_FIELD_20261104_v63';gate=read(source/'report.json');assert gate['all_passed']and gate['rank3_proved']and gate['normal42_adapter_SLP']
    for n,h in gate['input_sha256'].items():assert sha(ROOT/n)==h,n
    slp=read(source/'normal42-adapter-SLP.json');original=read(source/'plane-W-original.json.gz');rankproof=read(source/'rank-certificate.json.gz');pulled=read(pref/'pulled134-input.json.gz');piv=slp['pivot_q_rows'];normal=slp['normal_q_rows42'];assert len(normal)==42
    pool=Pool();memo={};leafbindings=dict(original_W_poly=str((source/'plane-W-original.json.gz').relative_to(ROOT)),rank_pivot_rf=str((source/'rank-certificate.json.gz').relative_to(ROOT)),selected_row_factor_poly=str((source/'rank-certificate.json.gz').relative_to(ROOT)))
    def imported(i):
        if i in memo:return memo[i]
        op,*a=slp['nodes'][i]
        if op=='qq':r=pool.const(a[0])
        elif op in['original_W_poly','rank_pivot_rf','selected_row_factor_poly']:
            if op=='original_W_poly':nonzero=bool(original['W45x3'][a[0]][a[1]])
            elif op=='rank_pivot_rf':nonzero=bool(rankproof['nonzero_pivot_factors'][a[0]]['numerator'])
            else:nonzero=bool(rankproof['selected_row_factors'][a[0]]['terms'])
            r=pool.node('external_parameter_'+op,a,0,nonzero)if nonzero else pool.z
        elif op=='add':r=pool.summ([imported(a[0]),imported(a[1])])
        elif op=='sub':r=pool.summ([imported(a[0]),pool.prod(pool.const(-1),imported(a[1]))])
        elif op=='mul':r=pool.prod(imported(a[0]),imported(a[1]))
        else:assert op=='inv';r=pool.inv(imported(a[0]))
        memo[i]=r;return r
    E=[[imported(i)for i in row]for row in slp['E45x3_nodes']]
    for a,i in enumerate(piv):
        for j in range(3):E[i][j]=pool.o if a==j else pool.z
    centers=list(map(nativepoly,original['center45']));C=[pool.atom(p)for p in centers];Z=[pool.node('z',[i],1,True)for i in range(3)];N=[pool.node('n',[i],1,True)for i in range(42)];qforward=[]
    for i in range(45):qforward.append(pool.summ([C[i],*[pool.prod(E[i][a],Z[a])for a in range(3)],*([N[normal.index(i)]]if i in normal else[])]))
    qinput=[pool.node('q_input',[i],1,True)for i in range(45)];zinverse=[pool.summ([qinput[i],pool.prod(pool.const(-1),C[i])])for i in piv];ninverse=[pool.summ([qinput[i],pool.prod(pool.const(-1),C[i]),*[pool.prod(pool.const(-1),E[i][a],zinverse[a])for a in range(3)]])for i in normal]
    # sigma uses the original N*T two directions and d as the third direction.
    commonP=read(tail/'shared-denominator-kernel-plane.json.gz')['common_denominator_P'];Pnode=pool.atom({tuple(t[1:]):Q(t[0])for t in commonP});gi=imported(slp['inverse_Gamma3_node']);sigmaCoefficients=[pool.prod(Pnode,imported(slp['adjugate3x3_nodes'][2][j]),gi)for j in range(3)];sigma=pool.summ([pool.prod(a,z)for a,z in zip(sigmaCoefficients,Z)]);sigma2=pool.prod(sigma,sigma);idempotents=[pool.summ([pool.o,sigma,pool.prod(pool.const(-2),sigma2)]),pool.prod(pool.const(Q(4,3)),pool.summ([sigma2,pool.prod(pool.const(-1),sigma)])),pool.prod(pool.const(Q(1,3)),pool.summ([pool.prod(pool.const(2),sigma2),sigma]))]
    equations=[]
    for encoded in pulled['equations']:
        row={}
        for m,e,c in encoded:row.setdefault(tuple(m),{})[tuple(e)]=Q(c)
        equations.append(row)
    @lru_cache(None)
    def qproduct(m):return pool.prod(*[qforward[i]for i in m])
    F=[]
    for row in equations:F.append(pool.summ([pool.prod(pool.atom(c),qproduct(m))for m,c in row.items()]))
    assert all(pool.deg[i]<=2 for i in F)
    proof=universal_checks(piv,normal);(out/'structural-proof.json').write_text(json.dumps(proof,indent=2)+'\n')
    def save(name,obj):
        b=json.dumps(obj,separators=(',',':')).encode()
        if len(b)>32*1024*1024:raise Stop(name+' uncompressed32MiB cap')
        (out/name).write_bytes(gzip.compress(b,mtime=0))
        if sum(p.stat().st_size for p in out.iterdir()if p.is_file())>32*1024*1024:raise Stop('total artifact32MiB cap')
    directn=len(pool.nodes);save('scalar-atoms-direct.json.gz',dict(variables=pulled['parameter_variables'],primitive_integer_atoms=pool.atoms));save('full134-direct-circuit.json.gz',dict(nodes=pool.nodes,scalar_atoms='scalar-atoms-direct.json.gz',external_leaf_files=leafbindings,E45x3=E,center45=C,forward_q45=qforward,inverse_z3=zinverse,inverse_n42=ninverse,F134=F,unknown_order=['z0','z1','z2',*[f'n{i}'for i in range(42)]],q_to_u_indices_zero_based=pulled['q_to_u_indices_zero_based'],pivot_q_rows=piv,normal_q_rows42=normal,sigma_coefficients=sigmaCoefficients,sigma_node=sigma,known5_idempotents_nodes=idempotents,known5_only_idempotents=True,all134_nonzero_by_inverse_affine_map=True,all_inverses_parameter_only=True,polynomial_unknown_degree_le2=True))
    error=None;blocks=[];pattern=None;shiftterms=0
    try:
        zz=list(itertools.combinations_with_replacement(range(3),2));roundtrips=0
        @lru_cache(None)
        def cc(i,j):return mul(centers[i],centers[j])
        for row in equations:
            c0=dict(row.get((),{}));L=[dict(row.get((i,),{}))for i in range(45)];Q2={m:c for m,c in row.items()if len(m)==2}
            for i in range(45):c0=add(c0,mul(row.get((i,),{}),centers[i]))
            for (i,j),c in Q2.items():c0=add(c0,mul(c,cc(i,j)));L[i]=add(L[i],mul(c,centers[j]));L[j]=add(L[j],mul(c,centers[i]))
            # Actual reverse-centering coefficient identity over QQ[thefiveparams].
            back0=dict(c0);backL=[dict(p)for p in L]
            for i in range(45):back0=add(back0,scale(mul(L[i],centers[i]),-1))
            for (i,j),c in Q2.items():back0=add(back0,mul(c,cc(i,j)));backL[i]=add(backL[i],scale(mul(c,centers[j]),-1));backL[j]=add(backL[j],scale(mul(c,centers[i]),-1))
            assert back0==row.get((),{})and all(p==row.get((i,),{})for i,p in enumerate(backL));roundtrips+=46
            assert max([len(c0),*[len(p)for p in L]])<=150000;shiftterms+=len(c0)+sum(map(len,L))
            linearZ=[[]for _ in range(3)];quadraticZ=[[]for _ in range(6)];cross=[[]for _ in range(126)];nn={};linearN=[pool.atom(L[i])for i in normal]
            for i in range(45):
                a=pool.atom(L[i])
                for j in range(3):linearZ[j].append(pool.prod(a,E[i][j]))
            for (i,j),c in Q2.items():
                a=pool.atom(c)
                for k,(h,l)in enumerate(zz):
                    quadraticZ[k].append(pool.prod(a,E[i][h],E[j][l]))
                    if h!=l:quadraticZ[k].append(pool.prod(a,E[i][l],E[j][h]))
                if i in normal:
                    for h in range(3):cross[42*h+normal.index(i)].append(pool.prod(a,E[j][h]))
                if j in normal:
                    for h in range(3):cross[42*h+normal.index(j)].append(pool.prod(a,E[i][h]))
                if i in normal and j in normal:nn[tuple(sorted((normal.index(i),normal.index(j))))]=a
            plane=[pool.atom(c0),*[pool.summ(p)for p in linearZ],*[pool.summ(p)for p in quadraticZ]];cross=[pool.summ(p)for p in cross];assert all(pool.deg[i]==0 for i in plane+cross+linearN+list(nn.values()))
            blocks.append(dict(plane_z_degree_le2=plane,linear_n42=linearN,cross_z_n126=cross,n2_sparse=[[i,j,a]for (i,j),a in sorted(nn.items())]))
        pattern=dict(n2=dict(entries=134*903,exact_nonzero=sum(len(b['n2_sparse'])for b in blocks),exact_zero=134*903-sum(len(b['n2_sparse'])for b in blocks)),linear_n=stats(pool,[i for b in blocks for i in b['linear_n42']]),cross_z_n=stats(pool,[i for b in blocks for i in b['cross_z_n126']]),plane=stats(pool,[i for b in blocks for i in b['plane_z_degree_le2']]),warning='Bounds use exact atom zeros/nonzeros and domain products; sums with unresolved field cancellations remain undetermined. No syntactic-nonzero conjecture.')
        save('scalar-atoms-full.json.gz',dict(variables=pulled['parameter_variables'],primitive_integer_atoms=pool.atoms));save('coefficient-block-circuit.json.gz',dict(prefix_source='full134-direct-circuit.json.gz',prefix_node_count=directn,appended_nodes=pool.nodes[directn:],scalar_atoms='scalar-atoms-full.json.gz',rows134=blocks,plane_order=['1','z0','z1','z2','z0^2','z0*z1','z0*z2','z1^2','z1*z2','z2^2'],cross_order='z_index*42+n_index',n2_order='i<=j',patterns=pattern,all_coefficients_parameter_only=True,actual_reverse_centering_QQ_identities=roundtrips,all1081_basis_assembly_identities=proof['all1081_qdegree_le2_basis_substitution_identities']))
    except Exception as exc:error=repr(exc)
    cost=dict(raw_cross_shape=[134,126],residual_cross_shape_after_plane_subtraction=[129,126],dense126_square_entries=126*126,plane_restriction_shape=[134,10],expected_generic_plane_restriction_rank5_from_existing_known3plus1plus1=True,
        fixed_plane5_and_cross126_minors_not_checked=True,construction='chooseanactualnonzero5minorofplanerestrictions;parameter-onlyrow subtraction gives129residuals; assemblecrossmatrixfromsaved134x126blocks',
        coordinate_invariance_scope='Forfixed full quadratic equation space and knownplane, changingaffineplane/normalframes or residualbasis multiplies crossmap by invertibleparameter matrices. This is not a new fixed126minor certificate on g=identity source.',
        known5_compatibility='sigma=P*row2(adj B_W)/GammaW*z; toric0,Q51=-1/2,Q52=1 by universal inverseidentity; idempotents are forknown5 quotient only, not fullunknownfiber.',
        old_sample_quintic_coefficients_NOT_reused=True,new_primitive_quintic_or_separator_not_computed=True,
        no_Macaulay_GB_Euler_J0_or_exhaustion=True)
    (out/'next-cost.json').write_text(json.dumps(cost,indent=2)+'\n')
    paths=[Path(__file__),ROOT/f'research/notes/{ID}.md',source/'report.json',source/'normal42-adapter-SLP.json',source/'plane-W-original.json.gz',source/'rank-certificate.json.gz',pref/'report.json',pref/'pulled134-input.json.gz',tail/'shared-denominator-kernel-plane.json.gz']
    report=dict(id=ID,all_passed=error is None,error=error,all134_direct_circuit_complete=True,coefficient_blocks_complete=len(blocks)==134,scalar_circuit_nodes=len(pool.nodes),scalar_atoms=len(pool.atoms),shifted_scalar_terms=shiftterms,
        all_inverses_parameter_only=True,all45_unknowns_retained=True,universal_inverse_and1081_assembly_proved=True,no_point_prime_GB_Macaulay_J0_or134_parameter_expansion=True,
        patterns=pattern,input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},evidence_sha256={p.name:sha(p)for p in out.iterdir()if p.is_file()},seconds=time.monotonic()-start)
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(dict(all_passed=report['all_passed'],blocks=report['coefficient_blocks_complete'],nodes=report['scalar_circuit_nodes'],atoms=report['scalar_atoms'],seconds=report['seconds'],error=error)))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args();main(a.output)
