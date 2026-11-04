#!/usr/bin/env python3
"""Direct H12(d=0) source and fixed DVR-table proof, no producer or GB."""
import hashlib,importlib.util,itertools,json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='AMBIENT_B18_NPAIR_RESTRICTION_DVR_20261104_v63_REVIEW'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main(out):
    started=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));out.mkdir(exist_ok=False,parents=True)
    helper=ROOT/'code/ambient_B18_proper_rational_elimination_review_20261104_v63.py';spec=importlib.util.spec_from_file_location('own_sparse_QQ',helper);ar=importlib.util.module_from_spec(spec);spec.loader.exec_module(ar)
    add,mul,scale,power,same=ar.add,ar.mul,ar.scale,ar.power,ar.same
    pb=ROOT/'reports/AMBIENT_B18_NPAIR_RESTRICTION_DVR_20261104_v63';pr=read(pb/'report.json');assert pr['all_passed']
    paths={Path(__file__),helper,ROOT/'research/notes'/(ID+'.md'),pb/'report.json'}
    for rel,h in pr['input_sha256'].items():assert sha(ROOT/rel)==h,rel;paths.add(ROOT/rel)
    for rel,h in pr['evidence_sha256'].items():assert sha(pb/rel)==h,rel;paths.add(pb/rel)
    e=read(pb/'exact-DVR-and-source-certificate.json')
    sp=ROOT/'reports/AMBIENT_B18_PROPER_RATIONAL_ELIMINATION_20261104_v63/actual12-flag-divided-polynomials.json';source=read(sp);paths.add(sp)
    cp=ROOT/'reports/AMBIENT_B18_PROPER_NPAIR_RESTRICTION_20261104_v63/actual-Npair-three-variable-polynomials.json';conn=read(cp);paths.add(cp)
    peer=ROOT/'reports/AMBIENT_B18_PROPER_RATIONAL_ELIMINATION_20261104_v63_REVIEW/report.json';assert read(peer)['all_passed'];paths.add(peer)
    assert conn['specialized8_variable_order']==['a','v2','v3','v4','v5','v6','v7','v8']
    # The primary used recorded H12(d0); independently bind it to literal H12.
    H=[]
    for i,record in enumerate(source['all12_actual_remaining_and_divided_records']):
        p=ar.decode(record['divided_primitive_polynomial']);restricted={m[:1]+m[2:]:c for m,c in p.items()if m[1]==0}
        same(restricted,ar.decode(conn['all12_exact_H_d0_polynomials'][i]));H.append(restricted)
    a,x,t=[ar.var(3,i)for i in range(3)];one=ar.one(3)
    roots3={3:scale(mul(a,add(a,x)),Q(2,5)),
      2:scale(add(add(scale(mul(a,a),5),scale(mul(a,x),3)),add(scale(mul(a,t),-6),scale(mul(x,t),2))),Q(1,4)),
      4:add(add(scale(power(a,3),Q(1,2)),scale(mul(mul(a,a),x),Q(3,10))),add(scale(mul(mul(a,a),t),Q(-1,5)),scale(mul(mul(a,x),t),Q(1,5)))),
      6:mul(t,t),7:scale(power(t,3),Q(1,3))}
    names=conn['specialized8_variable_order'];baseindices=[0,1,5]
    def embed(p):
        ans={}
        for m,c in p.items():
            n=[0]*8
            for i,k in zip(baseindices,m):n[i]=k
            ans[tuple(n)]=c
        return ans
    roots={i:embed(p)for i,p in roots3.items()}
    for i,p in roots.items():same(p,ar.decode(conn['all_polynomial_roots'][names[i]]))
    def substitute(p,values):
        out={};cache={}
        for m,c in p.items():
            n=list(m)
            for i in values:n[i]=0
            term={tuple(n):c}
            for i,value in values.items():
                if m[i]:
                    key=i,m[i]
                    if key not in cache:cache[key]=power(value,m[i])
                    term=mul(term,cache[key])
            out=add(out,term)
        return out
    prior={};stages=[(0,3,10),(4,2,8),(1,4,80),(2,6,16),(3,7,192)]
    for row,i,unit in stages:
        reduced=substitute(H[row],prior);same(reduced,scale(add(ar.var(8,i),roots[i],-1),unit));prior[i]=roots[i]
    reduced=[substitute(p,roots)for p in H]
    for row in range(5):same(reduced[row],{})
    actual=[]
    u=ar.var(3,0)
    for k,p in enumerate(reduced[5:]):
        assert all(all(not m[i]for i in range(8)if i not in baseindices)for m in p)
        q={tuple(m[i]for i in baseindices):c for m,c in p.items()}
        rec=conn['seven_remaining_actual_polynomials'][k];assert rec['source_H_index']==k+5
        same(q,ar.decode(rec['exact_polynomial3']))
        # a=u-x, with x and t unchanged.
        ans={}
        for m,c in q.items():
            ans=add(ans,mul({(0,m[1],m[2]):c},power(add(u,x,-1),m[0])))
        same(ans,ar.decode(e['seven_source_polynomials'][k]));actual.append(ans)
    F=add(mul(u,u),scale(mul(u,x),6),-1)
    G=add(add(mul(add(scale(x,23),scale(t,-6)),u),scale(mul(x,x),2)),scale(mul(x,t),4))
    Hgen=mul(u,mul(x,x));gen=[F,G,Hgen]
    for p,q in zip(gen,e['three_generators']):same(p,ar.decode(q))
    same(actual[0],F);same(actual[1],add(scale(F,10),scale(G,2)))
    same(actual[2],add(scale(mul(add(u,scale(x,3)),F),Q(-6,5)),scale(Hgen,-24)))
    for p,cofs in zip(actual,e['seven_membership_cofactors']):
        ans={}
        for c,g in zip(cofs,gen):ans=add(ans,mul(ar.decode(c),g))
        same(ans,p)
    qpoly=mul(u,x)
    for name,z in [('u',u),('x',x),('t',t)]:
        ans={}
        for c,g in zip(e['annihilator_cofactors'][name],gen):ans=add(ans,mul(ar.decode(c),g))
        same(ans,mul(z,qpoly))
    # Own R=QQ[t] algebra on R^3 + k*q, table checked literally.
    Z={};O=one
    table=[[[{}for _ in range(4)]for _ in range(4)]for _ in range(4)]
    for i in range(4):table[0][i][i]=O;table[i][0][i]=O
    table[1][1]=[Z,scale(t,-2),scale(t,3),scale(O,Q(-23,2))]
    table[1][2]=table[2][1]=[Z,Z,Z,O];table[2][2]=[Z,Z,Z,scale(O,6)]
    def scalar_string(s):
        if s.endswith('*t'):return scale(t,Q(s[:-2]))
        return scale(O,Q(s))
    for i in range(4):
        for j in range(4):
            for k in range(4):same(table[i][j][k],scalar_string(e['multiplication_table'][i][j][k]))
    def norm(v):return v[:3]+[{ar.zero(3):v[3][ar.zero(3)]}if v[3].get(ar.zero(3),Q())else{}]
    def multiply(v,w):
        ans=[{}for _ in range(4)]
        for i in range(4):
            for j in range(4):
                co=mul(v[i],w[j])
                for k in range(4):ans[k]=add(ans[k],mul(co,table[i][j][k]))
        return norm(ans)
    basis=[[O if i==j else{}for j in range(4)]for i in range(4)]
    for i,j,k in itertools.product(range(4),repeat=3):assert multiply(multiply(basis[i],basis[j]),basis[k])==multiply(basis[i],multiply(basis[j],basis[k]))
    def evaluate(p):
        ans=[{}for _ in range(4)]
        for (eu,ex,et),c in p.items():
            word=basis[0]
            for _ in range(eu):word=multiply(word,basis[2])
            for _ in range(ex):word=multiply(word,basis[1])
            coeff={(0,0,et):c}
            ans=norm([add(a,mul(coeff,v))for a,v in zip(ans,word)])
        return ans
    for p in gen+actual:assert evaluate(p)==[{}]*4
    assert evaluate(qpoly)==basis[3]
    # The reverse algebra map follows from these source table identities.
    same(add(mul(u,u),scale(qpoly,6),-1),F)
    same(add(add(add(mul(x,x),scale(mul(t,u),-3)),scale(mul(t,x),2)),scale(qpoly,Q(23,2))),scale(G,Q(1,2)))
    saturated=[mul(u,u),qpoly,add(add(mul(x,x),scale(mul(t,u),-3)),scale(mul(t,x),2))]
    # Literal components: u=2x/3 modulo x², and u=0,x=-2t.
    def substitute3(p,images):
        ans={}
        for m,c in p.items():
            term=scale(O,c)
            for i,image in enumerate(images):
                if m[i]:term=mul(term,power(image,m[i])if image else{})
            ans=add(ans,term)
        return ans
    for p in saturated:
        val=substitute3(p,[scale(x,Q(2,3)),x,t]);assert all(m[1]>=2 for m in val)
        same(substitute3(p,[{},scale(t,-2),t]),{})
    gluing=[[O,Z,Z],[Z,O,scale(O,Q(2,3))],[O,scale(t,-2),Z]]
    for i in range(3):
        for j in range(3):same(gluing[i][j],scalar_string(e['saturated_primary_gluing_matrix'][i][j]))
    determinant={}
    for perm in itertools.permutations(range(3)):
        term=O
        for i,j in enumerate(perm):term=mul(term,gluing[i][j])
        determinant=add(determinant,term,-1 if sum(perm[i]>perm[j]for i in range(3)for j in range(i+1,3))%2 else 1)
    same(determinant,scale(t,Q(4,3)))
    proof={'all_passed':True,'direct_source_H12_d0':'All12 restrictions were recomputed from the literal actual12-flag-divided-polynomials.json, then compared to the Conn8-variable input before using its roots.',
      'five_constant_unit_triangular_roots_and_all7_source_remainders':True,'all7_membership_and_three_reverse_source_generators':True,
      'all3_q_annihilator_identities':True,'all64_associativity_checks':True,'two_sided_algebra_model':True,
      'module':'QQ[t]^3 direct_sum QQ[t]/(t), ordered basis1,x,u,q=ux',
      'inverse_maps':'The table algebra is associative and generated by u,x (q=ux), and kills allF,G,H. Conversely its table entries follow in the source quotient from F,G and the three checked u*q,x*q,t*q identities; the module relation t*q=0 is included. Generator-preserving maps in both directions are inverse.',
      'H0':'A finite-length element must be killed by a power of t, hence belongs to k*q since the other three summands are QQ[t]-free. Conversely q is killed by u,x,t, and is nonzero in the table. Thus H0 at the origin has length1 and depth0.',
      'Hilbert_series':'(1+2z)/(1-z)+z^2 from basisdegrees0,1,1,2; dimension1,multiplicity3',
      'saturation_and_components':'The t-torsion is exactly k*q, so saturation removes q and yields (u²,ux,x²-3tu+2tx). The two displayed component maps land in QQ[t,x]/x² (u=2x/3) and QQ[t] (u=0,x=-2t). Their map on free basis(1,x,u) has determinant4t/3!=0, hence is injective over QQ[t]. Thus the saturated ideal equals (3u-2x,x²) intersect (u,x+2t).',
      'gluing_det':'4t/3','actual_core_quotient_by_d_only':True,'no_E4_quotient_by_d_D0_or_full_D_depth_inference':True,
      'new_GB_sourcejet_point_prime_or_CAS_jobs':0,'counts':ar.COUNT}
    (out/'independent-proof.json').write_text(json.dumps(proof,indent=2)+'\n')
    report={'id':ID,'all_passed':True,'direct_actual_H12_d0_five_roots_all7_source_checked':True,'two_sided_DVR_model_H0_length1':True,'actual_core_restriction_only':True,
      'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in sorted(paths)},'evidence_sha256':{'independent-proof.json':sha(out/'independent-proof.json')},
      'seconds':time.monotonic()-started,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'frozen':True}
    assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024;signal.alarm(0)
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2),flush=True)
if __name__=='__main__':
    import argparse
    pa=argparse.ArgumentParser();pa.add_argument('--output',type=Path,default=ROOT/'reports'/ID);main(pa.parse_args().output.resolve())
