#!/usr/bin/env python3
"""Independent Fraction(u) audit of all ten x^3y / x^2y^2 boundary sources.
Arithmetic reused only from the independent checker, never from producer/CAS.
"""
import classification_unimodular_logarithmic_review as arithmetic
from classification_unimodular_logarithmic_review import (
    F,RF,N,ZERO,PAIRS,atom,add,scale,mul,der,parse,scalar,matrix,transpose,
    rank,det,submatrix,annihilate,pf,cross,divergence,xycoeff,linearvec,roots_supported,
    combinations_with_replacement)
import hashlib,json,resource,time
from pathlib import Path

def main():
    resource.setrlimit(resource.RLIMIT_CPU,(110,111));resource.setrlimit(resource.RLIMIT_AS,(1024*1024**2,1024*1024**2))
    started=time.monotonic();root=Path(__file__).resolve().parents[1]
    path=root/'reports/CLASSIFICATION_UNIMODULAR_LOGARITHMIC_BOUNDARY_20261107.json'
    data=json.loads(path.read_text());records=[];coverage={}
    for c in data['cases']:
        begin=time.monotonic();arithmetic.U=RF((F(0),F(1))) if c['residue']=='generic u' else RF(F(c['residue']))
        U=RF(3 if c['name'].startswith('x3y') else 1)
        q=parse(c['quadric']);x,y=atom(1,0),atom(1,1);dq=[der(q,i) for i in range(4)]
        gamma=[add(scale(q,(U-1)*int(i==0 and j==1)),scale(mul(y,add(scale(dq[j],int(i==0)),scale(dq[i],-int(j==0)))),U),mul(x,add(scale(dq[j],int(i==1)),scale(dq[i],-int(j==1))))) for i,j in PAIRS]
        gamma[0]=add(gamma[0],mul(x,x if c['name'].startswith('x3y') else y))
        Q=[gamma[5],scale(gamma[4],-1),gamma[3],gamma[2],scale(gamma[1],-1),gamma[0]]
        assert Q==list(map(parse,c['highest_bivector'])) and not pf(Q) and not any(divergence(Q))
        columns=[];mon3=set()
        for i in range(6):
            for j in range(4):
                L=[atom(1,j) if k==i else {} for k in range(6)];cr=cross(Q,L);mon3.update(cr)
                columns.append(([v.get(ZERO,RF()) for v in divergence(L)],cr))
        mon3=sorted(mon3,reverse=True)
        M=transpose([d+[cr.get(e,RF()) for e in mon3] for d,cr in columns])
        assert M==matrix(c['full_linear_matrix'])
        K=matrix(c['full_kernel']);r=c['rank'];assert annihilate(M,K) and rank(K)==24-r
        pivot=det(submatrix(M,c['pivot_rows'],c['pivot_columns']));assert pivot and pivot==scalar(c['pivot_determinant']) and len(c['pivot_rows'])==r
        T=transpose([linearvec([der(p,j) for p in Q]) for j in range(4)])
        assert T==matrix(c['translation_basis']) and annihilate(M,T)
        tr=c['translation_rank'];assert rank(T)==tr
        trrows=c['translation_cut_rows'];trcols=[]
        for j in range(4):
            if rank(submatrix(T,trrows,trcols+[j]))>len(trcols):trcols.append(j)
        trminor=det(submatrix(T,trrows,trcols));assert trminor==scalar(c['translation_minor']) and trminor
        Ln=list(map(parse,c['normal_linear_bivector']));nd=c['normal_dimension']
        NK=transpose([linearvec([der(p,4+i) for p in Ln]) for i in range(nd)])
        assert nd+tr==24-r and annihilate(M,NK) and rank(NK)==nd
        assert all(not v for i in trrows for v in NK[i])
        assert rank([a+b for a,b in zip(T,NK)])==24-r
        mon2=[]
        for i,j in combinations_with_replacement(range(4),2):
            e=[0]*N;e[i]+=1;e[j]+=1;mon2.append(tuple(e))
        G=transpose([[cross(Q,[atom() if k==i else {} for k in range(6)]).get(e,RF()) for e in mon2] for i in range(6)])
        assert G==matrix(c['constant_equation_matrix'])
        gr,gc=c['constant_pivot_rows'],c['constant_pivot_columns'];gm=det(submatrix(G,gr,gc))
        assert gm and gm==scalar(c['constant_minor']) and rank(G)==len(gc)
        C=list(map(parse,c['constant_source']));free=[i for i in range(6) if i not in gc]
        assert all(C[i]==atom(1,9+j) for j,i in enumerate(free))
        Pf=pf([add(a,b,d) for a,b,d in zip(Q,Ln,C)]);coeff=xycoeff(Pf)
        assert all(mon2[i][:4] not in coeff for i in gr)
        actual={tuple(sorted(p.items())) for e,p in coeff.items() if sum(e)>0}
        expected={tuple(sorted(parse(p).items())) for p in c['remaining_full_Pfaffian_equations']}
        assert actual==expected and coeff.get((0,)*4,{})==parse(c['constant_Pfaffian'])
        if c['residue']=='generic u':
            allowed={F(d['residue']) for d in data['cases'] if d['name']==c['name'] and d['residue']!='generic u'}
            found=set()
            for v in [pivot,trminor,gm]:found |= roots_supported(v.p,allowed)|roots_supported(v.q,allowed)
            for row in K+NK:
                for v in row:found |= roots_supported(v.q,allowed)
            for p in C:
                for v in p.values():found |= roots_supported(v.q,allowed)
            assert found==allowed
            coverage[c['name']]=list(map(str,sorted(found)))
        records.append({'name':c['name'],'residue':c['residue'],'full_linear_rank':r,'translation_rank':tr,'normal_dimension':nd,'constant_rank':len(gc),'all_original_coefficient_matrices_and_complete_kernels':True,'true_translation_slice_and_inverse_constant_solve':True,'all_positive_degree_Pfaffian_and_constant_coefficients':True,'all_passed':True,'seconds':time.monotonic()-begin})
        print(c['name'],c['residue'],'PASS',round(time.monotonic()-begin,3),flush=True)
    assert len(records)==10 and len(coverage)==8
    result={'task_id':'G','all_passed':True,'method':'Independent standard-library Fraction(u) and sparse multivariate polynomial arithmetic; no producer or CAS imported. Whole generic identities and every exceptional residue, not point interpolation.','cases':records,'exceptional_coverage':coverage,'input_sha256':{str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [path,Path(arithmetic.__file__)]},'reviewer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'scope':'All ten boundary coefficient certificates. Partial-fraction geometric completeness and reduced-support implications are separate proof dependencies.','seconds':time.monotonic()-started}
    (root/'reports/CLASSIFICATION_UNIMODULAR_LOGARITHMIC_BOUNDARY_REVIEW_20261107.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'all_passed':True,'cases':10,'seconds':result['seconds']}))
if __name__=='__main__':main()
