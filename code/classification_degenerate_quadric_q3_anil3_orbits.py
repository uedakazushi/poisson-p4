#!/usr/bin/env python3
"""Full ambient tangent and all true normals of the three proved whole orbits."""
import hashlib,itertools,json,resource,time
from pathlib import Path
import sympy as s
from sympy.polys.matrices import DomainMatrix
from classification_three_hyperplanes_quadric_independent import ambient_matrices
from classification_rank3_mixed_component_review import COLUMNS,minor_certificate,jacobi,sparse
ROOT=Path(__file__).resolve().parents[1]
def dm(M):return DomainMatrix.from_Matrix(M).convert_to(s.QQ)
def encode(P):return {str(i)+str(j):str(s.expand(P[i,j])) for i,j in itertools.combinations(range(5),2) if P[i,j]!=0}
def tensor(v,xx):
    P=s.zeros(5)
    for (ij,e),c in zip(COLUMNS,v):P[ij]+=c*s.prod(x**k for x,k in zip(xx,e))
    return P-P.T

def main():
    resource.setrlimit(resource.RLIMIT_CPU,(80,81));resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2))
    start=time.monotonic();path=ROOT/'reports/CLASSIFICATION_DEGENERATE_QUADRIC_Q3_ANIL3_SUPPORT_20261107.json';data=json.loads(path.read_text());xx=s.symbols('x0:5');records={}
    weight_lists={'upper':[4,3,2,1,0],'Czero_nonzero_I':[3,2,1,0,2],'Czero_nonzero_II':[3,2,1,0,5]}
    for name,source in data['families'].items():
        begin=time.monotonic();P=s.zeros(5)
        for key,value in source['canonical_brackets'].items():i,j=map(int,key);P[i,j]=s.sympify(value)
        P-=P.T;den=s.ilcm(1,*[s.denom(v) for i,j in itertools.combinations(range(5),2) for v in s.Poly(P[i,j],*xx).coeffs()]);P*=den
        D,V,G,F=ambient_matrices(P,xx,[P])
        rr,piv=dm(D).rref();rr=rr.to_Matrix();g_rank=len(dm(G).rref()[1]);v_rank=len(dm(V).rref()[1])
        kernel=[]
        for j in range(150):
            if j in piv:continue
            vec=s.zeros(150,1);vec[j]=1
            for i,k in enumerate(piv):vec[k]=-rr[i,j]
            vec*=s.ilcm(1,*[s.denom(v) for v in vec]);kernel.append(vec)
        K=s.Matrix.hstack(*kernel);selected=dm(G.row_join(K)).rref()[1]
        normal_indices=[i-G.cols for i in selected if i>=G.cols]
        assert D*K==s.zeros(350,K.cols)
        weights=weight_lists[name]
        relative=[weights[i]+weights[j]-sum(w*k for w,k in zip(weights,e))+1 for (i,j),e in COLUMNS]
        normals=[tensor(K[:,i],xx) for i in normal_indices]
        normalweights=[sorted({relative[j] for j,v in enumerate(K[:,i]) if v}) for i in normal_indices]
        assert all(len(w)==1 for w in normalweights)
        assert all(weights[i]+weights[j]-sum(w*k for w,k in zip(weights,e))==-1 for (i,j),terms in sparse(P,xx).items() for e,v in terms.items() if v)
        record={'source_orbit':name,'integer_scale':str(den),'canonical_brackets':encode(P),'coordinate_weights':weights,'P_weight':-1,
          'full_D_rank':len(piv),'vertical_rank':v_rank,'actual_orbit_vertical_scale_rank':g_rank,'full_kernel_dimension':K.cols,
          'true_normal_dimension':len(normal_indices),'normal_kernel_indices':normal_indices,'normal_brackets':[encode(R) for R in normals],
          'normal_relative_weights':[w[0] for w in normalweights],'relative_ambient_weight_range':[min(relative),max(relative)],
          'normal_itself_Poisson':[jacobi(sparse(R,xx))==s.zeros(350,1) for R in normals],
          'full_D_and_actual_group_minors':{'D':minor_certificate(D,len(piv)),'G':minor_certificate(G,g_rank)},
          'full_D_kernel_zero':True,'actual_group_plus_normals_full_kernel':len(selected)==K.cols,
          'projective_orbit_dimension':g_rank-25,'projective_tangent_dimension':125-len(piv),'seconds':time.monotonic()-begin}
        assert record['actual_group_plus_normals_full_kernel'];records[name]=record
        print(json.dumps({k:record[k] for k in ['source_orbit','full_D_rank','actual_orbit_vertical_scale_rank','normal_relative_weights','normal_itself_Poisson','seconds']}),flush=True)
    result={'task_id':'E_reduced_nonnormal','all_passed':True,'source_support_report':str(path.relative_to(ROOT)),'input_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'cases':records,'producer_source':Path(__file__).read_text(),'seconds':time.monotonic()-start,'scope':'Full ambient tangent of the three complete source orbits, not componenthood or membership by itself.'}
    out=ROOT/'reports/CLASSIFICATION_DEGENERATE_QUADRIC_Q3_ANIL3_ORBITS_20261107.json'
    with out.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
if __name__=='__main__':main()
