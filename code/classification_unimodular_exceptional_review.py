#!/usr/bin/env python3
"""Independent Fraction reconstruction of the whole E3 principal-orbit exclusion."""
import hashlib,json,resource,time
from pathlib import Path
import classification_unimodular_logarithmic_review as r

def main():
    resource.setrlimit(resource.RLIMIT_CPU,(25,26));resource.setrlimit(resource.RLIMIT_AS,(768*1024**2,768*1024**2))
    start=time.monotonic();root=Path(__file__).resolve().parents[1]
    path=root/'reports/CLASSIFICATION_UNIMODULAR_EXCEPTIONAL_20261107.json';c=json.loads(path.read_text())
    X=[r.scale(r.atom(1,i),r.F(a,4)) for i,a in enumerate([-5,-1,3,7])]
    Y=[{},r.atom(4,0),r.atom(4,1),r.atom(4,2)]
    Q=[r.add(r.mul(Y[i],X[j]),r.scale(r.mul(X[i],Y[j]),-1)) for i,j in r.PAIRS]
    assert Q==list(map(r.parse,c['highest_bivector'])) and not r.pf(Q) and not any(r.divergence(Q))
    cols=[];mon3=set()
    for i in range(6):
        for j in range(4):
            L=[r.atom(1,j) if k==i else {} for k in range(6)];cross=r.cross(Q,L);mon3.update(cross)
            cols.append(([v.get(r.ZERO,r.RF()) for v in r.divergence(L)],cross))
    M=r.transpose([d+[p.get(e,r.RF()) for e in sorted(mon3,reverse=True)] for d,p in cols])
    assert M==r.matrix(c['full_linear_matrix'])
    md=r.det(r.submatrix(M,c['full_linear_pivot_rows'],c['full_linear_pivot_columns']))
    assert len(c['full_linear_pivot_rows'])==20 and md and md==r.scalar(c['full_linear_pivot_determinant'])
    T=r.transpose([r.linearvec([r.der(p,j) for p in Q]) for j in range(4)])
    assert T==r.matrix(c['translation_basis']) and r.annihilate(M,T)
    td=r.det(r.submatrix(T,c['translation_cut_rows'],c['translation_pivot_columns']))
    assert len(c['translation_cut_rows'])==4 and td and td==r.scalar(c['translation_pivot_determinant'])
    assert c['normal_dimension']==0 and not any(map(r.parse,c['normal_linear_bivector']))
    mon2=[]
    for i,j in r.combinations_with_replacement(range(4),2):
        e=[0]*r.N;e[i]+=1;e[j]+=1;mon2.append(tuple(e))
    G=r.transpose([[r.cross(Q,[r.atom() if k==i else {} for k in range(6)]).get(e,r.RF()) for e in mon2] for i in range(6)])
    assert G==r.matrix(c['constant_equation_matrix'])
    gd=r.det(r.submatrix(G,c['constant_pivot_rows'],c['constant_pivot_columns']))
    assert len(c['constant_pivot_rows'])==6 and gd and gd==r.scalar(c['constant_pivot_determinant'])
    assert not any(map(r.parse,c['constant_source'])) and not r.parse(c['constant_Pfaffian']) and c['all_remaining_positive_degree_Pfaffian_coefficients']==[]
    helper=Path(r.__file__)
    out={'task_id':'G','all_passed':True,'method':'Independent standard-library Fraction reconstruction from the two linear vector fields in Pym2014 Section3.6; no producer or CAS arithmetic imported.',
      'highest_Q_is_minus_X_wedge_Y':True,'highest_Pfaffian_and_divergence_zero':True,'full24_linear_matrix_rank':20,'full_kernel_equals_translation_image_dimension':4,'true_translation_normal_dimension':0,'full6_constant_matrix_rank':6,
      'linear_rank_minor':c['full_linear_pivot_determinant'],'translation_rank_minor':c['translation_pivot_determinant'],'constant_rank_minor':c['constant_pivot_determinant'],
      'conclusion':'Every affine quadratic unimodular extension over the principal E3 highest orbit with constant Pfaffian has Pfaffian zero. The prescribed nonzero-constant-Pfaffian support is empty.',
      'boundary_scope':'No orbit-closure extension is inferred; closed L112, pullback and common-factor boundary obligations remain separate.',
      'input_sha256':{str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [path,helper]},'reviewer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'seconds':time.monotonic()-start}
    (root/'reports/CLASSIFICATION_UNIMODULAR_EXCEPTIONAL_REVIEW_20261107.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'all_passed':True,'seconds':out['seconds'],'ranks':[20,4,6]}))
if __name__=='__main__':main()
