#!/usr/bin/env python3
"""One full-ambient witness for each of three proved surviving q4 nil sources."""
import argparse,hashlib,itertools,json,resource,time
from pathlib import Path
import sympy as s
from sympy.polys.matrices import DomainMatrix
from classification_three_hyperplanes_quadric_independent import ambient_matrices
from classification_rank3_mixed_component_review import minor_certificate

ROOT=Path(__file__).resolve().parents[1]
def decode(data,local):
    P=s.zeros(5)
    for ij in itertools.combinations(range(5),2):
        P[ij]=s.sympify(data.get(''.join(map(str,ij)),'0'),locals=local);P[ij[::-1]]=-P[ij]
    return P
def encode(A):return {str(ij):str(s.expand(A[ij])) for ij in itertools.combinations(range(5),2) if A[ij]!=0}
def qr(A):return len(DomainMatrix.from_Matrix(A).convert_to(s.QQ).rref()[1])

def main(output):
    resource.setrlimit(resource.RLIMIT_CPU,(55,56));resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2))
    begun=time.monotonic();xx=s.symbols('x0:5');nn=s.symbols('n0:15');local={str(v):v for v in (*xx,*nn)}
    support_path='reports/CLASSIFICATION_DEGENERATE_QUADRIC_Q4_NIL_SUPPORT_20261107.json'
    support=json.loads((ROOT/support_path).read_text());records=[]
    for label,kind,which in [('nil31_long_I','NIL31_LONG','I'),('nil31_long_II','NIL31_LONG','II'),('nil22_nonzero','NIL22_NONZERO',None)]:
        start=time.monotonic();path='reports/CLASSIFICATION_DEGENERATE_QUADRIC_Q4_'+kind+'_20261107.json'
        source=json.loads((ROOT/path).read_text());record=next(c for c in support['cases'] if c['case']==kind)
        if which:
            spec=next(c for c in record['surviving_families'] if c['family']==which and c['sign']==1)
            sub={local[k]:s.sympify(v,locals=local) for k,v in spec['substitution'].items()}
        else:sub={local[k]:s.sympify(v,locals=local) for k,v in record['surviving_source_substitution'].items()}
        P=decode(source['brackets'],local).subs(sub,simultaneous=True).applyfunc(s.expand)
        shear_certificate=None
        if which:
            # Prove the entire two-parameter source is one projective orbit before choosing the zero witness.
            k=s.Rational(5,2) if which=='I' else s.Integer(-5)
            G=s.eye(5);G[4,0]=-nn[0]/k;G[4,3]=-nn[1]/k
            inv=G.inv()*s.Matrix(xx);P0=P.subs({nn[0]:0,nn[1]:0})
            pushed=G*P0.subs(dict(zip(xx,inv)),simultaneous=True)*G.T
            assert (pushed-P).applyfunc(s.expand)==s.zeros(5)
            shear_certificate={'kappa':str(k),'coordinate_change_x4':'x4-n0*x0/kappa-n1*x3/kappa',
                'full_source_equals_pushforward_of_zero_parameter_member':True}
            # A constant complex linear change and a reparameterization give a rational full source.
            v=s.Symbol('v');T=s.diag(1,1,1,s.I,1)
            P=(T*P.subs({xx[3]:-s.I*xx[3],nn[1]:s.I*v},simultaneous=True)*T.T).applyfunc(s.expand)
            parameters=(nn[0],v);witness={nn[0]:0,v:0}
        else:
            parameters=tuple(nn[i] for i in [0,2,3,6]);witness={p:s.Integer(1) for p in parameters}
        scale=s.ilcm(1,*[s.denom(a) for ij in itertools.combinations(range(5),2) for a in s.Poly(P[ij],*xx,*parameters).coeffs()])
        P=P*scale
        assert all(c.is_Rational for ij in itertools.combinations(range(5),2) for c in s.Poly(P[ij],*xx,*parameters).coeffs())
        assert all(s.expand(sum(P[i,l]*s.diff(P[j,k],xx[l])+P[j,l]*s.diff(P[k,i],xx[l])+P[k,l]*s.diff(P[i,j],xx[l]) for l in range(5)))==0 for i,j,k in itertools.combinations(range(5),3))
        point=P.subs(witness);directions=[P.diff(p).subs(witness) for p in parameters]
        D,V,M,F=ambient_matrices(point,xx,directions)
        matrices={'full350_D':D,'vertical':V,'orbit_vertical_scale':M,'actual_full_source':F}
        ranks={key:qr(A) for key,A in matrices.items()}
        lower=ranks['actual_full_source']-25;upper=125-ranks['full350_D']
        result={'source':label,'source_report':path,'source_report_sha256':hashlib.sha256((ROOT/path).read_bytes()).hexdigest(),
            'original_source_substitution':{str(k):str(v) for k,v in sub.items()},
            'entire_source_shear_certificate':shear_certificate,
            'rational_coordinate_reparameterization':'x3_new=I*x3_old, n1=I*v' if which else None,
            'integer_scale':str(scale),'parameters':list(map(str,parameters)),
            'universal_rational_family':encode(P),'all_parameter_strict_Jacobi_zero':True,
            'witness_substitution':{str(k):str(v) for k,v in witness.items()},
            'point_brackets':encode(point),'all_actual_parameter_directions':[encode(R) for R in directions],
            'ranks_QQ':ranks,'certificates':{key:minor_certificate(A,ranks[key]) for key,A in matrices.items()},
            'exact_full_D_times_all_actual_directions_zero':True,
            'projective_family_dimension_lower':lower,'projective_tangent_dimension_upper':upper,
            'true_normal_dimension':upper-lower,'smooth_component_if_equal':lower==upper,
            'matrix_sha256':{key:hashlib.sha256(json.dumps(A.tolist(),default=str).encode()).hexdigest() for key,A in matrices.items()},
            'seconds':time.monotonic()-start}
        records.append(result);print(json.dumps({key:result[key] for key in ['source','ranks_QQ','projective_family_dimension_lower','projective_tangent_dimension_upper','seconds']}),flush=True)
    inputs=[support_path,'code/classification_three_hyperplanes_quadric_independent.py','code/classification_rank3_mixed_component_review.py']
    result={'task_id':'E_reduced_nonnormal','all_passed':True,'cases':records,
        'input_sha256':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in inputs},
        'producer_source':Path(__file__).read_text(),'seconds':time.monotonic()-begun,
        'scope':'One finite full-ambient witness for each proved surviving complete source. Equal tangent and family dimensions imply componenthood; no known-component membership or global exhaustion is inferred.',
        'global_exhaustion_proved':False}
    with output.open('x') as f:json.dump(result,f,indent=2);f.write('\n')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True)
    main(ap.parse_args().output)
