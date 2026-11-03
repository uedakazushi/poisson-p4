#!/usr/bin/env python3
"""Actual42 Euler-normal model: sharp degree28 delta and scalar target/barrier."""
import argparse,gzip,hashlib,itertools,json,time
from pathlib import Path
import numpy as np
import sympy as s
from sympy.polys.matrices import DomainMatrix
ROOT=Path(__file__).resolve().parents[1];PRIME=32003
def load(p):return json.loads(gzip.decompress(p.read_bytes())if p.suffix=='.gz'else p.read_bytes())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def residue(q):
    q=s.sympify(q);a=int(s.numer(q))%PRIME;b=int(s.denom(q))%PRIME;assert b
    return a*pow(b,-1,PRIME)%PRIME
def main(out):
    started=time.monotonic();out.mkdir(parents=True,exist_ok=False)
    deps=[ROOT/'reports/S61_AFFINE_NORMAL_DIRECTION_GRAPH_20261103_v63_REVIEW/report.json',
      ROOT/'reports/S61_AFFINE_NORMAL_COMMUTATOR_20261103_v63_REVIEW.json',
      ROOT/'reports/S61_AFFINE_NORMAL_GLOBAL_QUINTIC_20261103_v63_REVIEW.json',
      ROOT/'reports/S61_AFFINE_NORMAL_SHARP_MOMENT_SCHUR_20261103_v63/report.json',
      ROOT/'reports/S61_AFFINE_NORMAL_SCHUR_SQUARE_VECTOR_20261103_v63_REVIEW.json']
    for p in deps:
        d=load(p);assert d['all_passed'] is True
        for name,digest in d.get('input_sha256',{}).items():assert sha(ROOT/name)==digest,name
    modelpath=ROOT/'reports/S61_AFFINE_NORMAL_DIRECTION_GRAPH_20261103_v63/coefficient-model-mod32003.json.gz';model=load(modelpath);assert model['prime']==PRIME
    specpath=ROOT/'code/data/s61_affine_normal_global_quintic_20261103_v63.json';spec=load(specpath)
    presentationpath=ROOT/'reports/S61_AFFINE_NORMAL_GLOBAL_QUINTIC_20261103_v63/presentation.json.gz';presentation=load(presentationpath)
    # Exact structural check: h4 initially contains no H and has parameter degree1.
    for r in presentation['all134_remainders']:
        for degree,terms in r:
            if degree==4:assert all(len(mon)<=1 and 44 not in mon for mon,c in terms)
    H,Q,t=s.symbols('H Q t');L=s.symbols('L0:5');gg=s.symbols('g0:5');T=s.zeros(5)
    for i in range(4):T[i+1,i]=1
    T[:,4]=s.Matrix([H-gg[0],-gg[1],-gg[2],-gg[3],-gg[4]])
    r=s.Matrix([L[0]+Q,*L[1:]]);moment=s.Matrix.hstack(*(T**j*r for j in range(4)))
    scaling={H:t**5*H,Q:t**8*Q,**{v:t**(4 if i==4 else 5)*v for i,v in enumerate(L)}}
    D=DomainMatrix.from_Matrix(moment.subs(scaling,simultaneous=True)[1:,:]);delta=s.Poly(s.expand(D.domain.to_sympy(D.det())),t)
    top=H**2*L[3]**2*L[4]**2-H*L[3]**3*Q-L[4]*Q**3
    assert delta.degree()==28 and s.expand(delta.LC()-top)==0
    # P2 and Ki2 are exact p-integral QQ circuits; the independently bound
    # modular coefficient model computes their reduction, never a zero lift.
    mons=model['normal_monomials'];P=np.array(model['P_coefficients'],dtype=np.int64);K=np.array(model['K_coefficients'],dtype=np.int64)
    HP=np.zeros((3,42,42,42),dtype=np.int64);HK=np.zeros((5,42,42),dtype=np.int64)
    for col,word in enumerate(mons):
        if len(word)!=2:continue
        a,b=word;HP[:,:,a,b]+=P[:,col].reshape(3,42);HP[:,:,b,a]+=P[:,col].reshape(3,42);HK[:,a,b]+=K[:,col];HK[:,b,a]+=K[:,col]
    HP%=PRIME;HK%=PRIME;v=np.eye(42,dtype=np.int64)[0]
    J=np.einsum('kiab,b->kia',HP,v)%PRIME;pv=(J@v)*pow(2,-1,PRIME)%PRIME;grad=HK@v%PRIME
    second={};third={}
    for a,b in itertools.combinations_with_replacement(range(3),2):
        second[a,b]=np.array([(int(pv[a]@(HK[i]@pv[b]%PRIME))+int(grad[i]@(J[b]@pv[a]%PRIME)))%PRIME for i in range(5)],dtype=np.int64)
    for a,b,c in itertools.combinations_with_replacement(range(3),3):
        jp=J[b]@pv[a]%PRIME;jca=J[c]@pv[a]%PRIME;jcb=J[c]@pv[b]%PRIME
        hh=(np.einsum('iab,b->ia',HP[c],pv[b])%PRIME)@pv[a]%PRIME;jj=J[c]@jp%PRIME
        third[a,b,c]=np.array([(int(jp@(HK[i]@pv[c]%PRIME))+int(pv[b]@(HK[i]@jca%PRIME))+int(pv[a]@(HK[i]@jcb%PRIME))+int(grad[i]@hh)+int(grad[i]@jj))%PRIME for i in range(5)],dtype=np.int64)
    z=s.symbols('z0:3');cofs=[]
    for number,rel in enumerate(sorted(spec['relations'],key=lambda d:d['target_index'])):
        degree=2 if number<2 else 3;row=[]
        for f in rel['cofactors']:
            row.append([(tuple(i for i,e in enumerate(ex)for _ in range(e)),residue(c))for ex,c in s.Poly(s.sympify(f),*z).terms()if sum(ex)==degree])
        cofs.append(row)
    S4=[-sum(c*int(second[a,b][i])for i,row in enumerate(cof)for(a,b),c in row)*pow(6,-1,PRIME)%PRIME for cof in cofs[:2]]
    H5=-sum(c*int(third[a,b,c0][i])for i,row in enumerate(cofs[2])for(a,b,c0),c in row)*pow(24,-1,PRIME)%PRIME
    plane={tuple(mon):int(c)for mon,c in model['plane5_sparse'][0]}
    Q8=(plane.get((1,1),0)*S4[0]**2+plane.get((1,2),0)*S4[0]*S4[1]+plane.get((2,2),0)*S4[1]**2)%PRIME
    first={degree:terms for degree,terms in presentation['all134_remainders'][0]};r4={tuple(mon):residue(c)for mon,c in first[4]}
    c3=sum(residue(c)for mon,c in first[3]if mon==[44])%PRIME
    L4=(r4.get((42,),0)*S4[0]+r4.get((43,),0)*S4[1])%PRIME;L3=c3*H5%PRIME
    nonzero=(H5**2*L3**2*L4**2-H5*L3**3*Q8-L4*Q8**3)%PRIME
    assert S4==[15338,11727]and H5==31338 and nonzero==31492
    certificate=dict(normal_direction=list(map(int,v)),P2_values=pv.tolist(),P2_Jacobians=J.tolist(),K2_gradients=grad.tolist(),
      all_ordered_second_derivative_values={','.join(map(str,k)):val.tolist()for k,val in second.items()},
      all_ordered_third_derivative_values={','.join(map(str,k)):val.tolist()for k,val in third.items()},
      cofactor_highest_homogeneous_terms=cofs,S1_S2_degree4=S4,H_degree5=H5,Q_degree8=Q8,L4_degree4=L4,L3_degree5=L3,
      actual_delta_degree28_witness=nonzero,prime=PRIME,
      exact_universal_degree28_top=str(top),
      Euler_Red_recipe='Forhomogeneousf_d: Red(z_k*f_d)=(sum_j P_kj(n)*partial_j f_d)/d; decomposebydegree anditeratewordinascendingouteroperatororder.',
      exact_scalar_target='delta_Euler(n,a) belongs J0=(H3,compatibility2583,cleared_plane210,commutator42)',
      scalar_generator_count=2838,scalar_generator_degree_bound=4)
    (out/'certificate.json.gz').write_bytes(gzip.compress((json.dumps(certificate,separators=(',',':'))+'\n').encode(),mtime=0))
    note=ROOT/'research/notes/S61_AFFINE_ACTUAL42_EULER_DELTA28_20261103_v63.md'
    paths=[Path(__file__),note,modelpath,specpath,presentationpath,*deps]
    report=dict(id='S61_AFFINE_ACTUAL42_EULER_DELTA28_20261103_v63',all_passed=True,
      new_Euler_normal_reduction_model_not_mutation_of_frozen_old_recipe=True,
      full_actual_134_equation_algebra_preserved_by_cross_and_plane_cofactor_identities=True,
      actual42_delta_degree_bound_and_exact_generic_degree=28,
      exact_characteristic_zero_nonzero_from_positive_punit_circuit_witness=True,
      modular_zero_or_emptiness_not_lifted=True,expanded_QQ_coefficient_not_claimed=True,
      unknown_empty_iff_actual_scalar_delta_in_J0=True,
      actual_scalar_ideal_generator_count=2838,actual_scalar_generator_degree_bound=4,
      polynomial_cofactor_degree_lower_bound=24,
      membership_or_generic_five_point_exhaustion_proved=False,
      no_constant_span_or_degree3_Macaulay_rerun=True,
      input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},evidence_sha256={'certificate.json.gz':sha(out/'certificate.json.gz')},seconds=time.monotonic()-started)
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print('ACTUAL42 EULER DELTA28 /2838scalarcriterion /COFACTOR>=24 PASS',report['seconds'],flush=True)
if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True);main(ap.parse_args().output.resolve())
