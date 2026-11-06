#!/usr/bin/env python3
"""Independent Fraction reconstruction of three fixed-slice cubic support proofs.
Only the previously independent sparse arithmetic is reused. No producer,
SymPy, NumPy, claimed derivative matrix, or completed-ring assertion is imported.
"""
import hashlib,json,resource,time
from pathlib import Path
from itertools import product
import classification_two_lines_plane_cubic_component_review as q
ROOT=Path(__file__).resolve().parents[1]
NAMES=['node_two_flexes','node_flex_secant','cusp_tangent_cusp_third_secant']

def plus(A,B):return {ij:q.add(A.get(ij,{}),B.get(ij,{})) for ij in set(A)|set(B)}
def jvec(P):return q.jvec(q.jacobi(P))
def mixed(A,B,JA=None,JB=None):
    return [c-a-b for c,a,b in zip(jvec(plus(A,B)),jvec(A) if JA is None else JA,jvec(B) if JB is None else JB)]
def sparse_vector(v):return {str(i):str(c) for i,c in enumerate(v) if c}

def main():
    resource.setrlimit(resource.RLIMIT_CPU,(55,56));resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,)*2)
    started=time.monotonic()
    paths=[ROOT/'reports'/name for name in [
      'CLASSIFICATION_TWO_LINES_PLANE_CUBIC_20261107.json',
      'CLASSIFICATION_TWO_LINES_PLANE_CUBIC_THIRD_SECANT_20261107.json',
      'CLASSIFICATION_TWO_LINES_PLANE_CUBIC_NORMALS_20261107.json',
      'CLASSIFICATION_TWO_LINES_PLANE_CUBIC_TWO_NORMALS_20261107.json',
      'CLASSIFICATION_TWO_LINES_PLANE_CUBIC_CUBIC_SUPPORT_20261107.json']]
    sources={d['case']:d for path in paths[:2] for d in json.loads(path.read_text())['cases']}
    normals={d['case']:d for path in paths[2:4] for d in json.loads(path.read_text())['cases']}
    saved={d['case']:d for d in json.loads(paths[4].read_text())['cases']}
    records=[]
    for name in NAMES:
        begin=time.monotonic();src=sources[name];packet=saved[name]
        P,directions,source_check=q.full_source_check(src)
        R=[q.decode(v) for v in normals[name]['normals']];n=len(R);r=102-n
        assert n==packet['normal_count'] and r==packet['full_D_rank']
        assert all(not v for v in q.jacobi(P).values())
        Dcols=[]
        for ij,e in q.COLS:
            value={key:dict(poly) for key,poly in P.items()}
            value[ij]=q.add(value.get(ij,{}),{e+(0,):q.F(1)})
            value[ij[::-1]]=q.scale(value[ij],-1)
            Dcols.append(jvec(value))
        D=q.transpose(Dcols);VV=[];OO=[]
        for a,b in product(range(5),repeat=2):
            V={};O={}
            for i,j in q.PAIRS:
                V[i,j]=q.add(q.atom(1,i,b) if j==a else {},q.atom(-1,j,b) if i==a else {})
                O[i,j]=q.add(P.get((b,j),{}) if i==a else {},P.get((i,b),{}) if j==a else {},q.scale(q.mul(q.atom(1,b),q.der(P.get((i,j),{}),a)),-1))
            VV.append(q.vec(V));OO.append(q.vec(O))
        Gcols=VV+OO+[q.vec(P)]
        Fcols=Gcols+[q.vec(A) for A in directions]
        Kcols=Fcols+[q.vec(A) for A in R]
        assert q.rank(q.transpose(VV))[0]==24
        assert q.rank(q.transpose(Gcols))[0]==47
        frank,fcols=q.rank(q.transpose(Fcols));assert frank==48
        assert q.rank(q.transpose(Kcols))[0]==48+n
        assert all(sum(a*b for a,b in zip(row,col))==0 for row in D for col in Kcols)
        rows=packet['implicit_D_rows'];cols=packet['implicit_D_columns']
        assert len(rows)==len(cols)==r
        det=q.determinant([[D[i][j] for j in cols] for i in rows])
        assert det and det==q.F(packet['implicit_minor_determinant'])
        # These checks prove a literal direct-sum decomposition C150=F48+S102.
        # The coordinate pivot directions have injective D image, while F+N
        # is the full kernel of dimension48+n. No group slice theorem is used.
        pairs=[tuple(pair) for pair in packet['normal_pairs']]
        assert pairs==([(0,0)] if n==1 else [(0,0),(0,1),(1,1)])
        JR=[jvec(A) for A in R]
        Q=[JR[i] if i==j else mixed(R[i],R[j],JR[i],JR[j]) for i,j in pairs]
        T=[q.decode(data) for data in packet['quadratic_corrections']]
        assert len(T)==len(Q)
        tv=[q.vec(A) for A in T]
        assert all(not value for vector in tv for j,value in enumerate(vector) if j not in cols)
        residual=[]
        for vector,rawQ in zip(tv,Q):
            value=[sum(a*b for a,b in zip(row,vector))+c for row,c in zip(D,rawQ)]
            assert all(value[i]==0 for i in rows)
            residual.append(value)
        quadratic_rows=q.transpose(residual)
        expected_quadratic=([[q.F(1),q.F(0),q.F(-16)],[q.F(0),q.F(1),q.F(-4)]] if name=='node_two_flexes'
                            else [[q.F(1),q.F(0),q.F(9,4400)]] if name=='node_flex_secant' else [])
        qrank=q.rank(quadratic_rows)[0]
        assert qrank==len(expected_quadratic)
        if expected_quadratic:assert q.rank(quadratic_rows+expected_quadratic)[0]==qrank
        JT=[jvec(A) for A in T]
        B={(i,j):mixed(R[i],T[j],JR[i],JT[j]) for i in range(n) for j in range(len(T))}
        if n==1:C=[B[0,0]]
        else:C=[B[0,0],[a+b for a,b in zip(B[1,0],B[0,1])],
                [a+b for a,b in zip(B[1,1],B[0,2])],B[1,2]]
        assert [sparse_vector(v) for v in C]==packet['all_cubic_coefficient_vectors']
        pairings=[]
        for fun in packet['cubic_functionals']:
            left={int(i):q.F(v) for i,v in fun['nonzero_coefficients'].items()}
            assert all(sum(v*D[i][j] for i,v in left.items())==0 for j in range(150))
            assert all(sum(v*col[i] for i,v in left.items())==0 for col in Q)
            values=[sum(v*col[i] for i,v in left.items()) for col in C]
            assert values==[q.F(v) for v in fun['cubic_pairing']]
            pairings.append(values)
        if n==1:
            assert len(pairings)==1 and pairings[0][0]
            radical={'certificate':'nonzero literal cubic c*t0^3','coefficient':str(pairings[0][0])}
            expected_cubic=[[q.F(1)]]
        elif name=='node_two_flexes':
            line=[q.F(64),q.F(16),q.F(4),q.F(1)]
            value=sum(a*b for a,b in zip(pairings[3],line));assert value==230605056000000
            radical={'quadratic_support_line':'t0=4*t1','chosen_cubic_functional':3,
                     'cubic_on_line_coefficient':str(value),'conclusion':'t1^3=0 on the sole quadratic support line, hence reduced support is the origin'}
            expected_cubic=[[q.F(int(i==j)) for j in range(4)] for i in range(4)]
        else:
            c=q.F(9,4400);p=pairings[1]
            assert p[1]==p[3]==0
            rem=p[2]-c*p[0];assert rem==q.F(-174254976000,11)
            radical={'quadratic':'t0^2+(9/4400)*t1^2','chosen_cubic_functional':1,
                     'cubic_mod_quadratic':str(rem)+'*t0*t1^2',
                     'power_membership':'t1^2*(t0^2+c*t1^2)-t0*(t0*t1^2)=c*t1^4, c=9/4400!=0',
                     'conclusion':'t1^4 and t0^2 lie in the radical certificate; reduced support is the origin'}
            expected_cubic=[[q.F(1),q.F(0),q.F(0),q.F(0)],
                            [q.F(0),q.F(1),q.F(0),c],[q.F(0),q.F(0),q.F(1),q.F(0)]]
        cr=q.rank(pairings)[0]
        assert cr==len(expected_cubic) and q.rank(pairings+expected_cubic)[0]==cr
        record={'case':name,'all_passed':True,**source_check,
          'full_D_shape':[350,150],'full_D_exact_rank':r,'vertical_rank':24,'actual_group_scale_vertical_rank':47,
          'actual_full_source_rank':48,'actual_source_basis_columns':fcols,'normal_count':n,
          'full_kernel_rank':48+n,'full_D_times_all_actual_source_and_normal_directions_zero':True,
          'true_literal_slice_dimension':102,
          'true_literal_slice_basis':'The published normals and standard coordinate vectors indexed by implicit_D_columns; the verified kernel decomposition and nonzero D minor prove this complements F48.',
          'implicit_D_rows':rows,'implicit_D_columns':cols,'implicit_minor_determinant':str(det),
          'quadratic_corrections_in_literal_pivot_coordinate_span':True,
          'all_selected_IFT_rows_vanish_to_order_two':True,
          'complete_quadratic_initial_space_rank':qrank,'quadratic_initial_generators':packet['quadratic_initial_generators'],
          'all_published_cubic_coefficient_vectors_independently_reconstructed':True,
          'all_cubic_functionals_annihilate_all_D_columns_and_all_raw_normal_quadrics':True,
          'cubic_pairing_matrix':[[str(v) for v in row] for row in pairings],
          'verified_cubic_initial_subspace_rank':cr,'cubic_initial_generators':packet['cubic_initial_generators'],
          'explicit_reduced_origin_certificate':radical,
          'zero_dimensional_true_slice_reduced_support':True,'projective_source_closure_component_dimension':23,
          'component_logic':'The genuine102-dimensional affine linear slice is cut by48 linear equations. The selected r equations eliminate all pivot variables formally. The verified quadratic and cubic initial forms have maximal radical, so the cut has dimension0. Krull bounds ambient local dimension by48; the irreducible actual source has dimension at least48. Its closure is a component. No local uniqueness or dense orbit is inferred.',
          'seconds':time.monotonic()-begin}
        records.append(record);print(json.dumps({'case':name,'all_passed':True,'seconds':record['seconds'],'D':r,'source':48,'normal':n}),flush=True)
    result={'task_id':'E_reduced_nonnormal','all_passed':True,
      'method':'Independent standard-library Fraction sparse coefficient reconstruction using only the earlier independent helper. All published cubic pairings are recomputed. No producer arithmetic or symbolic algebra engine imported.',
      'input_sha256':{str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest() for path in paths+[Path(q.__file__)]},
      'reviewer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'reviewer_source':Path(__file__).read_text(),
      'cases':records,'independent_formal_review':'generic_structure independently accepted the fixed-slice IFT and Krull implication; full-kernel ambiguity invariance is unnecessary.',
      'scope':'Three whole six-magnetic-parameter source closures are23-dimensional components. Completeness of the geometric source list and distinctness from all previous components are separate dependencies. No exact nonreduced ring, local uniqueness, dense orbit, or exhaustive cubic ideal is claimed.',
      'seconds':time.monotonic()-started}
    out=ROOT/'reports/CLASSIFICATION_TWO_LINES_PLANE_CUBIC_CUBIC_REVIEW_20261107.json'
    with out.open('x') as f:json.dump(result,f,indent=2);f.write('\n')

if __name__=='__main__':main()
