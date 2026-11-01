#!/usr/bin/env python3
"""A true-core flag sample: finite-base rank5, with separate mod-p six-initial data."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import tempfile
import time
from itertools import combinations_with_replacement as cwr
import sympy as s
from s61_model import S61
from s61_minimal_infinity_local_data_20261101 import calculate
from s61_flag_transverse_obstruction_20261101 import specialize
from truncated_local_ideal_io_20261101 import write_input, m2_initial, monomials
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from weight1_structure import COORDS, sp_basis, infinitesimal_action

ROOT = Path(__file__).resolve().parents[1]
PRIME = 32003


def mod(value):
    value = s.Rational(value)
    assert int(value.q) % PRIME
    return int(value.p)*pow(int(value.q), -1, PRIME) % PRIME


def run(command, **kwargs):
    result = subprocess.run(list(map(str, command)), cwd=ROOT, text=True,
                            capture_output=True, **kwargs)
    if result.returncode:
        raise RuntimeError(result.stdout+'\n'+result.stderr)
    return result


def m2(program, directory):
    path = directory/'geometry.m2'
    path.write_text(program)
    env = os.environ.copy()
    env['DYLD_LIBRARY_PATH'] = str(ROOT/'work/m2-runtime/lib')+':/opt/homebrew/lib:/opt/homebrew/opt/givaro/lib'
    return run(['/opt/homebrew/Cellar/macaulay2/1.26.05/bin/M2-binary', '--script', path],
               env=env, timeout=90)


def actual_support(model, old, point, J, T):
    z, u, B = primitive_beta(model)
    vec = s.Matrix([B[i,j].diff(z[k]) for (i,j),k in COORDS])
    basis = vec.jacobian(u)
    rows = list(basis.T.rref()[1])
    inverse = basis[rows,:].inv()
    action = s.Matrix.hstack(*[inverse*infinitesimal_action(vec,A)[rows,:] for A in sp_basis()])
    one = [i for i,w in enumerate(model.reduced['parameter_weights']) if w == 1]
    action = action.subs(dict(zip(u, [point[i] for i in one])))
    H = s.Matrix(old['beta_kernel']).applyfunc(s.Rational)
    normal = s.Matrix.hstack(*H.T.nullspace()).T*action
    assert normal.rank() == 1
    allowed = s.Matrix.hstack(*normal.nullspace())
    O = s.Matrix(old['q_orbit_differential']).applyfunc(s.Rational).col_join(s.zeros(1,10))
    orbit = O*allowed
    selected = list(orbit.rref()[1])
    orbit = orbit[:,selected]
    assert orbit.rank() == 3
    free = [i for i in range(45) if i not in J.to_DM().rref()[1]]
    CT = orbit[free,:]
    assert T*CT == orbit
    parameters = [0,1,6]
    determinant = CT[parameters,:].det()
    assert determinant != 0
    hr = list(H.T.rref()[1])
    stabilizer = s.Matrix.hstack(*(O*allowed).nullspace())
    assert not any((H[hr,:].inv()*action[hr,:]*allowed*stabilizer)[1,:])
    b_action = H[hr,:].inv()*action[hr,:]*allowed[:,selected]
    differential = -b_action[1,:]
    coordinate_differential = differential*CT[parameters,:].inv()
    assert any(coordinate_differential)
    ann = s.zeros(7,10)
    for i,j in enumerate([9,8,7,5,4]):
        ann[i,j] = 1
    ann[5,3],ann[5,6],ann[6,2],ann[6,6] = 1,2,1,4
    assert ann*CT == s.zeros(7,3) and ann.rank() == 7
    return dict(tangent=[[str(c) for c in row] for row in CT.tolist()],
                parameter_coordinates=parameters, parameter_minor=str(determinant),
                beta_normal=list(map(str,normal)),
                compensated_b1_derivative_in_parameters=list(map(str,coordinate_differential)),
                tangent_annihilator=[[str(c) for c in row] for row in ann.tolist()],
                same_smooth_support_generic_length_lower_bound=5)


def verify():
    start = time.time()
    model = S61()
    source = 'reports/S61_FLAG_UNIFORM_OPEN_20261101.json'
    old = json.loads((ROOT/source).read_text())['flag']
    H = s.Matrix(old['beta_kernel']).applyfunc(s.Rational)
    b = s.zeros(15,1)
    b[0] = b[12] = 1
    b[2] = -s.Rational(4,25)
    for i,value in zip([4,7,8,9,10,11,13],[1,2,3,5,7,11,13]):
        b[i] = value
    point = s.Matrix(old['point61']).applyfunc(s.Rational)
    one = [i for i,w in enumerate(model.reduced['parameter_weights']) if w == 1]
    for i,value in zip(one,H*b):
        point[i] = value
    local = calculate(model,True,point)
    J = s.Matrix(local['jacobian']).applyfunc(s.Rational)
    T = s.Matrix(local['tangent_basis']).applyfunc(s.Rational)
    n, rank = T.cols, J.rank()
    assert (rank,n) == (35,10) and not any(T[44,:])
    support = actual_support(model,old,point,J,T)
    columns = local['chart_columns']
    fixed = next(i for i in range(45) if i not in columns)
    _, equations = specialize(model,{i:point[i] for i in one},True)
    terms = []
    for equation in equations:
        row = []
        for monomial,c in equation.items():
            if fixed not in monomial:
                a,b0 = [columns.index(k) for k in monomial]+[44]*(2-len(monomial))
                mod(c)
                row.append((a,b0,c))
        terms.append(row)
    pivots = list(J.to_DM().rref()[1])
    selected = list(J[:,pivots].T.to_DM().rref()[1])
    U = s.zeros(45,rank)
    inverse = J[selected,pivots].inv(method='DM')
    for i,col in enumerate(pivots):
        U[col,:] = inverse[i,:]
    assert J*T == s.zeros(135,n) and J[selected,:]*U == s.eye(rank)
    factors = [s.ilcm(*[s.denom(c) for c in list(J[i,:])+[c for _,_,c in terms[i]]]) for i in range(135)]
    assert all(int(c) % PRIME for c in factors)
    Jn = s.diag(*factors)*J
    Un = U*s.diag(*[1/s.Integer(factors[i]) for i in selected])
    denominator = s.ilcm(*[s.denom(c) for c in Un])
    td = s.ilcm(*[s.denom(c) for c in T])
    assert int(denominator) % PRIME and int(td) % PRIME
    Un = (denominator*Un).applyfunc(s.cancel)
    Tn = td*T
    depth = 6
    mons = {d:list(cwr(range(n),d)) for d in range(depth+1)}
    indices = {d:{mon:j for j,mon in enumerate(mons[d])} for d in mons}
    integer_native = 'code/integer_implicit_jets_20261101.cpp'
    closure_native = 'code/truncated_local_ideal_filtered_trace_20261101.cpp'
    code = (ROOT/integer_native).read_text().split(' ofstream out(argv[2]);')[0]+r'''
 const long prime=32003;
 auto power=[&](long a,long b){long z=1;while(b){if(b&1)z=z*a%prime;a=a*a%prime;b>>=1;}return z;};
 vector<long>iv(D+1);for(int d=2;d<=D;d++){if(mpz_divisible_ui_p(scale[d].get_mpz_t(),prime))return 8;iv[d]=power(mpz_fdiv_ui(scale[d].get_mpz_t(),prime),prime-2);}
 ofstream out(argv[2]);if(!out)return 9;
 for(int i=0;i<F;i++){for(int d=2;d<=D;d++)for(auto&a:A[d][i])out<<mpz_fdiv_ui(a.get_mpz_t(),prime)*iv[d]%prime<<" ";out<<"\n";}
 for(int d=2;d<=D;d++){vector<long>r(sizes[d]);for(auto&[k,c]:W[d][44])r[k]=mpz_fdiv_ui(c.get_mpz_t(),prime)*iv[d]%prime;for(auto c:r)out<<c<<" ";}out<<"\n";
 }
'''
    tail = r'''vector<int>ob(N);for(int&i:ob)in>>i;for(int i=0;i<N;i++)if(ob[i]&&!base[i].empty()){int a=ob[i];for(auto[k,c]:base[i]){ob[k]=(ob[k]-(long long)a*c)%p;if(ob[k]<0)ob[k]+=p;}}cout<<"rank "<<rank<<" length "<<N-rank<<"\n";for(int d=1;d<=cutoff+1;d++){int ct=0;for(int i=0;i<N;i++)if(deg[i]<d&&bid[i]<0)ct++;cout<<"modm "<<d<<" length "<<ct<<"\n";}cout<<"t";for(int i=0;i<N;i++)if(ob[i])cout<<" "<<i<<":"<<ob[i];cout<<"\n";'''
    closure_code = (ROOT/closure_native).read_text().replace('cout<<rank<<"\\n";',tail)
    assert 'vector<int>ob(N)' in closure_code
    with tempfile.TemporaryDirectory(prefix='true-core-rank5-') as temp:
        directory = Path(temp)
        cpp = directory/'jets.cpp'; cpp.write_text(code)
        native = directory/'jets'
        run(['c++','-O3','-std=c++17','-I/opt/homebrew/include',cpp,'-L/opt/homebrew/lib','-lgmpxx','-lgmp','-o',native])
        inp = directory/'jets-input'
        with inp.open('w') as output:
            output.write(f'135 {n} {depth} {rank}\n{denominator} {td}\n')
            output.write(' '.join(str(len(mons[d])) for d in mons)+'\n')
            output.write(' '.join(map(str,selected))+'\n')
            for matrix in [Jn,Tn,Un]:
                for row in matrix.tolist():
                    output.write(' '.join(map(str,row))+'\n')
            for i,row in enumerate(terms):
                output.write(str(len(row))+' '+' '.join(f'{a} {b0} {c*factors[i]}' for a,b0,c in row)+'\n')
            for a in range(1,depth):
                for b0 in range(1,depth+1-a):
                    output.write(' '.join(str(indices[a+b0][tuple(sorted(x+y))]) for x in mons[a] for y in mons[b0])+'\n')
        coefficients = directory/'coefficients'
        native_result = run([native,inp,coefficients],timeout=120)
        raw = [list(map(int,row.split())) for row in coefficients.read_text().splitlines()]
        assert len(raw) == 136
        A, observable, offset = {}, {}, 0
        for d in range(2,depth+1):
            A[d] = [row[offset:offset+len(mons[d])] for row in raw[:135]]
            observable[d] = raw[135][offset:offset+len(mons[d])]
            offset += len(mons[d])
        assert all(len(row)==offset for row in raw)
        cp = directory/'closure.cpp'; cp.write_text(closure_code)
        closure = directory/'closure'
        run(['c++','-O3','-std=c++17',cp,'-o',closure])
        def close(jets,obs,variables,label):
            input_path = directory/(label+'-input')
            write_input(input_path,PRIME,variables,jets,obs,depth)
            trace,order,initial = [directory/(label+x) for x in ['-trace','-order','-initial']]
            result = run([closure,input_path,trace,order,initial],timeout=180)
            lengths = [int(line.split()[-1]) for line in result.stdout.splitlines() if line.startswith('modm ')]
            tline = next(line for line in result.stdout.splitlines() if line.startswith('t'))
            terms0 = [int(x.split(':')[0]) for x in tline.split()[1:]]
            mm = monomials(variables,depth)
            order_t = min((len(mm[j]) for j in terms0),default=None)
            return initial, lengths, order_t, hashlib.sha256(initial.read_bytes()).hexdigest()
        initial, lengths, order_t, initial_hash = close(A,observable,n,'full')
        assert lengths == [1,11,51,152,357,689,1163] and order_t == 4
        normal = [2,3,4,5,7,8,9]
        cut,cut_t = {},{}
        for d in range(2,depth+1):
            cc = [indices[d][tuple(normal[i] for i in mon)] for mon in cwr(range(7),d)]
            cut[d] = [[row[j] for j in cc] for row in A[d]]
            cut_t[d] = [observable[d][j] for j in cc]
        _, transverse_lengths, transverse_order_t, _ = close(cut,cut_t,7,'transverse')
        assert transverse_lengths == [1,8,21,35,43,49,49] and transverse_order_t == 4
        minimal_path = directory/'minimal.txt'
        program = m2_initial(initial,n,depth,PRIME)+r'''
C=ideal(v_9,v_8,v_7,v_5,v_4,v_3+2*v_6,v_2+4*v_6);
assert(radical J==C);assert(dim(R/J)==3 and degree(R/J)==5);
h=v_1-18*v_6/25;B=saturate(J,ideal h);E=(module B)/(module J);
assert(dim(R/B)==3 and degree(R/B)==5);
assert(dim(R/(B+ideal(v_0,v_1,v_6)))==0 and degree(R/(B+ideal(v_0,v_1,v_6)))==5);
pn=poincare(R/B);z=(ring pn)_0;
assert(pn==(1+4*z)*(1-z)^7);
assert(dim E==2 and degree E==111);
pn=poincare E;z=(ring pn)_0;num=3*z+16*z^2+26*z^3+38*z^4+18*z^5+10*z^6;
assert(pn==num*(1-z)^8);
Q=E/(ideal(v_0,v_6)*E);assert(dim Q==0 and degree Q==111);
assert(apply(toList(0..7),i->hilbertFunction(i,Q))=={0,3,16,26,38,18,10,0});
print "MODULAR INITIAL AND EMBEDDED MODULE PASS";
'''
        program += 'G=mingens J;out=openOut "'+str(minimal_path)+'";scan(numColumns G,i->out<<toString(G_(0,i))<<endl);close out;\n'
        geometry = m2(program,directory)
        assert 'MODULAR INITIAL AND EMBEDDED MODULE PASS' in geometry.stdout
        minimal = minimal_path.read_text().splitlines()
        assert len(minimal)==71
    paths = [str(Path(__file__).relative_to(ROOT)),integer_native,closure_native,
             'code/truncated_local_ideal_io_20261101.py','code/s61_model.py',
             'code/s61_minimal_infinity_local_data_20261101.py','code/s61_flag_transverse_obstruction_20261101.py',
             'code/primitive_jet_elliptic_curve_20261101.py','code/weight1_structure.py',source,
             'research/notes/S61_FLAG_ZERO12_RANK_FIVE_20261101.md',
             'research/notes/S61_FLAG_TRUE_CORE_RANK_FIVE_20261101.md']
    return dict(id='S61_FLAG_TRUE_CORE_RANK_FIVE_20261101',all_passed=True,point61=list(map(str,point)),
                core='-4/25',nuisance=[1,2,3,5,7,11,13],jacobian_rank=35,tangent_dimension=10,
                integral_model_prime=PRIME,all_denominators_and_selected_inverse_prime_units=True,
                actual_support=support,
                characteristic_zero_conclusions=dict(finite_over_three_parameter_formal_base=True,
                    module_generators_at_most=49,generic_module_rank=5,homogenizer_is_base_torsion=True),
                modular_actual_sixth_hilbert=[1,10,40,101,205,332,474],modular_quotient_by_m7=1163,
                modular_t_order_through_six=order_t,modular_transverse_hilbert=[1,7,13,14,8,6],
                modular_transverse_length=49,modular_transverse_m6_zero=True,
                modular_six_initial=dict(dimension=3,degree=5,primary_part_free_rank=5,
                    embedded_difference_free_rank=111,embedded_generator_degrees=[3,16,26,38,18,10],
                    initial_sha256=initial_hash,minimal_generators=minimal),
                native_log=native_result.stderr,m2_log=geometry.stdout+geometry.stderr,
                limitations=['The six-initial ideal and embedded-module computations are over F32003 only.',
                    'No characteristic-zero transverse-length equality or t-order equality is inferred.',
                    'No full initial ideal, ambient t-nilpotence or whole-core saturation removal.'],
                proof_note='research/notes/S61_FLAG_TRUE_CORE_RANK_FIVE_20261101.md',
                input_sha256={path:hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in paths},
                seconds=time.time()-start)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();result=verify()
    with args.output.open('x') as stream:
        json.dump(result,stream,indent=2);stream.write('\n')
    print('TRUE FLAG CORE FINITE RANK5 PASS',result['seconds'])
