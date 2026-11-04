#!/usr/bin/env python3
"""ONE row-leading test, then bounded primitive-row rank3, never big3adj."""
import argparse,gzip,itertools,json,resource,shlex,subprocess,time
from fractions import Fraction as Q
from pathlib import Path
import sympy as s
from q_new_family_v_zero_degree26_support_preflight_20261104_v63 import read,sha
from s61_sparse_generic5_input_preflight_20261104_v63 import mul
ROOT=Path(__file__).resolve().parents[1]
ID='S61_SPARSE_GENERIC5_PLANE_RANK3_CIRCUIT_20261104_v63'
def dec(p):return {tuple(t[1:]):Q(t[0])for t in p}
def enc(p):return[[str(c),*e]for e,c in sorted(p.items())]
def circuit(rows,rankproof=None):
    nodes=[]
    def new(op,*a):nodes.append([op,*a]);return len(nodes)-1
    zero=new('qq','0');one=new('qq','1');minus=new('qq','-1');W=[[new('original_W_poly',i,j)for j in range(3)]for i in range(45)];B=[W[i]for i in rows]
    def prod(a):
        q=one
        for i in a:q=new('mul',q,i)
        return q
    def total(a):
        q=zero
        for i in a:q=new('add',q,i)
        return q
    terms=[]
    for p in itertools.permutations(range(3)):
        v=prod([B[i][p[i]]for i in range(3)]);terms.append(v if sum(p[i]>p[j]for i in range(3)for j in range(i+1,3))%2==0 else new('mul',minus,v))
    gamma=total(terms);factor_identity=None
    if rankproof:
        # Keep inverse determinant factored; do not multiply its large polynomials.
        ff=[new('selected_row_factor_poly',i)for i in range(3)]+[new('rank_pivot_rf',i)for i in range(3)]
        sign=new('qq',str(rankproof['column_permutation_sign']));inv=prod([sign,*[new('inv',f)for f in ff]])
        factor_identity='Gamma(original W block)=column_permutation_sign*product(selected row gcd factors)*product(nonzero field pivots); original row order is selected order.'
    else:inv=new('inv',gamma)
    adj=[[None]*3 for _ in range(3)]
    for i in range(3):
        for j in range(3):
            rr=[k for k in range(3)if k!=j];cc=[k for k in range(3)if k!=i];z=new('sub',new('mul',B[rr[0]][cc[0]],B[rr[1]][cc[1]]),new('mul',B[rr[0]][cc[1]],B[rr[1]][cc[0]]));adj[i][j]=z if (i+j)%2==0 else new('mul',minus,z)
    E=[[new('mul',inv,total([new('mul',W[i][k],adj[k][j])for k in range(3)]))for j in range(3)]for i in range(45)];normal=[i for i in range(45)if i not in rows];zz=[new('z',j)for j in range(3)];nn=[new('n',j)for j in range(42)];q=[]
    for i in range(45):
        v=total([new('center',i),*[new('mul',E[i][j],zz[j])for j in range(3)]])
        q.append(new('add',v,nn[normal.index(i)])if i in normal else v)
    return dict(nodes=nodes,source_parameters=['H','L','s0','q0','t0'],leaf_source='plane-W-original.json.gz',factor_leaf_source='rank-certificate.json.gz'if rankproof else None,
        leaf_bindings={'original_W_poly':'plane-W-original.json.gz/W45x3[row][column]','center':'plane-W-original.json.gz/center45[row]','rank_pivot_rf':'rank-certificate.json.gz/nonzero_pivot_factors[index]','selected_row_factor_poly':'rank-certificate.json.gz/selected_row_factors[index]'},
        Gamma3_node=gamma,inverse_Gamma3_node=inv,factored_Gamma_identity=factor_identity,adjugate3x3_nodes=adj,E45x3_nodes=E,q45_nodes=q,pivot_q_rows=rows,normal_q_rows42=normal,
        inverse='z=q_pivot-center_pivot; n_i=q_nonpivot_i-center_i-sum_j E_ij*z_j',
        proof='Use ORIGINAL W, not row-primitive W. B*adj(B)=adj(B)*B=Gamma*I, Gamma!=0; henceE_pivotI3. W=P*V globally, so E equals original plane embedding.',
        all45_unknowns_retained=True,no_unknown_q_inverted=True,no_expanded3det_or3adjugate=True,no134_orJ0_substitution=True)
def main(out,work):
    start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False);work.mkdir(parents=True,exist_ok=False)
    source=ROOT/'reports/S61_SPARSE_GENERIC5_SMALL_TAIL_FIELD_20261104_v63';gate=read(source/'report.json');assert gate['all_passed']and gate['kernel_complete']and gate['old_span_verified']
    for n,h in gate['input_sha256'].items():assert sha(ROOT/n)==h,n
    sh=read(source/'shared-denominator-kernel-plane.json.gz');P=dec(sh['common_denominator_P']);assert P
    mixed=[[dec(p)for p in row]for row in sh['plane45x3_mixed_numerators']];W=[[row[0],row[1],mul(P,row[2])]for row in mixed];center=read(source/'plane-rational45x3.json')['center45']
    original=dict(variables=['H','L','s0','q0','t0'],W45x3=[[enc(p)for p in row]for row in W],center45=center,common_P=sh['common_denominator_P'],relation='W=P*[N*T0,N*T1,Q52-center]',row_primitive_never_used_as_actual_plane=True)
    with gzip.open(out/'plane-W-original.json.gz','wt')as f:json.dump(original,f,separators=(',',':'))
    lead=[max((m for p in row for m in p),default=None)for row in W];M=s.Matrix([[s.Rational(p.get(m,0))if m is not None else 0 for p in row]for row,m in zip(W,lead)]);rr=list(M.T.rref()[1]);cc=list(M[rr,:].rref()[1])if rr else[];det=M.extract(rr,cc).det()if rr else s.Integer(1);assert det
    leading=dict(rank=len(rr),row_wise_common_leading_monomials=[list(m)if m is not None else None for m in lead],matrix45x3=[[str(c)for c in row]for row in M.tolist()],rows=rr,columns=cc,nonzero_leading_coefficient_minor=str(det),method='lex H,L,s0,q0,t0; rowwise maximum common toall3cols, NOT flattened coefficients')
    (out/'row-leading-certificate.json').write_text(json.dumps(leading,indent=2)+'\n');rank=len(rr);selected=rr if rank==3 else None;result=None;error=None
    native=ROOT/'code/s61_sparse_generic5_plane_rank3_circuit_20261104_v63.cpp';helper=ROOT/'code/q_new_family_v_zero_degree26_restriction_20261104_v63.cpp'
    try:
        if rank<3:
            inp=work/'input.txt'
            with inp.open('w')as f:
                f.write('45 3\n')
                for row in W:
                    for p in row:
                        f.write(str(len(p))+'\n')
                        for e,c in p.items():f.write(str(c)+' '+' '.join(map(str,e))+'\n')
            with gzip.open(out/'native-input.txt.gz','wb')as f:f.write(inp.read_bytes())
            flags=shlex.split(subprocess.check_output(['pkg-config','--cflags','--libs','flint'],text=True));binary=work/'plane-rank3';subprocess.run(['c++','-O2','-std=c++17',str(native),*flags,'-lgmp','-o',str(binary)],check=True,timeout=max(1,60-(time.monotonic()-start)))
            with(out/'native.log').open('w')as f:subprocess.run([str(binary),str(inp),str(out)],stdout=f,stderr=subprocess.STDOUT,check=True,timeout=max(1,60-(time.monotonic()-start)))
            result=read(out/'rank-certificate.json');assert result['all_passed'];rank=result['rank'];selected=result['selected_original_q_rows']if rank==3 else None
        if rank==3:
            a=circuit(selected,result)
            qids=[i for i,w in enumerate(read(ROOT/'legacy/v02/data/symplectic_slice_reduced_equations.json')['parameter_weights'])if w==2]
            a.update(q_to_u_indices_zero_based=qids,pivot_u_numbers_one_based=[qids[i]+1 for i in selected],leaf_format='Wterms:[QQcoefficient,H,L,s0,q0,t0 exponents]; center/rowfactor:QQcontent times integer terms; rankpivot:integer numerator/denominator polynomial pair')
            (out/'normal42-adapter-SLP.json').write_text(json.dumps(a,indent=2)+'\n')
        assert sum(p.stat().st_size for p in out.iterdir()if p.is_file())<32*1024*1024
        # Native already enforced the raw32MiB cap; archive large literal certificates.
        for name in['row-primitive45.json','field-op-0.json','field-op-1.json','field-op-2.json','rank-certificate.json']:
            p=out/name
            if p.exists():
                (out/(name+'.gz')).write_bytes(gzip.compress(p.read_bytes(),mtime=0));p.unlink()
    except Exception as e:error=repr(e)
    paths=[Path(__file__),native,helper,ROOT/f'research/notes/{ID}.md',source/'report.json',source/'shared-denominator-kernel-plane.json.gz',source/'plane-rational45x3.json',ROOT/'legacy/v02/data/symplectic_slice_reduced_equations.json']
    report=dict(id=ID,all_passed=error is None,error=error,leading_rank=len(rr),rank3_proved=rank==3,selected_q_rows=selected,normal42_adapter_SLP=(out/'normal42-adapter-SLP.json').exists(),
        original_plane_used_for_adapter=True,row_primitive_rank_scaling_only=True,
        generic_field_only=True,no_point_prime_GB_oldCramer_or134_J0_substitution=True,no_big3det_adjugate_expansion=True,
        input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},evidence_sha256={p.name:sha(p)for p in out.iterdir()if p.is_file()},seconds=time.monotonic()-start)
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(dict(all_passed=report['all_passed'],leading_rank=len(rr),rank3=report['rank3_proved'],adapter=report['normal42_adapter_SLP'],rows=selected,seconds=report['seconds'],error=error)))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--work',type=Path,required=True);a=p.parse_args();main(a.output,a.work)
