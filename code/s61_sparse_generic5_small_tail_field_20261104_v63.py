#!/usr/bin/env python3
"""ONE row-primitive/canonical generic-field rank7 kernel method."""
import argparse,gzip,itertools,json,resource,shlex,subprocess,time
from fractions import Fraction as Q
from pathlib import Path
from q_new_family_v_zero_degree26_support_preflight_20261104_v63 import read,sha
ROOT=Path(__file__).resolve().parents[1]
ID='S61_SPARSE_GENERIC5_SMALL_TAIL_FIELD_20261104_v63'
def repack_shared(kernel,plane):
    rf=[p for row in kernel['T20x2_rational']for p in row if len(p['denominator'])>1];assert rf
    D=rf[0]['denominator'];P={tuple(t[1:]):Q(int(t[0]),int(D[0][0]))for t in D}
    def num(p,poly=False):
        d=p['denominator'];dd={tuple(t[1:]):int(t[0])for t in d};nn={tuple(t[1:]):Q(int(t[0]))for t in p['numerator']}
        if not nn:return[]
        if len(dd)==1 and next(iter(dd))==(0,0,0,0,0):
            c=next(iter(dd.values()));nn={e:v/c for e,v in nn.items()}
            if poly:return[[str(v),*e]for e,v in sorted(nn.items())]
            assert len(nn)==1 and next(iter(nn))==(0,0,0,0,0);factor=next(iter(nn.values()));return[[str(factor*v),*e]for e,v in sorted(P.items())]
        c=int(d[0][0]);assert {tuple(t[1:]):Q(int(t[0]),c)for t in d}==P
        return[[str(v/c),*e]for e,v in sorted(nn.items())]
    return dict(common_denominator_P=[[str(v),*e]for e,v in sorted(P.items())],kernel20x2_numerators=[[num(p)for p in row]for row in kernel['T20x2_rational']],
        plane45x3_mixed_numerators=[[num(p)if j<2 else num(p,True)for j,p in enumerate(row)]for row in plane['V45x3']],
        plane_column_denominators=['P','P','1'],all_denominator_QQ_scalar_comparisons_verified=True,no_parameter_gcd_or_new_source_computation=True)
def make_adapter(plane,rank):
    nodes=[]
    def node(op,*args):nodes.append([op,*args]);return len(nodes)-1
    zero=node('qq','0');minus=node('qq','-1');one=node('qq','1');V=[[node('plane_rf',i,j)for j in range(3)]for i in range(45)];rows=rank['selected_original_q_rows'];B=[V[i]for i in rows]
    def times(xs):
        z=one
        for x in xs:z=node('mul',z,x)
        return z
    def plus(xs):
        z=zero
        for x in xs:z=node('add',z,x)
        return z
    parts=[]
    for p in itertools.permutations(range(3)):
        sign=(-1)**sum(p[i]>p[j]for i in range(3)for j in range(i+1,3));v=times([B[i][p[i]]for i in range(3)]);parts.append(v if sign==1 else node('mul',minus,v))
    gamma=plus(parts);gi=node('inv',gamma);adj=[[None]*3 for _ in range(3)]
    for i in range(3):
        for j in range(3):
            rr=[k for k in range(3)if k!=j];cc=[k for k in range(3)if k!=i];z=node('sub',node('mul',B[rr[0]][cc[0]],B[rr[1]][cc[1]]),node('mul',B[rr[0]][cc[1]],B[rr[1]][cc[0]]));adj[i][j]=z if (i+j)%2==0 else node('mul',minus,z)
    E=[[node('mul',gi,plus([node('mul',V[i][k],adj[k][j])for k in range(3)]))for j in range(3)]for i in range(45)];normal=[i for i in range(45)if i not in rows];zz=[node('z',j)for j in range(3)];nn=[node('n',j)for j in range(42)];center=[node('center',i)for i in range(45)];q=[]
    for i in range(45):
        v=plus([center[i],*[node('mul',E[i][j],zz[j])for j in range(3)]])
        if i in normal:v=node('add',v,nn[normal.index(i)])
        q.append(v)
    return dict(nodes=nodes,source_variables=['H','L','s0','q0','t0'],leaf_source='plane-rational45x3.json',Gamma3_node=gamma,adjugate3x3_nodes=adj,embedding45x3_nodes=E,q45_forward_nodes=q,pivot_q_rows=rows,normal_q_rows=normal,
        inverse_rule='z_j=q_pivot_j-center_pivot_j; n_i=q_nonpivot_i-center_i-sum_j E_ij*z_j',
        nonzero_Gamma_proof='rank3 factors are exact nonzero coefficient-field pivots; determinant is their product up to recorded column permutation sign',
        E_pivot_I3_proof='universal B*adj(B)=adj(B)*B=det(B)*I over a commutative field; Gamma!=0',
        no_expanded_determinant_or_adjugate=True,no_unknown_q_inverted=True,all45_unknowns_retained=True,no134_orJ0_substitution=True)
def main(out,work):
    start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False);work.mkdir(parents=True,exist_ok=False)
    pref=ROOT/'reports/S61_SPARSE_GENERIC5_INPUT_PREFLIGHT_20261104_v63';prefix=ROOT/'reports/S61_SPARSE_GENERIC5_CONSTANT_UNIT_PREFIX_20261104_v63';old=ROOT/'reports/S61_SPARSE_GENERIC5_PLANE_ADAPTER_20261104_v63';base=read(prefix/'constant-prefix.json');assert base['constant_QQ_pivots']==11;data=read(pref/'known-five-input.json.gz');oldk=read(old/'kernel20x2.json')
    for path in[pref/'report.json',prefix/'report.json',old/'report.json']:
        for n,h in read(path)['input_sha256'].items():assert sha(ROOT/n)==h,n
    inp=work/'input.txt'
    def primitive(p):return [[t[1:],str(Q(p['content'])*int(t[0]))]for t in p['terms']]
    def write(f,p):
        f.write(str(len(p))+'\n')
        for e,c in p:f.write(c+' '+' '.join(map(str,e))+'\n')
    with inp.open('w')as f:
        f.write('5\n')
        for row in base['reduced51x20']:
            for p in row:write(f,primitive(p))
        for A in[data['L51x20'],data['N45x20']]:
            for row in A:
                for p in row:write(f,p)
        for key in['q_center','Q52_minus_center']:
            for p in data[key]:write(f,p)
        f.write(' '.join(map(str,base['column_order']))+'\n');write(f,primitive(oldk['Delta']))
        for row in oldk['T_numerators20x2']:
            for p in row:write(f,primitive(p))
        f.write(' '.join(map(str,oldk['free_columns']))+'\n')
    with gzip.open(out/'native-input.txt.gz','wb')as f:f.write(inp.read_bytes())
    native=ROOT/'code/s61_sparse_generic5_small_tail_field_20261104_v63.cpp';helper=ROOT/'code/q_new_family_v_zero_degree26_restriction_20261104_v63.cpp';binary=work/'small-tail-field';flags=shlex.split(subprocess.check_output(['pkg-config','--cflags','--libs','flint'],text=True));subprocess.run(['c++','-O2','-std=c++17',str(native),*flags,'-lgmp','-o',str(binary)],check=True,timeout=max(1,60-(time.monotonic()-start)))
    result=None;error=None;adapter=False
    try:
        with(out/'native.log').open('w')as f:subprocess.run([str(binary),str(inp),str(out)],stdout=f,stderr=subprocess.STDOUT,check=True,timeout=max(1,60-(time.monotonic()-start)))
        result=read(out/'native-certificate.json');assert result['all_passed']and result['tail_rank7']and result['all51_kernel']and result['old_span_verified']
        shared=repack_shared(read(out/'generic-field-kernel.json'),read(out/'plane-rational45x3.json'))
        with gzip.open(out/'shared-denominator-kernel-plane.json.gz','wt')as f:json.dump(shared,f,separators=(',',':'))
        if result['plane_rank3_factored']:
            a=make_adapter(read(out/'plane-rational45x3.json'),read(out/'plane-rank3-factors.json'));(out/'normal42-adapter-SLP.json').write_text(json.dumps(a,indent=2)+'\n');adapter=True
        assert sum(p.stat().st_size for p in out.iterdir()if p.is_file())<32*1024*1024
    except Exception as e:error=repr(e)
    paths=[Path(__file__),native,helper,ROOT/f'research/notes/{ID}.md',pref/'report.json',pref/'known-five-input.json.gz',prefix/'report.json',prefix/'constant-prefix.json',old/'report.json',old/'kernel20x2.json',ROOT/'reports/TORIC_DEGREE5_COVARIANT_20261101.json',ROOT/'reports/TORIC_C9_POLYNOMIAL_20261101.json',ROOT/'reports/TORIC_JACOBIAN_TWO_TORSION_20261101.json']
    report=dict(id=ID,all_passed=error is None,error=error,kernel_complete=(out/'generic-field-kernel.json').exists(),old_span_verified=(out/'old-kernel-span.json').exists(),plane_rank3_factored=(out/'plane-rank3-factors.json').exists(),normal42_adapter_SLP=adapter,
        generic_coefficient_field_only_not_FULLU=True,no_oldCramer_retry_big3adjugate_expansion_point_prime_GB_or134_substitution=True,
        completed_blocks=[n for n in['row-primitive.json','field-prefix-checkpoint.json','generic-field-kernel.json','old-kernel-span.json','plane-rational45x3.json','plane-rank3-factors.json','normal42-adapter-SLP.json']if(out/n).exists()],
        input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},evidence_sha256={p.name:sha(p)for p in out.iterdir()if p.is_file()},seconds=time.monotonic()-start)
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(dict(all_passed=report['all_passed'],kernel=report['kernel_complete'],oldspan=report['old_span_verified'],plane=report['plane_rank3_factored'],adapter=adapter,seconds=report['seconds'],error=error)))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--work',type=Path,required=True);a=p.parse_args();main(a.output,a.work)
