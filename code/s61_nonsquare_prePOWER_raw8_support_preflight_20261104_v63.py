#!/usr/bin/env python3
"""Source-only total-degree/lattice/lifetime estimate for ONE fixed raw8 minor."""
import argparse,hashlib,itertools,json,math,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_PREPOWER_RAW8_SUPPORT_PREFLIGHT_20261104_v63'
def read(p):return json.loads(p.read_text())
def sha(p):
    with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def dump(p,d):p.write_text(json.dumps(d,indent=2)+'\n')
def main(out):
    started=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(parents=True,exist_ok=False);pins={}
    dump(out/'phase.json',dict(phase='source-only fixed raw8 preflight; no scalar job',fixed_current15_rows=[0,1,3,6,2,4,5,8],fixed_original32_columns=[22,11,13,25,15,19,27,31]))
    def load(st,name='certificate.json'):
        p=ROOT/'reports'/st/name;pins[str(p.relative_to(ROOT))]=sha(p);return read(p)
    q=load('S61_NONSQUARE_CUBIC_Q2_L_UNIT_COMPLETION_20261104_v63');cp=load('S61_NONSQUARE_CUBIC_ZW25_UNITS_20261104_v63','checkpoint.json');pd=load('S61_NONSQUARE_CUBIC_ZW25_PAIR_DUAL_PREFLIGHT_20261104_v63');h=load('S61_NONSQUARE_PURE_KERNEL12_SOURCE_PREFLIGHT_20261104_v63','prePOWER-h9-source.json');fx=load('S61_NONSQUARE_PREPOWER_H22_PURE3_FIXED_MINOR_20261104_v63','native.json');st=load('S61_NONSQUARE_PREPOWER_H22_PURE4_NAMED_STRIP_20261104_v63','native.json');direct=load('S61_NONSQUARE_PREPOWER_DIRECT_P5_H5_SCHUR_20261104_v63','operator.json')
    assert direct['context']['no_R434_D4_h22_inverse'];oldN=q['scalar_DAG']['nodes'];cols=q['new22_target_original32_columns'];rows=[0,1,3,6,2,4,5,8];labels=[22,11,13,25,15,19,27,31]
    cachepoly={x['old_DAG_root']:x['Laurent_polynomial']for x in fx['node_values']};DAG=[];atomkeys={};cache={};source_atom_terms=0
    def lattice(box,lo,hi):return sum(max(0,min(box[2][1],hi-a-e)-max(box[2][0],lo-a-e)+1)for a in range(box[0][0],box[0][1]+1)for e in range(box[1][0],box[1][1]+1))
    def put(op,args,bound,den,norm,**kw):
        i=len(DAG);DAG.append(dict(id=i,op=op,args=args,bound=bound,QQ_denominator=str(den),integer_l1_upper=str(norm),**kw));return i
    def atom(poly,ref,unit=Q(1),ashift=0):
        nonlocal source_atom_terms
        p=[([m[0]+ashift,m[1],m[2]],Q(c)/unit)for m,c in poly if Q(c)]
        key=tuple((tuple(m),str(c))for m,c in p)
        if key in atomkeys:return atomkeys[key]
        D=math.lcm(*(c.denominator for m,c in p))if p else 1;N=sum(abs(c.numerator)*(D//c.denominator)for m,c in p)
        if not p:b=None
        else:b=dict(support=len(p),degree=[min(sum(m)for m,c in p),max(sum(m)for m,c in p)],box=[[min(m[j]for m,c in p),max(m[j]for m,c in p)]for j in range(3)])
        i=put('atom',[],b,D,N,reference=ref,QQ_unit=str(unit),a_shift=ashift);atomkeys[key]=i;source_atom_terms+=len(p);return i
    def const(x):return atom([[[0,0,0],str(Q(x))]],dict(constant=str(x)))
    z=const(0);one=const(1)
    def bound(i):return DAG[i]['bound']
    def dn(i):return int(DAG[i]['QQ_denominator']),int(DAG[i]['integer_l1_upper'])
    def add(i,j,sign=1):
        if i==z:return j if sign==1 else mul(const(-1),j)
        if j==z:return i
        if sign==-1 and i==j:return z
        A,B=bound(i),bound(j);a,n=dn(i);b,m=dn(j);D=math.lcm(a,b);N=n*(D//a)+m*(D//b)
        lo=min(A['degree'][0],B['degree'][0]);hi=max(A['degree'][1],B['degree'][1]);box=[[min(A['box'][k][0],B['box'][k][0]),max(A['box'][k][1],B['box'][k][1])]for k in range(3)]
        support=min(A['support']+B['support'],lattice(box,lo,hi),N)
        return put('add'if sign==1 else'sub',[i,j],dict(support=support,degree=[lo,hi],box=box),D,N)
    def mul(i,j):
        if i==z or j==z:return z
        if i==one:return j
        if j==one:return i
        A,B=bound(i),bound(j);a,n=dn(i);b,m=dn(j);lo=A['degree'][0]+B['degree'][0];hi=A['degree'][1]+B['degree'][1];box=[[A['box'][k][0]+B['box'][k][0],A['box'][k][1]+B['box'][k][1]]for k in range(3)]
        support=min(A['support']*B['support'],lattice(box,lo,hi),n*m)
        return put('mul',[i,j],dict(support=support,degree=[lo,hi],box=box),a*b,n*m)
    def shift(i,k):
        if i==z or k==0:return i
        A=bound(i);box=[x[:]for x in A['box']];box[0]=[x+k for x in box[0]];D,N=dn(i)
        return put('known_a_shift',[i],dict(support=A['support'],degree=[x+k for x in A['degree']],box=box),D,N,a_exponent_shift=k)
    def imp(i):
        if i in cache:return cache[i]
        n=oldN[i];op=n[0]
        if i in cachepoly:x=atom(cachepoly[i],dict(path='reports/S61_NONSQUARE_PREPOWER_H22_PURE3_FIXED_MINOR_20261104_v63/native.json',old_DAG_root=i))
        elif op=='z':x=z
        elif op=='one':x=one
        elif op in ['a','s','m']:
            a,b=imp(n[1]),imp(n[2]);x=mul(a,b)if op=='m'else add(a,b,-1 if op=='s'else 1)
        else:
            num,den=cp['matrix25x32'][n[1]][n[2]]if op=='checkpoint'else pd['Lframe_dual6x2'][n[1]][1]if op=='accepted_pair_second'else q['scalar_DAG']['new_unit_atoms_num_den'][n[1]]
            assert len(den)==1 and den[0][0][1:]==[0,0],('Unexpected B-pole in SAMEselectedraw8 closure',i,n)
            x=atom(num,dict(path='accepted_Q2L_source_atom',old_DAG_root=i,definition=n),Q(den[0][1]),-den[0][0][0])
        cache[i]=x;return x
    matrix=[]
    for ii,r in enumerate(rows):
        row=[]
        for j in labels:
            if ii<4 and j in [27,31]:x=z
            elif r==0 and j in h['original_target_labels']:x=atom(h['source_h9'][h['original_target_labels'].index(j)],dict(path='reports/S61_NONSQUARE_PURE_KERNEL12_SOURCE_PREFLIGHT_20261104_v63/prePOWER-h9-source.json',h9_target=j))
            else:x=imp(q['remaining15x22_scalar_roots'][r][cols.index(j)])
            row.append(x)
        matrix.append(row)
    # Exact known a factors only, no C/p/Delta/R factor assumptions.
    rowa=[min(bound(x)['box'][0][0]for x in row if x!=z)for row in matrix]
    matrix=[[shift(x,-a)for x in row]for row,a in zip(matrix,rowa)]
    assert all(bound(x)['box'][0][0]>=0 for row in matrix for x in row if x!=z)
    def det2(r1,r2,c1,c2):return add(mul(matrix[r1][c1],matrix[r2][c2]),mul(matrix[r1][c2],matrix[r2][c1]),-1)
    heavy={rr:det2(rr[0],rr[1],6,7)for rr in itertools.combinations(range(4,8),2)}
    def det4(chosen):
        if tuple(chosen)==(0,1,2,3):return atom(st['exact_unscaled_pure4_minor'],dict(path='reports/S61_NONSQUARE_PREPOWER_H22_PURE4_NAMED_STRIP_20261104_v63/native.json',key='exact_unscaled_pure4_minor'),ashift=-sum(rowa[:4]))
        s=z
        for pp in itertools.permutations(range(4)):
            t=one
            for r,c in enumerate(pp):t=mul(t,matrix[r][chosen[c]])
            parity=sum(pp[i]>pp[j]for i in range(4)for j in range(i+1,4));s=add(s,t,-1 if parity%2 else 1)
        return s
    output=z;laplace=[]
    for chosen in itertools.combinations(range(6),4):
        rem=[j for j in range(6)if j not in chosen];bottom=z
        for rr,H in heavy.items():
            pure_rows=[i for i in range(4,8)if i not in rr];v=det2(*pure_rows,*rem);sgn=-1 if(5+sum(i-4 for i in rr))%2 else 1;bottom=add(bottom,mul(H,v),sgn)
        top=det4(chosen);term=mul(top,bottom);output=add(output,term,-1 if sum(chosen)%2 else 1)
        laplace.append(dict(pure_top4_columns=list(chosen),pure_bottom2_columns=rem,top4_root=top,bottom4_root=bottom,term_root=term))
    needed=set();todo=[output]
    while todo:
        i=todo.pop()
        if i in needed:continue
        needed.add(i);todo+=DAG[i]['args']
    use={i:0 for i in needed}
    for i in needed:
        for j in DAG[i]['args']:use[j]+=1
    use[output]+=1;live=peak=0;peakid=0;records=[]
    for i in sorted(needed):
        n=DAG[i];s=n['bound']['support']if n['bound']else 0;live+=s;scratch=s if n['op']in['mul','add','sub']else 0
        if live+scratch>peak:peak=live+scratch;peakid=i
        records.append(dict(node=i,live_upper=live,scratch_upper=scratch))
        for j in n['args']:
            use[j]-=1
            if not use[j]:live-=DAG[j]['bound']['support']if DAG[j]['bound']else 0
    actual_atoms={i for i in needed if DAG[i]['op']=='atom'};atom_reserve=sum(DAG[i]['bound']['support']if DAG[i]['bound']else 0 for i in actual_atoms)
    D,N=dn(output);ob=bound(output);maxnum=max([N]+[int(DAG[i]['integer_l1_upper'])for i in actual_atoms]);maxden=max([D]+[int(DAG[i]['QQ_denominator'])for i in actual_atoms])
    bytes_upper=(ob['support']+atom_reserve)*(len(str(maxnum))+len(str(maxden))+110)+1000000
    total=peak+atom_reserve+15000
    recipe=dict(source_current15_rows=rows,source_original32_columns=labels,raw_input_roots=[[q['remaining15x22_scalar_roots'][r][cols.index(j)]for j in labels]for r in rows],
      accepted_top4_heavy_literal_zero=True,row_known_a_factors=rowa,represented_raw8='detraw8=a^(sum rowa)*output',
      fixed_border6_to_raw8_multiplier='rhoLCM4*rhoLCM7*delta2/(rho1*rho3*rho4*rho7*gamma)',
      no_p_inverse_in_polynomial_identity=True,all8_other_Q9_tails_retained_in_future_complement=True,
      future_pivot11=dict(rows=rows+[7,13,10],original32_columns=labels+[28,30,29],upper_right8x3_literal_zero=True,determinant='raw8det*(rho6*rho12*rho9)*Delta407; no p/R/Saux/h22inverse',remaining_current15=[9,11,12,14],remaining_original32=[1,2,3,4,6,7,8,9,14,17,18]),
      scalar_nodes=DAG,output=output,laplace15=laplace,heavy2_roots=[dict(rows=list(rr),root=v)for rr,v in heavy.items()],
      all_selected_B_denominator_powers_zero=True,all_source_denominators_a_and_QQ_only=True,
      exact_parameter_products_or_coefficients_NOT_evaluated=True)
    dump(out/'operator-and-cost.json',recipe)
    cert=dict(totaldegree_and_box_lattice_count=True,rigorous_l1_zero_only=True,all_8_row_a_factors_by_lower_exponent_only=rowa,
      all_source_shapes_and_indices_fixed=True,live_records=records,peak_node=peakid,last_use_peak_upper=peak,literal_atom_reserve=atom_reserve,work_reserve=15000,total_live_upper=total,
      output_bound=ob,output_integer_l1_bits=N.bit_length(),QQ_denominator_bits=D.bit_length(),raw_output_and_atom_artifact_bytes_upper=bytes_upper,
      resource_gate_passes=total<=150000 and bytes_upper<32*1024*1024,
      no_actual_minor_nonzero_or_4x11rank_claim=True)
    dump(out/'certificate.json',cert)
    note=ROOT/f'research/notes/{ID}.md';pins[str(note.relative_to(ROOT))]=sha(note);pins[str(Path(__file__).relative_to(ROOT))]=sha(Path(__file__))
    size=sum(p.stat().st_size for p in out.iterdir());assert size<32*1024*1024
    report=dict(id=ID,all_passed=True,frozen=True,source_only=True,resource_gate_passes=cert['resource_gate_passes'],typed_scalar_nodes=len(DAG),live_nodes=len(needed),
      row_known_a_factors=rowa,output_support_upper=ob['support'],output_degree_range=ob['degree'],output_integer_l1_bits=N.bit_length(),QQ_denominator_bits=D.bit_length(),
      last_use_plus_source_scratch_reserve_upper=total,raw_output_and_atom_artifact_bytes_upper=bytes_upper,
      no_new_actual_value_minor_rank_point_prime_GB=True,input_sha256=pins,evidence_sha256={p.name:sha(p)for p in out.iterdir()},seconds=time.monotonic()-started)
    dump(out/'report.json',report);signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in ['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
