#!/usr/bin/env python3
"""Bounded exact four-jets and true quartics from frozen constant cancellations.

This first stage uses only the 158 already certified zero cubic rows.
Every nonzero quartic is carried through the SAME constant QQ row operations;
raw quartic tails alone are not called initial equations.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
import gzip
import hashlib
import json
from pathlib import Path
import resource
import time

from q_ordinary_quadratic_cone_20261103_v63 import TruncatedQQ
from q_ordinary_cubic_initial_20261103_v63 import subtract_scaled,encode
from q_weighted_finite_model_validation_20261103_v63 import physical_equations

ROOT=Path(__file__).resolve().parents[1]


class Truncated4(TruncatedQQ):
    inverse_constants=[]

    @staticmethod
    def mul(left,right):
        if not left or not right:return {}
        a=[[]for _ in range(5)];b=[[]for _ in range(5)]
        for m,c in left.items():a[len(m)].append((m,c))
        for m,c in right.items():b[len(m)].append((m,c))
        out={}
        for i in range(5):
            for j in range(5-i):
                for m,c in a[i]:
                    for n,d in b[j]:
                        mon=tuple(sorted(m+n));value=out.get(mon,F(0))+c*d
                        if value:out[mon]=value
                        else:out.pop(mon,None)
        return out

    @classmethod
    def inv(cls,value):
        c=value.get((),F(0));assert c
        cls.inverse_constants.append(str(c))
        positive={m:-a/c for m,a in value.items()if m}
        out=cls.one;p=cls.one
        for _ in range(4):p=cls.mul(p,positive);out=cls.add(out,p)
        out={m:a/c for m,a in out.items()}
        assert cls.mul(value,out)==cls.one
        return out


def evaluate_free(nodes,inputs,keep):
    """Release every DAG value after its last use, except requested outputs."""
    refs=Counter(keep)
    for node in nodes:
        if node[0]in('+','*'):refs.update(node[1])
        elif node[0]=='inv':refs[node[1]]+=1
    values=[]
    for index,node in enumerate(nodes):
        op,arg=node
        if op=='q':value=Truncated4.q(arg);used=[]
        elif op=='v':value=inputs[arg];used=[]
        elif op in('+','*'):
            value=Truncated4.zero if op=='+'else Truncated4.one
            operation=Truncated4.add if op=='+'else Truncated4.mul
            for k in arg:assert values[k]is not None;value=operation(value,values[k])
            used=arg
        elif op=='inv':value=Truncated4.inv(values[arg]);used=[arg]
        else:raise AssertionError(op)
        values.append(value if refs[index]else None)
        for k in used:
            refs[k]-=1
            if refs[k]==0:values[k]=None
        if index%1000==0:print('DAG',index,'RSS_KB',resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,flush=True)
    return values


def unpack(data):return {tuple(m):F(c)for m,c in data}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--cpu-limit',type=int,default=240)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=False);started=time.monotonic()
    resource.setrlimit(resource.RLIMIT_AS,(3000000000,3000000000))
    resource.setrlimit(resource.RLIMIT_CPU,(args.cpu_limit,args.cpu_limit+1))
    mpath=ROOT/'reports/Q_WEIGHTED_FINITE_MODEL_20261103_v63/model.json.gz'
    qpath=ROOT/'reports/Q_ORDINARY_QUADRATIC_CONE_20261103_v63/quadrics.json.gz'
    cpath=ROOT/'reports/Q_ORDINARY_CUBIC_INITIAL_20261103_v63/cubic-initials.json.gz'
    model=json.loads(gzip.decompress(mpath.read_bytes()));oldq=json.loads(gzip.decompress(qpath.read_bytes()))
    oldc=json.loads(gzip.decompress(cpath.read_bytes()))
    inputs=[{(i,):F(1)}for i in range(35)];inputs[0][()]=F(2,9)
    keep=set(model['raw150_reconstruction_nodes']+model['all350_residual_nodes'])
    values=evaluate_free(model['nodes'],inputs,keep)
    raw=[values[k]for k in model['raw150_reconstruction_nodes']]
    residual=[values[k]for k in model['all350_residual_nodes']]
    assert all(len(m)in(2,3,4)for p in residual for m in p)
    parts=[[{m:c for m,c in p.items()if len(m)==degree}for p in residual]for degree in(2,3,4)]
    assert parts[0]==[{(i,j):F(c)for i,j,c in row}for row in oldq['all350_quadratic_rows']]
    assert parts[1]==[unpack(row)for row in oldc['raw350_cubic_terms']]
    print('ALL350 degree4 DAG complete, Q/C match frozen',time.monotonic()-started,flush=True)
    normalized={};normalized_cf={};quartic_rows=[];relations={entry['source']:unpack(entry['terms'])for entry in oldc['cubic_initial_relations']}
    zero_cubic=[]
    for op in oldc['exact_constant_row_elimination']:
        source=op['source'];d=parts[2][source].copy();cf={(source,):F(1)}
        for pivot,c in op['reductions']:
            pivot=tuple(pivot);c=F(c)
            subtract_scaled(d,normalized[pivot],c)
            subtract_scaled(cf,normalized_cf[pivot],c)
        if 'new_pivot'in op:
            factor=F(op['normalizer']);pivot=tuple(op['new_pivot'])
            normalized[pivot]={m:c/factor for m,c in d.items()}
            normalized_cf[pivot]={m:c/factor for m,c in cf.items()}
        else:
            quartic_rows.append({'source':source,'cubic':relations[source],'quartic':d,'cofactors':cf})
            if not relations[source]:
                zero_cubic.append(source)
                # Fresh explicit equality in each homogeneous degree, including
                # all lower cancellations, is the true-initial certificate.
                for degree in range(3):
                    total={}
                    for (k,),c in cf.items():subtract_scaled(total,parts[degree][k],-c)
                    assert total==(d if degree==2 else {})
    assert len(quartic_rows)==293 and len(zero_cubic)==158
    checkpoint={'variables':oldq['variables'],'old_weights':model['input_weights'],
                'all350_quartic_terms':list(map(encode,parts[2])),
                'all293_constant_quadratic_cancellation_rows':[
                    {'source':row['source'],'cubic_terms':encode(row['cubic']),
                     'quartic_terms':encode(row['quartic']),
                     'constant_full350_cofactors':[[m[0],str(c)]for m,c in sorted(row['cofactors'].items())]}
                    for row in quartic_rows]}
    cp=args.output/'quartic-tails.json.gz'
    with gzip.open(cp,'wt')as handle:json.dump(checkpoint,handle,separators=(',',':'))
    print('CHECKPOINT quartic tails; true zeroC quartics',sum(bool(row['quartic'])for row in quartic_rows if not row['cubic']),
          'terms',sum(len(row['quartic'])for row in quartic_rows if not row['cubic']),time.monotonic()-started,flush=True)
    # Full physical differentiation/multiplication uses no saved equation rows.
    divergence,physical=physical_equations(raw,Truncated4)
    assert not any(divergence)and physical==residual
    print('ALL150 physical divergence/ALL350 degree4 Jacobi replay PASS',time.monotonic()-started,flush=True)
    inputs=[Path(__file__),mpath,qpath,cpath,ROOT/'code/q_ordinary_quadratic_cone_20261103_v63.py',
            ROOT/'code/q_ordinary_cubic_initial_20261103_v63.py',ROOT/'code/q_weighted_finite_model_validation_20261103_v63.py']
    report={'id':'Q_ORDINARY_QUARTIC_INITIAL_20261103_V63','all_passed':True,
            'coefficient_field':'QQ','ordinary_truncation':4,'refcount_DAG_freeing':True,
            'all350_Q_C_match_frozen':True,'all350_degree4_physical_replay':True,
            'all_inverse_Taylor_products_verified_to4':True,'constant_cubic_rows':293,
            'identically_zero_cubic_rows':158,
            'true_quartics_from_zero_cubic_rows':sum(bool(row['quartic'])for row in quartic_rows if not row['cubic']),
            'true_quartic_term_count':sum(len(row['quartic'])for row in quartic_rows if not row['cubic']),
            'all_saved_zeroC_rows_explicit_lower2_lower3_cancellations_verified':True,
            'full350_constant_cofactors_saved':True,'all293_quartic_carries_saved_for_followup_kernels':True,
            'full_quartic_initial_space_claimed':False,'local_dimension_upper12_claimed':False,
            'memory_limit_bytes':3000000000,'cpu_limit_seconds':args.cpu_limit,
            'peak_rss_kb':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'seconds':time.monotonic()-started,
            'quartic_tails_sha256':hashlib.sha256(cp.read_bytes()).hexdigest(),
            'input_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in inputs}}
    (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
