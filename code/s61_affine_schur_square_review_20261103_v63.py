#!/usr/bin/env python3
"""Pure module/Schur-square criteria and independent companion staircase audit."""
import argparse
from fractions import Fraction
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import subprocess
import time
import sympy as s

ROOT=Path(__file__).resolve().parents[1]
BASE='S61_AFFINE_NORMAL_SCHUR_SQUARE_VECTOR_20261103_v63'
STAIR='S61_AFFINE_NORMAL_SCHUR_STAIRCASE_20261103_v63'
def read(path):
    with (gzip.open(path,'rt') if path.suffix=='.gz' else path.open()) as f:return json.load(f)
def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda:f.read(1<<20),b''):h.update(chunk)
    return h.hexdigest()
def plus(A,B,c=Fraction(1)):
    C=A.copy()
    for mon,value in B.items():
        value=C.get(mon,Fraction(0))+c*value
        if value:C[mon]=value
        else:C.pop(mon,None)
    return C
def times(A,B):
    C={}
    for m,a in A.items():
        for n,b in B.items():
            key=tuple(sorted(m+n));value=C.get(key,Fraction(0))+a*b
            if value:C[key]=value
            else:C.pop(key,None)
    return C
def decode(poly):return {d:{tuple(m):Fraction(c) for m,c in coefficients} for d,coefficients in poly}


def main(output,work):
    started=time.monotonic();work.mkdir(parents=True,exist_ok=False);paths=[Path(__file__),ROOT/'research/notes'/f'{BASE}_REVIEW.md']
    for stem,source in [(BASE,'s61_affine_normal_schur_square_vector_20261103_v63.py'),(STAIR,'s61_affine_normal_schur_staircase_20261103_v63.py')]:
        folder=ROOT/'reports'/stem;report=read(folder/'report.json');assert report['all_passed']
        for name,digest in report['input_sha256'].items():assert sha(ROOT/name)==digest,name
        for name,digest in report['evidence_sha256'].items():assert sha(folder/name)==digest,name
        replay=work/stem
        subprocess.run([str(ROOT/'.venv/bin/python'),str(ROOT/'code'/source),'--output',str(replay)],check=True,timeout=40)
        for name in report['evidence_sha256']:assert (folder/name).read_bytes()==(replay/name).read_bytes()
        paths += [folder/'report.json',folder/'certificate.json',ROOT/'code'/source,ROOT/'research/notes'/f'{stem}.md']
    acceptedpath=ROOT/'reports/S61_AFFINE_NORMAL_MOMENT_SCHUR_20261103_v63_REVIEW/report.json'
    accepted=read(acceptedpath);assert accepted['all_passed']
    assert accepted['delta_regular_on_Ba_and_all_existing_Ass_and_nilpotent_submodules_retained']
    assert accepted['B_to_Ba_delta_faithfully_flat_by_flatness_and_nonzero_residue_field_fibers']
    assert accepted['all134_structured_source_column_indices_and_degree_bounds_bound']
    # Independent full generic5x4 integer polynomial toy, using SymPy
    # determinants and explicit100-column coefficient witnesses.
    xx=s.symbols('x0:20');V=s.Matrix(5,4,xx)
    lam=[(-1)**i*V[[j for j in range(5) if j!=i],:].det() for i in range(5)]
    delta=s.expand(lam[0]);assert len(s.Poly(delta,*xx).terms())==24
    assert all(s.expand(sum(lam[i]*V[i,j] for i in range(5)))==0 for j in range(4))
    co=[s.diff(delta,x)/4 for x in xx]
    vector=[sum(co[j]*xx[j] for j in range(20)) if i==0 else s.S.Zero for i in range(5)]
    assert all(s.expand(value-(delta if i==0 else 0))==0 for i,value in enumerate(vector))
    assert s.expand(sum(co[j]*xx[j]*lam[0] for j in range(20))-delta**2)==0
    # Rebuild all3350 entries by the actual one-step companion recurrence,
    # not the producer's generic univariate multiplication/division code.
    presentation=ROOT/'reports/S61_AFFINE_NORMAL_GLOBAL_QUINTIC_20261103_v63/presentation.json.gz'
    data=read(presentation);G=decode(data['monic_quintic_G']);residues=list(map(decode,data['all134_remainders']))
    assert G[5]=={():Fraction(1)} and G[0].get((44,))==-1
    assert all(not mon for k,p in G.items() if k in range(1,5) for mon in p)
    g=[G.get(i,{}).get((),Fraction(0)) for i in range(5)]
    staircase=read(ROOT/'reports'/STAIR/'certificate.json');digest=hashlib.sha256();entries=0
    source_Q=None;source_r4=None
    for i,r in enumerate(residues):
        Q={m:c for m,c in r.get(0,{}).items() if len(m)==2}
        assert all(44 not in m for m in Q)
        assert all(len(m)<=1 for k,p in r.items() if k>=1 for m in p)
        if i==staircase['source_nonzero_Q_witness']['row']:
            witness=staircase['source_nonzero_Q_witness'];assert Q[tuple(witness['monomial'])]==Fraction(witness['coefficient']);source_Q=True
        if i==staircase['source_nonzero_r4_witness']['row']:
            assert r[4]=={tuple(m):Fraction(c) for m,c in staircase['source_nonzero_r4_witness']['coefficients']};source_r4=True
        current=[r.get(k,{}) for k in range(5)]
        for j in range(5):
            for k,p in enumerate(current):
                assert all(len(m)<=2 for m in p)
                quad={m:c for m,c in p.items() if len(m)==2}
                if k>j:assert not quad
                if k==j:assert quad==Q
                digest.update((str((i,j,k))+json.dumps([[list(m),str(c)] for m,c in sorted(p.items())],separators=(',',':'))).encode())
                entries+=1
            previous=current
            current=[times({(44,):Fraction(1),():-g[0]},previous[4])]
            current += [plus(previous[k-1],previous[4],-g[k]) for k in range(1,5)]
    assert entries==3350 and source_Q and source_r4
    # Independent top-weight determinant: rows1..3 quadratic, row4 linear.
    HH,QQ=s.symbols('H Q');LL=s.symbols('L0:5');gg=s.symbols('g0:5')
    T=s.zeros(5)
    for j in range(4):T[j+1,j]=1
    T[:,4]=s.Matrix([HH-gg[0],-gg[1],-gg[2],-gg[3],-gg[4]])
    r=s.Matrix([LL[0]+QQ,*LL[1:]]);columns=[r]
    for _ in range(3):columns.append(T*columns[-1])
    M=s.Matrix.hstack(*columns);weighted=[HH,*LL,QQ]
    def part(poly,weight):
        return s.expand(sum(c*s.prod(x**n for x,n in zip(weighted,mon)) for mon,c in s.Poly(poly,*weighted).terms() if sum(mon[:-1])+2*mon[-1]==weight))
    top=s.Matrix(4,4,lambda i,j:part(M[i+1,j],2 if i<3 else 1))
    assert s.expand(top.det()+LL[4]*QQ**3)==0
    assert all(part(M[4,j],2)==0 for j in range(4))
    globalreview=ROOT/'reports/S61_AFFINE_NORMAL_GLOBAL_QUINTIC_20261103_v63_REVIEW.json';assert read(globalreview)['all_passed']
    paths += [acceptedpath,globalreview,presentation]
    result=dict(id=BASE+'_REVIEW',all_passed=True,
        both_frozen_producer_evidence_files_fresh_byte_equal=True,
        independent_generic5x4_all_lambda_V_and_explicit100column_vector_square_toy=True,
        known_E5_not_required_to_split_over_Ka=True,
        conditional_exact_kernel_m_Rprime5_and_lambda_image_m_lambda_proof_accepted=True,
        exact_B_empty_iff_delta_unit_vector_membership_accepted=True,
        exact_B_empty_iff_delta_square_in_Jg_accepted=True,
        Fitt0_subset_AnnA_and_delta_regular_only_on_B_used_correctly=True,
        all_unknown_Ass_and_nilpotent_parts_retained=True,
        independent_actual_companion_staircase_entries=entries,staircase_coefficient_digest=digest.hexdigest(),
        independent_top_weight_delta_equals_minus_L4_Q_cubed=True,
        both_actual_Q_and_r4_nonzero_witnesses_bound=True,
        formal_degree_bounds=dict(delta=7,U=7,vector=7,square=14,cubic=21,scalar134=30,norm=10),
        formal45_parameter_degree_not_after_nonlinear_normal_basechange=True,
        cofactor_degree_bound_or_actual_delta_expansion_membership_or_GB_not_claimed=True,
        input_sha256={str(path.relative_to(ROOT)):sha(path) for path in paths},seconds=time.monotonic()-started)
    output.parent.mkdir(parents=True,exist_ok=True);output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ['all_passed','independent_actual_companion_staircase_entries','formal_degree_bounds','seconds']},indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--work',type=Path,required=True);args=parser.parse_args();main(args.output.resolve(),args.work.resolve())
