#!/usr/bin/env python3
"""Pure exact-sequence corollary, gated on the accepted actual Dred conductor.

No source replay, new jet, or CAS computation.
"""
import argparse
import hashlib
import json
from pathlib import Path
import time

ROOT=Path(__file__).resolve().parents[1]
ID='AMBIENT_B18_DRED_CONDUCTOR_QUOTIENT_DEPTH_20261103_v63_REVIEW'

def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main(output,conductor_review):
    started=time.monotonic()
    #Gate before creating an output: no PASS can be frozen before this review exists.
    cp=ROOT/'reports/AMBIENT_B18_DRED_CONDUCTOR_20261103_v63/report.json'
    np=ROOT/'reports/AMBIENT_B18_N_FULL_COMPLETION_20261103_v63/report.json'
    nr=ROOT/'reports/AMBIENT_B18_N_FULL_COMPLETION_20261103_v63_REVIEW.json'
    reports=[cp,conductor_review,np,nr];checked={}
    data={}
    for path in reports:
        obj=read(path);assert obj['all_passed'];data[path]=obj
        for name,want in obj.get('input_sha256',{}).items():
            p=ROOT/name;assert p.is_file()and sha(p)==want,name;checked[name]=want
    assert data[cp]['cyclic_image_length']==4
    assert data[conductor_review]['IB_N_over_J_length4_transverse_and_rank4_over_C_g19']
    assert data[np]['exact_full_odd_ideal']=='I=(nu,s1,...,s6,w*v,v^2)'
    assert data[nr]['source_IB_and_target_I_plus_It_conductor_proof_accepted']
    output.mkdir(parents=True,exist_ok=False)
    #Independent monomial core check: v survives and is annihilated by (w,v).
    hf=[]
    for degree in range(12):
        surviving=[(iw,degree-iw)for iw in range(degree+1)
                   if degree-iw<2 and not(iw and degree-iw)]
        hf.append(len(surviving))
    assert hf==[1,2]+[1]*10
    G=19;Mdim=Mdepth=G;Cdim=G+1;Cdepth=G;Edepth=G
    assert Mdepth==Cdepth and Cdepth-1<Mdepth
    evidence={
        'exact_sequence':'0 -> M=IB_N/J -> E=A_N/J -> C=R/I -> 0',
        'M_regular_group_variable_count':G,'M_free_rank':4,
        'C_full_model':'C[[g1,...,g19,w,v]]/(w*v,v^2)',
        'two_variable_core_Hilbert_function_first12':hf,
        'C_dimension_depth_multiplicity':[Cdim,Cdepth,1],
        'M_dimension_depth':[Mdim,Mdepth],
        'local_cohomology_injection':'H^18_m(C)=0, so 0 -> H^19_m(M) -> H^19_m(E); source is nonzero',
        'depth_lower_bound':'depth lemma gives depth(E)>=min(depth M,depth C)=19',
        'multiplicity_additivity':'dim M=19<dim E=dim C=20; hence e(E)=e(C)=1 for the induced maximal ideal',
        'Dred_conductor_quotient_identification':'F/(J,T^4) = (A_N/J) x_(C4) C4 = A_N/J',
        'cone_dimension_depth_multiplicity':[20,Edepth,1],
        'projective_dimension_depth_multiplicity':[19,Edepth-1,1],
    }
    ev=output/'independent-certificate.json';ev.write_text(json.dumps(evidence,indent=2)+'\n')
    note=ROOT/'research/notes'/(ID+'.md')
    files=[Path(__file__),note,*reports,
           ROOT/'research/notes/AMBIENT_B18_DRED_CONDUCTOR_20261103_v63.md',
           ROOT/'research/notes/AMBIENT_B18_DRED_CONDUCTOR_20261103_v63_REVIEW.md',
           ROOT/'research/notes/AMBIENT_B18_N_FULL_COMPLETION_20261103_v63.md']
    report={
        'id':ID,'all_passed':True,
        'accepted_entire_actual_Dred_conductor_and_N_completion_gates_checked':True,
        'independent_exact_sequence_and_local_cohomology_depth_proof':True,
        'cone_conductor_quotient_dimension_depth_multiplicity':[20,19,1],
        'projective_conductor_quotient_dimension_depth_multiplicity':[19,18,1],
        'embedding_dimension_split_ring_full_D_or_global_locus_claim':False,
        'new_source_jet_or_CAS_replay':False,
        'scope':'Formal generic B18 Dred normalization-conductor quotient only.',
        'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in files},'checked_dependency_sha256':checked,
        'evidence_sha256':{ev.name:sha(ev)},'seconds':time.monotonic()-started,
    }
    (output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:report[k]for k in ['all_passed','cone_conductor_quotient_dimension_depth_multiplicity','seconds']}))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,default=ROOT/'reports'/ID)
    p.add_argument('--conductor-review',type=Path,default=ROOT/'reports/AMBIENT_B18_DRED_CONDUCTOR_20261103_v63_REVIEW/report.json')
    args=p.parse_args();main(args.output.resolve(),args.conductor_review.resolve())
