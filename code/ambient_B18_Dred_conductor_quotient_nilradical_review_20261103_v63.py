#!/usr/bin/env python3
"""Pure nilradical corollary of the independently accepted conductor sequence."""
import hashlib
import json
from pathlib import Path
import time

ROOT=Path(__file__).resolve().parents[1]
ID='AMBIENT_B18_DRED_CONDUCTOR_QUOTIENT_NILRADICAL_20261103_v63_REVIEW'

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    started=time.monotonic()
    depth=ROOT/'reports/AMBIENT_B18_DRED_CONDUCTOR_QUOTIENT_DEPTH_20261103_v63_REVIEW/report.json'
    conductor=ROOT/'reports/AMBIENT_B18_DRED_CONDUCTOR_20261103_v63_REVIEW/report.json'
    N=ROOT/'reports/AMBIENT_B18_N_FULL_COMPLETION_20261103_v63_REVIEW.json'
    checked={};reports={}
    for path in [depth,conductor,N]:
        report=json.loads(path.read_text());assert report['all_passed'];reports[path]=report
        for name,want in report.get('input_sha256',{}).items():
            assert sha(ROOT/name)==want,name;checked[name]=want
    assert reports[conductor]['IB_N_over_J_length4_transverse_and_rank4_over_C_g19']
    assert reports[N]['exact_I']=='(nu,s1,...,s6,w*v,v^2)'
    assert reports[depth]['cone_conductor_quotient_dimension_depth_multiplicity']==[20,19,1]
    note=ROOT/'research/notes'/(ID+'.md');output=ROOT/'reports'/ID
    output.mkdir(parents=True,exist_ok=False)
    files=[Path(__file__),note,depth,conductor,N]
    report={
        'id':ID,'all_passed':True,
        'E':'Dred/conductor = A_N/J, G=C[[g1,...,g19]]',
        'exact_sequence':'0 -> M=IB_N/J (G-free rank4) -> Nil(E) -> G*v -> 0',
        'reduced_conductor_quotient':'G[[w]]',
        'reduced_cone_dimension':20,'reduced_projective_dimension':19,
        'nilradical_as_regular_G_module_free_rank':5,
        'nilradical_reduced_support_dimension_cone':19,
        'nilradical_reduced_support_dimension_projective':18,
        'transverse_1D_core_H0_length':5,
        'transverse_core_mod_H0':'C[[w]]',
        'proof':'M is a finite-length ideal in the 1D transverse local ring, hence nilpotent. Its quotient has nilradical(v), so the preimage is exactly Nil(E). The G-free rank1 quotient splits the module sequence, giving rank5. The quotient by Nil(E) is a DVR transversely, so H0=Nil(E).',
        'ring_split_nilpotency_index_embedding_dimension_full_D_or_global_extension_claim':False,
        'new_jet_CAS_or_source_replay':False,
        'scope':'Formal generic B18 normalization-conductor quotient of Dred only.',
        'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in files},
        'checked_dependency_sha256':checked,'seconds':time.monotonic()-started,
    }
    (output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'all_passed':True,'nilradical_G_free_rank':5,'seconds':report['seconds']}))

if __name__=='__main__':main()
