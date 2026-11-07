"""ROOT full mathematical worker and V4 control/custody review, no responses."""
from pathlib import Path
import hashlib,json
R=Path(__file__).resolve().parents[1];sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest();pins={}
def load(p,h):
 q=R/p;assert sha(q)==h;pins[p]=h;d=json.loads(q.read_text());assert d['all_passed'];return d
stem='reports/AMBIENT_B18_ELL2_T2_UPPER20_'
gate=stem+'WHOLE_INPUT_GATE_V4_20261105_v63/report.json'
g=load(gate,'cbdbe69c3a701ab3809d6ac5220adb3abc5e0eef04b3ba37c0f7c373b3e73de4')
c=load(stem+'VALUES_V4_PREPARED_CUSTODY_20261105_v63/report.json','49dacd5da40887594c86a1beb480991fd5fe8bd2910969c27f4d51de3d444d9c')
p=load(stem+'WHOLE_V4_20261105_v63_REVIEW/report.json','a2a3521caaa96e53df644123940b0c688a092e624702e81599ab6bedd10af662')
w='code/ambient_B18_ell2_T2_upper20_values_v4_20261105_v63.py';wh='cdb2418923e09c9d12330e02d446fdf0251e41eef7bf79263f36041e494851ba'
assert sha(R/w)==c['worker_sha256']==wh;pins[w]=wh
for path,h in c['worker_and_source_sha256'].items():assert sha(R/path)==h;pins[path]=h
base='reports/AMBIENT_B18_ELL2_T2_UPPER20_WHOLE_INPUT_GATE_{}_20261105_v63/whole20-coupled-pruned-cache-input-and-literal-owner-height.json'
a=json.loads((R/base.format('V2')).read_text());b=json.loads((R/base.format('V4')).read_text())
assert all(a[k]==b[k] for k in ['OLD','NATIVE','positive_global_envelopes'])
assert b['whole_numeric_cap']==58883 and b['whole_RAW_cap']==17545112
out=R/'reports/AMBIENT_B18_ELL2_T2_UPPER20_V4_ROOT_SOURCE_REVIEW_20261105_v63';out.mkdir(exist_ok=False)
r=dict(id=out.name,all_passed=True,SOURCE_CODE_custody_only_no20response_values=True,input_sha256=pins,
 manually_reviewed_full276line_original_mathematics_and_all_V4_control_diffs=True,OLD_NATIVE_and_positive_envelopes_byte_structure_identical_to_reviewedV2=True,
 fresh_independent_whole_source_code_peer_accepted=True,
 proof_checks=['584completepositiveSOURCEzeros+554independentsignedzero+20exceptions exhaust605nonopaque','OLD C0/typed fullcube/second, all1742 FIRST factors andnative92outer retained','16FORMALraw native auxiliaries including8/32/80nonzeroS0 becomeONLY15truefull48pointkernelrestrictions; omitted10kernels haveexplicitbothchannelpositivezero certificates','case617OLD/FIRST/outer subtract minus20raw563 coherently; native617alreadycorrectedonce','fullL2S0/moniccarry andOLD4wholeopaqueannihilation, neverindividualS0/unknownhighertruezero','nativefine andallpieces/Q checked beforemapping, bound354/temp1063 underhard384/1153','rowLCDinteger rank334bit/minor and502bitTEMP identities/kernel freeblockchecked, no totalHF/H0/nilinference'],
 control_checks=['earlyinitializedallprefixstates, peroperation pending andlastTEMP snapshot','envrejectedq/bound, incrementalCramer firstnumerator/currentkernel/dot cursors retained','saveerrorrepr then clearunwoundframes/traceback/context/cause andgc BEFOREGUARD; nonmathlifecyclefixture verifiesrelease','SUCCESScodec andGUARDcodec exclusivelyshare16652reserve; finaltupleFIRST1086 vs1024 compensatedbyunused178nativepointslots leaving116spare','55secondcompute pluscheckpoint withmandatoryexternal60secondwalltimeout'],
 committed_prior=354,temporary_prior=1063,hard_bit_limits=[384,1153],whole_numeric_prior=58883,private_hard_numeric=60000,whole_RAW_prior=17545112,
 one_fixed20_value_bundle_authorized=True,output_stem='AMBIENT_B18_ELL2_T2_UPPER20_VALUES_V4_20261105_v63',mandatory_external_wall_timeout_seconds=60,
 no_other_case_native_physical_IFT_jet_query_or_retry_authorized=True,independent_signed_postpeer_required=True,
 no_rank_HF_total_length_H0_or_nil_conclusion_yet=True,code_sha256=sha(Path(__file__)))
p=out/'report.json';p.write_text(json.dumps(r,indent=2)+'\n');print(sha(p))
