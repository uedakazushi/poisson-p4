#!/usr/bin/env python3
"""SOURCE/positive only: two gamma dots on already validated37 scalar values."""
import hashlib,json,math,signal,time
from fractions import Fraction as Q
from pathlib import Path
import ambient_B18_C_genuine_projection7_source_v2_20261105_v63 as S
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_C_KBAR_GAMMA12_VALUE_GATE_20261105_v63'
PATHS={
 'gamma':'reports/AMBIENT_B18_C_KBAR_GAMMA12_SAME_CACHE_SOURCE_20261105_v63/two-unpaired-filtered-gamma-functionals-on-same-thirtyseven-cache-values.json',
 'gamma_report':'reports/AMBIENT_B18_C_KBAR_GAMMA12_SAME_CACHE_SOURCE_20261105_v63/report.json',
 'initial_peer':'reports/AMBIENT_B18_C_KBAR_R2_INITIAL_SOURCE_20261105_v63_REVIEW/report.json',
 'mu_source_peer':'reports/AMBIENT_B18_C_KBAR_MU23_SOURCE_CODE_20261105_v63_REVIEW/report.json',
 'validated_t':'reports/AMBIENT_B18_C_KBAR_MU23_CACHED_VALUES_20261105_v63/phase1-all37-actual-Pg-second-values-and-lower-readbacks.json',
 'mu_value_report':'reports/AMBIENT_B18_C_KBAR_MU23_CACHED_VALUES_20261105_v63/report.json',
 'integer_peer':'reports/AMBIENT_B18_C_KBAR_MU23_ROOT_INTEGER_POSTPEER_20261105_v63/report.json',
}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def env(a):
 a=list(map(Q,a));d=math.lcm(*(x.denominator for x in a));u=max(abs(x.numerator)*(d//x.denominator)for x in a);return[u,d]
def main():
 t=time.monotonic();signal.alarm(60);pins={p:sha(ROOT/p)for p in PATHS.values()};rec=[]
 fields={'gamma':['gamma1_gamma2_weights37','same_actual37_namespace'],'validated_t':['all37_actual_t_head4r2']}
 d={k:S.selective(ROOT/p,fields[k],rec)if k in fields else json.loads((ROOT/p).read_text())for k,p in PATHS.items()}
 assert d['initial_peer']['all_passed']and d['mu_source_peer']['all_passed']and d['integer_peer']['all_passed']
 assert d['integer_peer']['input_sha256'][PATHS['validated_t']]==pins[PATHS['validated_t']]
 assert d['integer_peer']['all37_rr_values_match']and d['mu_value_report']['evidence_sha256'][Path(PATHS['validated_t']).name]==pins[PATHS['validated_t']]
 vals=d['validated_t']['all37_actual_t_head4r2'];ws=d['gamma']['gamma1_gamma2_weights37'];assert len(vals)==37 and len(ws)==2 and all(len(w)==37 for w in ws)
 ev=env(vals);ew=[env(w)for w in ws];bounds=[[37*u*ev[0],v*ev[1]]for u,v in ew]
 H=max(max(u.bit_length(),v.bit_length())for u,v in[ev]+ew+bounds);assert H<=136
 evidence={'scope':'SOURCE positiveprofile; onlyTWOuncomputedgamma dots onVALIDATED37t. No21rootreevaluation, noC/F/Delta/H value.',
 'input_sha256':pins,'gamma_weights2_by37':ws,'same_actual37_namespace':d['gamma']['same_actual37_namespace'],
 'validated_t_foreign_reference':{'path':PATHS['validated_t'],'sha256':pins[PATHS['validated_t']],'field':'all37_actual_t_head4r2','count':37,'ROOT_integer_postpeer':{'path':PATHS['integer_peer'],'sha256':pins[PATHS['integer_peer']]}},
 'mathematical_source':'gamma1=-pi_q2 M2u(rh1-alpha1v),gamma2=-pi_q2 M2u(rh2-rh3/4-alpha2v); exactactual37lifts,oldM5x825/target179. Notnativehomogeneousquadratic.',
 'positive_before_any_gamma_pairing':{'given37_U_D':ev,'weight_rows_U_D':ew,'gamma_partial_sum_U_D_upper':bounds,'calculated_H':H,'chosen_committed':136,'chosen_temporary':273,'hard_committed':384,'hard_temporary':1153},
 'future_numeric':{'static':8000,'hard':12000,'global':150000,'SCI_including_all74termrecords_upper':1000,'complete_failure_saved_pluslive_upper':7000},
 'future_RAW':{'input_source_peer_report_upper':60000,'all_input_and_outputs_prefix_upper':400000,'global_cap':33554432},
 'required_control':'Initializedscientificcursor/pendingparse/mul/add+candidate,completepartial37inputparse/twooutputs/termrecords;SIGALRMexceptionprefix. Phaseinputreceiptbeforetwodots. Noretry/replay/capraise.',
 'future_seconds':60,'future_memory_bytes':2000000000,'no_signed_gamma_values_computed':True,'selective_receipts':rec}
 owned=S.count(d)+S.count(evidence)+12000;assert owned<150000;evidence['whole_SOURCE_owned_upper']=owned
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);p=out/'two-gamma-positive-source-and-validated37-input-only-gate.json';p.write_text(json.dumps(evidence,indent=2)+'\n')
 r={'id':STEM,'all_passed':True,'status':'TWO_GAMMA_SOURCE_POSITIVE_ONLY_UNPAIRED','input_sha256':pins,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/(STEM+'.md')),'evidence_sha256':{p.name:sha(p)},'calculated_H':H,'static_value_upper':8000,'whole_SOURCE_owned_upper':owned,'seconds':time.monotonic()-t};(out/'report.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r));signal.alarm(0)
if __name__=='__main__':main()
