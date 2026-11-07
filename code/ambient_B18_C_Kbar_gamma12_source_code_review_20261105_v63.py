#!/usr/bin/env python3
"""Independent SOURCE/code review of exactly two unpaired gamma functionals."""
import ast,hashlib,importlib.util,json,math,pathlib,resource,signal,time
from fractions import Fraction as Q
from ambient_B18_C3_alpha4_cached_source_review_20261105_v63 import field
from ambient_B18_C_genuine_projection7_source_review_20261105_v63 import count,sha
ROOT=pathlib.Path(__file__).resolve().parents[1]
GATE='AMBIENT_B18_C_KBAR_GAMMA12_VALUE_GATE_20261105_v63';ID=GATE+'_REVIEW'
WORKER='code/ambient_B18_C_Kbar_gamma12_cached_values_20261105_v63.py'
PIN='79f227eb16b031b352644cbca915b7f4dcf33c908de71627c2e4f5e6782ea4fe'
PRE=0
def qa(a,b):
 global PRE
 if a and b:PRE=max(PRE,2*max(abs(a.numerator).bit_length(),a.denominator.bit_length(),abs(b.numerator).bit_length(),b.denominator.bit_length())+1)
 assert PRE<=1153;v=a+b;assert max(abs(v.numerator).bit_length(),v.denominator.bit_length())<=384;return v
def qm(a,b):
 global PRE
 if a and b:PRE=max(PRE,max(abs(a.numerator).bit_length()+abs(b.numerator).bit_length(),a.denominator.bit_length()+b.denominator.bit_length()))
 assert PRE<=1153;v=a*b;assert max(abs(v.numerator).bit_length(),v.denominator.bit_length())<=384;return v
def total(xs):
 z=Q(0)
 for x in xs:z=qa(z,x)
 return z
def put(d,k,c):
 if c:
  v=qa(d.get(k,Q(0)),c)
  if v:d[k]=v
  elif k in d:del d[k]
def envelope(xs):
 xs=list(xs);D=math.lcm(*(q.denominator for q in xs));return [max(abs(q.numerator)*(D//q.denominator) for q in xs),D]
def main():
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));out=ROOT/'reports'/ID;out.mkdir(exist_ok=False)
 pins={};raw=0
 def pin(p,h=None,read=False):
  nonlocal raw
  p=pathlib.Path(p);pp=p if p.is_absolute() else ROOT/p;k=str(pp.relative_to(ROOT));z=sha(pp)
  if h:assert h==z,(k,z,h)
  if k not in pins:raw+=pp.stat().st_size
  pins[k]=z;return json.loads(pp.read_text()) if read else pp
 gp=f'reports/{GATE}/report.json';g=pin(gp,'8b243861ed0252f3425329468661e58c92ea9a9549e7620fa603dc973b0578bf',True);assert g['all_passed']
 for p,h in g['input_sha256'].items():pin(p,h)
 ep=f'reports/{GATE}/two-gamma-positive-source-and-validated37-input-only-gate.json';D=pin(ep,'18a9e953eaaf8f7966b2ae4455a60644214996e15fef3ba2a9f538318f83abfc',True)
 pin(WORKER,PIN)
 cr=pin('reports/AMBIENT_B18_C_KBAR_GAMMA12_VALUE_CUSTODY_20261105_v63/report.json',read=True);assert cr['all_passed'] and cr['code_sha256']==PIN
 ce='reports/AMBIENT_B18_C_KBAR_GAMMA12_VALUE_CUSTODY_20261105_v63/prepared-two-gamma-only-worker-full-prefix-and-whole-custody.json'
 C=pin(ce,'320070ded953c0dde40f70f0233ef526d30509d08e1c697742b92c3dee5f7ded',True)
 for p,h in C['input_sha256'].items():pin(p,h)
 gamma=ROOT/'reports/AMBIENT_B18_C_KBAR_GAMMA12_SAME_CACHE_SOURCE_20261105_v63/two-unpaired-filtered-gamma-functionals-on-same-thirtyseven-cache-values.json'
 G=json.loads(gamma.read_text())
 for p,h in G['input_sha256'].items():pin(p,h)
 owninitial=pin('reports/AMBIENT_B18_C_KBAR_R2_INITIAL_SOURCE_20261105_v63_REVIEW/report.json','662785338116da297cfe47ea8b8e117fee0b4ee2834d808458600f59f95384c7',True)
 assert owninitial['all_passed']
 initial=ROOT/'reports/AMBIENT_B18_C_KBAR_R2_INITIAL_SOURCE_20261105_v63_REVIEW/independent-linear-system-alphas-cofactor-multiplybacks-and-projected-r2.json'
 pin(initial,owninitial['evidence_sha256'][initial.name])
 cofs=[[Q(c) for c in row] for row in field(initial,'two_independent_source_cofactor_vectors')]
 mu=ROOT/'reports/AMBIENT_B18_C_KBAR_MU23_SOURCE_INVENTORY_20261105_v63/fourteen-z-actions-common37-M2-source-functionals-and-minimal-generator-criterion.json'
 T=[[Q(c) for c in row] for row in field(mu,'normalized37_actual37_constant_lift')]
 namespace=field(mu,'actual_selected37_sources')
 assert G['same_actual37_namespace']==D['same_actual37_namespace']==namespace
 orig=ROOT/'reports/AMBIENT_B18_C_KBAR_R2_INITIAL_SOURCE_20261105_v63/two-cached-alpha-reductions-and-projected-r2-seven-source-certificates.json'
 records=field(orig,'two_exact_linear_reductions_and_membership_cofactors')
 assert len(records)==2 and len(T)==37 and all(len(row)==37 for row in T)
 # Our prior independent bottom-pivot solve supplied these unique source
 # coordinates; compare all normalized source coefficients, not gamma values.
 for row,b in zip(records,cofs):
  fromtrace=[Q(0)]*37
  for j,shift,c in row['exact_linear_source_cofactors']:
   assert shift==[0]*8;fromtrace[j]=qa(fromtrace[j],Q(c))
  assert fromtrace==b
 weights=[]
 for b in cofs:weights.append([-total(qm(b[k],T[k][j]) for k in reversed(range(37))) for j in range(37)])
 assert [[str(c) for c in row] for row in weights]==G['gamma1_gamma2_weights37']==D['gamma_weights2_by37']
 # Optional lineage check involves SOURCE weights only, never root values.
 mcache=ROOT/'reports/AMBIENT_B18_C_KBAR_MU23_CACHE_SOURCE_20261105_v63/complete-old-Pg-second-cache-source-and-fourteen-action-positive-gate.json'
 r37=field(mcache,'all37_Pg_source_weights');rw=[]
 for w in weights:
  v={}
  for a,row in zip(w,r37):
   for root,c in row:put(v,root,qm(a,Q(c)))
  rw.append([[j,str(c)] for j,c in sorted(v.items())])
 assert rw==G['gamma1_gamma2_direct_given21_root_weights']
 # Existing validated values are inspected ONLY for a common absolute
 # envelope. No weight is paired with a requested scalar value.
 ref=D['validated_t_foreign_reference'];assert ref['count']==37 and ref['field']=='all37_actual_t_head4r2'
 tp=pin(ref['path'],ref['sha256']);integer=pin(ref['ROOT_integer_postpeer']['path'],ref['ROOT_integer_postpeer']['sha256'],True)
 assert integer['all_passed'] and integer['all37_rr_values_match'] and integer['input_sha256'][ref['path']]==ref['sha256']
 vals=list(map(Q,field(tp,ref['field'])));assert len(vals)==37
 E=envelope(vals);W=[envelope(row) for row in weights];bounds=[[37*a*E[0],d*E[1]] for a,d in W]
 H=max(max(a.bit_length(),d.bit_length()) for a,d in [E]+W+bounds);assert H==72
 P=D['positive_before_any_gamma_pairing'];assert P['given37_U_D']==E and P['weight_rows_U_D']==W and P['gamma_partial_sum_U_D_upper']==bounds
 assert P['chosen_committed']==136 and P['chosen_temporary']==273 and H<136
 # Full worker/code/control review. No call to main or a requested dot.
 path=ROOT/WORKER;text=path.read_text();tree=ast.parse(text)
 assert {n.module for n in ast.walk(tree) if isinstance(n,ast.ImportFrom)}=={'pathlib','fractions'}
 for token in ["original=json.loads(tp.read_text())","given=original[src['field']]","del original","SCI['pending']['computed_candidate']=str(q)","terms.append([j,str(a),str(v),str(value)]);z=add(z,value)","except BaseException as exc:preserve(exc);raise","integer['all37_rr_values_match']","assert peer['input_sha256'][GATE]==GATE_SHA","assert count(SCI)<=1000"]:
  assert token in text,token
 assert text.index("save('phase0-validated-existing37-input-reference.json'")<text.index("SCI['phase']='only_two_gamma_dots'")
 assert 'rootweights' not in text and 'C_F_Delta_H_not_computed' in text
 spec=importlib.util.spec_from_file_location('gamma_guard_fixture_ONLY',path);worker=importlib.util.module_from_spec(spec);spec.loader.exec_module(worker)
 assert worker.GATE_SHA==pins[ep]
 assert worker.add(Q(1,7),Q(2,11))==Q(25,77)
 worker.mark('unrelated_candidate_fixture')
 try:worker.checked(Q(1<<136));raise RuntimeError('missing guard')
 except AssertionError:assert worker.SCI['pending']['computed_candidate']==str(1<<136) and not worker.SCI['pending']['completed']
 # Count a complete maximum scientific prefix with unrelated136bit text.
 huge=str((1<<136)-1)+'/'+str((1<<135)-1)
 SCI={'phase':'only_two_gamma_dots','row':1,'term':36,'pending':{'operation':'add','completed':False,'operands':[huge,huge],'computed_candidate':huge},'parsed_t37':[huge]*37,'values':[huge]*2,'terms':[[[j,str(weights[k][j]),huge,huge] for j in range(37)] for k in range(2)],'accumulator':huge}
 sci=worker.count(SCI);compiled=worker.count(D);peer=1000;iv=worker.count(integer)
 # Charge complete foreign artifact while extracting37, current SCI,
 # saved output and its serialized prefix again, plus full controller.
 phase_artifact_numeric=541;body={'gamma1_gamma2':[huge]*2,'all74_termwise_products':SCI['terms'],'same37_inputs_reused_without_root_or_tdot_replay':True,'C_F_Delta_H_not_computed':True}
 savemax=compiled+peer+iv+37+sci+worker.count(body)*2+2000
 failure=compiled+peer+iv+37+2*sci+worker.count(body)+2200
 assert sci<=1000 and max(savemax,failure)<8000<12000
 rawprefix=len(json.dumps({'scientific_state':SCI,'traceback':'x'*10000},indent=2).encode())
 actualREAD=sum((ROOT/p).stat().st_size for p in [ep,ref['path'],ref['ROOT_integer_postpeer']['path']])+60000
 assert actualREAD+rawprefix+100000<400000
 assert not (ROOT/'reports/AMBIENT_B18_C_KBAR_GAMMA12_CACHED_VALUES_20261105_v63').exists()
 owned=count((g,D,C,G,cr,owninitial,cofs,T,namespace,records,weights,r37,rw,pins))+18000
 assert owned<150000 and raw+400000<32*1024**2
 evidence={'all_passed':True,'SOURCE_CODE_ONLY_no_gamma_scalar_paired':True,'independent_previous_fixed_linear_source_cofactors_used':True,
  'all74_true_actual37_gamma_weights_rebuilt':True,'direct21root_SOURCE_weights_only_match':True,'same37_true_OLD825_168plus657_namespace':True,
  'true_filtered_gamma_proof':'rh1-alpha1v=M1u1 and rh2-rh3/4-alpha2v=M1u2. Using true whole relations M, not just M1, their q2 coefficients modulo oldm3 are -pi_q2*M2u. Complete pi_q2 is literal179 and oldm3 tails are killed for gamma. Negative sign and normalized37->actual37 lift are rechecked; no C/F/Delta/H or action/type/H0 value follows.',
  'ROOT_validated_foreign37_reference':ref,'NO37dot_or21root_recomposition':True,'positive_envelopes_only':{'provided_U_D':E,'weight_U_D':W,'partial_bounds':bounds,'calculated_H':H,'committed136_pre273':True},
  'full_SCI_shape_numeric':sci,'complete_save_conservative':savemax,'complete_failure_conservative':failure,'foreign541_parse_in8000':compiled+peer+iv+phase_artifact_numeric+2000<8000,
  'worst_prefix_RAW_with10000traceback':rawprefix,'READ_plus_prefix_and100000_outputs_upper':actualREAD+rawprefix+100000,'hard_RAW400000':True,
  'pending_operands_candidate_full74terms_and_accumulator_preserved':True,'phase0_reference_saved_before_two_dots':True,'SIGALRM_exception_prefix_no_retry':True,'external_AS2GB_CPU60_required':True,
  'whole_SOURCE_numeric_upper':owned,'SOURCE_actual_pre_bits':PRE,'RAW_SOURCE_withreserve':raw+400000,'synthetic_candidate_guard_only':True,'actual_gamma_values_UNCOMPUTED':True}
 p=out/'independent-two-gamma-true37-source-code-and-whole-prefix-proof.json';save=lambda p,v:p.write_text(json.dumps(v,indent=2)+'\n');save(p,evidence)
 report={'id':ID,'all_passed':True,'SOURCE_code_custody_ONLY':True,'input_sha256':pins,'code_sha256':sha(pathlib.Path(__file__)),'evidence_sha256':{p.name:sha(p)},'worker_sha256':PIN,'calculated_prior_H72_coarse136_pre273':True,'static8000_hard12000':True,'no_gamma_pairing':True,'seconds':time.monotonic()-start}
 save(out/'report.json',report);signal.alarm(0);print(json.dumps({'all_passed':True,'report_sha256':sha(out/'report.json'),'source_owned':owned,'seconds':report['seconds']}))
if __name__=='__main__':main()
