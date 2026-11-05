#!/usr/bin/env python3
"""Explicit accepted-custody reuse; no hiddenfullmatrix read in nextjob."""
from pathlib import Path
import hashlib,json,time
ROOT=Path(__file__).resolve().parents[1];ID='S61_SPARSE_GENERIC5_CONSTRAINED83_FIFTH_SOURCE_CUSTODY_PLAN_20261105_v64'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def enc(v):return(json.dumps(v,separators=(',',':'))+'\n').encode()
begin=time.monotonic();out=ROOT/'reports'/ID;out.mkdir(exist_ok=False);inputs={Path(__file__):sha(Path(__file__)),ROOT/f'research/notes/{ID}.md':sha(ROOT/f'research/notes/{ID}.md')}
sp=ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_FIFTH_UNIT_COL81_SOURCE_20261105_v64/report.json';assert sha(sp)=='d668d98214a0f755afe96990eb143f355c90ecaf8df7a9c9a988eeb4509b11fe';s=read(sp);inputs[sp]=sha(sp);np=sp.parent/'next-col81-fixed-source-operands.json';assert sha(np)==s['evidence_sha256'][np.name];inputs[np]=sha(np);new=read(np)['source_operands']
op=ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_FOURTH_UNIT_COL76_SOURCE_20261105_v64/col76-fixed-source-operands.json';assert sha(op)==s['input_sha256'][str(op.relative_to(ROOT))];inputs[op]=sha(op);old=read(op)['source_operands']
gp=ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_FOURTH_COL76_IMPLEMENTATION_GATE_20261105_v64/report.json';assert sha(gp)=='1c5f79ea5b373e2c8b7f95aff5ccefe9a058c7d4d7bd20df6bc7429bd6c14b23';g=read(gp);assert g['all_passed']and g['all63_leaf_custody'];inputs[gp]=sha(gp);ep=gp.parent/'implementation-custody-storage-bit-addendum.json';assert sha(ep)==g['evidence_sha256'][ep.name];inputs[ep]=sha(ep);accepted={z['label']:z for z in read(ep)['all63_leaf_custody_checks']};reused=[];cached=[];quot=[];aux=[]
for label,z in new.items():
 if label=='Q'or label.startswith('U')or label.startswith('V'):
  assert label in accepted and accepted[label]['actual_original_source_equal']and z==old[label];reused.append(label)
 elif label=='T'or label=='R'or label.startswith('W'):
  assert z['source']['kind']=='actual_stage2_col81';cached.append({'label':label,'source':z['source'],'future_check':'Matchactualsource row label/path/SHA to accepted22cacheAPI, then decodeitsliteralbinary and comparevalue.'})
 elif label.startswith('C'):
  assert z['source']['kind']=='known_Q_quotient_existing_col76';quot.append({'label':label,'source':z['source'],'future_check':'Readoriginalcurrentcol76numeratorbinary; verifyQ*newliteral equalsoriginalFi EXACT. This normalizesexistingvalue, notnewfuturecol.'})
 else:
  assert label=='G'and z['source']['kind']=='fixed_source_row_sum';aux.append({'label':label,'source':z['source'],'future_check':'Readthetwoacceptedstage2cachebinaries14/35 and verifynewliteral=M2_1481+6M2_3581.'})
assert(len(new),len(reused),len(cached),len(quot),len(aux))==(80,39,21,19,1)
v={'all_passed':True,'source_only':True,'all80_new_operand_custody_classes':{'exact_accepted_old_Q_U_V_reuse39':reused,'actual_stage2_cache_bindings21':cached,'knownQ_quotient_bindings19':quot,'fixed_row_sum_binding1':aux},'accepted_physical_old_custody_gate':str(gp.relative_to(ROOT)),'accepted_gate_sha256':sha(gp),'proof':'EveryreusedQ/U/V object is exactlyequaltoitsold63bundleobject, whosephysicalsourcebindingswereacceptedunderthissamehash. ReuseONLYthese39labels. NewR/T/W labels are NOToldR/W andmustcheckthe22stage2actualcache.','no_unlisted_oldfullmatrix_reads_in_next_RAW_budget':'Nexthelpermayreuseacceptedphysicalgate ratherthan rehash/rereadcurrent-matrix.bin. Ifitchoosesfresholdmatrixrevalidation, anewactualREAD/RAW gateisrequired; doNOTclaimthe2.9MBboundthen.','all_other_leaf_checks_literal_smallbinaries':True,'no_new_family_coefficients':True,'future_valueworker_still_NOT_authorized':True}
(out/'accepted-custody-reuse-and-all80-source-interfaces.json').write_bytes(enc(v));report={'id':ID,'all_passed':True,'source_only':True,'all80_operand_sources_classified':True,'accepted_exact_old_reuse':39,'literal_small_cache_Qquotient_row_sum':41,'no_old_matrix_bytes_read':True,'no_future_values_computed':True,'input_sha256':{str(p.relative_to(ROOT)):h for p,h in inputs.items()},'evidence_sha256':{p.name:sha(p)for p in out.iterdir()},'seconds':time.monotonic()-begin};(out/'report.json').write_bytes(enc(report));print(json.dumps({k:v for k,v in report.items()if k not in('input_sha256','evidence_sha256')}));print('report_sha256',sha(out/'report.json'))
