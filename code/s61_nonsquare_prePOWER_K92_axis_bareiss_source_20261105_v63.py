#!/usr/bin/env python3
"""SOURCE integerpolynomial Bareiss bounds for thewhole finitea-axis matrix."""
import argparse,hashlib,json,math,resource,signal,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_PREPOWER_K92_WHOLE_AXIS_BAREISS_SOURCE_20261105_v63'
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def read(p):return json.loads(p.read_text())
def dump(p,x):p.write_text(json.dumps(x,indent=2,default=str)+'\n')
def main(out):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));out.mkdir(parents=True,exist_ok=False);pins={}
 def src(rel):p=ROOT/rel;pins[rel]=sha(p);return read(p)
 base='reports/S61_NONSQUARE_PREPOWER_K92_WHOLE_AXIS_SOURCE_V2_20261105_v63';r=src(base+'/report.json');assert r['all_passed'];a=src(base+'/whole-axis-literal-source.json');kr=src('reports/S61_NONSQUARE_PREPOWER_K92_WHOLE_AXIS_KERNEL_ENTRY_20261105_v63/report.json');assert kr['all_passed']and not kr['right55_ALL_ZERO']and not kr['left55_ALL_ZERO']
 A={(e['row'],e['column']):{int(k):int(v)for k,v in e['axis_integer_polynomial']}for e in a['axis_coefficients']};assert len(A)==121 and sum(map(len,A.values()))==579;rowM=[max(sum(abs(x)for x in A[i,j].values())for j in range(11))for i in range(11)];assert all(rowM);H=math.factorial(11)*math.prod(max(1,x)for x in rowM);bits=H.bit_length();maxD=max(max(p,default=0)for p in A.values());D=11*maxD;rawbits=2*bits+10
 cells=121*(D+1)+5*(2*D+1)+2*579+15000;patchcells=sum((10-k)**2*((k+1)*maxD+1)for k in range(10));checkpointcells=121*(D+1)+5*(2*D+1)
 mainRAW=(patchcells+579+D+1)*((bits+2)//3+100)+1000000;stateRAW=checkpointcells*((rawbits+2)//3+100)+500000;totalRAW=mainRAW+stateRAW
 pairbound=3*385*(D+1)**2
 dump(out/'integer-polynomial-source-proof.json',dict(input_literal_source=base+'/whole-axis-literal-source.json',input_polynomial_entries=121,input_integer_coefficients=579,max_entry_degree=maxD,
  row_entry_L1max=list(map(str,rowM)),every_minor_L1_upper=str(H),minor_coefficient_bits=bits,minor_degree_upper=D,
  exact_proof='eachk-minor L1<=k!*productselectedrowentryL1max<=11!*productallmax(1,rowentryL1max)=H; degree<=k*8<=88',
  raw_product_degree_upper=2*D,raw_product_coefficient_bound='(D+1)*H²',raw_Bareiss_difference_bound='2*(D+1)*H²',manual_integer_longdivision_residual_bound='3*(D+1)*H²',raw_integer_bit_upper=rawbits,
  division='ownleadingcoefficient exactZZ division and subtract knownquotientterm; quotientis a minor soeveryquotientcoef<=H; remainderZERO required; NOQQ(a)RFcanonicalization or polynomialgcd',
  rank='lexnonzeroremainingpolynomial pivot withrow/colswaps; allremainingZERO certifiesupperrank; lastnonzeropivot isselectedinputminor certifyinglower',
  determinant='ifrank11 signedlastpivot; otherwiseaxisdetZERO',original_detM_axis='a^5*detCaxis/prod(axisrowLCD)',
  axisd0_OUTSIDEFULLU=True,scope='genericrankoverQQ(a) only; notallapointwise, notfullFULLUcubic57'))
 dump(out/'whole-worker-source-plan.json',dict(algorithm='completepivot fractionfreeINTEGERPOLYNOMIAL Bareiss, n=11degree8; no fractionfieldGaussian retry',input_source=a,
  fixed_dense_minor_capacity_perentry=D+1,fixed_dense_raw_capacity_perbuffer=2*D+1,actualmatrix_coefficient_capacity=121*(D+1),five_buffers_coefficient_capacity=5*(2*D+1),two_literal_inputcopies=2*579,
  all_owned_with_metadata15000=cells,coefficient_bit_bound=bits,all_raw_temporary_bit_bound=rawbits,naive_integer_mult_upper=pairbound,
  proof_storage='initialmatrixsource immutable; eachstage recordsrow/colswaps and ONLYnewremainingblock (max10065coeffs total). Independentchecker reconstructsoldmatrix and verifieseveryBareissdivision byremultiplication; no repeatedwholematrix snapshots.',
  all_patch_coefficient_capacity=patchcells,main_RAW_upper=mainRAW,full_guard_checkpoint_RAW_upper=stateRAW,total_possible_new_RAW_upper=totalRAW,
  guard='45s insidechecked multiply/exactdivide, ownedslots beforecommit; preservefullcurrentmatrix+5buffers+operationcursor pluscompletedpatches; NOautoresume',
  proposed_one_job='exactwholeaxisrank/determinant withALLpolynomialdivisioncertificates; boundFIRST, no newa-value/point/prime orcoefwindow',
  hard_seconds=60,hard_address_bytes=2000000000,hard_scalar_terms=150000,hard_total_RAW=33554432))
 pins[str(Path(__file__).relative_to(ROOT))]=sha(Path(__file__));note=ROOT/'research/notes'/f'{ID}.md';pins[str(note.relative_to(ROOT))]=sha(note);new=sum(p.stat().st_size for p in out.iterdir());external=sum((ROOT/p).stat().st_size for p in pins);assert new+external<33554432
 report=dict(id=ID,all_passed=True,frozen=True,source_only=True,no_new_polynomial_division_or_rank_values=True,finite_wholeaxis_matrix=True,polynomial_minor_degree=D,minor_coefficient_bits=bits,raw_temporary_bits=rawbits,
  whole_owned_slots_upper=cells,main_RAW_upper=mainRAW,checkpoint_RAW_upper=stateRAW,total_new_RAW_upper=totalRAW,naive_integer_product_upper=pairbound,
  preflight_fits=cells<150000 and totalRAW+external<33554432,no_generic_QQ3variable_fieldGaussian_or_inverse_search=True,no_fullU_rank_claim=True,
  total_source_new_RAW_before_report=new,external_read_RAW=external,input_sha256=pins,evidence_sha256={p.name:sha(p)for p in out.iterdir()},seconds=time.monotonic()-start)
 dump(out/'report.json',report);signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
