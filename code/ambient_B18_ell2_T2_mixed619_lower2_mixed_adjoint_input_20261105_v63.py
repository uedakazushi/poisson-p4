#!/usr/bin/env python3
"""Source-only136 mixed outputs on210 first slots; no rank/adjoint solve."""
import json,math,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
from ambient_B18_ell2_T2_nonactive171_source_review_20261105_v63 import sha,sp,count,add
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_ELL2_T2_MIXED619_LOWER2_MIXED_ADJOINT_INPUT_20261105_v63'
def main():
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024));inputs={}
 def read(p):inputs[p]={'path':p,'sha256':sha(ROOT/p)};return json.loads((ROOT/p).read_text())
 kp='reports/AMBIENT_B18_ELL2_T2_MIXED619_LOWER4_KNOWN_INPUT_GATE_20261105_v63/pruned-known619-lower4-input-and-prior-resource-gate.json';d=read(kp)
 pre='reports/AMBIENT_B18_ELL2_T2_MIXED619_LOWER4_COVECTORS_20261105_v63';r=read(pre+'/report.json');assert r['all_passed']
 cp=pre+'/actual-four-full918-source-covectors-and52-plane-adjoint-input.json';c=read(cp);assert sha(ROOT/cp)==r['evidence_sha256'][Path(cp).name]
 G2=list(map(sp,d['G2linear_source102_columns']));alpha=list(map(sp,c['full918_covectors4']));assert not alpha[0] and not alpha[1]
 mono=[(7,7),(7,8),(0,7),(0,8),(4,7),(4,8),(5,7),(5,8)]
 idx={tuple(sorted(m)):i for i,m in enumerate(mono)};cols=[{} for _ in range(17*8)]
 for j,col in enumerate(G2):
  for k,v in col.items():
   target,var=divmod(k,9)
   if target>=17:continue
   for first in range(9):
    m=tuple(sorted((var,first)))
    if m in idx:
     out=8*target+idx[m];cols[out]=add(cols[out],{9*j+first:v})
 domain=[9*j+var for j in range(42) for var in (0,4,5,7,8)];slot={v:i for i,v in enumerate(domain)}
 assert len(domain)==210
 assert all(set(v)<=set(domain) for v in cols+alpha[2:])
 nonzero=[i for i,v in enumerate(cols) if v];A=[{slot[k]:v for k,v in cols[i].items()} for i in nonzero]
 target=[{slot[k]:v for k,v in a.items()} for a in alpha[2:]]
 # Bound formation before any solve: coefficients copied from exact G2linear,
 # each monomial output coefficient has at most2 contributing entries.
 D=math.lcm(*(x.denominator for v in G2 for x in v.values()));U=max([0]+[abs(x.numerator)*(D//x.denominator) for v in G2 for x in v.values()])
 assert all(D%x.denominator==0 and abs(x.numerator)*(D//x.denominator)<=2*U for v in A for x in v.values())
 inputnum=count([A,target,domain,nonzero,mono]);rankmax=min(len(A),210)
 # Solver frontier stores one original/normalized/repr basis, input sparse
 # columns, two target reductions, fixed receipts and serialized overlaps.
 owned=2*inputnum+2*rankmax*210+rankmax*rankmax+4*rankmax+2*210+10000
 rawbound=owned*(2*math.ceil(384*math.log10(2))+40)+1_000_000
 ready=owned<150000 and rawbound<32*1024*1024
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);ep=out/'actual136-mixed-output210-slot-adjoint-source-input.json'
 ev={'monomial_output_order8':mono,'active_target_order17':list(range(17)),'all136_output_indices':list(range(136)),'nonzero_output_indices':nonzero,
 'full918_domain_slots210':domain,'matrix_columns_on210_domain':[ [[k,str(x)] for k,x in sorted(v.items())] for v in A],
 'target_covectors2_on210_domain':[[[k,str(x)] for k,x in sorted(v.items())] for v in target],
 'all_other708_domain_slots_identically_zero':True,'all136_source_outputs_present_before_pruning':True,'no_unknown_value_or_Gaussian':True,
 'matrix_integer_clear_D':str(D),'matrix_abs_clear_entry_bound':str(2*U),'coefficient_formation_prior_bits':max(D.bit_length(),(2*U).bit_length()),
 'whole_solver_numeric_bound':owned,'whole_solver_RAW_bound_at384':rawbound,'prior_ownership_RAW_fit':ready,'rank_height_not_claimed_a_priori':True,
 'future_guard_protocol':'Only with parent authorization: fixedQQ elimination dynamic384 commit guard, uncommitted1153 maximum, complete frontier and pivot history saved before stop;60s2GB32MiB150k. No implicit/physical/gauge evaluation.',
 'mathematical_scope':'Seekrho withrho(G2linear*Y1)=alpha onALL918slots, usingall136listedmixedquadraticoutputs. Existencewouldreducefourlower619channels(twoalreadyzero)toX2rr/rellplusz-planeoutputcoefficients; failureisrestricted-familyobstruction,notneedfornewphysicalsource.'}
 ep.write_text(json.dumps(ev,separators=(',',':'))+'\n')
 report={'id':STEM,'all_passed':True,'no_rank_or_adjoint_solve':True,'nonzero_output_count':len(A),'domain_count':210,'whole_solver_numeric_bound':owned,'whole_solver_RAW_bound_at384':rawbound,'prior_ownership_RAW_fit':ready,'input_sha256':inputs,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-t};(out/'report.json').write_text(json.dumps(report,separators=(',',':'))+'\n');print(json.dumps({'all_passed':True,'nonzero_outputs':len(A),'domain':210,'owned':owned,'RAW':rawbound,'ready':ready,'report_sha256':sha(out/'report.json')}))
if __name__=='__main__':main()
