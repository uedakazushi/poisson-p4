#!/usr/bin/env python3
"""SOURCE terminal degreeface of SAME finiteQQ[t] matrix, no newweight/rank."""
import argparse,hashlib,importlib.util,json,math,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_PREPOWER_K92_ORIGINAL_TERMINAL_FACE_SOURCE_20261105_v63';W=[1,2244,68];NEG=-10**15
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def read(p):return json.loads(p.read_text())
def dump(p,x):p.write_text(json.dumps(x,indent=2,default=str)+'\n')
def main(out):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));out.mkdir(parents=True,exist_ok=False);pins={}
 def src(rel):p=ROOT/rel;pins[rel]=sha(p);return read(p)
 ar=src('reports/S61_NONSQUARE_PREPOWER_K92_EXACT_ABC72_20261105_v63/report.json');dr=src('reports/S61_NONSQUARE_PREPOWER_K92_EXACT_D49_20261105_v63/report.json');assert ar['all_passed']and dr['all_passed'];src('reports/S61_NONSQUARE_PREPOWER_K92_EXACT_ABC_SOURCE_PREFLIGHT_20261105_v63/all72-linear-operator.json');frame=src('reports/S61_NONSQUARE_PREPOWER_K92_ONE_ROW_PRECONDITION_20261105_v63/retropicalized-source.json')
 hp=ROOT/'code/s61_nonsquare_prePOWER_K92_tropical_source_preflight_v2_20261105_v63.py';pins[str(hp.relative_to(ROOT))]=sha(hp);sp=importlib.util.spec_from_file_location('match',hp);h=importlib.util.module_from_spec(sp);sp.loader.exec_module(h)
 binding=src('reports/S61_NONSQUARE_PREPOWER_K92_DET11_GUARDED_ATTEMPT_20261104_v63/source-binding.json');state=ROOT/'reports/S61_NONSQUARE_PREPOWER_K92_DET11_GUARDED_ATTEMPT_20261104_v63/state.bin';pins[str(state.relative_to(ROOT))]=sha(state);assert pins[str(state.relative_to(ROOT))]=='cd5fb2431c2f65649153d82890d52b52ce026f9d0e478bc30592687f17b20f02'
 M,inventory,kept,skipped=h.selected_inputs(state,binding['actual_pure_matrix_roots']);row=src('reports/S61_NONSQUARE_PREPOWER_K92_ONE_ROW_PRECONDITION_20261105_v63/actual-row11.json')
 for p in row['columns']:M[1][p['column']]=dict(state_id=-1001-p['column'],c=Q(p['QQ_common_factor']),co=[(tuple(m),int(v))for m,v in p['primitive_integer_polynomial']])
 source_terms=sum(len(p['co'])for rr in M for p in rr);V=[[NEG]*11 for i in range(11)];coeff=[[Q(0)]*11 for i in range(11)];records=[]
 for i in range(11):
  for j in range(11):
   p=M[i][j];groups={};members={}
   for mon,z in p['co']:
    v=-sum(a*b for a,b in zip(mon,W));groups[v]=groups.get(v,0)+z;members.setdefault(v,[]).append(dict(monomial=mon,integer=str(z)))
   groups={g:z for g,z in groups.items()if z}
   if groups:
    g=max(groups);V[i][j]=g;coeff[i][j]=p['c']*groups[g]
   records.append(dict(row=i,column=j,state_polynomial_id=p['state_id'],source_coefficient_terms=len(p['co']),terminal_degree=None if not groups else g,terminal_primitive_integer=None if not groups else str(groups[g]),source_QQ_content=str(p['c']),terminal_source_monomials=[]if not groups else members[g]))
 degree,count=h.assignment(V);rp,cp,matching=h.dual_potentials(V);assert sum(rp)+sum(cp)==degree
 tight=[[V[i][j]!=NEG and V[i][j]==rp[i]+cp[j]for j in range(11)]for i in range(11)];face=[[coeff[i][j]if tight[i][j]else Q(0)for j in range(11)]for i in range(11)];D=[math.lcm(*(v.denominator for v in r))for r in face];C=[[int(d*x)for x in rr]for d,rr in zip(D,face)]
 wb=src('reports/S61_NONSQUARE_PREPOWER_K92_T18_G_FACTORED_GUARDED_PREFLIGHT_20261104_v63/weighted-output.json');box=wb['output_metadata']['box'];weights=wb['weights'];interval=wb['output_metadata']['linear'];targetweight=-degree;z=targetweight-5;mon=(5+z%68,(z//68)//33,(z//68)%33)
 inbox=all(lo<=x<=hi for x,(lo,hi)in zip(mon,box));inall=inbox and all(lo<=sum(a*b for a,b in zip(mon,w))<=hi for w,(lo,hi)in zip(weights,interval));minimum=None;maximum=None;minmon=maxmon=None
 for a in range(box[0][0],box[0][1]+1):
  for e in range(box[1][0],box[1][1]+1):
   lo,hi=box[2]
   for w,(lower,upper)in zip(weights,interval):
    b=w[0]*a+w[1]*e
    if w[2]:lo=max(lo,(lower-b+w[2]-1)//w[2]);hi=min(hi,(upper-b)//w[2])
    elif not lower<=b<=upper:hi=lo-1;break
   if lo<=hi:
    low=a+2244*e+68*lo;high=a+2244*e+68*hi
    if minimum is None or low<minimum:minimum=low;minmon=(a,e,lo)
    if maximum is None or high>maximum:maximum=high;maxmon=(a,e,hi)
 bits=max(abs(x).bit_length()for r in C for x in r);detbits=11*bits+math.factorial(11).bit_length();structural_nz=sum(bool(x)for r in C for x in r)
 linp=ROOT/'code/s61_nonsquare_prePOWER_K92_leading_structure_20261105_v63.py';pins[str(linp.relative_to(ROOT))]=sha(linp);sp=importlib.util.spec_from_file_location('lin',linp);lin=importlib.util.module_from_spec(sp);sp.loader.exec_module(lin)
 R1=[[Q(x)for x in r]for r in frame['row_transform']];Q1=[[Q(x)for x in r]for r in frame['column_transform']];qfactor=Q(1,math.prod(D))
 dump(out/'terminal-face-and-normalization.json',dict(same_fixed_weight=W,formalvariable='same t, maxdegree endpoint; NOT newweight',full_input_rows=records,degree_matrix=V,terminal_max_degree=degree,max_degree_matching_count=count,row_dual=rp,column_dual=cp,one_matching=matching,
  exact_dual_conditions='r_i+c_j>=degreeMij; equalityonmatching; degree=sumr+sumc; everymaxdegreepermutation usesonlytightentries',terminal_tight_pattern=tight,terminal_QQ_face_matrix=face,row_integer_denominator=D,terminal_integer_face_matrix=C,
  originalM_relation='Mprime(a=t^-1,e=t^-2244,d=t^-68) finiteLaurentpolys; Mprime=row1+2row2 onoriginalM, determinant1. NO R1/Q1 orpotential mixing',
  terminal_integer_det_to_original_coefficient_factor=str(qfactor),source_all121_original_exact_polynomial_atoms=True,NO_terminal_face_determinant_or_rank_computed=True))
 dump(out/'terminal-Newton-filter.json',dict(original_target_weight=targetweight,radix_candidate=mon,inside_original_box=inbox,inside_all_Newton_halfspaces=inall,accepted_box=box,accepted_weights=weights,accepted_intervals=interval,
  tight_original_Newton_minimum=minimum,minimum_monomial=minmon,tight_original_Newton_maximum=maximum,maximum_monomial=maxmon,
  terminal_face_forced_singular_ifoutside=not inall,proof='a nonzero face determinant would yield a nonzero originalM coefficient at thissameweight. InjectiveNewton box permits atmostone and forbidsit ifoutside. Thus face det=0 ifoutside; no newscalar arithmetic.',
  no_detM_ZERO_or_NONZERO_claim=True))
 live=121*12+11*11*5+15000;raw=(11*11*8+200)*((detbits+2)//3+250)+500000
 dump(out/'resource-gate.json',dict(job='ONE constant11face rank/kernel/unitminor profile; detnonzero attempt onlyifinsideoriginalNewton',terminal_face_can_be_NONZERO=bool(inall),structural_nonzero_constants=structural_nz,constant_integer_entry_bit_upper=bits,integer_det_bit_upper=detbits,
  finite_matching_count=count,unique_matching=count==1,owned_slots_upper=live,new_RAW_upper=raw,hard_seconds=60,hard_address_bytes=2000000000,hard_scalar_terms=150000,hard_total_RAW=33554432))
 pins[str(Path(__file__).relative_to(ROOT))]=sha(Path(__file__));note=ROOT/'research/notes'/f'{ID}.md';pins[str(note.relative_to(ROOT))]=sha(note);new=sum(p.stat().st_size for p in out.iterdir());external=sum((ROOT/p).stat().st_size for p in pins);assert new+external<33554432
 report=dict(id=ID,all_passed=True,frozen=True,source_only=True,same_fixed_weight_and_complete_original121_atoms=True,terminal_max_degree=degree,maximum_degree_matching_count=count,terminal_target_original_weight=targetweight,target_radix_monomial=mon,
  inside_original_Newton=inall,terminal_face_forced_singular=not inall,no_new_face_det_or_rank_scalar_job=True,proposed_constant_profile_owned=live,proposed_constant_profile_RAW=raw,no_full_detM_zero_nonzero_or_global57_claim=True,
  total_new_RAW_before_report=new,external_read_RAW=external,input_sha256=pins,evidence_sha256={p.name:sha(p)for p in out.iterdir()},seconds=time.monotonic()-start)
 dump(out/'report.json',report);signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
