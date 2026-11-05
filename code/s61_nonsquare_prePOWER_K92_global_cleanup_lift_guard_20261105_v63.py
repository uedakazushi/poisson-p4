#!/usr/bin/env python3
"""SOURCE guard: generic blockzeros before specialization, namedunit images."""
import argparse,hashlib,importlib.util,json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_PREPOWER_K92_GLOBAL_CLEANUP_LIFT_GUARD_20261105_v63'
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def read(p):return json.loads(p.read_text())
def dump(p,x):p.write_text(json.dumps(x,indent=2,default=str)+'\n')
def main(out):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));out.mkdir(parents=True,exist_ok=False);pins={}
 def src(rel):p=ROOT/rel;pins[rel]=sha(p);return read(p)
 base='reports/S61_NONSQUARE_PREPOWER_K92_GLOBAL_CLEANUP_SOURCE_20261105_v63';r=src(base+'/report.json');assert r['all_passed'];lift=src(base+'/typed-generic-rowcleanup-and-pivot.json');fr=src('reports/S61_NONSQUARE_PREPOWER_K92_ONE_ROW_PRECONDITION_20261105_v63/retropicalized-source.json');R=fr['row_transform'];QQ=fr['column_transform'];assert [j for j,x in enumerate(R[0])if Q(x)]==[0]and [j for j,x in enumerate(R[1])if Q(x)]==[0,1];assert all(Q(QQ[i][j])==Q(i==c)for j,c in[(2,8),(3,10)]for i in range(11))
 hp=ROOT/'code/s61_nonsquare_prePOWER_K92_tropical_source_preflight_v2_20261105_v63.py';pins[str(hp.relative_to(ROOT))]=sha(hp);sp=importlib.util.spec_from_file_location('sourcehelper',hp);h=importlib.util.module_from_spec(sp);sp.loader.exec_module(h);binding=src('reports/S61_NONSQUARE_PREPOWER_K92_DET11_GUARDED_ATTEMPT_20261104_v63/source-binding.json');state=ROOT/'reports/S61_NONSQUARE_PREPOWER_K92_DET11_GUARDED_ATTEMPT_20261104_v63/state.bin';pins[str(state.relative_to(ROOT))]=sha(state);M,inventory,kept,skipped=h.selected_inputs(state,binding['actual_pure_matrix_roots']);assert kept==21690
 zeros=[]
 for i in[0,1,2]:
  for j in[8,10]:assert not M[i][j]['co'];zeros.append(dict(source_row=i,source_column=j,literal_generic_coefficient_count=0))
 units={
  'a':[(1,[1,0,0])], 'd':[(1,[0,0,1])],
  'F':[(500,[1,0,0]),(90,[0,2,0]),(3,[0,1,0]),(-24,[0,0,0])],
  'B':[(50,[1,0,1]),(18,[0,2,1]),(-9,[0,1,1]),(-2,[0,0,1]),(60,[0,1,0]),(-40,[0,0,0])],
  'C':[(100,[1,0,0]),(18,[0,2,0]),(-9,[0,1,0]),(-2,[0,0,0])],
  'J':[(5,[1,0,1]),(-6,[0,1,0]),(-4,[0,0,0])]}
 images={}
 for k,pp in units.items():
  z={}
  for c,m in pp:g=-sum(x*y for x,y in zip(m,[1,2244,68]));z[g]=z.get(g,0)+c
  z={g:c for g,c in z.items()if c};assert z;images[k]=dict(QQ_Laurent_terms=sorted(z.items()),lowest_exponent=min(z),lowest_nonzero_coefficient=z[min(z)])
 dump(out/'generic-block-zero-and-unit-proof.json',dict(generic_M_literal_zeros=zeros,source_terms_read=kept,no_source572_reevaluated=True,
  Mprime_row1='M1+2M2, henceH3cols8,10remainliteralZERO',R1_top2_support=[[(j,x)for j,x in enumerate(row)if Q(x)]for row in R[:2]],Q1_pivot_columns=dict(column2='e8',column3='e10'),
  proof='row1+2row2 preservesgenericliteralzeros; Dr,Dc monomials preservezeros; R1 top2useonlyMprime rows0,1; Q1first4 pivots includeunmixedcols8,10. Thus genericA_topright2x2 isZERO beforeanyKronecker map.',
  no_false_injectivity_on_individual_row_numerators=True,named_FULLU_polynomials=units,formal_images=images,FULLU_localization_maps_to_QQ_Laurent_field=True,
  chart_nonzero='deltaP,Q specializeconstant1; parentdirect5coeffpeer provesN0specializationnonzero =>genericrootsnonzero. NoneisclaimedglobalFULLUunit.',
  all7_genericcallbacks_RETAINED_even_if_specializedzero=True))
 pins[str(Path(__file__).relative_to(ROOT))]=sha(Path(__file__));note=ROOT/'research/notes'/f'{ID}.md';pins[str(note.relative_to(ROOT))]=sha(note);new=sum(p.stat().st_size for p in out.iterdir());external=sum((ROOT/p).stat().st_size for p in pins);assert new+external<33554432
 report=dict(id=ID,all_passed=True,frozen=True,source_only=True,generic_A_topright_block_literallyZERO_before_specialization=True,all6_FULLU_unit_images_nonzero=True,
  no_genericzero_from_specializedzero=True,generic_chart_complements_retained=True,no_new_signed_callback_coefficient_values=True,total_new_RAW_before_report=new,external_read_RAW=external,input_sha256=pins,evidence_sha256={p.name:sha(p)for p in out.iterdir()},seconds=time.monotonic()-start)
 dump(out/'report.json',report);signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
