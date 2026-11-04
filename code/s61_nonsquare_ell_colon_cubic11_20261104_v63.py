#!/usr/bin/env python3
"""One old frozen QQ cubic matrix test and exact ell-colon criterion recipe."""
import argparse,hashlib,itertools,json,math,resource,signal,subprocess,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_ELL_COLON_CUBIC11_20261104_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def add(d,k,c):
 z=d.get(k,Q())+c
 if z:d[k]=z
 else:d.pop(k,None)
def main(out):
 start=time.monotonic();deadline=start+60;resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(parents=True,exist_ok=False)
 prior=ROOT/'reports/S61_NONSQUARE_CUBIC_NINE_ROW_CRITERION_20261103_v63';prp=prior/'report.json';pr=json.loads(prp.read_text());assert pr['all_passed']and pr['chosen_chart_nonempty_by_exact_QQ_point']
 pointfile=prior/'point-matrix.txt';assert sha(pointfile)==pr['evidence_sha256']['point-matrix.txt'];pc=json.loads((prior/'circuit.json').read_text());assert pc['point']==[1,2,3]
 sourcep=ROOT/'reports/S61_BINARY_NONSQUARE_CUBIC_REPLAY_20261103_v62/source-report.json';source=json.loads(sourcep.read_text());assert sha(sourcep)==pr['input_sha256'][str(sourcep.relative_to(ROOT))]
 fullsource=ROOT/pc['source'];assert sha(fullsource)==pr['input_sha256'][str(fullsource.relative_to(ROOT))]
 mons=source['normal_monomials'];assert len(mons)==91 and all(len(m)==8 for m in mons)
 lines=pointfile.read_text().splitlines();assert lines[0]=='57 91';M=[[Q(x)for x in line.split()]for line in lines[2:]];assert len(M)==57 and all(len(row)==91 for row in M)
 Qp=ROOT/'reports/S61_BINARY_HYPERPLANE_QUADRATIC_CODIMENSION_20261102.json';qd=json.loads(Qp.read_text())
 zwp=ROOT/'reports/S61_NONSQUARE_ZW_PRODUCT_RANK19_20261103_v63/report.json';zw=json.loads(zwp.read_text());assert zw['all_exact_QQ_checks_passed']and sha(Qp)==zw['source_sha256']
 # Target monomials after Z=3v0-v1+v6=W=v3=0.
 puremons=[]
 for c in itertools.combinations_with_replacement(range(5),3):puremons.append(tuple(c.count(j)for j in range(5)))
 pureidx={m:i for i,m in enumerate(puremons)}
 def restrict(terms):
  result={}
  for mon,c in terms:
   if mon[3]:continue
   assert mon[5]<=1
   for k in range(mon[6]+1):
    ex=(mon[0]+mon[6]-k,mon[1]+k,mon[2],mon[4],mon[7]);add(result,(mon[5],ex),c*math.comb(mon[6],k)*(-3)**(mon[6]-k))
  return result
 quads=[]
 for i in[0,1,3,5]:
  terms=[]
  for mon,ts in qd['quadratic_rows_sparse'][i]:
   c=sum((Q(c)*2**ex[1]*3**ex[2]for ex,c in ts),Q());terms.append((mon,c))
  q=restrict(terms);assert all(h==0 and sum(m)==2 for h,m in q);quads.append(q)
 Qrows=[]
 for q in quads:
  for k in range(5):
   row=[Q()]*35
   for(h,m),c in q.items():ex=list(m);ex[k]+=1;row[pureidx[tuple(ex)]]+=c
   Qrows.append(row)
 G=[]
 for row in M[41:]:
  poly=restrict(zip(mons,row));v=[Q()]*40
  for(h,m),c in poly.items():
   if not h:v[pureidx[m]]+=c
   else:assert sum(m)==1;v[35+m.index(1)]+=c
  G.append(v)
 oldlight=[0,1,2,3,4,6,7];hcols=[]
 for v in oldlight:
  ex=[0]*8;ex[5]=ex[v]=1;hcols.append(mons.index(ex))
 H7=[[row[j]for j in hcols]for row in M[41:]]
 inp=out/'oldpoint-native-input.txt'
 with inp.open('w')as f:
  f.write('ELLCOLON11_FROZEN_QQ_POINT\n')
  for matrix in[Qrows,G,H7]:f.write(str(len(matrix))+' '+str(len(matrix[0]))+'\n');f.write('\n'.join(' '.join(map(str,row))for row in matrix)+'\n')
 cpp=ROOT/'code/s61_nonsquare_ell_cubic11_20261104_v63.cpp';exe=out/'numeric-cubic11'
 flags=subprocess.check_output(['pkg-config','--cflags','--libs','flint'],text=True).split()
 with(out/'compile.log').open('w')as f:subprocess.run(['c++','-O2','-std=c++17',str(cpp),*flags,'-o',str(exe)],stdout=f,stderr=subprocess.STDOUT,check=True,timeout=max(.1,deadline-time.monotonic()))
 with(out/'native.log').open('w')as f:subprocess.run([str(exe),str(inp),str(out/'native.json')],stdout=f,stderr=subprocess.STDOUT,check=True,timeout=max(.1,deadline-time.monotonic()))
 exe.unlink();n=json.loads((out/'native.json').read_text());positive=n['cubic11_rank']==11
 assert n['quad_product_rank']==19 and n['restricted_heavy_rank']==5 and n['original_heavy_rank']==7
 if positive:assert Q(n['full35_minor_determinant'])==Q(n['quad19_minor'])*Q(n['heavy5_minor'])*Q(n['nonzero11_minor'])and Q(n['nonzero11_minor'])
 # Symbolic recipe uses actual polynomial source entries, not point interpolation.
 circuit={'parametric_source':str(fullsource.relative_to(ROOT)),'source_uncompressed_sha256':pc['source_uncompressed_sha256'],'genuine_basis_rows':list(range(41,57)),
 'normal_monomials_source':str(sourcep.relative_to(ROOT)),'quad_source':str(Qp.relative_to(ROOT)),'quad_source_rows':[0,1,3,5],
 'restriction_v8':['x0','x1','x2','0','x3','H','x1-3*x0','x4'],'restricted40_column_order':{'pure35_exponents':puremons,'heavy_last5':['H*x0','H*x1','H*x2','H*x3','H*x4']},
 'restricted36_row_order':'first20 qfour timesfive x; last16 actualgenuinecubicbasis rows41..56',
 'quad19_pivot_rows':n['quad_pivot_rows'],'quad19_pivot_columns':n['quad_pivot_columns'],'genuine_heavy5_pivot_rows':n['restricted_heavy5_pivot_rows'],'remaining11_genuine_rows':n['remaining11_genuine_rows'],'remaining16_pure_columns':n['quotient16_columns'],
 'entry_definition':'R11ij equals bordered25minor of restricted36x40 using fixed24pivot dividedbyDelta24; pivot rows=quad19 plusselectedheavy5, pivotcols=quadpure19 plusheavy5',
 'source_chart':'D(DeltaQuad19*DeltaHeavy5), extra chart factors only; FULLU units not enlarged globally',
 'symbolic_rank_upper':11,'symbolic_rank11_proved_by_saved_oldpoint_nonzero35minor':positive,
 'point':[1,2,3],'no_new_point_evaluation_of57_source':True,
 'colon_hypotheses':{'FULLU':True,'raw_quadratic_cofactor84x126_rank78':True,'heavy16x7_rank7':True,'ZW_fivevariable_products_rank19_actual':True,'full_CI_or_global_ell_nonzerodivisor_unmixedness_not_assumed':True},
 'rank16_source':'actual35−19, not naive4quadCI','exact_bundle_sequence':'0->V16 --ell->B64->V48->0 onU78',
 'canonical11_source':'ell*r9 andH*ellZ,H*ellW; imagesV16=(r9|ZW,pZ|ZW,pW|ZW)',
 'quartic_exact_remainders':{'if_projection122_toV48_surjective':'afterknown15kernel,59source->16target','if_additionally_cubic11_image_rank11':'after11source,targeteliminations and48projection,48source->5target'},
 'all32_genuine_quartic_residuals_retained':True,'remaining43kernel_not_constructed':True,'full239_or_allpoint57_not_proved':True,'rankdrop_stop_no_symbolicfield_or_newpoint_retry':not positive}
 (out/'circuit.json').write_text(json.dumps(circuit,indent=2)+'\n')
 deps=[ROOT/'reports'/z/'report.json'for z in['S61_NONSQUARE_VERTEX_SCROLL_FINITE_20261103_v63','S61_NONSQUARE_SOURCE282_ATLAS_20261104_v63','S61_NONSQUARE_KERNEL22_EXTRACTOR_20261104_v63','S61_NONSQUARE_QUARTIC_122_CIRCUIT_20261103_v63']]
 paths=[Path(__file__),cpp,ROOT/'research/notes'/(ID+'.md'),pointfile,prp,prior/'circuit.json',sourcep,fullsource,Qp,zwp,*deps]
 report={'id':ID,'all_saved_source_and_exact_QQ_rank_checks_passed':True,'same_old_core_point_only':[1,2,3],'quad_product_rank':19,'restricted_heavy_rank':5,'cubic11_rank_at_saved_point':n['cubic11_rank'],'symbolic_QQ_a_e_d_rank11_proved':positive,'new_exact_ell_colon_filtration_and59x16_or48x5_criteria':True,'no_allpoint_rank_or_full239_claim':True,'no_new_point_prime_GB_field_coefficient_expansion_or_quartic_expansion':True,'rankdrop_no_retry':not positive,'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},'evidence_sha256':{p.name:sha(p)for p in out.iterdir()if p.is_file()},'seconds':time.monotonic()-start,'frozen':True}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=ROOT/'reports'/ID);a=p.parse_args();main(a.output.resolve())
