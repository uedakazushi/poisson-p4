#!/usr/bin/env python3
"""Source-only input/bounds for G18,G39,G68 times seven lights; no values."""
import argparse,gzip,hashlib,json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
from s61_nonsquare_P100_pure98_source_preflight_20261104_v63 import addbound,mb,width,bdict,ZERO
ROOT=Path(__file__).resolve().parents[1]
ID='S61_NONSQUARE_P100_PURE21_SOURCE_PREFLIGHT_20261104_v63'
def read(p):return json.loads(p.read_text())
def sha(p):
 h=hashlib.sha256()
 with p.open('rb')as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return h.hexdigest()
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(57);out.mkdir(parents=True,exist_ok=False)
 base=ROOT/'reports';rowbase=base/'S61_BINARY_NONSQUARE_CUBIC_ROW_BASIS_20261103_v62'
 matrix=rowbase/'remaining.txt.gz';basispath=rowbase/'report.json';sp=base/'S61_BINARY_NONSQUARE_CUBIC_REPLAY_20261103_v62/source-report.json'
 ap=base/'S61_NONSQUARE_ORACLE_FRAME_ADAPTER_20261103_v63/adapter.json';hp=base/'S61_NONSQUARE_H_SQUARE_H1_PREFLIGHT_20261104_v63/certificate.json'
 rp=base/'S61_NONSQUARE_RAW_H_LINEAR_UNITS_20261104_v63/preflight.json';fp=base/'S61_NONSQUARE_CUBIC_PURE3_FUNCTIONAL_INPUT_20261104_v63/certificate.json'
 fdpath=base/'S61_NONSQUARE_JOINT_LELL_FREE164_20261104_v63/certificate.json';up=base/'S61_NONSQUARE_H1_ATOMIC_UNIT_DENOMINATORS_20261104_v63/certificate.json'
 pm=base/'S61_NONSQUARE_P100_PURE98_SOURCE_PREFLIGHT_20261104_v63/manifest.json'
 paths=[Path(__file__),matrix,basispath,sp,ap,hp,rp,fp,fdpath,up,pm]
 basis,source,ad,h,raw,old,fd,units,prior=map(read,[basispath,sp,ap,hp,rp,fp,fdpath,up,pm])
 assert basis['all_passed']and basis['remaining_matrix_shape']==[57,91]
 mons=source['normal_monomials'];assert len(mons)==91 and all(len(m)==8 for m in mons)
 pure=[j for j,m in enumerate(mons)if m[5]==0];heavy=[j for j,m in enumerate(mons)if m[5]]
 assert len(pure)==84 and len(heavy)==7
 selected=[41,44,48];genuine=[0,3,7];actual=[18,39,68]
 rows={};H=hashlib.sha256();linecount=0
 with gzip.open(matrix,'rb')as f:
  first=f.readline();H.update(first);assert first.strip()==b'57 91 57'
  for row in range(57):
   rr=[]
   for col in range(91):
    line=f.readline();assert line;H.update(line);linecount+=1
    if row in selected:
     t=line.decode().split();n=int(t[0]);assert len(t)==1+4*n
     rr.append([[list(map(int,t[2+4*k:5+4*k])),str(Q(t[1+4*k]))]for k in range(n)])
   if row in selected:rows[row]=rr
  assert not f.read()
 assert linecount==5187 and H.hexdigest()==basis['remaining_matrix_sha256']
 defs=h['all112_cubic_source319_definitions'];records=[];products=[]
 for pos,i,eq in zip(selected,genuine,actual):
  assert ad['basis57_source_definitions'][pos]=={'kind':'genuine_cubic','original_equation':eq}
  assert basis['remaining_original_row_indices'][pos]==defs[7*i][1]
  assert all(not rows[pos][j]for j in heavy)and all(not p for p in raw['raw16x7'][i])
  assert old['relative16_genuine_indices']==genuine and old['raw_genuine_source135_indices']==actual
  rr={'genuine_index':i,'basis57_position':pos,'original135_equation':eq,'original63_normalized_row':basis['remaining_original_row_indices'][pos],
      'oracle_to_saved_row_unit':ad['basis57_encoded_source_units'][pos],'full_pure84_coefficients':[rows[pos][j]for j in pure],
      'seven_full_Hlinear_coefficients_literal_zero':True}
  records.append(rr)
  for j,v in enumerate([0,1,2,3,4,6,7]):
   label=defs[7*i+j];assert label==['cubic_normalized',basis['remaining_original_row_indices'][pos],eq,v]
   products.append({'source319_index':174+7*i+j,'definition':label,'physical_light':v,'saved_source_record':len(records)-1})
 assert old['actual_G11_source_subtraction_coefficients']==['0','0','0']
 # Bounds only: no coefficient is specialized at a curve point or expanded in z.
 def bits(p):return max([max(abs(Q(c).numerator).bit_length(),Q(c).denominator.bit_length())for m,c in p]+[0])
 def kpoly(p):return tuple(sorted((tuple(m),str(Q(c)))for m,c in p))
 denproof={kpoly(p):units['literal_constant_and_powers_by_denominator'][i]['a_d_F_B_C_J_powers']for i,p in enumerate(units['unique_atomic_denominators'])}
 def corebound(num,den):
  if not num:return ZERO
  key=kpoly(den)
  if key not in denproof:
   assert len(den)==1 and den[0][0]==[0,0,0]and Q(den[0][1]);powers=[0]*6
  else:powers=denproof[key]
  aa,dd,ff,bb,cc,jj=powers;clear=max(2*m[0]+m[1]for m,c in num)
  nd=clear+max(2*m[2]-m[0]for m,c in num);zp=2*aa-dd+ff-bb+cc-clear
  return nd+max(zp,0),(max(-zp,0),aa+dd,bb+cc,ff,jj)
 pb=[corebound(p,[[[0,0,0],'1']])for rr in records for p in rr['full_pure84_coefficients']]
 pi=addbound(pb);bb=prior['bounds']['physical_to_joint_linear_common'];bb=(bb['numerator_degree_upper'],tuple(bb['known_linear_denominator_powers_upper']))
 productbound=mb(pi,(4*bb[0],tuple(4*x for x in bb[1])))
 # Reuse the accepted literal normal-form operator and its frozen scalar closure.
 LN=fd['linear_nodes'];SN=fd['scalar_circuit']['nodes'];SA=fd['scalar_circuit']['atoms_num_den']
 purecols=h['target122_Hdegree_classes']['0'];op=fd['full_target_projection239_to122'];roots=[op['output122_linear_roots'][j]for j in purecols]
 need=set(prior['NF_needed_linear_ids']);sneed=set(prior['NF_needed_scalar_ids']);satom=set(prior['NF_needed_atom_ids'])
 sb={0:ZERO};sab={i:corebound(*SA[i])for i in satom}
 for i in sorted(sneed):
  n=SN[i];sb[i]=ZERO if n[0]=='z'else sab[n[1]]if n[0]=='c'else mb(sb[n[1]],sb[n[2]])if n[0]=='m'else addbound([sb[n[1]],sb[n[2]]])
 lb={0:ZERO}
 for i in sorted(need):
  n=LN[i];lb[i]=ZERO if n[0]=='z'else productbound if n[0]=='v'else addbound([mb(sb[c],lb[j])for c,j in n[1]])
 totalout=sum(width(lb[i])for i in roots)
 # Exact typed dependency audit that pure inputs do not create H1 or H2 outputs.
 inputmap={v:j for j,v in enumerate(op['input_variables'])};targetmons=h['parent_joint_target239_monomials'];deps={0:set()}
 allneed=set();todo=op['output122_linear_roots'][:]
 while todo:
  i=todo.pop()
  if i in allneed:continue
  allneed.add(i);n=LN[i]
  if n[0]=='l':todo.extend(j for c,j in n[1])
 for i in sorted(allneed):
  n=LN[i]
  if n[0]=='v':assert n[1]in inputmap
  deps[i]={targetmons[inputmap[n[1]]][7]}if n[0]=='v'else set().union(*(deps[j]for c,j in n[1]))if n[0]=='l'and n[1]else set()
 assert all(not deps[op['output122_linear_roots'][j]]or deps[op['output122_linear_roots'][j]]=={int(deg)}for deg,cols in h['target122_Hdegree_classes'].items()for j in cols)
 counts={'full_saved_pure_coefficients':252,'nonzero_coefficients':sum(bool(p)for r in records for p in r['full_pure84_coefficients']),
         'QQ_core_terms':sum(len(p)for r in records for p in r['full_pure84_coefficients']),
         'max_entry_terms':max(len(p)for r in records for p in r['full_pure84_coefficients']),
         'max_coefficient_bits':max(bits(p)for r in records for p in r['full_pure84_coefficients']),
         'max_core_total_degree':max(sum(m)for r in records for p in r['full_pure84_coefficients']for m,c in p)}
 data={'physical_light_order':[0,1,2,3,4,6,7],'physical_pure_cubic_monomials84':[mons[j]for j in pure],
       'three_full_sources':records,'twentyone_actual_source319_products':products,'profiles':counts,
       'normalization':'Already normalized FULL Tn 57-by-91 matrix; identical to FrameAdapter basis57_coefficient. No extra row-unit multiplication or division is now applied.',
       'old_restricted_PURE3_not_used_for_full_polynomials':True,'source_value_or_oracle_query_count':0}
 (out/'source-input.json').write_text(json.dumps(data,separators=(',',':'))+'\n')
 cert={'full_source_and_all_Hzeros_verified':True,'all21_sources_actual_pure_quartics_without_H1_lift':True,
       'NF117_preserves_Hdegree_by_literal_linear_dependency':True,'NF_pure98_target122_labels':purecols,
       'NF_input_joint_not_original_physical':True,'required_future_coordinate_transport':'Transport each full normalized physical cubic and each physical light through the accepted light inverse, assemble joint210, then apply the fixed NF117 pure98 linear operator.',
       'source_ideal_transport':'No restricted-ZW representative is lifted. Each source is the literal full polynomial Gi*vj. The normal form differs from that actual source by the already certified quadratic source rows, so its image belongs to the original quartic row span.',
       'bounds':{'original_252_coefficients_live_upper':sum(width(b)for b in pb),'common_P100_pure_cubic_coefficient':bdict(pi),'accepted_physical_to_joint_linear_common':bdict(bb),'one_joint210_row_common':bdict(productbound),'one_pure98_output_terms_upper':totalout,'all21_output_terms_upper_if_retained':21*totalout,'all7_of_one_cubic_output_terms_upper_if_retained':7*totalout},
       'NF_nodes_reused_from_frozen_preflight':{'linear':len(need),'scalar':len(sneed),'atoms':len(satom)},
       'no_curve_evaluation_new_oracle_query_rank_Gaussian_GB_point_or_minor':True,'next_coefficient_job_NOT_launched':True}
 (out/'certificate.json').write_text(json.dumps(cert,indent=2)+'\n')
 report={'id':ID,'full_source_input_bound':True,'profiles':counts,'no_new_source_values':True,'no_actual_curve_coefficients_or_normal_forms_evaluated':True,
         'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},'uncompressed_full57_matrix_sha256':H.hexdigest(),'source_input_sha256':sha(out/'source-input.json'),'certificate_sha256':sha(out/'certificate.json'),
         'seconds':time.monotonic()-start,'frozen':True,'arithmetic_job_requires_separate_authorization':True}
 assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024;signal.alarm(0)
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
