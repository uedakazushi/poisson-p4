#!/usr/bin/env python3
"""ONLY support/degree/live forecasts for10 actual pure cubics times7 lights."""
import argparse,itertools,json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
from s61_nonsquare_P100_H1_right_inverse_20261104_v63 import read,sha
from s61_nonsquare_P100_pure98_source_preflight_20261104_v63 import ZERO,addbound,mb,width,bdict
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_P100_PURE70_SOURCE_PREFLIGHT_20261104_v63'
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(57);out.mkdir(parents=True,exist_ok=False);base=ROOT/'reports'
 kb=base/'S61_NONSQUARE_P100_PURE_CUBIC_KERNEL_20261104_v63';kr=read(kb/'report.json');binding=read(kb/'source-binding.json');assert kr['complete']
 gp=base/'S61_NONSQUARE_P100_G18_PURE7_20261104_v63/joint-cubic-frame.json';bp=base/'S61_NONSQUARE_P100_PURE98_SOURCE_PREFLIGHT_20261104_v63/manifest.json'
 fdpath=base/'S61_NONSQUARE_JOINT_LELL_FREE164_20261104_v63/certificate.json';hp=base/'S61_NONSQUARE_H_SQUARE_H1_PREFLIGHT_20261104_v63/certificate.json';up=base/'S61_NONSQUARE_H1_ATOMIC_UNIT_DENOMINATORS_20261104_v63/certificate.json'
 ids=[0,2,3,4,5,7,10,12,13,15];kpaths=[kb/f'pure-cubic84-{i}.json'for i in ids];paths=[Path(__file__),kb/'report.json',kb/'source-binding.json',gp,bp,fdpath,hp,up]+kpaths
 for p in kpaths:assert sha(p)==kr['evidence_sha256'][p.name]
 B=read(gp)['physical_light_inverse7x7'];m=read(bp);fd,h,u=map(read,[fdpath,hp,up]);oldmons=binding['pure_monomials84'];physical=[0,1,2,3,4,6,7]
 nn=[tuple(mon[i]for i in physical)for mon in oldmons];j3=sorted(nn);idx3={v:i for i,v in enumerate(j3)}
 joint=h['parent_joint_target239_monomials'];j4=[tuple(v[:7])for v in joint if v[7]==0];assert len(j4)==210
 units=[(Q(1),Q(0)),(Q(1),Q(20)),(Q(1),Q(60)),(Q(203),Q(2700)),(Q(1),Q(12))];dcache={}
 def dv(poly):
  key=tuple((m[0],c)for m,c in poly)
  if key in dcache:return dcache[key]
  p={m[0]:Q(c)for m,c in poly};v=[]
  for a,b in units:
   e=0
   while p and max(p)>=1:
    r=p.copy();q={}
    while r and max(r)>=1:
     k=max(r);c=r[k]/a;q[k-1]=q.get(k-1,Q())+c
     for j,t in[(k,c*a),(k-1,c*b)]:
      r[j]=r.get(j,Q())-t
      if not r[j]:del r[j]
    if r:break
    p=q;e+=1
   v.append(e)
  assert p and max(p)==0;dcache[key]=tuple(v);return tuple(v)
 def rb(v):return ZERO if not v['numerator']else(max(m[0]for m,c in v['numerator']),dv(v['denominator']))
 bb=[[rb(v)for v in row]for row in B];Bterms=sum(len(v['numerator'])+len(v['denominator'])for row in B for v in row);Bbound=sum(width(v)or 1 for row in bb for v in row)
 def key(p):return tuple(sorted((tuple(m),str(Q(c)))for m,c in p))
 dproof={key(p):u['literal_constant_and_powers_by_denominator'][i]['a_d_F_B_C_J_powers']for i,p in enumerate(u['unique_atomic_denominators'])}
 def corebound(n,d):
  if not n:return ZERO
  k=key(d)
  if k in dproof:a,dd,f,b,c,j=dproof[k]
  else:assert len(d)==1 and d[0][0]==[0,0,0]and Q(d[0][1]);a=dd=f=b=c=j=0
  clear=max(2*m[0]+m[1]for m,q in n);nd=clear+max(2*m[2]-m[0]for m,q in n);zp=2*a-dd+f-b+c-clear
  return nd+max(zp,0),(max(-zp,0),a+dd,b+c,f,j)
 SN=fd['scalar_circuit']['nodes'];SA=fd['scalar_circuit']['atoms_num_den'];LN=fd['linear_nodes'];nf=fd['full_target_projection239_to122'];roots=m['NF_pure98_linear_roots']
 sb={0:ZERO}
 for i in m['NF_needed_scalar_ids']:
  n=SN[i];sb[i]=ZERO if n[0]=='z'else corebound(*SA[n[1]])if n[0]=='c'else mb(sb[n[1]],sb[n[2]])if n[0]=='m'else addbound([sb[n[1]],sb[n[2]]])
 scalar_retained=sum(width(sb[i])or 1 for i in m['NF_needed_scalar_ids']);scalar_atom_terms=sum(len(SA[i][0])+len(SA[i][1])for i in m['NF_needed_atom_ids'])
 # Symbolic Sym^3(B) coefficient bounds; no RF coefficient is multiplied.
 image=[]
 for mon in nn:
  d={(0,)*7:(0,(0,)*5)}
  for j,e in enumerate(mon):
   for repeat in range(e):
    terms={}
    for p,b in d.items():
     for k,c in enumerate(bb[j]):
      if c[0]>=0:
       q=list(p);q[k]+=1;q=tuple(q);terms.setdefault(q,[]).append(mb(b,c))
    d={q:addbound(v)for q,v in terms.items()}
  image.append(d)
 # Flatten only coefficient BOUNDS for the sparse98x210 normal-form operator.
 column={v:i for i,v in enumerate(nf['input_variables'])};coeff={0:{}}
 for i in m['NF_needed_linear_ids']:
  n=LN[i]
  if n[0]=='v':coeff[i]={column[n[1]]:(0,(0,)*5)}
  elif n[0]=='z':coeff[i]={}
  else:
   terms={}
   for c,j in n[1]:
    for t,b in coeff[j].items():terms.setdefault(t,[]).append(mb(sb[c],b))
   coeff[i]={t:addbound(v)for t,v in terms.items()}
 operator_profiles=[]
 for light in range(7):
  entries=0;rfupper=0;maxrow=0
  for root in roots:
   terms={}
   for t,b in coeff[root].items():
    v=joint[t];assert v[7]==0
    for k,c in enumerate(bb[light]):
     if c[0]>=0 and v[k]:
      q=list(v[:7]);q[k]-=1;terms.setdefault(idx3[tuple(q)],[]).append(mb(b,c))
   vv=[addbound(v)for v in terms.values()];entries+=len(vv);maxrow=max(maxrow,len(vv));rfupper+=sum(width(v)for v in vv)
  operator_profiles.append({'physical_light':physical[light],'sparse98x84_nonzero_support_upper':entries,'max_output_row_support_upper':maxrow,'all_numeric_operator_RF_terms_upper':rfupper})
 # Last-use value bound for the existing547-node operator, not a matrix expansion.
 def nflive(bounds):
  uses={i:0 for i in m['NF_needed_linear_ids']}
  for i in uses:
   n=LN[i]
   if n[0]=='l':
    for c,j in n[1]:uses[j]+=1
  for r in roots:uses[r]+=1
  alive=peak=0
  for i in sorted(uses):
   alive+=width(bounds[i])or 1;peak=max(peak,alive);n=LN[i]
   if n[0]=='l':
    for c,j in n[1]:
     uses[j]-=1
     if not uses[j]:alive-=width(bounds[j])or 1
  return peak
 records=[];all_output=0;maxlive=0;allqq_p_integral=True
 for kid,p in zip(ids,kpaths):
  d=read(p);vals=d['pure84_values'];vb=[rb(v)for v in vals];jp=[ZERO]*84
  for target in range(84):jp[target]=addbound([mb(vb[k],image[k][j3[target]])for k in range(84)if j3[target]in image[k]])
  jpterms=sum(width(b)or 1 for b in jp);inputdense=sum(width(b)or 1 for b in vb);inputliteral=sum(len(v['numerator'])+len(v['denominator'])for v in vals);rows=[]
  for light in range(7):
   pb=[]
   for v in j4:
    summands=[]
    for k,c in enumerate(bb[light]):
     if c[0]>=0 and v[k]:
      q=list(v);q[k]-=1;summands.append(mb(jp[idx3[tuple(q)]],c))
    pb.append(addbound(summands))
   pindex={tuple(v):i for i,v in enumerate(j4)};lb={0:ZERO}
   for i in m['NF_needed_linear_ids']:
    n=LN[i]
    if n[0]=='v':v=joint[column[n[1]]];assert not v[7];lb[i]=pb[pindex[tuple(v[:7])]]
    elif n[0]=='z':lb[i]=ZERO
    else:lb[i]=addbound([mb(sb[c],lb[j])for c,j in n[1]])
   ot=sum(width(lb[r])or 1 for r in roots);pt=sum(width(v)or 1 for v in pb);lp=nflive(lb)
   reserve=inputliteral+Bterms+scalar_atom_terms
   live=reserve+inputdense+Bbound+jpterms+scalar_retained+pt+lp+ot+10000
   maxlive=max(maxlive,live);all_output+=ot;rows.append({'physical_light':physical[light],'joint210_retained_RF_upper':pt,'NF_linear_value_peak_upper':lp,'pure98_output_RF_upper':ot,'total_live_upper_including_output_copy_and10000temps':live})
  terms=7*0+sum(r['pure98_output_RF_upper']for r in rows)
  records.append({'kernel_nonpivot_genuine':kid,'source_actual135':next(r['actual135']for r in binding['nonpivot_sources10']if r['genuine_index']==kid),'source_cached_RF_terms':inputliteral,'source84_dense_upper':inputdense,'joint84_dense_upper':jpterms,'seven_output_rows':rows,'all7_retained_output_RF_upper':terms,'stream_each_output_and_release_required':True,'max_live_upper':max(r['total_live_upper_including_output_copy_and10000temps']for r in rows)})
  for v in vals:
   for part in['numerator','denominator']:
    allqq_p_integral&=all(Q(c).denominator%32003!=0 for _,c in v[part])
 quartic={'complete_H1_cache_source_indices':list(range(113,134)),'H2free_source_labels21':h['all_remaining21_quartic_H2free_labels'],
          'source_definition':'psi=-b131*Qhat77+b77*Qhat131; every other retainedQhat_l becomes Qhat_l-(b_l/U)*chi, chi=uQhat77+vQhat131,U=(50aBC)^2. Then subtract the actual112cubic H1 lift to get pure98; keep all pure tails.',
          'frame_rule':'Qhat_l=(50/B)^2 F_l(T34*diag(Tden)*v), all33; cache113..133 already has joint/NF/H2 correction. Pure98 coefficients remain lazy, no quartic query in this preflight.',
          'source117_projection_and_chi_root_bound_by':str(hp.relative_to(ROOT)),'future_pure_quartic_coefficient_queries_NOT_prepared_or_evaluated':True}
 cert={'one_cubic7rows_per_future_job':True,'blocks10':records,'all70_outputs_if_retained_RF_upper':all_output,'max_one_row_live_upper':maxlive,'all_one_row_live_bounds_pass150k':maxlive<=150000,
       'operator_support':{'NF98x210_support_upper':sum(len(coeff[r])for r in roots),'Sym3B84x84_support_upper':sum(len(v)for v in image),'seven_numeric98x84_operator_forecasts':operator_profiles,'recommendation':'Reuse fixed244scalar/547linear DAG. Materialize no full7 numericoperator family unless its own live/output bound passes; operator forecasts are bounds, not expanded coefficients.'},
       'source_quartics21_needed_later':quartic,'source_H0_and_source_kernel_proof':'The accepted whole-curve source kernel gives literal pure cubics. Multiplicationby7physical lights preservesH0; NF117 H-degree is strictly preserved. No independence of70images is assumed.',
       'Laurent_rank_proposal_only':{'formal_variable':'z at0 is an excluded pole, not a new physical source point','old_prime':32003,'saved_kernel_QQ_coefficients_p_integral':allqq_p_integral,'required_before_positive_rank_gate':['Prove leading denominators after known row-pole clearing are32003-units for every circuit input.','Use exact truncatedLaurent arithmetic withglobalprecision9.','Recordfixed selectedminor/LU andtotalprecisionloss or determinantvaluation<=8; not8perpivot.','A nonzero leading minor coefficient then provesQQ(z)lower only. Zero/insufficientprecision gives no rank upper or geometric negative.'],'series_or_rank_calculation_executed':False},
       'coefficient_values_or_70NF_rows_computed':False,'no_new_oracle_point_prime_GB_or_automatic_coefficient_job':True}
 (out/'certificate.json').write_text(json.dumps(cert,indent=2)+'\n')
 report={'id':ID,'source_only':True,'blocks10_profiled':10,'max_one_row_live_upper':maxlive,'all70_retained_RF_upper':all_output,'operator_support_only_no_numeric_expansion':True,'no70NF_value_or_rank_job':True,'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},'certificate_sha256':sha(out/'certificate.json'),'seconds':time.monotonic()-start,'frozen':True}
 assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024;signal.alarm(0);(out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items()if k!='input_sha256'},indent=2),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
