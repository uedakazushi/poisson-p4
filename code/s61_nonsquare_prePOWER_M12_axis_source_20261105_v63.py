#!/usr/bin/env python3
"""SOURCE ONLY added G134/target17 axis supports, reusing frozen579 coefficients."""
import argparse,functools,hashlib,json,math,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_PREPOWER_M12_AXIS_SOURCE_20261105_v63'
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def dump(p,x):p.write_text(json.dumps(x,indent=2)+'\n')
def main(out):
 started=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(parents=True,exist_ok=False);pins={}
 def load(rel):p=ROOT/rel;pins[rel]=sha(p);return json.loads(p.read_text())
 def text(rel):p=ROOT/rel;pins[rel]=sha(p);return p.read_text()
 def subtree(s,key):
  tag=json.dumps(key)+':';assert s.count(tag)==1;i=s.index(tag)+len(tag)
  while s[i].isspace():i+=1
  return json.JSONDecoder().raw_decode(s,i)[0]
 def skip(s,i):
  while s[i].isspace():i+=1
  if s[i]in'[{':
   depth=0;quoted=False;escape=False
   for j in range(i,len(s)):
    c=s[j]
    if quoted:
     if escape:escape=False
     elif c=='\\':escape=True
     elif c=='"':quoted=False
    elif c=='"':quoted=True
    elif c in'[{':depth+=1
    elif c in']}':
     depth-=1
     if not depth:return j+1
  return json.JSONDecoder().raw_decode(s,i)[1]
 def ws(s,i):
  while s[i].isspace()or s[i]==',':i+=1
  return i
 qr='reports/S61_NONSQUARE_CUBIC_Q2_L_UNIT_COMPLETION_20261104_v63/certificate.json';q=load(qr);nd=q['scalar_DAG']['nodes'];cols=q['new22_target_original32_columns']
 cp_rel='reports/S61_NONSQUARE_CUBIC_ZW25_UNITS_20261104_v63/checkpoint.json';ct=text(cp_rel);tag='"matrix25x32":';i=ct.index(tag)+len(tag);i=ws(ct,i);assert ct[i]=='[';i+=1;spans={}
 for r in range(25):
  i=ws(ct,i);assert ct[i]=='[';i+=1
  for c in range(32):i=ws(ct,i);j=skip(ct,i);spans[r,c]=(i,j);i=j
  i=ws(ct,i);assert ct[i]==']';i+=1
 pd_rel='reports/S61_NONSQUARE_CUBIC_ZW25_PAIR_DUAL_PREFLIGHT_20261104_v63/certificate.json';pdtext=text(pd_rel);pd=subtree(pdtext,'Lframe_dual6x2');del pdtext
 hc_rel='reports/S61_NONSQUARE_CUBIC_GENUINE14_H5_PREFLIGHT_20261104_v63/certificate.json';ht=text(hc_rel);hc=subtree(ht,'raw_compact14x5');rho=subtree(ht,'known_row_units_rho14');del ht
 assert rho[13]['a_d_F_B_C_J_powers']==[1,0,1,0,2,0];c13=Q(rho[13]['constant']);assert c13==Q(8000000,16119)
 hr='reports/S61_NONSQUARE_PURE_KERNEL12_SOURCE_PREFLIGHT_20261104_v63/prePOWER-h9-source.json';h9=load(hr)
 bp={(1,0,1):Q(50),(0,2,1):Q(18),(0,1,1):Q(-9),(0,0,1):Q(-2),(0,1,0):Q(60),(0,0,0):Q(-40)}
 atoms=[];profiles={};visits=[]
 def make(num,den,ref):
  dd={tuple(m):Q(c)for m,c in den if Q(c)};assert dd;aa=bb=0
  if len(dd)==1:
   m,u=next(iter(dd.items()));assert m[1:]==(0,0);aa=m[0]
  else:
   assert set(dd)==set(bp);u=dd[next(iter(bp))]/bp[next(iter(bp))];assert all(dd[m]==u*bp[m]for m in bp);bb=1
  assert u and all(m[1]>=0 and m[2]>=0 for m,c in num)
  # Only absolute L1 metadata and surviving exponent intervals, never a signed
  # new axis coefficient vector.
  ax=[(m[0]-aa,abs(Q(c)/(u*((-40)**bb))))for m,c in num if m[1]==0 and m[2]==0 and Q(c)]
  if ax:
   D=math.lcm(*(c.denominator for e,c in ax));p=dict(lo=min(e for e,c in ax),hi=max(e for e,c in ax),terms=len(ax),D=D,L1=sum(c.numerator*(D//c.denominator)for e,c in ax))
  else:p=None
  b=dict(axis=p,Bpole=bb,generic_a_lower=min([m[0]-aa for m,c in num if Q(c)]+[0]),generic_nonnegative_e_d_after_aBclear=True)
  atoms.append(dict(reference=ref,literal_numerator_terms=len(num),axis_surviving_literal_terms=len(ax),denominator_QQ_unit=str(u),denominator_a_power=aa,denominator_B_power=bb,bound=b));return b
 def scale(x,c=Q(1),a=0):
  z=dict(x);z['generic_a_lower']+=a
  if x['axis']:
   p=dict(x['axis']);p['lo']+=a;p['hi']+=a;p['D']*=c.denominator;p['L1']*=abs(c.numerator);z['axis']=p
  return z
 def combine(x,y,mul=False):
  a=x['axis'];b=y['axis'];pole=x['Bpole']+y['Bpole']if mul else max(x['Bpole'],y['Bpole']);lower=x['generic_a_lower']+y['generic_a_lower']if mul else min(x['generic_a_lower'],y['generic_a_lower'])
  if mul and(a is None or b is None):p=None
  elif a is None:p=b
  elif b is None:p=a
  elif mul:
   lo=a['lo']+b['lo'];hi=a['hi']+b['hi'];p=dict(lo=lo,hi=hi,terms=min(a['terms']*b['terms'],hi-lo+1),D=a['D']*b['D'],L1=a['L1']*b['L1'])
  else:
   lo=min(a['lo'],b['lo']);hi=max(a['hi'],b['hi']);D=math.lcm(a['D'],b['D']);p=dict(lo=lo,hi=hi,terms=min(a['terms']+b['terms'],hi-lo+1),D=D,L1=a['L1']*(D//a['D'])+b['L1']*(D//b['D']))
  return dict(axis=p,Bpole=pole,generic_a_lower=lower,generic_nonnegative_e_d_after_aBclear=True)
 @functools.lru_cache(None)
 def imp(i):
  a=nd[i];op=a[0]
  if op=='z':z=dict(axis=None,Bpole=0,generic_a_lower=0,generic_nonnegative_e_d_after_aBclear=True)
  elif op=='one':z=make([[[0,0,0],'1']],[[[0,0,0],'1']],dict(constant=1))
  elif op in['a','s','m']:z=combine(imp(a[1]),imp(a[2]),op=='m')
  else:
   if op=='checkpoint':i0,i1=spans[a[1],a[2]];v=json.JSONDecoder().raw_decode(ct,i0)[0];ref=dict(path=cp_rel,row=a[1],column=a[2],text_span=[i0,i1]);assert skip(ct,i0)==i1
   elif op=='accepted_pair_second':v=pd[a[1]][1];ref=dict(path=pd_rel,json_pointer=f'/Lframe_dual6x2/{a[1]}/1')
   else:assert op=='new_unit_atom';v=q['scalar_DAG']['new_unit_atoms_num_den'][a[1]];ref=dict(path=qr,json_pointer=f'/scalar_DAG/new_unit_atoms_num_den/{a[1]}')
   z=make(*v,ref)
  profiles[i]=z;visits.append(i);return z
 def src(r,t):return imp(q['remaining15x22_scalar_roots'][r][cols.index(t)])
 h17=make(h9['source_h9'][h9['original_target_labels'].index(17)],[[[0,0,0],'1']],dict(path=hr,target=17))
 c={j:Q(rho[j]['constant'])for j in[6,8,9,10,11,12]};assert rho[11]==rho[9]
 addedcol=[h17,scale(src(1,17),a=-1),scale(src(3,17),a=-1),src(6,17),src(2,17),scale(src(11,17),1/c[10]),scale(combine(src(12,17),src(10,17)),1/c[11]),scale(src(9,17),1/c[8]),scale(src(7,17),1/c[6]),scale(src(13,17),1/c[12],a=1),scale(src(10,17),1/c[9])]
 # The B row is a subtraction. Its unsigned support/L1 bound is the same sum.
 addedcol[1]=combine(addedcol[1],scale(addedcol[2],Q(2)))
 pure=[22,11,13,25,15,19,14,18];bottom=[src(14,t)for t in pure]
 for j in[1,3,2]:bottom.append(make(hc[13][j],[[[0,0,0],str(1/c13)]],dict(path=hc_rel,compact=13,heavy_column=j,known_exact_rho13_over_L=str(c13))))
 corner=src(14,17)
 axrel='reports/S61_NONSQUARE_PREPOWER_K92_WHOLE_AXIS_SOURCE_V2_20261105_v63/whole-axis-literal-source.json';axis=load(axrel);assert axis['axis_total_integer_coefficients']==579
 termrel='reports/S61_NONSQUARE_PREPOWER_K92_ORIGINAL_TERMINAL_FACE_SOURCE_20261105_v63/terminal-face-and-normalization.json';term=load(termrel);rs=term['row_dual'];cs=term['column_dual'];oldD=list(map(int,axis['axis_row_denominators']));assert sum(rs)+sum(cs)==-5
 rb=max([0]+[-(b['axis']['lo']+cs[j])for j,b in enumerate(bottom)if b['axis']]);cb=max([0]+[-(b['axis']['lo']+rs[j])for j,b in enumerate(addedcol)if b['axis']]+([-corner['axis']['lo']-rb]if corner['axis']else[]))
 Br=max([b['Bpole']for b in bottom]+[0]);Bc=max([b['Bpole']for b in addedcol]+[corner['Bpole']]);extraD=[math.lcm(oldD[i],b['axis']['D']if b['axis']else 1)//oldD[i]for i,b in enumerate(addedcol)]
 D12=math.lcm(*[b['axis']['D']if b['axis']else 1 for b in bottom+[corner]])
 rowL=[0]*12;maxdegree=0;entries=[]
 for rec in axis['axis_coefficients']:
  i=rec['row'];poly=rec['axis_integer_polynomial'];rowL[i]+=extraD[i]*sum(abs(int(z))for e,z in poly);maxdegree=max(maxdegree,max([e for e,z in poly]+[0]))
 for i,b in enumerate(addedcol):
  p=b['axis'];shift=rs[i]+cb;L=0 if p is None else p['L1']*(oldD[i]*extraD[i]//p['D'])*(40**Bc)
  if p:assert p['lo']+shift>=0;maxdegree=max(maxdegree,p['hi']+shift)
  rowL[i]+=L;entries.append(dict(row=i,column=11,source_axis_bound=p,a_shift=shift,known_B_axis_power=Bc,integer_L1_bound=str(L)))
 for j,b in enumerate(bottom+[corner]):
  p=b['axis'];shift=rb+(cs[j]if j<11 else cb);Bpower=Br+(Bc if j==11 else 0);L=0 if p is None else p['L1']*(D12//p['D'])*(40**Bpower)
  if p:assert p['lo']+shift>=0;maxdegree=max(maxdegree,p['hi']+shift)
  rowL[11]+=L;entries.append(dict(row=11,column=j,source_axis_bound=p,a_shift=shift,known_B_axis_power=Bpower,integer_L1_bound=str(L)))
 H=math.prod(max(1,z)for z in rowL);DB=12*maxdegree;HB=H.bit_length();TB=(DB+4).bit_length()+2*HB+3
 axisinput=579+sum(0 if x['source_axis_bound']is None else x['source_axis_bound']['terms']for x in entries)
 owned=144*(DB+1)+5*(2*DB+1)+2*axisinput+15000
 patchslots=sum((11-j)**2*(min(DB,(j+2)*maxdegree)+1)for j in range(11));raw_upper=2*1024*1024+(patchslots+owned)*((TB*30103+99999)//100000+48)
 source=dict(added_column11=addedcol,added_G134_row11=bottom,corner=corner,all_atom_profiles=atoms,finite_DAG_profiles=[dict(node=i,definition=nd[i],bound=profiles[i])for i in sorted(profiles)],
  old_M11_precondition='row1+2row2 applied to BOTH old11 block and new column; G134 row unchanged',
  exact_H3_carry='rho13=(8000000/16119)*aFC², normalized H3 order[28,30,29] is compact13 indices[1,3,2] times that QQconstant',
  checkpoint_selective_decode=True,no_old572_replay=True,no_signed_new_axis_coefficient_vector=True,
  generic_nonnegative_e_d_after_literal_a_B_clearing=True,B_axis_image='-40 known QQ unit',no_kappa3_or_unknown_inverse_in_input=True)
 dump(out/'added-source-profiles.json',source)
 forecast=dict(new_row_a_shift=rb,new_column_a_shift=cb,new_row_B_clear=Br,new_column_B_clear=Bc,old11_row_extra_QQ_integer_scales=extraD,new_row_LCD=str(D12),new_entries=entries,
  determinant_identity='detC12=(prod(oldD_i*extraD_i))*D12*(-40)^(Brow+Bcol)*a^(-5+rnew+cnew)*detM12_axis',
  old579_coefficient_source_reference=axrel,known_all_row_integer_L1_upper=[str(z)for z in rowL],all_minor_L1_upper=str(H),input_degree_upper=maxdegree,minor_degree_upper=DB,all_minor_coefficient_bits_upper=HB,
  uncombined_integer_product_and_longdivision_bits_upper=TB,owned_cells_upper=owned,all_patch_coefficient_slots_upper=patchslots,new_RAW_upper=raw_upper,
  prospective_slot_gate_passes=owned<=150000,prospective_RAW_gate_passes=raw_upper<32*1024*1024,
  fixed_axis_det12_NOT_computed=True,axis_zero_would_be_inconclusive=True,axis_not_a_FULLU_point=True)
 dump(out/'whole-axis-source-worker-forecast.json',forecast)
 peer=load('reports/S61_NONSQUARE_PREPOWER_TARGET17_THREE_CRAMER_ENTRY_20261105_v63_REVIEW/report.json');assert peer['all_passed']
 pins[str(Path(__file__).relative_to(ROOT))]=sha(Path(__file__));note=ROOT/f'research/notes/{ID}.md';pins[str(note.relative_to(ROOT))]=sha(note)
 total=sum(p.stat().st_size for p in out.iterdir());assert total<32*1024*1024
 rep=dict(id=ID,all_passed=True,frozen=True,source_only=True,known579_axis_coefficients_reused=True,added_source_DAG_nodes=len(profiles),no_axis_values_or_determinant_rank_point_prime_GCD=True,
  generic_axis_legal_after_named_unit_clearing=True,exact_H3_carry_checked=True,input_degree_upper=maxdegree,minor_degree_upper=DB,minor_bits_upper=HB,temporary_bits_upper=TB,owned_upper=owned,new_RAW_upper=raw_upper,
  total_new_RAW_before_report=total,input_sha256=pins,evidence_sha256={p.name:sha(p)for p in out.iterdir()},seconds=time.monotonic()-started)
 dump(out/'report.json',rep);signal.alarm(0);print(json.dumps({k:v for k,v in rep.items()if k not in['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
