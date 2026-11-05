#!/usr/bin/env python3
"""Whole ONE554 positive/source/height gate; no signed response values."""
import argparse,json,math,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
import ambient_B18_ell2_T2_mu_20261105_v63 as H
from ambient_B18_ell2_T2_G3_active4_source8_fold_20261105_v63 import count
ROOT=H.ROOT;STEM='AMBIENT_B18_ELL2_T2_UPPER554_WHOLE_INPUT_GATE_20261105_v63'
def main(out):
 start=time.monotonic();signal.alarm(60)
 resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024))
 out.mkdir(exist_ok=False,parents=True);pins={};raw=0;cursor={};fmax=0;imax=0
 def read(p):
  nonlocal raw
  pins[p]=H.sha(ROOT/p);raw+=2*(ROOT/p).stat().st_size;assert raw<32*1024*1024
  return json.loads((ROOT/p).read_text())
 def field(p,k):
  nonlocal raw
  pins[p]=H.sha(ROOT/p);raw+=2*(ROOT/p).stat().st_size;assert raw<32*1024*1024
  return H.selected(p,k)
 def fck(q):
  nonlocal fmax
  fmax=max(fmax,abs(q.numerator).bit_length(),q.denominator.bit_length())
  if fmax>152:
   H.dump(out/'BEFORE_RESPONSE_POSITIVE_F_GUARD.json',{'cursor':cursor,'uncommitted_positive_QQ':str(q),'no_signed_response_values':True});raise RuntimeError('positiveF152STOP')
  return q
 def pa(a,b):return fck(a+b)
 def pm(a,b):return fck(a*b)
 def ick(n):
  nonlocal imax
  imax=max(imax,abs(n).bit_length());assert abs(n).bit_length()<=384;return n
 Z=(1,0)
 def ex(q):q=F(q);return (ick(q.denominator),ick(abs(q.numerator))) if q else Z
 def ea(a,b):
  if not a[1]:return b
  if not b[1]:return a
  d=ick(math.lcm(a[0],b[0]));u=ick(a[1]*(d//a[0]));v=ick(b[1]*(d//b[0]));return d,ick(u+v)
 def em(a,b):return (ick(a[0]*b[0]),ick(a[1]*b[1])) if a[1] and b[1] else Z
 def decode(q):return int(q[0]),int(q[1])
 def esum(xs):
  z=Z
  for q in xs:z=ea(z,decode(q))
  return z
 np='reports/AMBIENT_B18_NATIVE_ONE554_POSITIVE_INPUT_20261105_v63/'
 op='reports/AMBIENT_B18_ELL2_T2_UPPER554_OLD_CARRY_INPUT_GATE_20261105_v63/'
 nr=read(np+'report.json');orr=read(op+'report.json');assert nr['all_passed'] and orr['all_passed']
 assert pins[np+'report.json']=='e9dbeecb12be941d83a97b701ce7b723c6db06200bdc698d8233300d5b6e11bd' and pins[op+'report.json']=='8373bcb85e3fb3f2bff7751e6f451d69db5850fd2747537102b82af69bc29d5d'
 n=read(np+'ONE554-positive-four-native-and-nonzero644-outer-input.json');o=read(op+'ONE554-OLD-two-columns-and57-FIRST-positive-component-input.json');p=n['full_primitive_input_J18_N_G_nu_GammaD']
 assert n['selected_case619']==o['selected_case619']==554 and p['Gamma_two56']==o['Gamma_two'] and p['D_two_by10_by56']==o['D_two_by10_by56'] and n['all29_a1_kept_UNEVALUATED']==o['full_a1_29']
 sr='reports/AMBIENT_B18_ELL2_T2_UPPER554_OLD_CARRY_SOURCE_20261105_v63/ONE554-six-c0-OLD-and-full29-a1-known-FIRST-source.json'
 c0=field(sr,'original768_source_all6');assert c0==sorted(o['OLD_source2']+[[r['original768'],r['constant_source_coefficient']] for r in n['parameterized_four_native_source_records']])
 assert n['plane_a1_pivot5_original644']==[[52,'-175/3'],[53,'-377/60']]
 assert read('reports/AMBIENT_B18_NATIVE92_ROOT_DIRECT_REVIEW_20261105_v63/report.json')['all_passed']
 hp='reports/AMBIENT_B18_NATIVE_ONE554_COMMON_DEN_HEIGHT_20261105_v63/'
 hr=read(hp+'report.json');assert hr['all_passed'];height=read(hp+'native-C72-known-input-common-denominator-and-direct-path-height-proof.json')
 assert height['committed_prior_bits']==339 and height['temporary_prior_bits']==1018
 qpath='reports/AMBIENT_B18_ELL2_T2_UPPER617_FULL_COORDINATE_SOURCE_20261105_v63/known-Psi-point-correction-and-full-b1-source-annihilation.json';Q=field(qpath,'Q_two6')
 source_before_positive=count([n,o,p,nr,orr,height,Q,c0,pins])+30000;assert source_before_positive<60000
 basis=[tuple(v) for v in p['basis18']];bi={m:i for i,m in enumerate(basis)};known=p['N_known_positions'];assert known==[0,1,2,3,4,5,6,15,16] and len(basis)==18
 active=[0,1,4,5,7,8];vf={v:bi[tuple(int(i==v) for i in range(9))] for v in active};assert vf=={0:1,1:2,4:3,5:4,7:5,8:6}
 pairs=[(0,7),(0,8),(1,7),(1,8),(4,7),(4,8),(5,7),(5,8),(7,7),(7,8)]
 assert all(bi[tuple(int(i==a)+int(i==b) for i in range(9))]==7+j for j,(a,b) in enumerate(pairs))
 assert bi[(0,0,0,0,0,0,0,2,1)]==17
 N=[[[[fck(abs(F(q))) for q in v] for v in cell] for cell in row] for row in p['N_known12_by12_by9_by4']]
 G=[[fck(abs(F(q))) for q in row] for row in p['G8_C72']];nu=[[fck(abs(F(q))) for q in row] for row in p['nu_two_J18']]
 Gamma=[{i:fck(abs(F(q))) for i,q in row if i<48} for row in p['Gamma_two56']]
 D=[[{i:fck(abs(F(q))) for i,q in enumerate(row) if q!='0' and i<48} for row in ch] for ch in p['D_two_by10_by56']]
 assert all(not row[0] and not row[vf[7]] and not row[vf[8]] for row in nu)
 for a in nu:
  for b in nu:
   assert all(tuple(x+y for x,y in zip(basis[i],basis[j])) not in bi for i,u in enumerate(a) if u for j,v in enumerate(b) if v)
 expectedR={tuple(k):{i:(int(d),int(u)) for i,d,u in row} for k,row in n['reused68_positive_covectorL_rows']};assert len(expectedR)==68
 rows={};checked_R=0
 # Independent normalized nonnegative-Fraction replay, grouped first by
 # physical output/trace entry rather than the producer's D/U operator code.
 for key,eb in sorted(expectedR.items()):
  kind,t,pair,target=key;cov=Gamma[t] if kind==0 else D[t][pair];ans={}
  for source in range(12):
   for im in range(4):
    cursor={'stage':'independent_positive68_Lrow_numeric','key':list(key),'inputcol':4*source+im};z=F(0);den,u=eb.get(4*source+im,(1,0))
    for physical in range(12):
     for kp,nb in enumerate(known):
      for mask,x in enumerate(N[source][physical][kp]):
       a=cov.get(4*physical+(mask^im),F(0))
       if not a or not x:continue
       collision=mask&im
       if not collision:
        if nb!=target:continue
        term=pm(a,x);assert den%term.denominator==0;z=pa(z,term)
       elif collision in [1,2]:
        for vb,y in enumerate(nu[collision-1]):
         if not y or tuple(a+b for a,b in zip(basis[nb],basis[vb]))!=basis[target]:continue
         term=pm(pm(a,x),y);assert den%term.denominator==0;z=pa(z,term)
       else:assert collision==3
    assert z==F(u,den);checked_R+=1
    if z:ans[4*source+im]=z
  rows[key]=ans
 assert checked_R==3264
 expectedS=[{i:(int(d),int(u)) for i,d,u in row} for row in n['positive_weighted_intrinsic_S18']];S=[{} for _ in basis]
 def sadd(base,col,term):
  den,u=expectedS[base].get(col,(1,0));assert den%term.denominator==0;S[base][col]=pa(S[base].get(col,F(0)),term)
 for rec in n['parameterized_four_native_source_records']:
  im=rec['input_trace_mask'];weight=fck(abs(F(rec['constant_source_coefficient'])))
  for kind,row,g,sign in rec['intrinsic_source_entries']:
   assert kind=='G8' and sign in [-1,1]
   for gb in range(18):
    for mask in range(4):
     x=G[g][4*gb+mask]
     if not x:continue
     collision=mask&im;col=4*row+(mask^im);cursor={'stage':'independent_positive_weightedS','source':rec['native128'],'base':gb,'mask':mask}
     if not collision:sadd(gb,col,pm(weight,x))
     elif collision in [1,2]:
      for vb,y in enumerate(nu[collision-1]):
       m=tuple(a+b for a,b in zip(basis[gb],basis[vb]))
       if y and m in bi:sadd(bi[m],col,pm(pm(weight,x),y))
     else:assert collision==3
 assert not S[0]
 for base in range(18):
  for col in range(48):den,u=expectedS[base].get(col,(1,0));assert S[base].get(col,F(0))==F(u,den)
 def dot(a,b):
  z=F(0)
  for i,x in a.items():
   if i in b:z=pa(z,pm(x,b[i]))
  return z
 cubic=[(0,17),(5,16),(6,15),(15,6),(16,5)];checked_fine=0
 for t in range(2):
  for j,(a,b) in enumerate(cubic):
   q=dot(rows[0,t,-1,a],S[b]);d,u=decode(n['positive_Gamma_weighted_cubic_two_by5'][t][j]);assert q==F(u,d);checked_fine+=1
  for k,(a,b) in enumerate(pairs):
   vals=[dot(rows[1,t,k,0],S[7+k])]
   for x,y in ([(a,b)] if a==b else [(a,b),(b,a)]):vals.append(dot(rows[1,t,k,vf[x]],S[vf[y]]))
   assert len(vals)==len(n['positive_D_weighted_second_two_by10_unique_first_decomp'][t][k])
   for q,v in zip(vals,n['positive_D_weighted_second_two_by10_unique_first_decomp'][t][k]):d,u=decode(v);assert q==F(u,d);checked_fine+=1
 assert checked_fine==68
 # Rebuild the native alternate algorithm common-denominator proof from the
 # literal signed primitive cache, without a signed response calculation.
 def values(v):
  if isinstance(v,(list,tuple)):
   for x in v:yield from values(x)
  else:yield F(v)
 A=[q for row in p['Gamma_two56'] for i,q in row]+[q for ch in p['D_two_by10_by56'] for row in ch for q in row]
 classes={'known_partial_N9':p['N_known12_by12_by9_by4'],'canonical_G8':p['G8_C72'],'monics':p['nu_two_J18'],'Gamma_and_D_scalar_coefficients':A,'native_c0_constants':[v['constant_source_coefficient'] for v in n['parameterized_four_native_source_records']]}
 rebuilt={}
 for key,v in classes.items():
  xs=list(values(v));den=1
  for x in xs:den=ick(math.lcm(den,x.denominator))
  u=max([den]+[ick(abs(x.numerator)*(den//x.denominator)) for x in xs]);rebuilt[key]=(den,u)
  old=height['scalar_operand_classes'][key];assert (str(den),str(u))==(old['D'],old['U_including_scalar_one'])
 den=1;u=2**28
 for d,v in rebuilt.values():den=ick(den*d);u=ick(u*v)
 assert str(den)==height['common_D'] and str(u)==height['common_absolute_scaled_U_with_path2pow28'] and max(den.bit_length(),u.bit_length())==339
 imax=0 # Separate the fine assembly metric from the coarse339 class book.
 # Final target scalars use the verified FINE component envelopes. The coarse
 #339 only bounds native intermediates; it is not multiplied by OLD orQ.
 ng=[];nd=[];outer=[];upper=[]
 for t in range(2):
  ng.append(esum(n['positive_Gamma_weighted_cubic_two_by5'][t]));nd.append(esum([q for row in n['positive_D_weighted_second_two_by10_unique_first_decomp'][t] for q in row]));outer.append(esum(n['positive_outer644_two_by2'][t]))
  z=Z
  for q in [decode(o['positive_OLD_Gamma_cubic_two'][t]),decode(o['positive_OLD_D_second_two'][t]),decode(o['positive_FIRSTcarry_total_two'][t]),ng[t],nd[t],outer[t]]:z=ea(z,q)
  upper.append(z)
 six=[]
 for i in range(6):
  q=Z
  for t in range(2):q=ea(q,em(ex(dict(Q[t]).get(i,'0')),upper[t]))
  six.append(q)
 assembly_bits=imax;global_bits=max(339,assembly_bits,152,o['source_committed_bits']);assert global_bits<=384 and 3*global_bits+1<=1153
 enc=lambda q:[str(q[0]),str(q[1])]
 data={'selected_case619':554,'literal_source_all6':c0,'OLD_FIRST_component':o,'NATIVE_outer_component':n,'fixed_Q_embedding_two6':Q,'native_common_denominator_height':height,'positive_fine_native_Gamma_two':list(map(enc,ng)),'positive_fine_native_D_two':list(map(enc,nd)),'positive_outer_total_two':list(map(enc,outer)),'positive_WHOLE_upper_two':list(map(enc,upper)),'positive_WHOLE_full_six':list(map(enc,six)),'independent_normalized_positive_Rrows_all3264_checked':True,'all864_weightedS_point_and_higher_bound_cells_checked':True,'all68_native_component_fine_bounds_checked':True,'all_nu_i_nu_j_products_supportzero_BEFORE_coefficient_arithmetic':True,'whole_sign_recipe':'upper=OLD_Gamma+OLD_D+native_Gamma+native_D-FIRST57-outer644; full6=Q*upper','value_native_algorithm':'Direct C4/J18 monomial algebra. BuildweightedSc; assertfulltracepoint0. Useonlyknown9 partialNbasecells. Unsupportedbase exponents mustbeprunedBEFOREcoefficientproducts; νAνBcollision3dropsbysource-supportproof. ComputeN^tSc onlyasintermediateforactualΓcubic+Dnonzz2; unspecifiedN2mixed/N3areformal,andprojectedindependencefollowsS0=0/zfactor. VerifyFINE nativeΓ/D boundsbeforeWHOLEassembly.','native_coarse_intermediate_bits339':True,'whole_assembly_bits':assembly_bits,'global_value_committed_bits':global_bits,'global_value_temporary_bits':3*global_bits+1,'hard_guards':[384,1153],'all_actual_response_values_UNRUN':True,'fresh_native_frame_and_OLDpositive_independent_peer_required_before_GO':True}
 basecount=count(data);firstsupport=sum(len(v[2]) for v in o['positive_Xp_FIRST_columns28']);ndest=len({i//4 for ch in p['D_two_by10_by56'] for row in ch for i,q in enumerate(row) if q!='0'}|{i//4 for row in p['Gamma_two56'] for i,q in row if i<48})
 assert ndest<=12
 workspace={'native_knownN_required4rows_by_all12targets_9base4mask':4*12*9*4,'G8dense576_nu36_Gamma18_D1120':576+36+18+1120,'OLD_physical7x21':147,'OLD_program11x21_sparse_index_value':11*21*2,'native_Sc4targets18base4mask_sparse_index_value':4*18*4*2,'native_NSc_targets_by18base4mask_sparse_index_value':ndest*18*4*2,'point_FIRST_cache_possible_support_index_value':2*firstsupport+4*28,'OLD_selected_L0_Le_and_OTHER_parse_margin':2500,'two_monomial_index_tables_and_route_metadata':1200}
 derived=sum(workspace.values());success=2*(basecount+1000)+derived+4000;guard=2*(basecount+1000)+derived+12000+4000;cap=max(success,guard)
 rawcap=2_000_000+cap*(2*math.ceil(global_bits*math.log10(2))+42)
 assert cap<60000 and cap<150000 and rawcap<32*1024*1024
 data.update({'input_numeric_before_finalmeta':basecount,'metadata_reserve_per_input_copy':1000,'explicit_derived_workspace':workspace,'success_live_upper':success,'guard_live_upper':guard,'future_whole_numeric_cap':cap,'future_whole_RAW_cap':rawcap,'guard_prefix_rule':'SaveallNEWSc/NSc/OLDjets/FIRSTcaches/componentprefix+pendingoperands/uncommittedvalue/cursor; knownprimitiveinputsrepresentedbyimmutablefilepins, no duplicateddenseN or72cubetensor.'})
 ep=out/'ONE554-whole-positive-input-fine-assembly-and-native339-height.json';H.dump(ep,data)
 source_owned=count([n,o,p,nr,orr,height,Q,c0,N,G,nu,Gamma,D,expectedR,rows,expectedS,S,classes,rebuilt,data,pins])+15000;assert source_owned<150000
 report={'id':STEM,'all_passed':True,'whole_positive_only_no_signed_responses':True,'fresh_numeric68_Rrows_verified_all3264slots':True,'fresh_weightedS_all864_and68fine_components_verified':True,'native_common_den_classes_rebuilt':True,'positive_F_stage_bits':fmax,'whole_assembly_stage_bits':assembly_bits,'global_value_bits':global_bits,'global_temporary_bits':3*global_bits+1,'source_before_positive_upper':source_before_positive,'source_owned_upper':source_owned,'success_live_upper':success,'guard_live_upper':guard,'future_whole_numeric_cap':cap,'future_whole_RAW_cap':rawcap,'whole_SOURCE_RAW_upper':raw,'fresh_SOURCE_code_peer_and_immutable_VALUEcustody_GO_pending':True,'input_sha256':pins,'evidence_sha256':{ep.name:H.sha(ep)},'code_sha256':H.sha(Path(__file__)),'note_sha256':H.sha(ROOT/'research/notes'/f'{STEM}.md'),'seconds':time.monotonic()-start};H.dump(out/'report.json',report);print(json.dumps({k:report[k] for k in ['all_passed','positive_F_stage_bits','whole_assembly_stage_bits','global_value_bits','source_owned_upper','future_whole_numeric_cap','future_whole_RAW_cap','seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output)
