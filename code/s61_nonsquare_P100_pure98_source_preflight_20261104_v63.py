#!/usr/bin/env python3
"""Source/DAG/degree budgets for ONE q5H-cleared pure98 row. No core values."""
import argparse,gzip,itertools,json,math,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
from s61_nonsquare_oracle_frame_adapter_20261103_v63 import FrameAdapter
from s61_nonsquare_P100_H1_right_inverse_20261104_v63 import sha,read,native_matrix
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_P100_PURE98_SOURCE_PREFLIGHT_20261104_v63'
ZERO=(-1,(0,0,0,0,0))
def addbound(bs):
 bs=[b for b in bs if b[0]>=0]
 if not bs:return ZERO
 e=tuple(max(b[1][j]for b in bs)for j in range(5))
 return max(n+sum(x-y for x,y in zip(e,d))for n,d in bs),e
def mb(a,b):return ZERO if min(a[0],b[0])<0 else(a[0]+b[0],tuple(x+y for x,y in zip(a[1],b[1])))
def width(b):return 0 if b[0]<0 else b[0]+sum(b[1])+2
def bdict(b):return{'numerator_degree_upper':b[0],'known_linear_denominator_powers_upper':b[1],'expanded_numden_terms_upper':width(b)}
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(57);out.mkdir(parents=True,exist_ok=False)
 base=ROOT/'reports';ad=base/'S61_NONSQUARE_ORACLE_FRAME_ADAPTER_20261103_v63/adapter.json';oraclepath=base/'S61_NONSQUARE_LAZY_QUARTIC_ORACLE_20261103_v63/oracle.json';fdpath=base/'S61_NONSQUARE_JOINT_LELL_FREE164_20261104_v63/certificate.json';hp=base/'S61_NONSQUARE_H_SQUARE_H1_PREFLIGHT_20261104_v63/certificate.json';up=base/'S61_NONSQUARE_H1_ATOMIC_UNIT_DENOMINATORS_20261104_v63/certificate.json';rip=base/'S61_NONSQUARE_P100_H1_RIGHT_INVERSE_20261104_v63/native.json';cp=base/'S61_NONSQUARE_H1_KNOWN_UNIT_CONTINUATION_20261104_v63/checkpoint.json.gz';inv=base/'S61_NONSQUARE_JOINT_LELL_PREFLIGHT_20261104_v63/coordinate-inverse-numerators7x7.txt'
 paths=[Path(__file__),ROOT/'code/s61_nonsquare_oracle_frame_adapter_20261103_v63.py',ROOT/'code/s61_nonsquare_lazy_quartic_oracle_20261103_v63.py',ad,oraclepath,fdpath,hp,up,rip,cp,inv]
 f=FrameAdapter.load(ad);o=f.oracle;initial_nodes=len(o.nodes);initial_atoms=len(o.atoms)
 def forbidden(*a,**k):raise AssertionError('coefficient_value_evaluation_forbidden_in_source_preflight')
 o.core_value=forbidden;f.core_value=forbidden
 light=[0,1,2,3,4,6,7];mons=[tuple(p.count(i)for i in range(8))for p in itertools.combinations_with_replacement(light,3)];assert len(mons)==84
 h=read(hp);defs=h['all112_cubic_source319_definitions'];genuine=[1,6,8,9,11,14];rows=[]
 for i in genuine:
  actual=defs[7*i][2];assert f.frame['basis57_source_definitions'][41+i]['original_equation']==actual
  roots=[f.basis57_coefficient(41+i,m)for m in mons];rows.append({'genuine_index':i,'basis57_position':41+i,'actual135_row':actual,'roots84':roots})
  assert time.monotonic()-start<45
 roots=[v for r in rows for v in r['roots84']];assert len(roots)==504 and all(o.weights[v]==0 or v==0 for v in roots)
 needed=set();atoms=set();todo=roots[:]
 while todo:
  i=todo.pop()
  if i in needed:continue
  needed.add(i);nd=o.nodes[i]
  if nd[0]=='c':atoms.add(nd[1])
  elif nd[0]=='m':todo+=nd[1:]
  elif nd[0]=='l':
   for c,j in nd[1]:atoms.add(c);todo.append(j)
  else:assert nd[0]=='z'
 R=o.R;a,e,d=o.a,o.e,o.d;F=500*a+90*e**2+3*e-24;B=(50*a+18*e**2-9*e-2)*d+60*e-40;C=100*a+18*e**2-9*e-2;J=5*a*d-6*e-4;known=[a,d,F,B,C,J]
 def enc(p):return[[list(m),str(v)]for m,v in p.items()]
 def key(p):return tuple(sorted((tuple(m),str(Q(v)))for m,v in p))
 units=read(up);fd=read(fdpath);assert units['source_input_sha256']==sha(hp)and h['scalar_circuit']['atoms_num_den'][:len(fd['scalar_circuit']['atoms_num_den'])]==fd['scalar_circuit']['atoms_num_den']
 denproof={key(p):units['literal_constant_and_powers_by_denominator'][i]for i,p in enumerate(units['unique_atomic_denominators'])};factcache={};fallback=0
 def factor_den(den):
  nonlocal fallback
  k=key(den)
  if k in factcache:return factcache[k]
  if k in denproof:v=denproof[k]['a_d_F_B_C_J_powers']
  else:
   fallback+=1;p=R.from_dict({tuple(m):o.R.domain.convert(Q(c))for m,c in den});v=[]
   for q in known:
    n=0
    while not p.rem(q):p=p.exquo(q);n+=1;assert n<100
    v.append(n)
   assert p.is_ground and p
  factcache[k]=tuple(v);return tuple(v)
 def corebound(num,den):
  if not num:return ZERO
  aa,dd,ff,bb,cc,jj=factor_den(den);clear=max(2*m[0]+m[1]for m,c in num)
  nd=clear+max(2*m[2]-m[0]for m,c in num);zpower=2*aa-dd+ff-bb+cc-clear
  return nd+max(zpower,0),(max(-zpower,0),aa+dd,bb+cc,ff,jj)
 def bits(p):return max([max(abs(Q(c).numerator).bit_length(),Q(c).denominator.bit_length())for m,c in p]+[0])
 ob={};obstats=[]
 for ai in sorted(atoms):
  v=o.atoms[ai];n,q=enc(v.numer),enc(v.denom);ob[ai]=corebound(n,q);obstats.append({'atom':ai,'numterms':len(n),'denterms':len(q),'coefficient_bits':max(bits(n),bits(q)),'P100_bound':bdict(ob[ai])})
 nb={0:ZERO}
 for i in sorted(needed):
  nd=o.nodes[i]
  if nd[0]=='z':v=ZERO
  elif nd[0]=='c':v=ob[nd[1]]
  elif nd[0]=='m':v=mb(nb[nd[1]],nb[nd[2]])
  else:v=addbound([mb(ob[c],nb[j])for c,j in nd[1]])
  nb[i]=v
 pi_bound=addbound([nb[i]for i in roots])
 # Exact source closure of the saved pure98 linear normal form, no evaluation.
 purecols=h['target122_Hdegree_classes']['0'];NF=fd['full_target_projection239_to122'];Lnodes=fd['linear_nodes'];Lroots=[NF['output122_linear_roots'][j]for j in purecols];Lneed=set();Sroots=set();todo=Lroots[:]
 while todo:
  i=todo.pop()
  if i in Lneed:continue
  Lneed.add(i);nd=Lnodes[i]
  if nd[0]=='l':
   for c,j in nd[1]:Sroots.add(c);todo.append(j)
  else:assert nd[0]in['z','v']
 Snodes=fd['scalar_circuit']['nodes'];Satoms=fd['scalar_circuit']['atoms_num_den'];Sneed=set();SA=set();todo=list(Sroots)
 while todo:
  i=todo.pop()
  if i in Sneed:continue
  Sneed.add(i);nd=Snodes[i]
  if nd[0]=='c':SA.add(nd[1])
  elif nd[0]in['a','s','m']:todo+=nd[1:]
  else:assert nd[0]=='z'
 sb={0:ZERO};sab={ai:corebound(*Satoms[ai])for ai in SA}
 for i in sorted(Sneed):
  nd=Snodes[i];sb[i]=ZERO if nd[0]=='z'else sab[nd[1]]if nd[0]=='c'else mb(sb[nd[1]],sb[nd[2]])if nd[0]=='m'else addbound([sb[nd[1]],sb[nd[2]]])
 # Bound known already-computed right-inverse cofactors, not their new q5H sum.
 with gzip.open(cp,'rt')as zf:cache=json.load(zf)
 q5=next(row for i,row in cache['materialized_raw_rows']if i==0);del cache
 qb=[corebound(n,q)for n,q in q5];ri=read(rip)['right_inverse112x23'];linear=[(1,0),(1,20),(1,60),(203,2700),(1,12)];ucache={}
 def univden(poly):
  k=key(poly)
  if k in ucache:return ucache[k]
  p={m[0]:Q(c)for m,c in poly if Q(c)};es=[]
  for aa,bb in linear:
   n=0
   while p and max(p)>=1:
    q={};r=p.copy()
    while r and max(r)>=1:
     j=max(r);v=r[j]/aa;q[j-1]=q.get(j-1,Q())+v
     for k0,c0 in[(j,v*aa),(j-1,v*bb)]:
      r[k0]=r.get(k0,Q())-c0
      if not r[k0]:del r[k0]
    if r:break
    p=q;n+=1
   es.append(n)
  assert p and max(p)==0;ucache[k]=tuple(es);return tuple(es)
 rib=[]
 for row in ri:
  rr=[]
  for v in row:
   rr.append(ZERO if not v['numerator']else(max(m[0]for m,c in v['numerator']),univden(v['denominator'])))
  rib.append(rr)
 kappa=[addbound([mb(b,q)for b,q in zip(row,qb)])for row in rib]
 # Undo delta7 to the six ORIGINAL pure cubic tails.
 delta_b=corebound(enc(-64000*C**2),[[[0,0,0],'1']]);oldk=kappa[:]
 for j in range(7):oldk[6*7+j]=mb(delta_b,kappa[6*7+j]);oldk[1*7+j]=addbound([kappa[1*7+j],oldk[6*7+j]])
 active=[i for i,b in enumerate(oldk)if b[0]>=0];assert sorted({i//7 for i in active})==genuine and len(active)==42;k_bound=addbound([oldk[i]for i in active])
 # Common bound for the accepted old-physical -> joint linear inverse.
 _,_,n49=native_matrix(inv);jacden=enc(R.domain.convert(Q(3200,81))*a*C*d**2);bbounds=[corebound(p,jacden)for p in n49];Bphys=[]
 for i in range(7):
  for j in range(7):
   if i==5:Bphys.append(addbound([bbounds[5*7+j],bbounds[1*7+j],bbounds[0*7+j]]))
   else:Bphys.append(bbounds[[0,1,2,6,3,5,4][i]*7+j])
 Bbound=addbound(Bphys)
 psib=(pi_bound[0]+k_bound[0]+4*Bbound[0],tuple(pi_bound[1][j]+k_bound[1][j]+4*Bbound[1][j]for j in range(5)))
 lb={0:ZERO};inputset=set(NF['input_variables']);vidcol={v:i for i,v in enumerate(NF['input_variables'])};newmons=h['parent_joint_target239_monomials']
 for i in sorted(Lneed):
  nd=Lnodes[i]
  if nd[0]=='z':v=ZERO
  elif nd[0]=='v':assert nd[1]in inputset and newmons[vidcol[nd[1]]][7]==0;v=psib
  else:v=addbound([mb(sb[c],lb[j])for c,j in nd[1]])
  lb[i]=v
 outbound=sum(width(lb[i])for i in Lroots)
 # A live last-use budget for a native coefficient-first implementation.
 def live_budget(nodes,needed,roots,bounds,coeffs=False):
  uses={i:0 for i in needed}
  for i in needed:
   n=nodes[i]
   if n[0]=='m' or n[0]in['a','s']:
    for j in n[1:]:uses[j]+=1
   elif n[0]=='l':
    for c,j in n[1]:uses[j]+=1
  for i in roots:uses[i]+=1
  alive=0;peak=0
  for i in sorted(needed):
   alive+=width(bounds[i]);peak=max(peak,alive);n=nodes[i]
   ds=n[1:]if n[0]in['m','a','s']else[j for c,j in n[1]]if n[0]=='l'else[]
   for j in ds:
    uses[j]-=1
    if not uses[j]:alive-=width(bounds[j])
  return{'node_value_peak_terms_upper':peak,'retained_output_terms_upper':alive}
 op_live=live_budget(o.nodes,needed,set(roots),nb);nf_live=live_budget(Lnodes,Lneed,Lroots,lb)
 atomterms=sum(v['numterms']+v['denterms']for v in obstats);nf_atomterms=sum(len(Satoms[i][0])+len(Satoms[i][1])for i in SA)
 coefficient_record={'normalization':'FrameAdapter.basis57_coefficient(41+i,mon): F_actual(D*v) divided by the frozen row-unit, not multiplied.',
  'old_physical_monomials84':[list(m)for m in mons],'sources6':rows,'needed_nodes':[[i,o.nodes[i]]for i in sorted(needed)],
  'needed_atoms':[[ai,[enc(o.atoms[ai].numer),enc(o.atoms[ai].denom)]]for ai in sorted(atoms)],
  'needed_atom_profiles':obstats,'core_value_calls':0,'coefficient_values_evaluated':False}
 (out/'cubic-coefficient-input.json').write_text(json.dumps(coefficient_record,separators=(',',':'))+'\n')
 manifest={'one_curve':'P100','one_full_source_row':'q5H minus its actual cubic112 H1 rightinverse lift; originalpurepart q5H=0',
  'source42_cofactor_transport':'k_old6=delta*k_new6, k_old1=k_new1-delta*k_new6; allotherssame',
  'active_original_sources':genuine,'coefficient_roots504':True,'oracle_nodes_before_after':[initial_nodes,len(o.nodes)],'oracle_atoms_before_after':[initial_atoms,len(o.atoms)],
  'needed_oracle_nodes':len(needed),'needed_oracle_atoms':len(atoms),'needed_oracle_atom_terms':atomterms,'max_oracle_atom_coefficient_bits':max(x['coefficient_bits']for x in obstats),
  'needed_NF98_linear_nodes':len(Lneed),'needed_NF98_scalar_nodes':len(Sneed),'needed_NF98_scalar_atoms':len(SA),'needed_NF_atom_terms':nf_atomterms,
  'pure98_target122_columns':purecols,'pure98_parent239_indices':[h['parent_target122_monomial_indices'][j]for j in purecols],
  'NF_pure98_linear_roots':Lroots,'NF_needed_linear_ids':sorted(Lneed),'NF_needed_scalar_ids':sorted(Sneed),'NF_needed_atom_ids':sorted(SA),
  'known_curve_factor_order':['z','z+20','z+60','203z+2700','z+12'],
  'bounds':{'six_cubic_coeff_common':bdict(pi_bound),'42_clearing_cofactors_common':bdict(k_bound),'physical_to_joint_linear_common':bdict(Bbound),'joint_pure210_input_common':bdict(psib),
   'all98_expanded_output_terms_upper':outbound,'oracle_coefficient_live':op_live,'NF_linear_live':nf_live,
   'all504_coefficient_outputs_upper':sum(width(nb[i])for i in roots),'all_required_NF_scalar_values_if_retained_upper':sum(width(sb[i])for i in Sneed)},
  'source_only_no_actual_curve_specialization_no_new_field_values':True,
  'structured89_assessment':{'not_constructed_or_claimed_independent_in_this_packet':True,'potential_sectors':[{'kind':'ten_Hfree_cubics_times7','count':70},{'kind':'Koszul_of_six_linear_Hrows','count':15},{'kind':'five_quadH_intersection_with_evalkernel_using_q1_unit','count':4}],
   'conditional_module_argument':'If the actual invertible16-source split and6-row linear frame are recorded, multiplication kernel is15, quadraticNFkernel K5 has evaluation split by q1unit andintersection rank4; lifting4gives an exact split extension. Together70 gives89. Actualcoefficientvectors/sourcebasis identities still need construction/verification.',
   'potential_redundant_pure98_generators':{'pure_quad_rows':43,'cube_kernel_candidate':89,'q5H_clear':1,'quartic_H2free_H1clear':21,'total':154},'no_rank98_full239_claim':True},
  'next_job_not_launched':True}
 (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
 result={'id':ID,'all504_normalized_cubic_coefficient_roots_prepared':True,'no_coefficient_value_evaluation':True,
  'source_atoms_term_total':atomterms+nf_atomterms,'oracle_needed_nodes':len(needed),'NF98_needed_linear_scalar_nodes':[len(Lneed),len(Sneed)],
  'output98_terms_upper':outbound,'future_full98_output_guard_passes':outbound<=150000,'future_live_NF_bound':nf_live,
  'no_new_point_GB_native_coefficient_job_or_rank':True,'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},
  'cubic_input_sha256':sha(out/'cubic-coefficient-input.json'),'manifest_sha256':sha(out/'manifest.json'),'seconds':time.monotonic()-start,'frozen':True}
 assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024;signal.alarm(0);(out/'report.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2),flush=True)
if __name__=='__main__':
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True);main(ap.parse_args().output.resolve())
