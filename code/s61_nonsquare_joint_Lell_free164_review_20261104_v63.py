#!/usr/bin/env python3
"""Independent inverse-transpose117 restoration and free164 circuit proof audit."""
import argparse,json,resource,signal,subprocess,time
from pathlib import Path
from s61_nonsquare_free281_circuit_review_20261104_v63 import ProofExpressions,polynomial,sha,read
from s61_nonsquare_joint_Lell_free183_review_20261104_v63 import LinearProof
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_JOINT_LELL_FREE164_20261104_v63_REVIEW'
def main(out):
 began=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(parents=True,exist_ok=False)
 b=ROOT/'reports/S61_NONSQUARE_JOINT_LELL_FREE164_20261104_v63';r=read(b/'report.json');c=read(b/'certificate.json')
 assert r['frozen']and r['actual_free_source_target_shape']==[164,122]and sha(b/'certificate.json')==r['certificate_sha256']
 for p,h in r['input_sha256'].items():assert sha(ROOT/p)==h,p
 jb=ROOT/'reports/S61_NONSQUARE_JOINT_LELL_FREE183_20261104_v63';jd=read(jb/'certificate.json');assert c['initial_joint183_certificate_sha256']==sha(jb/'certificate.json')
 ob=ROOT/'reports/S61_NONSQUARE_JOINT_LELL_OVERLAP_UNIT_20261104_v63';od=read(ob/'overlap-unit-certificate.json');pair=read(ob/'remaining-pair-block.json');qb=ROOT/'reports/S61_NONSQUARE_JOINT_LELL_QPAIR_UNITS_20261104_v63';qs=read(qb/'checkpoint.json')
 sg=c['scalar_circuit'];ln=c['linear_nodes'];assert len(sg['nodes'])==22244 and len(sg['atoms_num_den'])==1140 and len(ln)==4616
 assert sg['nodes'][:len(jd['scalar_circuit']['nodes'])]==jd['scalar_circuit']['nodes']and sg['atoms_num_den'][:len(jd['scalar_circuit']['atoms_num_den'])]==jd['scalar_circuit']['atoms_num_den']and ln[:len(jd['linear_nodes'])]==jd['linear_nodes']
 sc=ProofExpressions();sm=[]
 for i,nd in enumerate(sg['nodes']):
  if nd[0]=='z':assert i==0;v=0
  elif nd[0]=='c':v=sc.atom(*sg['atoms_num_den'][nd[1]])
  else:assert nd[0]in['a','s','m']and all(0<=j<i for j in nd[1:]);v={'a':sc.plus,'s':sc.minus,'m':sc.times}[nd[0]](sm[nd[1]],sm[nd[2]])
  sm.append(v)
 lin=LinearProof(sc);lm=[];allowed=set(range(685))|set(range(1000,1057))|set(range(2000,2141))|set(range(2200,2336))
 for i,nd in enumerate(ln):
  if nd[0]=='z':assert i==0;v=0
  elif nd[0]=='v':assert nd[1]in allowed;v=lin.variable(nd[1])
  else:assert nd[0]=='l'and all(0<=k<len(sm)and 0<=j<i for k,j in nd[1]);v=lin.combine((sm[k],lm[j])for k,j in nd[1])
  lm.append(v)
 def S(key,expected):assert[sm[j]for j in c[key]]==expected,key
 def L(key,expected):assert[lm[j]for j in c[key]]==expected,key
 def rat(v):return sc.atom(*v)
 def inverse(v):return sc.atom(v[1],v[0])
 def compose(roots,values):
  cache={0:0}
  def rec(i):
   if i in cache:return cache[i]
   nd=lin.data[i]
   if nd[0]=='variable':v=values[nd[1]]
   else:v=lin.combine((co,rec(j))for co,j in nd[1])
   cache[i]=v;return v
  return[rec(i)for i in roots]
 oldE=jd['retained76_quadratic_source_rows'];overlap=od['source_original_overlap_rows'];qrows=[{j:sm[k]for j,k in row}for row in jd['quadric_product_coefficient_sparse_rows174']]
 newm2=[m for _,m in jd['new_prime_quadric_source174_definitions'][:29]];mi={tuple(m):j for j,m in enumerate(newm2)};plabels=[i*29+mi[tuple(m)]for i,m in pair['source_rows']]
 assert oldE[:18]==overlap and oldE[18:]==plabels
 ost=od['unit_pivots'];qst=qs['all_completed_unit_steps'];sp=[v['pivot_original_overlap_row']for v in ost];sq=[18+v['pivot_row']for v in qst];newE=[i for i in range(76)if i not in sp+sq]
 assert len(c['reverse14_steps'])==14 and len(c['reverse5_steps'])==5
 assert len(newE)==57 and len(set(sp+sq))==19 and c['additional_pivots_in_old76_source_coordinates']==sp+sq and c['remaining57_in_old76_source_coordinates']==newE
 var57=list(range(1000,1057));assert c['restoration57_input_variables']==var57;coeff=[0]*76
 for i,j in enumerate(newE):coeff[j]=lin.variable(var57[i])
 for step,record in zip(reversed(qst),c['reverse14_steps']):
  i=18+step['pivot_row'];ip=inverse(step['pivot']);assert record['source76_pivot']==i and sm[record['inverse_pivot_root']]==ip
  assert [[t,sm[j]]for t,j in record['elimination_coefficient_roots']]==[[18+t,rat(v)]for t,v in step['eliminations']]
  terms=[(sc.one,coeff[i])]+[(sc.minus(0,rat(v)),coeff[18+t])for t,v in step['eliminations']];coeff[i]=lin.combine([(ip,lin.combine(terms))])
 L('coefficient76_after_reverse14_linear_roots',coeff)
 nf=jd['target_projection_operator'];nfroots=[lm[j]for j in nf['output141_linear_roots']];tvars=nf['input_variables'];orows=[[rat(v)for v in od['after_rows'][step['pivot_original_overlap_row']]]for step in ost];ocols=[step['pivot_target141_column']for step in ost]
 # Reconstruct all58x5 missing source coefficients from the INDEPENDENTLY
 # accepted full NF98 operator and actual174 quadric products.
 paircofs=[]
 for source in plabels:
  vals=lin.substitute([nfroots[j]for j in ocols],{tvars[j]:qrows[source].get(j,0)for j in range(239)});cs=[]
  for k,v in enumerate(vals):
   for i in range(k):v=sc.minus(v,sc.times(cs[i],orows[i][ocols[k]]))
   cs.append(v)
  paircofs.append(cs)
 assert [[sm[j]for j in row]for row in c['qpair58_at_overlap5_pivot_coefficient_scalar_roots']]==paircofs
 for k,record in zip(range(4,-1,-1),c['reverse5_steps']):
  step=ost[k];i=step['pivot_original_overlap_row'];terms=[(sc.one,coeff[i])];entries=[]
  for t,v in step['eliminations']:co=rat(v);terms.append((sc.minus(0,co),coeff[t]));entries.append([t,co])
  for t in range(58):
   co=paircofs[t][k]
   if co:terms.append((sc.minus(0,co),coeff[18+t]));entries.append([18+t,co])
  ip=inverse(step['pivot_value']);assert record['source76_pivot']==i and sm[record['inverse_pivot_root']]==ip and [[t,sm[j]]for t,j in record['elimination_coefficient_roots']]==entries
  coeff[i]=lin.combine([(ip,lin.combine(terms))])
 L('restored76_after_reverse19_linear_roots',coeff)
 oldrest=jd['quad_restoration_operator'];values=dict(zip(oldrest['input_variables'],coeff))
 L('restored_prime174_linear_roots',compose([lm[j]for j in oldrest['prime174_linear_roots']],values));L('restored_old174_linear_roots',compose([lm[j]for j in oldrest['old174_linear_roots']],values))
 P57=compose([lm[j]for j in jd['descended_P22_operator']['output22_linear_roots']],values);L('P57_output22_linear_roots',P57)
 k76=jd['descended_K22_product_operator'];K57=[lm[k76['output76_linear_roots'][i]]for i in newE];L('K57_output57_linear_roots',K57);assert c['K57_input22_variables']==k76['input_variables']
 x57=[sm[jd['new_quad76_kernel_after_old15_roots'][i]]for i in newE];u=[sm[j]for j in jd['new_B22_kernel_after_old15_roots']];h=[sm[j]for j in jd['new_C108_kernel_after_old15_roots']]
 S('actual_new_kernel_quad57_roots',x57);S('actual_new_kernel_B22_roots',u);S('actual_new_kernel_C108_roots',h)
 Px=lin.substitute(P57,dict(zip(var57,x57)));S('P57_on_actualx_roots',Px);kk=[sc.minus(Px[j],u[j])for j in range(22)];S('kernel22_coordinate_roots',kk)
 Kk=lin.substitute(K57,dict(zip(k76['input_variables'],kk)));S('K57_times_actualk_roots',Kk);v165=[sc.minus(x57[i],Kk[i])for i in range(57)]+h;S('actual_free165_kernel_roots',v165)
 assert c['Bnew22_original319_labels']==jd['Bnew22_original319_labels']and c['Cnew108_original319_labels']==jd['Cnew108_original319_labels']
 B=c['Bnew22_original319_labels'];C=c['Cnew108_original319_labels'];assert all(i>=174 for i in B+C)and not set(B)&set(C)and all(i not in B for i in[182,186])
 p,q=[57+C.index(i)for i in[182,186]];assert[p,q]==c['separator165_coordinates']==[57,58];theta=sm[c['known_unit_separator_root']]
 assert c['known_unit_separator_root']==jd['separator_known_unit_root']and c['separator_identity_zero_root']==jd['separator_zero_check_root']and sc.minus(sc.plus(v165[p],v165[q]),theta)==sm[c['separator_identity_zero_root']]
 changed=v165[:];changed[p]=theta;S('changed165_kernel_roots',changed)
 # The reciprocal is the exact pinned free183 unit atom, so no new factor inversion.
 norm183=jd['normalized_unit_kernel_roots'];first=next(j for j in range(len(norm183))if jd['changed184_kernel_roots'][j]and j!=76);oldnode=jd['scalar_circuit']['nodes'][norm183[first]]
 assert oldnode[0]=='m';inverseTheta=sm[oldnode[2]]
 normalized=[sc.times(v,inverseTheta)for v in changed];normalized[p]=sc.one;S('unit_normalized_kernel_roots',normalized)
 assert c['actual_free164_shape']==[164,122]and c['retained_source165_indices_after_constant_pair_change']==[i for i in range(165)if i!=p]
 # Exact target projection through the same19 unit image pivots.
 vars141=list(range(2000,2141));vars136=list(range(2200,2336));a141=[lin.variable(i)for i in vars141]
 for step,row in zip(ost,orows):
  j=step['pivot_target141_column'];v=a141[j];assert row[j]==sc.one
  for k,co in enumerate(row):
   if k!=j and co:a141[k]=lin.combine([(sc.one,a141[k]),(sc.minus(0,co),v)])
  a141[j]=0
 keep136=pair['target141_remaining_columns'];pr136=[a141[j]for j in keep136];assert c['target_projection141_to136']['input_variables']==vars141 and [lm[j]for j in c['target_projection141_to136']['output136_linear_roots']]==pr136
 a136=[lin.variable(i)for i in vars136];used=[]
 for step in qst:
  row=[rat(v)for v in qs['matrix58x136'][step['pivot_row']]];j=step['pivot_column'];v=a136[j];assert row[j]==sc.one and all(row[k]==0 for k in used)
  for k,co in enumerate(row):
   if k!=j and co:a136[k]=lin.combine([(sc.one,a136[k]),(sc.minus(0,co),v)])
  a136[j]=0;used.append(j)
 keep122=qs['remaining_columns'];pr122=[a136[j]for j in keep122];assert c['target_projection136_to122']['input_variables']==vars136 and [lm[j]for j in c['target_projection136_to122']['output122_linear_roots']]==pr122
 firstproj=compose(pr136,dict(zip(vars141,nfroots)));fullNF=compose(pr122,dict(zip(vars136,firstproj)));assert c['full_target_projection239_to122']['input_variables']==tvars and [lm[j]for j in c['full_target_projection239_to122']['output122_linear_roots']]==fullNF
 targetD=jd['retained141_target_columns'];source117=jd['unit98_source_rows']+[oldE[i]for i in sp+sq];target117=jd['unit98_target_columns']+[targetD[j]for j in ocols]+[targetD[keep136[s['pivot_column']]]for s in qst]
 source57=[oldE[i]for i in newE];target122=[targetD[keep136[j]]for j in keep122]
 assert c['unit117_original_prime_quad_source_rows']==source117 and c['unit117_original_joint_target_columns']==target117 and len(set(source117))==len(set(target117))==117
 assert c['remaining57_original_prime_quad_source_rows']==source57 and c['remaining122_original_joint_target_columns']==target122 and sorted(source117+source57)==list(range(174))and sorted(target117+target122)==list(range(239))
 expecteddefs=[['prime_quadric_product',i,jd['new_prime_quadric_source174_definitions'][i]]for i in source57]+[['old_source319',i,jd['old_source319_definitions'][i]]for i in range(174,319)];assert c['source202_definitions']==expecteddefs
 with(out/'scalar-proof-input.txt').open('w')as f:
  f.write(f'{len(sg["atoms_num_den"])} {len(sg["nodes"])} {c["separator_identity_zero_root"]} {c["known_unit_separator_root"]}\n')
  for num,den in sg['atoms_num_den']:
   for terms in[num,den]:
    po=polynomial(terms);f.write(str(len(po))+''.join(' '+str(co)+' '+' '.join(map(str,m))for m,co in po)+'\n')
  for nd in sg['nodes']:
   val=[0,0,0]if nd[0]=='z'else[1,nd[1],0]if nd[0]=='c'else[{'a':2,'s':3,'m':4}[nd[0]],nd[1],nd[2]];f.write(' '.join(map(str,val))+'\n')
 cpp=ROOT/'code/s61_nonsquare_joint_Lell_free164_review_20261104_v63.cpp';exe=out/'review-native';flags=subprocess.check_output(['pkg-config','--cflags','--libs','flint'],text=True).split();rem=lambda:max(.1,60-(time.monotonic()-began))
 with(out/'compile.log').open('w')as f:subprocess.run(['c++','-O2','-std=c++17',str(cpp),*flags,'-o',str(exe)],check=True,stdout=f,stderr=subprocess.STDOUT,timeout=rem())
 run=subprocess.run([str(exe),str(out/'scalar-proof-input.txt')],capture_output=True,text=True,timeout=rem());exe.unlink();(out/'native.log').write_text(run.stdout+run.stderr);assert run.returncode==0,run.stderr;native=json.loads(run.stdout);assert native['all_passed']
 deps=[ROOT/'reports'/s/'report.json'for s in['S61_NONSQUARE_JOINT_LELL_FREE183_20261104_v63_REVIEW','S61_NONSQUARE_JOINT_LELL_OVERLAP_UNIT_20261104_v63_REVIEW','S61_NONSQUARE_JOINT_LELL_QPAIR_UNITS_20261104_v63_REVIEW']];assert all(read(p)['all_passed']for p in deps)
 note=ROOT/'research/notes'/(ID+'.md');utils=[ROOT/'code'/n for n in['s61_nonsquare_free281_circuit_review_20261104_v63.py','s61_nonsquare_joint_Lell_free183_review_20261104_v63.py','s61_nonsquare_ell_Lframe_review_20261104_v63.cpp']]
 result=dict(id=ID,all_passed=True,all117_source_target_pivots_and57_122_complement_indices_checked=True,
  all14_plus5_inverse_transpose_coefficient_steps_and_missing58x5_cofs_checked=True,
  restored_prime174_old174_P57_K57_and_all_actual165_roots_match=True,
  pure_quad_kernel117_restoration_proves_P57K57_without_rank_flatness_assumption=True,
  only_pivot_generators_scaled_so_never_pivot_coefficients_unchanged=True,
  old15_nonquad145_and_Bnew22_C108_separator_preserved=True,
  all22244_scalar4616_linear_nodes_topological_and1140_denominators_unit=True,
  exact_actual319_cokernel_equals_global_free164_to122_accepted=True,
  no_expanded164matrix_rank122_full239_newpoint_GB_solver_or_author_rerun=True,
  known38_missing42_unchanged=True,native=native,
  input_sha256={str(p.relative_to(ROOT)):sha(p)for p in[Path(__file__),cpp,note,b/'report.json',b/'certificate.json',jb/'certificate.json',ob/'overlap-unit-certificate.json',ob/'remaining-pair-block.json',qb/'checkpoint.json',*utils,*deps]},
  evidence_sha256={p.name:sha(p)for p in out.iterdir()if p.is_file()},seconds=time.monotonic()-began,frozen=True)
 (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');signal.alarm(0);print(json.dumps({k:v for k,v in result.items()if k not in['input_sha256','evidence_sha256']},indent=2),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);main(p.parse_args().output.resolve())
