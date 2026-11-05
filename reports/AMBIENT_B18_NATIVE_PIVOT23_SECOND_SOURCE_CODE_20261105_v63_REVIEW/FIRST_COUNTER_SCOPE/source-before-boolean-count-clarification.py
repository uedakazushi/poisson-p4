#!/usr/bin/env python3
"""Independent source/frame, ordinary-product and unrun native92 code audit."""
import ast,itertools,json,math,pathlib,resource,signal,time
from fractions import Fraction as Q
from ambient_B18_C3_alpha4_cached_source_review_20261105_v63 import field
from ambient_B18_C_genuine_projection7_source_review_20261105_v63 import count,sha
ROOT=pathlib.Path(__file__).resolve().parents[1]
PRI='AMBIENT_B18_ELL2_NATIVE_PIVOT23_SECOND_SOURCE_20261105_v63_V2'
CUS='AMBIENT_B18_NATIVE_PIVOT23_SECOND_VALUES_PREPARED_CUSTODY_20261105_v63'
ID='AMBIENT_B18_NATIVE_PIVOT23_SECOND_SOURCE_CODE_20261105_v63_REVIEW'
GATE_SHA='ab2e8d901a32da8e2f0e410f86e5bf46512b74fcc3342ecc12edbb52ccb316a0'
WORKER_SHA='910b6d8feabdeb5e3156ce4aa9c561221ef88fa332dd932c0f0577b051b0262d'
MON=[(0,0),(1,0),(0,1),(2,0),(1,1)];MI={m:i for i,m in enumerate(MON)}
def save(p,d):p.write_text(json.dumps(d,indent=2)+'\n')
def main():
 t=time.monotonic();signal.alarm(55);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30))
 out=ROOT/'reports'/ID;out.mkdir(exist_ok=False);pins={};raw=0
 def pin(p,h=None,parse=True):
  nonlocal raw
  p=pathlib.Path(p);k=str(p.relative_to(ROOT));s=sha(p)
  if h:assert s==h,(k,s,h)
  if k not in pins:raw+=p.stat().st_size
  pins[k]=s
  return json.loads(p.read_text()) if parse else None
 pr=pin(ROOT/'reports'/PRI/'report.json',GATE_SHA)
 for p,h in pr['input_sha256'].items():
  pin(ROOT/p,h,False)
  if p.endswith('/report.json'):assert json.loads((ROOT/p).read_text())['all_passed']
 en=next(iter(pr['evidence_sha256']));d=pin(ROOT/'reports'/PRI/en,pr['evidence_sha256'][en])
 custody=pin(ROOT/'reports'/CUS/'report.json')
 assert custody['all_passed'] and custody['gate_sha256']==GATE_SHA and custody['worker_sha256']==WORKER_SHA
 for p,h in custody['literal_input_sha256'].items():pin(ROOT/p,h,False)
 led=pin(ROOT/'reports'/CUS/'literal-input-and-whole-ownership.json',custody['ownership_ledger_sha256'])
 pin(ROOT/'code/ambient_B18_ell2_native_pivot23_second_source_v2_20261105_v63.py',pr['code_sha256'],False)
 worker=ROOT/'code/ambient_B18_native_pivot23_second_values_20261105_v63.py';pin(worker,WORKER_SHA,False)
 def src(s):return ROOT/next(p for p in pr['input_sha256'] if p.endswith(s))
 def fields(p,ks):return {k:field(p,k) for k in ks}
 op=fields(src('actual-complete-Res4-syzygy-operator.json'),['basis','native_S12_by32','actual_d2','chain_identity1','matrix_Res4','all48_actual_Qcolumns_independently_bound_to_mu_F'])
 chart=fields(src('actual27-plus99-unit-complements-and-source-combinations.json'),['original27_X_source_columns','free30_V_indices','remaining_free3_originalV_rows'])
 frame=fields(src('pruned-known-Cfull224-input-and-positive-source-envelopes.json'),['Gamma_two','pivot27_source_indices','pivot27_ell_recipes','P27_columns56','Pi3_columns56','native_pivot_ell_is_SUM_L0Se_LeS0','full_coordinate_not_pureU'])
 origin=fields(src('phase2-native-adapter.json'),['N','G8_C0','all_source_H12_to_G8_identities','d_not_inverted'])
 nr=field(src('phase2-actual-N-B11.json'),'N_C_B44')
 ne=fields(src('phase2-actual-T-N-G8-C16.json'),['N','actualG8','all_N0_rfirst_matches','all_G8_firstell_matches'])
 ca=fields(src('actual-G8-C68-and-monic-second-values.json'),['B17_basis','G8_order','full_G8_C68','monic_second10','all352_B11_and128_C16_overlap_coefficients_verified'])
 gamma=field(src('actual-two-rell-source-folds-and-retained-native-first-operators.json'),'two_Gamma_on56_actual_X2')
 assert op['basis']==['1','tauA','tauB','tauA*tauB'] and op['all48_actual_Qcolumns_independently_bound_to_mu_F']
 assert op['actual_d2']=={'shape':[48,128],'value':'Res4(N^t)*Res4(S12_by32)','exact_block_recipe':'block(i,j)=mu(sum_k N[k,i]*S[k,j]); native12rowblock then4basis columns'}
 assert op['chain_identity1']=='Q*Res4(N^t)*Res4(S)=Res4((dG,0)^t_as_row*S)=0'
 assert 'never transpose within mu' in op['matrix_Res4']
 piv=[604,613,622,631]+[640+4*i for i in range(7)]+list(range(752,768))
 assert piv==chart['original27_X_source_columns']==frame['pivot27_source_indices']==d['original27_X_source_columns']
 native=[j-640 for j in piv[4:]];assert native==[0,4,8,12,16,20,24]+list(range(112,128))
 assert frame['pivot27_ell_recipes'][4:]==[['native',i] for i in native]
 assert frame['Gamma_two']==gamma==d['Gamma_two'] and any(i==55 for row in gamma for i,q in row)
 assert d['P27_columns56']==frame['P27_columns56'] and d['Pi3_columns56']==frame['Pi3_columns56']
 assert frame['native_pivot_ell_is_SUM_L0Se_LeS0'] and frame['full_coordinate_not_pureU']
 assert origin['all_source_H12_to_G8_identities'] and origin['d_not_inverted'] and ne['all_N0_rfirst_matches'] and ne['all_G8_firstell_matches']
 N=[]
 for i in range(12):
  row=[]
  for j in range(12):
   a=nr[i][j];b=ne['N'][i][j];assert len(a)==44 and len(b)==16 and a[:8]==b[:8] and a[:4]==origin['N'][i][j]
   row.append([a[:4],a[4:8],b[8:12],a[8:12],b[12:16]])
  N.append(row)
 assert N==d['cached_N_12_by12_plane5_free4_UNEVALUATED_IN_PIVOT_RESPONSE']
 pos=[0,5,6,15,16];basis=ca['B17_basis']
 assert [basis[k][-2:] for k in pos]==[[0,0],[1,0],[0,1],[2,0],[1,1]] and all(not any(basis[k][:7]) for k in pos)
 G=[[[g[4*k+i] for i in range(4)] for k in pos] for g in ca['full_G8_C68']]
 assert G==d['canonical_G8_plane5_free4_UNEVALUATED_IN_PIVOT_RESPONSE'] and ca['G8_order']==d['canonical_G8_order']
 assert ca['all352_B11_and128_C16_overlap_coefficients_verified']
 for i in range(8):assert G[i][0]==origin['G8_C0'][i] and sum(G[i][:3],[])==ne['actualG8'][i][:12]
 nu=[['0','0','0']+a[8:10] for a in ca['monic_second10']]
 assert nu==d['canonical_nu_two_plane5']==[['0','0','0','-25/192','0'],['0','0','0','725/36','41/9']]
 # All28 Koszul columns, including signs, are reconstructed independently.
 sparse=[];by={j:[] for j in range(32)}
 for j,(i,k) in enumerate(itertools.combinations(range(8),2)):
  sparse += [[i,j,k,1],[k,j,i,-1]];by[j]=[['G8',i,k,1],['G8',k,i,-1]]
 assert sparse==op['native_S12_by32']['koszul28_in_first8_rows']
 assert op['native_S12_by32']['free4_zero_source_columns']==[[i,20+i,1] for i in range(8,12)]
 for i in range(8,12):by[20+i]=[['constant',i,1]]
 for k,(j,r) in enumerate(zip(native,d['native23_pivot_recipes'])):
  assert (r['pivot27_position'],r['original768_source'],r['native128_source'],r['S32_source'],r['input_trace_mask'])==(k+4,j+640,j,j//4,j%4)
  assert r['S_source_entries']==by[j//4] and r['individual_S0_is_NOT_assumed_zero']
  assert r['target_Gamma_Xp_rr']==['Gamma','L0*Srr','Lr*Sr','Lrr*S0']
  assert r['target_Gamma_Xp_rell']==['Gamma','L0*Sre','Lr*Se','Le*Sr','Lre*S0']
 trace=[]
 for a,c,b in itertools.product(range(4),repeat=3):
  pa=sum((x&1) for x in [a,c,b])//2;pb=sum((x>>1) for x in [a,c,b])//2
  row=[a,c,b,a^c^b,pa,pb];trace.append(row)
  # Both associations of the generic tau monomials give exactly these powers.
  for x,y,z in [(a,c,b),(c,b,a)]:
   assert (int(bool(x&y&1))+int(bool((x^y)&z&1)),int(bool(x&y&2))+int(bool((x^y)&z&2)))==(pa,pb)
 assert trace==d['universal_trace_product_a_c_inputmask_output_nuA_nuB']
 assert all(pa<=1 and pb<=1 for a,c,b,o,pa,pb in trace)
 def flatten(v):
  if isinstance(v,(list,tuple)):
   for x in v:yield from flatten(x)
  else:yield Q(v)
 def du(v):
  a=list(flatten(v));den=math.lcm(*(x.denominator for x in a));u=max([den]+[abs(x.numerator)*(den//x.denominator) for x in a])
  return {'D':str(den),'U':str(u),'D_bits':den.bit_length(),'U_bits':u.bit_length(),'entries':len(a)}
 classes={'Gamma':du([q for row in gamma for i,q in row]),'N_plane5':du(N),'G8_plane5_with_one':du([G,'1']),'nu_plane5_with_one':du([nu,'1'])}
 assert classes==d['source_positive_classes_D_U']
 den=math.prod(int(v['D']) for v in classes.values());U=math.prod(int(v['U']) for v in classes.values())*(1<<20)
 assert max(den.bit_length(),U.bit_length())==280 and 3*280+1==841
 # Per individual coefficient there are at most two S terms,12 physical rows,
 # 4x4 tau masks,5x5 ordinary bases and two nu2 directions. The row selector
 # selects one output mask, so its support is NOT another multiplicity factor.
 path_upper=2*12*4*4*5*5*2;assert path_upper<1<<20
 def routes(base,mask,im):
  collision=mask&im;outmask=mask^im
  if not collision:return [(base,outmask,1)]
  if collision==3:return []
  ans=[]
  for b in [3,4]:
   e=tuple(MON[base][j]+MON[b][j] for j in range(2))
   if e in MI and Q(nu[collision-1][b]):ans.append((MI[e],outmask,1))
  return ans
 As=[[set() for _ in MON] for _ in range(2)]
 for i,im,physical,base,mask in itertools.product(range(12),range(4),range(12),range(5),range(4)):
  if not Q(N[i][physical][base][mask]):continue
  for b,om,c in routes(base,mask,im):
   for t0,row in enumerate(gamma):
    if any(j==4*physical+om and Q(q) for j,q in row):As[t0][b].add(4*i+im)
 Acharge=2*sum(len(s) for row in As for s in row);assert Acharge==312
 maxS=0
 for recipe in d['native23_pivot_recipes']:
  ss=[set() for _ in MON];im=recipe['input_trace_mask']
  for entry in recipe['S_source_entries']:
   row=entry[1];f=G[entry[2]] if entry[0]=='G8' else [['1','0','0','0']]+[['0']*4 for _ in range(4)]
   for base,mask in itertools.product(range(5),range(4)):
    if not Q(f[base][mask]):continue
    for b,om,c in routes(base,mask,im):ss[b].add(4*row+om)
  maxS=max(maxS,2*sum(map(len,ss)))
 assert maxS==36 and count(d)==3883 and list(map(count,[N,G,nu]))==[2880,160,10]
 assert led['parsed_N_G_nu_Gamma']==[2880,160,10,32] and led['Gamma_L_support_numeric_upper']==Acharge and led['current_intrinsic_S_numeric_upper']==maxS
 success=2*count(d)+2880+160+10+32+Acharge+maxS+1000+2*(92+322)+1000
 assert success==led['explicit_success_upper']==14024
 assert success+Acharge+maxS==led['explicit_guard_upper']==14372<d['prospective_value_whole_numeric_bound']==21034
 # Fresh literal AST/code audit. No main or actual input products are executed.
 text=worker.read_text();tree=ast.parse(text)
 assert {a.name for n in tree.body if isinstance(n,ast.Import) for a in n.names}<={'argparse','hashlib','json','math','resource','signal','time'}
 assert {n.module for n in tree.body if isinstance(n,ast.ImportFrom)}=={'fractions','pathlib'}
 mainnode=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='main')
 funcs={n.name:n for n in mainnode.body if isinstance(n,ast.FunctionDef)}
 routefn=funcs['trace_routes'];fixture_nu=[[Q(0),Q(0),Q(0),Q(2),Q(3)],[Q(0),Q(0),Q(0),Q(5),Q(7)]]
 env={'F':Q,'nu':fixture_nu,'mon':MON,'idx':MI};exec(compile(ast.Module(body=[routefn],type_ignores=[]),str(worker),'exec'),env)
 for base,mask,im in itertools.product(range(5),range(4),range(4)):
  coll=mask&im;want=[]
  if coll==0:want=[(base,mask^im,Q(1))]
  elif coll!=3:
   for b in [3,4]:
    ex=tuple(MON[base][j]+MON[b][j] for j in range(2))
    if ex in MI:want.append((MI[ex],mask^im,fixture_nu[coll-1][b]))
  assert env['trace_routes'](base,mask,im)==want
 required=[
  "Gamma=[{i:ck(F(q)) for i,q in row if i<48 and F(q)} for row in d['Gamma_two']]",
  'v in enumerate(N[source][physical][base])','row=4*physical+om','madd(A[t][b],col,mul(mul(Gamma[t][row],v),carry))',
  'rr=[dot(A[t][0],S[3]),dot(A[t][1],S[1]),dot(A[t][3],S[0])]',
  're=[dot(A[t][0],S[4]),dot(A[t][1],S[2]),dot(A[t][2],S[1]),dot(A[t][4],S[0])]',
  'for v in rr:vr=add(vr,v)','for v in re:ve=add(ve,v)',
  "assert owned(state)<cap;dump('guard-frontier.json',state)",
  "'completed_target_prefix':targets,'completed_piece_prefix':pieces",'owner=owned(ev);assert owner<cap',
  "c['worker_sha256']==sha(Path(__file__)) and c['gate_sha256']==GATE_SHA",
  "for path,h in c['literal_input_sha256'].items():assert sha(ROOT/path)==h"]
 for s in required:assert s in text,s
 assert 'def add(a,b):return ck(a+b)' in text and 'def mul(a,b):return ck(a*b)' in text
 # Committed280 gives rational add/mul temporary <=561. The source's841 also
 # covers unfused triple products and the final integer envelope comparisons.
 assert max(2*280+1,3*280+1,den.bit_length()+280,U.bit_length()+280)<=841<1153
 assert len(d['native23_pivot_recipes'])*2*2==92 and len(d['native23_pivot_recipes'])*2*(3+4)==322
 owned=count([pr,custody,led,d,op,chart,frame,origin,nr,ne,ca,gamma,N,G,nu,classes,pins])+25000
 assert owned<150000 and raw+3000000<32<<20
 assert not (ROOT/'reports/AMBIENT_B18_NATIVE_PIVOT23_SECOND_VALUES_20261105_v63').exists()
 assert all(sha(ROOT/p)==h for p,h in pins.items())
 proof={'literal_original27_native23_and_all_S32_signs_verified':True,'full_Gamma56_preserved_native_bottom8_zero_only':True,
  'same_original_Nrr_Nre_canonical_C68_G8_and_monic_frame_verified':True,'normalizer_order':'N[source12][physical12]; transpose12only; never transpose mu4',
  'all64_formal_tau_triple_associations_and80_unrelated_worker_route_fixtures':True,
  'ordinary_rr_three_and_re_four_pieces_L2S0_kept':True,
  'free16_individual_pivots_have_constant_S_and_nonvanishing_S0_not_zero_filled':True,
  'classes_D_U_independently_rebuilt':classes,'per_output_path_upper':path_upper,'source_path_reserve':1<<20,
  'committed_prior280_temporary841_hard384_1153_unchanged':True,'ordinary_rational_before_commit_temporary561':True,
  'normal_owned_upper':success,'failure_owned_upper':success+Acharge+maxS,'hard_owner':21034,
  'worker_only92_values_and322_pieces_UNRUN':True,'no_signed_native_pivot_products_evaluated':True,
  'SOURCE_review_owned':owned,'RAW_input_bytes':raw,'CODE_STANDARD_LIBRARY_ONLY':True,
  'no_full_response_honest_upper_HF_length_or_original_H0_claim':True,
  'scope':'Source/frame/operator/code/custody readiness for the same ONE bounded native23 second92 job only. No value job executed by reviewer.'}
 ep=out/'independent-native23-source-order-trace-product-and-unrun-worker-custody.json';save(ep,proof)
 report={'id':ID,'all_passed':True,'frozen':True,'SOURCE_CODE_ONLY':True,'input_sha256':pins,'source_gate_sha256':GATE_SHA,'worker_sha256':WORKER_SHA,'code_sha256':sha(__file__),'evidence_sha256':{ep.name:sha(ep)},'worker_values_UNRUN':True,'seconds':time.monotonic()-t}
 save(out/'report.json',report);signal.alarm(0);print(json.dumps({'all_passed':True,'sha256':sha(out/'report.json'),'SOURCE_owned':owned,'success_owned':success,'failure_owned':success+Acharge+maxS,'seconds':report['seconds']}))
if __name__=='__main__':main()
