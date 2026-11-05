#!/usr/bin/env python3
"""Independent literal IAIB and integer action-based 660 NF descent audit."""
import hashlib,itertools,json,math,pathlib,resource,signal,time
from fractions import Fraction as Q
from ambient_B18_C3_alpha4_cached_source_review_20261105_v63 import field
from ambient_B18_C_genuine_projection7_source_review_20261105_v63 import count,sha
ROOT=pathlib.Path(__file__).resolve().parents[1]
PRI='AMBIENT_B18_ELL2_NATIVE_NONZZ10_DESCENT_SOURCE_20261105_v63'
HEIGHT='AMBIENT_B18_ELL2_NATIVE_NONZZ10_DESCENT_HEIGHT_20261105_v63'
ID=PRI+'_REVIEW';ZERO=(0,)*9
def save(p,d):p.write_text(json.dumps(d,indent=2)+'\n')
def sp(v):return {int(i):Q(x) for i,x in v if Q(x)}
def fastfield(p,k):
 text=p.read_text();m=list(__import__('re').finditer(r'"'+__import__('re').escape(k)+r'"\s*:',text));assert len(m)==1
 pos=m[0].end()
 while text[pos].isspace():pos+=1
 return json.JSONDecoder().raw_decode(text,pos)[0]
def exp(slots):
 m=[0]*9
 for i in slots:assert 0<=i<9;m[i]+=1
 return tuple(m)
def main():
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30))
 out=ROOT/'reports'/ID;out.mkdir(exist_ok=False);pins={};RAW=0
 def pin(p,h=None,parse=True):
  nonlocal RAW
  p=pathlib.Path(p);s=sha(p)
  if h:assert s==h,(str(p),s,h)
  name=str(p.relative_to(ROOT))
  if name not in pins:RAW+=p.stat().st_size
  pins[name]=s
  return json.loads(p.read_text()) if parse else None
 pr=pin(ROOT/'reports'/PRI/'report.json','d490c810468a9b5fc3011e762406e7342fbc8dbdc1978f9837e235664eb00c3b')
 hr=pin(ROOT/'reports'/HEIGHT/'report.json','fe9a9594f49cc6d97ddf0ebde9da5e954bea9c77ea4922a604702be0804d257a')
 for report in [pr,hr]:
  for p,h in report['input_sha256'].items():
   pin(ROOT/p,h,False)
   if p.endswith('_REVIEW/report.json'):assert json.loads((ROOT/p).read_text())['all_passed']
 pin(ROOT/'code/ambient_B18_ell2_native_nonzz10_descent_source_20261105_v63.py',pr['code_sha256'],False)
 pin(ROOT/'code/ambient_B18_ell2_native_nonzz10_descent_height_20261105_v63.py',hr['code_sha256'],False)
 pin(ROOT/'research/notes'/f'{PRI}.md',pr['note_sha256'],False)
 en=next(iter(pr['evidence_sha256']));D=pin(ROOT/'reports'/PRI/en,pr['evidence_sha256'][en])
 def locate(s):return ROOT/next(p for p in pr['input_sha256'] if p.endswith(s))
 special=locate('actual-special-fiber20-coefficients.json')
 keys=['actual_physical_rows9','source_F9_ids','YB_in_YA_linear','A0','actual_generator_normalization','source_conductor_generators_rows_and_ids','ell_exact']
 S={k:field(special,k) for k in keys}
 ca=[7148,7199,7454,7505,7688,7943,8045];cb=[27897,27948,28203,28254,28437,28692,28794]
 assert S['actual_physical_rows9']==[27,108,109,112,113,117,118,119,61]
 assert S['ell_exact']=='5*(s_CA-r_CA)'
 assert S['source_F9_ids']=={'CA':ca+[17493,19106],'CB':cb+[38273,39886]}
 IAids=[i for row,i in S['source_conductor_generators_rows_and_ids']['CA']]
 IBids=[i for row,i in S['source_conductor_generators_rows_and_ids']['CB']]
 assert IAids==ca+[7790] and IBids==cb+[28539]
 full=locate('full9-firstjet-certificate.json');slotmons=field(full,'source2_monomials')
 source2=field(full,'both57_source2');source2_peak=count(source2)
 mons2=[exp(slots) for slots in slotmons]
 assert len(mons2)==55 and len(set(mons2))==55
 unitmons=[exp([i]) for i in range(9)];firstix=[mons2.index(ZERO)]+[mons2.index(e) for e in unitmons]
 sourceproof=[]
 for fam,ids,f63 in [('CA',ca+[17493,19106],7790),('CB',cb+[38273,39886],28539)]:
  x=source2[fam]
  assert x['F9_nodes']==ids and x['f63']==f63 and x['all57_equations_all55_coefficients_verified']
  values=x['pure_function2jets'][str(f63)]
  assert len(values)==55 and all(Q(values[j])==0 for j in firstix)
  sourceproof.append([fam,f63,firstix])
 del source2
 B=[[Q(x) for x in row] for row in S['YB_in_YA_linear']]
 assert len(B)==9 and all(len(row)==9 for row in B)
 assert [[str(x) for x in row[:7]] for row in B[:7]]==S['A0']
 assert all(not B[i][j] for i in range(7) for j in [7,8])
 IA=[{unitmons[i]:Q(1)} for i in range(7)]+[{}]
 IB=[{unitmons[j]:B[i][j] for j in range(7) if B[i][j]} for i in range(7)]+[{}]
 c0=locate('honest-C0-D0-D4-DAG.json')
 assert fastfield(c0,'actual_C0_basis8_scalar_nodes')==ca+[505938]
 mg=locate('actual-T2-mu-input-source-and-cache-gate.json')
 transport=locate('actual-active-unit-bases-and26-frame-transport.json')
 M=field(mg,'Mnew0_columns');assert field(transport,'actual_Mnew0_columns')==M
 nfpath=locate('actual660-monomial-NFs-and-combined-rho-raw33-factor-input.json')
 mons=[tuple(x) for x in field(nfpath,'monomials660_degree_le4')]
 NF=[sp(x) for x in field(nfpath,'all660_actual_raw33_normal_forms')]
 basisidx=field(nfpath,'actual33_basis_monomial_indices')
 rawpath=locate('actual-ell2-product33-and-raw45-actions.json')
 b74=field(rawpath,'basis74');free=field(rawpath,'free_column_indices');basis=[]
 for k in free:
  sector,r,e=b74[k];v=[0]*9;v[7]=r;v[8]=e
  if sector:v[sector-1]=1
  basis.append(tuple(v))
 assert len(basis)==len(set(basis))==33 and basis[0]==ZERO
 assert basis==[mons[j] for j in basisidx]
 expected=[exp(slots) for n in range(5) for slots in itertools.combinations_with_replacement(range(9),n) if slots.count(8)<2]
 assert mons==expected==[tuple(x) for x in D['all660_monomial_inputs']] and len(mons)==660
 pairs=[(0,7),(0,8),(1,7),(1,8),(4,7),(4,8),(5,7),(5,8),(7,7),(7,8)]
 targets=[exp(pair) for pair in pairs];selectors=[basis.index(e) for e in targets]
 assert selectors==[22,23,20,21,14,15,12,13,24,25]
 assert D['pair_order10']==[list(x) for x in pairs]
 assert [tuple(x) for x in D['ordinary_target_monomials10']]==targets
 assert D['lambda10_by33_sparse_rows']==[[[j,'1']] for j in selectors]
 assert D['lambda10_by33_dense_rows']==[['1' if j==k else '0' for j in range(33)] for k in selectors]
 # Independent nonhomogeneous NF regeneration from the ACTUAL raw33
 # coordinate actions, using integer matrices with one common denominator.
 actions=[[sp(col) for col in cols] for cols in field(rawpath,'coordinate_actions9_sparse_columns')]
 assert len(actions)==9 and all(len(cols)==33 for cols in actions)
 den=math.lcm(*(c.denominator for cols in actions for col in cols for c in col.values()))
 integer=[[[ (i,c.numerator*(den//c.denominator)) for i,c in col.items()] for col in cols] for cols in actions]
 L=max(sum(abs(v) for i,v in col) for cols in integer for col in cols)
 assert den==11828436480 and den.bit_length()==34 and L.bit_length()==38
 numerator_bound=L**4;denominator_bound=den**4
 nfheight=max(max(abs(x.numerator).bit_length(),x.denominator.bit_length()) for nf in NF for x in nf.values())
 comparison=max(numerator_bound.bit_length()+nfheight,denominator_bound.bit_length()+nfheight)+1
 literal=max(max(abs(x.numerator).bit_length(),x.denominator.bit_length()) for row in B[:7] for x in row[:7])
 assert literal==7 and numerator_bound.bit_length()==150 and denominator_bound.bit_length()==134 and comparison<2048
 sourceowned=count([pr,hr,D,S,slotmons,mons2,firstix,sourceproof,B,IA,IB,M,mons,NF,basisidx,b74,free,basis,pairs,targets,selectors,actions,integer,pins])+source2_peak+20000
 # Include future660 integer sparse states, all9 controller labels perstate,
 # transient33 accumulators, evidence copies, and all positive metadata.
 future=sourceowned+660*(9+2*33)+3000
 assert future<150000 and RAW+2000000<32<<20
 prior={'source_only_before_product64_and_new_integer660NF_checks':True,
  'literal_copy_bits':literal,'primary_selector_NF_prior_scope_bits':pr['positive_prior_bits'],
  'integer_common_denominator':str(den),'integer_action_column_L1_max':str(L),
  'four_action_integer_numerator_bound':str(numerator_bound),'four_action_denominator_bound':str(denominator_bound),
  'integer_numerator_bits':150,'denominator_bits':134,'comparison_preoperation_bound':comparison,
  'source_owned_including_FULL9_transient':sourceowned,'whole_owned_future_with_dense_state_reserve':future,
  'RAW_unique_inputs_plus_output_reserve':RAW+2000000,'hard_numeric':150000,'hard_bits':2048}
 save(out/'positive-before-independent-source-arithmetic.json',prior)
 # Rebuild all64 ORIGINAL products; since IA is one coordinate orzero,
 # each coefficient is copied, so no arithmetic height is concealed.
 prod=[]
 for i,aa in enumerate(IA):
  for j,bb in enumerate(IB):
   row={}
   for x,one in aa.items():
    assert one==1
    for y,c in bb.items():
     e=tuple(a+b for a,b in zip(x,y));assert e not in row;row[e]=c
   assert all(sum(e)==sum(e[:7])==2 for e in row)
   assert not any(row.get(e,Q(0)) for e in targets)
   saved=D['all64_literal_product_degree2_initials'][8*i+j]
   assert saved['source_nodes']==[IAids[i],IBids[j]]
   assert row=={tuple(e):Q(c) for e,c in saved['homogeneous_degree2_initial']}
   prod.append([i,j,len(row)])
 # Each monomial is recomputed from a lower monomial through a different
 # source path than the saved NF producer. Nonhomogeneous matrix entries
 # are retained exactly; no degree truncation of a resulting raw33 vector.
 states={ZERO:{0:1}};actualbits=1;nfhash=hashlib.sha256();lambda_hash=hashlib.sha256();checks=0
 for e,nf in zip(mons,NF):
  degree=sum(e)
  if degree:
   var=max(i for i,x in enumerate(e) if x);lo=list(e);lo[var]-=1;lo=tuple(lo);acc={}
   for j,c in states[lo].items():
    for i,v in integer[var][j]:
     # The positive column-L1 bound controls EACH product/partialsum,
     # even before later exact cancellations.
     assert abs(c*v)<=L**degree
     acc[i]=acc.get(i,0)+c*v;assert abs(acc[i])<=L**degree
   states[e]={i:x for i,x in acc.items() if x}
  current=states[e];dn=den**degree
  for i in set(current)|set(nf):
   num=current.get(i,0);q=nf.get(i,Q(0))
   left=num*q.denominator;right=q.numerator*dn
   actualbits=max(actualbits,abs(left).bit_length(),abs(right).bit_length())
   assert actualbits<=comparison and left==right,(e,i)
  nfhash.update(json.dumps([[i,str(x)] for i,x in sorted(current.items())],separators=(',',':')).encode())
  outputs=[]
  for k,(t,sel) in enumerate(zip(targets,selectors)):
   assert current.get(sel,0)==dn*int(e==t),(e,k)
   assert nf.get(sel,Q(0))==int(e==t);checks+=1
   if e==t:outputs.append([k,'1'])
  assert outputs==D['all660_lambda_values10_sparse'][len(states)-1] if degree else outputs==D['all660_lambda_values10_sparse'][0]
  lambda_hash.update(json.dumps(outputs,separators=(',',':')).encode())
 assert checks==6600 and len(states)==660
 assert all(NF[j]=={i:Q(1)} for i,j in enumerate(basisidx))
 # Pin the original slot-multiset reader failure and explicit7-vs1 scope.
 oldphase=pin(ROOT/'reports'/PRI/'phase-positive-before-source-checks.json',hr['input_sha256'][f'reports/{PRI}/phase-positive-before-source-checks.json'])
 assert 'copiesofliteralBlinearcoefficients' in oldphase['positive_literal_product_bound']
 assert hr['literal_product_copy_numerator_denominator_maxbits']==literal==7
 assert hr['primary_positive_prior1_scope']=='lambda10timesNF660dotsonly'
 assert pr['positive_prior_bits']==1 and hr['metadata_only_no_product_or_NF_rerun']
 for p,h in pins.items():assert sha(ROOT/p)==h
 proof={'all64_actual_literal_IA8_times_IB8_quadratic_initials_rebuilt':True,
  'all660_physical_monomials_and_all660_integer_action_NFs_independently_rebuilt':True,
  'all6600_lambda_NF_identities_and_nonhomogeneous_carry_checked':True,
  'selector_indices10':selectors,'lambda_shape':[10,33],
  'ordinary_coefficient_no_factorial':True,
  'exact_source_descent_proof':[
   'The actual IA generators are CAz7 andrawf637790; IB generators areCBz7 andrawf6328539. Source2 multisets were independently converted to9exponents. Both rawf63 have exactzero constant/first in their OWN centeredF9 charts. Centered common pullback preserves that order.',
   'Actual common CBlinearfirst7 have no r/s entries. Thus IAdegree1=zi/0 andIBdegree1=A0zi/0, and ALL64 literal productdegree2 initials arezz. Arbitrary formal multipliers contribute onlytheir constants atdegree2; allten nonzz degree2 ordinarycoefficient functionals annihilate the entire actualproductideal,ell2,andm3.',
   'Accepted source-equality and m5zero identify actualR/(IAIB,ell2) withraw33. Lambda is specified on33LITERAL monomialbasis elements. We independently rebuilt everydegree<=4,ellpower<2 monomialNF fromnine actualnonhomogeneous raw33coordinate actions, then checkedall6600 ordinarycoefficient identities. Higher carry is retained exactly, includingdegree changes in the actiontables.',
   'OLD Q9 uses actualC0 CAz7+gA505938, notrawf637790. The samephysicalF9 rows,ell=5(s-r),andcurrentM0 transport were matched literally. No M1 orcomparison gauge coefficient was assignedzero.'
  ],
  'future_scope':'Componentwise application to actualcanonical nativekappa/G8requires a separate samebranch/free4/monic/overlap/resource sourcegate. No response was read. Honest regular g/hquotients remain opaque; raw33 supplies no nonunitdivision orhighermixedfunctional certificate.',
  'height_clarification_checked':'Primaryprior1 covers selector*NF dots only. All64 literal products arecopieswithsourceheight7, as recorded BEFORE checks. The separate HEIGHT receipt expands that scope without changing primarybytes or budgets. Our DIFFERENT independentintegerNF method has its own prior150/134bit and bounded comparison.',
  'independent_integer_NF_hash':nfhash.hexdigest(),'all660_lambda_stream_hash':lambda_hash.hexdigest(),
  'actual_integer_comparison_maxbits':actualbits,'positive_budget':prior,
  'no_kappa_G8_response_honest_jet_or_upper_values':True}
 ep=out/'independent-64-original-products-660-integer-NFs-and6600-descent-checks.json';save(ep,proof)
 r={'id':ID,'all_passed':True,'frozen':True,'SOURCE_only_no_functionresponses':True,
  'input_sha256':pins,'code_sha256':sha(__file__),'evidence_sha256':{ep.name:sha(ep)},
  'all64_all660_and6600_nonhomogeneous_source_identities_PASS':True,
  'source_owned_upper':sourceowned,'whole_future_owned_upper':future,'seconds':time.monotonic()-start}
 save(out/'report.json',r);signal.alarm(0)
 print(json.dumps({'all_passed':True,'report_sha256':sha(out/'report.json'),'owned':future,'source':sourceowned,'actualbits':actualbits,'seconds':r['seconds']}))
if __name__=='__main__':main()
