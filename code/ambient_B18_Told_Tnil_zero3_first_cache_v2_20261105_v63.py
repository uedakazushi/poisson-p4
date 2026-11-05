#!/usr/bin/env python3
"""Cache-only old canonical m1, full81 cleanup, and FORMAL second callbacks.
No new source/gauge jet or assignment of an unavailable second coefficient.
"""
import collections,hashlib,itertools,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
import ambient_B18_honest_nil_T_first_values_review_20261104_v63 as Q
ROOT=Path(__file__).resolve().parents[1]
ID='AMBIENT_B18_TOLD_TNIL_ZERO3_FIRST_CACHE_V2_20261105_v63';OUT=ROOT/'reports'/ID
MON=[()]+[(i,)for i in range(8)]+list(itertools.combinations_with_replacement(range(8),2));IX={m:i for i,m in enumerate(MON)}
EXP=[]
for mon in MON:
 e=[0]*8
 for i in mon:e[i]+=1
 EXP.append(e)
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def enc(v):
 if isinstance(v,F):return str(v)
 if isinstance(v,dict):return {str(k):enc(x)for k,x in v.items()}
 if isinstance(v,(list,tuple)):return [enc(x)for x in v]
 return v
def mul(a,b):
 z={}
 for i,x in a.items():
  for j,y in b.items():
   m=tuple(sorted(MON[i]+MON[j]))
   if len(m)<=2:Q.put(z,IX[m],x*y)
 return z
def inverse(a):
 c=a[0];v=Q.plus(Q.scaled(a,-1/c),{0:F(1)})
 return Q.scaled(Q.plus({0:F(1)},v,mul(v,v)),1/c)
def main():
 st=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));OUT.mkdir(exist_ok=True);pins={}
 def load(s,n='report.json'):
  p=ROOT/'reports'/s/n;pins[str(p.relative_to(ROOT))]=sha(p);return json.loads(p.read_text())
 def save(n,v):(OUT/n).write_text(json.dumps(enc(v),indent=2)+'\n')
 assert load('AMBIENT_B18_TOLD_TNIL_ZERO3_FIRST_SOURCE_PREFLIGHT_20261104_v63')['all_passed']
 gate=load('AMBIENT_B18_TOLD_TNIL_ZERO3_FIRST_SOURCE_PREFLIGHT_20261104_v63','three-oldm1-cache-and-cleanup-callback-gate.json')
 save('phase.json',dict(phase='SOURCE gate pinned; cache recovery not yet run',no_new_source_or_jet=True))
 J=load('AMBIENT_B18_T3_FULL9_FIRSTJET_20261104_v63','full9-firstjet-certificate.json');src=load('AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63','actual-source-coefficients.json')
 phys=load('AMBIENT_B18_HONEST_NIL_T_FIRST_LINEAR_20261104_v63','phase1-cached-source-monics-and-physical-adapter.json')
 O=load('AMBIENT_B18_FULLV56_COMPARISON_UNIT_20261104_v63','actual-fullV56-m0-and-typed-comparison.json');N=load('AMBIENT_B18_TOLD_TNIL_NATURAL_MAP_20261104_v63','actual-natural-Told-to-Tnil-source-map-and-origin.json')['actual_origin_data']
 X=load('AMBIENT_B18_HONEST_NIL_T_FIRST_VALUES_20261104_v63','phase2-actual-radical-X-firsts.json');W=load('AMBIENT_B18_HONEST_NIL_T_FIRST_VALUES_20261104_v63','phase3-actual-W84-rank-and-mu.json')
 M=load('AMBIENT_B18_HONEST_NIL_T_MINIMAL_OPERATOR_20261104_v63','actual-minimal-Tnil5-849-and-d3-linear-operator.json');mD=load('AMBIENT_B18_HONEST_NIL_T_FIRST_VALUES_20261104_v63','phase1-actual-m0-and-section.json')
 first=load('AMBIENT_B18_HONEST_NIL_T_FIRST_LINEAR_20261104_v63','phase2-row20-surrogates-and-true-m1row7.json');lin=load('AMBIENT_B18_HONEST_NIL_T_FIRST_LINEAR_20261104_v63','phase3-all849-minimal-first-relations-and-rank.json')
 whole=load('AMBIENT_B18_HONEST_NIL_T_LINEAR_SOURCE_WHOLE_PREFLIGHT_20261104_v63','actual-final681-source-sensitivities-and-cached-directions.json')
 B=load('AMBIENT_B18_COMPARISON_X_FIRST640_LINEAR_20261104_v63','actual-first640-X-mod-m2.json');At=load('AMBIENT_B18_COMPARISON_X_FIRST640_LINEAR_20261104_v63','phase1-actual-atom-firsts.json')
 R=load('AMBIENT_B18_RADICAL_COMPARISON_OPERATOR_20261104_v63','actual-radical-normalizer-comparison-operator.json')
 # All values here are in the saved physical ell=0 coordinate frame.
 def fold(v):
  z={}
  for mon,c in zip(J['source2_monomials'],v):
   if F(c):Q.put(z,IX[tuple(sorted(7 if i==8 else i for i in mon))],F(c))
  return z
 Y=[Q.sparse(v)for v in phys['CB_F9_in_CA_F9_ell0_degree2']]
 def compose(v):
  z={}
  for mon,c in zip(J['source2_monomials'],v):
   if not F(c):continue
   p={0:F(c)}
   for i in mon:p=mul(p,Y[i])
   z=Q.plus(z,p)
  return z
 leaves={f:{i:(fold(v)if f=='CA'else compose(v))for i,v in J['both57_source2'][f]['all57_source2jets']}for f in['CA','CB']}
 def source(root,family,closure):
  nodes=src['DAG_scalar_nodes'];uses=collections.Counter([root]);seen=set(closure)
  for i in seen:
   if nodes[i][0]in['a','m','n','i']:uses.update(nodes[i][1:])
  vals={}
  for i in sorted(seen):
   n=nodes[i];op=n[0]
   if op=='c':v={0:F(n[1],n[2])}if n[1]else{}
   elif op=='v':v=leaves[family][i]
   elif op=='a':v=Q.plus(vals[n[1]],vals[n[2]])
   elif op=='m':v=mul(vals[n[1]],vals[n[2]])
   elif op=='n':v=Q.scaled(vals[n[1]],-1)
   elif op=='i':v=inverse(vals[n[1]])
   else:raise AssertionError(op)
   vals[i]=v
   if op in['a','m','n','i']:
    for j in n[1:]:
     uses[j]-=1
     if not uses[j]:del vals[j]
  return vals[root]
 points=[source(p['root'],p['family'],p['closure'])for p in gate['four_cached_point_function_profiles']]
 delta5=Q.plus(points[0],Q.scaled(points[1],-1));delta7=Q.plus(points[2],Q.scaled(points[3],-1))
 fA=fold(J['both57_source2']['CA']['pure_function2jets']['7790']);fB=compose(J['both57_source2']['CB']['pure_function2jets']['28539'])
 assert all(not p.get(i,F(0))for p in[delta5,delta7,fA,fB]for i in range(9))
 # Cross-check CA f63 against the independently saved seven ordered firsts.
 for p in gate['CA_rawf63_first_prefixes']:
  i=p['prefix_index'];have=Q.pairs(p['first9_bar'])
  expect={q+1:fA.get(IX[tuple(sorted((i,q)))],F(0))for q in range(8)if q<=i or q==7}
  expect={k:v for k,v in expect.items()if v};assert have==expect
 nu=[Q.sparse(v)for v in phys['actual_nu_A_B_degree2']];atom={i:Q.pairs(v)for i,v in At['all79_atom_values']}
 bA5=atom[7943].get(0,F(0));bB5=atom[28692].get(0,F(0));assert bA5==bB5==0
 functions=[[Q.scaled(mul(atom[28692],nu[1]),-1),mul(atom[7943],nu[0]),delta5],[delta7,{},Q.scaled(fB,-1)],[{},delta7,fA]]
 old_m1=[]
 for col,curves in zip(functions,gate['nine_r_first_values_cached']):
  v=[{}for _ in range(24)]
  for sec,p in enumerate(col):
   for i in range(7):v[8*sec+i]={q:p[IX[tuple(sorted((i,q)))]]for q in range(8)if(q<=i or q==7)and p.get(IX[tuple(sorted((i,q)))],F(0))}
   if F(curves[sec]):v[8*sec+7]={7:F(curves[sec])}
  old_m1.append(v)
 oldn=[Q.sparse(v)for v in N['old_n0_56_by24']];cold1=[Q.mprod(oldn,v)for v in old_m1];cold1=[[Q.scaled(v,-1)for v in col]for col in cold1]
 # Reconstruct full Z0/Z1, not only the three initially selected derivatives.
 C36=Q.FiniteC(EXP[:9],[[v.get(i,F(0))for i in range(9)]for v in nu]);C180=Q.FiniteC(EXP,[[v.get(i,F(0))for i in range(45)]for v in nu])
 Z0=[{}for _ in range(44)];Z1=[[{}for _ in range(56)]for _ in range(44)];direct20=[]
 for i,row in enumerate(X['new_A12_by9_C36']):
  for gen,a in enumerate(row):
   for mask in range(4):
    t=C36.mul(Q.sparse(a),{mask:F(1)})
    for out in range(4):
     rr=4*gen+out;cc=4*i+mask
     if t.get(out,F(0)):Z0[rr][cc]=t[out]
     Z1[rr][cc]={q:t[4*(q+1)+out]for q in range(8)if t.get(4*(q+1)+out,F(0))}
 for i in range(8):Z0[36+i][48+i]=F(1)
 assert Z0==[Q.sparse(v)for v in N['actual_Z0_44_by56']]
 for oldcol in[23,29,30]:
  i,mask=divmod(oldcol,4);v=C180.mul(Q.sparse(X['new_A12_by9_C36'][i][5]),{mask:F(1)})
  direct20.append(C180.component(v,0))
 C0=[Q.sparse(N['all_five_new_cycle_origins44'][j])for j in[1,2,3]];C1=[];C20=[]
 for q,(oldcol,cold)in enumerate(zip([23,29,30],cold1)):
  c1=[Q.plus(*(Q.scaled(cold[k],a)for k,a in row.items()),Z1[i][oldcol])for i,row in enumerate(Z0)];C1.append(c1)
  z=Q.plus(direct20[q],*(mul({p+1:a for p,a in Z1[20][k].items()},{p+1:a for p,a in v.items()})for k,v in enumerate(cold)if v and Z1[20][k]));C20.append(z)
  assert [z.get(i,F(0))for i in range(9)]==[C0[q].get(20,F(0))]+[c1[20].get(i,F(0))for i in range(8)]
 save('phase.json',dict(phase='actual canonical oldm1 and new c1 recovered from caches; no A2 assigned',seconds=time.monotonic()-st))
 # Only fixed source eliminations; arbitrary second callbacks are kept separate.
 m0=[Q.sparse(v)for v in mD['m0_24_by44']];xp=X['X0_pivot_columns'];Xcol=[{int(r):Q.pairs(v)for r,v in c}for c in X['all780_source_columns44_sparse_dual9']]
 free=[i for i in range(44)if i not in mD['selected_source_columns24']];I=[free[i]for i in M['fixed_minor_data']['compatible_X0_top18_rows_in_C20']];Ji=[Q.sparse(v)for v in M['fixed_minor_data']['inverse']]
 Orel=W['residual_original_X_columns'];wp=W['W_pivot_residual_columns'];wr=W['W_pivot_rows'];Wi=[Q.sparse(v)for v in W['W_fixed_minor_inverse']];wcols=[Q.pairs(v)for v in W['all_W84_columns_sparse']];xc=[Q.sparse(v)for v in X['all780_X0_factorization_coefficients']]
 def cres(j):return Q.plus({j:F(1)},*({xp[k]:-a}for k,a in xc[j].items()))
 def combine(vals,c):return Q.plus(*(Q.scaled(vals[k],a)for k,a in c.items()))
 raw15=[[F(v)for v in row]for row in first['all780_row20_second_surrogates15']];needed=[tuple(m)for m in whole['all15_needed_row20_X2_monomials']]
 bold0=[]
 for col in B['source_columns640_rows56_sparse_dual9']:
  bold0.append({row:Q.pairs(sp)[0]for row,sp in col if row<48 and Q.pairs(sp).get(0,F(0))})
 bold0 += [{}for _ in range(140)]
 sensitiv=[{}for _ in range(780)]
 for j,col in enumerate(R['native_tensor_Koszul6_HB2_resolution'][1]['columns']):
  for row,pol in col:
   if row==5:
    for e,c in pol:
     if sum(e)==1 and any(e[:5]):Q.put(sensitiv[640+4*j],e.index(1),F(c))
 PAIRS=list(itertools.combinations(range(8),2));records=[]
 for j in range(3):
  a0=Q.mvec(Ji,{k:C0[j][i]for k,i in enumerate(I)if C0[j].get(i,F(0))});co={xp[k]:-a for k,a in a0.items()}
  d1=[Q.plus(C1[j][i],*({v-1:coef*a for v,coef in Xcol[k].get(i,{}).items()if v}for k,a in co.items()))for i in range(44)]
  v=Q.mprod(m0,d1);h0={}
  for sec in range(3):
   for k,(a,b)in enumerate(PAIRS):
    c=v[8*sec+7].get(a,F(0))if b==7 else(v[8*sec+b].get(a,F(0))-v[8*sec+a].get(b,F(0)))/2
    if c:h0[28*sec+k]=c
  beta=Q.mvec(Wi,{k:h0[r]for k,r in enumerate(wr)if h0.get(r,F(0))});assert combine([wcols[wp[k]]for k in range(81)],beta)==h0
  co=Q.plus(co,*(Q.scaled(cres(Orel[wp[k]]),-a)for k,a in beta.items()))
  final1=[Q.plus(C1[j][i],*({v-1:coef*a for v,coef in Xcol[k].get(i,{}).items()if v}for k,a in co.items()))for i in range(44)]
  final20=[C20[j].get(IX[m],F(0))+sum((a*raw15[k][t]for k,a in co.items()),F(0))for t,m in enumerate(needed)]
  Aweight=combine(bold0,co);Uweight=combine(sensitiv,co)
  a1=[Q.plus(*(Q.scaled(final1[I[k]],a)for k,a in row.items()))for row in Ji];m1r=list(map(F,first['true_actual_m1_row7_r_coefficients44']));hi=lin['all54_actual_pivot_h0_readback'];out=[]
  for t in[42,43]:out.append(Q.plus(final1[t],*(Q.scaled(a1[k],-Xcol[xp[k]].get(t,{}).get(0,F(0)))for k in range(18))))
  for hidx,i in enumerate([0,4,5]):
   z={}
   for q in range(8):
    c=F(5,2)*final20[needed.index(tuple(sorted((i,q))))]if q<=i or q==7 else F(0)
    if q==7:c+=sum((m1r[k]*final1[k].get(i,F(0))for k in range(44)),F(0))
    c-=sum((F(hi[hidx][k])*a1[k].get(q,F(0))for k in range(18)),F(0))
    if c:z[q]=c
   out.append(z)
  records.append(dict(old_free56=[23,29,30][j],old_m1_24x8=old_m1[j],new_c1_44x8=C1[j],first18=a0,full81=beta,total_actual_X780_constant_coefficients=co,final_known40={8*i+q:a for i,p in enumerate(out)for q,a in p.items()},remaining_actual_adapter_A2_weights=Aweight,remaining_native_even2_weights=Uweight,NOT_actual_first_column_until_callbacks_closed=True))
 checks=0
 for j,r in enumerate(records):
  co=r['total_actual_X780_constant_coefficients'];c0=Q.plus(C0[j],*(Q.scaled({i:v[0]for i,v in Xcol[k].items()if v.get(0,F(0))},a)for k,a in co.items()));assert not c0
  final1=[Q.plus(C1[j][i],*({v-1:coef*a for v,coef in Xcol[k].get(i,{}).items()if v}for k,a in co.items()))for i in range(44)]
  assert all(not v for v in Q.mprod(m0,final1)),('FULL first Koszul coefficients',j)
  checks+=192
 assert checks==576
 assert sum(len(v)for r in records for v in r['old_m1_24x8'])<=324
 save('actual-three-oldm1-and-full81-formal-callbacks.json',dict(records=records,base_monomials=MON,required_second_monomials15=needed,old_n1_and_old_c2_cancel=True,all576_first_Koszul_coefficients_verified=True,full_b_times_nu_including_b1_nu1=True,no_unavailable_second_callback_assigned_zero=True,new_source_or_jet=False))
 for p in[Path(__file__),Path(Q.__file__),ROOT/f'research/notes/{ID}.md']:pins[str(p.relative_to(ROOT))]=sha(p)
 rep=dict(id=ID,all_passed=True,frozen=True,cache_only_old_m1_and_full81_cleanup=True,actual_zero3_first_or_HF_NOT_yet_claimed=True,input_sha256=pins,evidence_sha256={p.name:sha(p)for p in OUT.iterdir()if p.is_file()and p.name!='report.json'},seconds=time.monotonic()-st,max_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
 assert sum(p.stat().st_size for p in OUT.iterdir()if p.is_file())<32<<20
 save('report.json',rep);signal.alarm(0);print(json.dumps({k:v for k,v in rep.items()if k not in['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':main()
