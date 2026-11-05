#!/usr/bin/env python3
"""Independent homogeneous600 closure, exact24 quotient, fixed rank21.
No rank search, new cubic values, source solve or Groebner computation.
"""
import hashlib,itertools,json,math,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
import ambient_B18_honest_nil_T_first_values_review_20261104_v63 as Q
ROOT=Path(__file__).resolve().parents[1];BASE='AMBIENT_B18_HONEST_NIL_T_CUBIC_CLOSURE_20261105_v63';ID=BASE+'_REVIEW';OUT=ROOT/'reports'/ID
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def enc(v):
 if isinstance(v,F):return str(v)
 if isinstance(v,dict):return {str(k):enc(x)for k,x in v.items()}
 if isinstance(v,(list,tuple)):return [enc(x)for x in v]
 return v
def main():
 st=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));OUT.mkdir(exist_ok=True);pins={}
 def load(s,n='report.json'):
  p=ROOT/'reports'/s/n;pins[str(p.relative_to(ROOT))]=sha(p);return json.loads(p.read_text())
 def save(n,v):(OUT/n).write_text(json.dumps(enc(v),indent=2)+'\n')
 r=load(BASE);assert pins[f'reports/{BASE}/report.json']=='34788c870152e53c0a7417543d9b4450fdca4c9dfc1fcf0fcf6e6e4805b8972b'
 for p,h in r['input_sha256'].items():assert sha(ROOT/p)==h
 for n,h in r['evidence_sha256'].items():assert sha(ROOT/'reports'/BASE/n)==h
 pf='AMBIENT_B18_HONEST_NIL_T_CUBIC_CLOSURE_PREFLIGHT_20261105_v63';pre=load(pf);assert pre['all_passed']
 for p,h in pre['input_sha256'].items():assert sha(ROOT/p)==h
 I=load(pf,'literal-known-cubic-input-and-exact24-row-reduction.json');R=load(BASE,'fixed-known-cubic-rank-and-length-completion.json')
 assert load('AMBIENT_B18_HONEST_NIL_T_QUADRATIC_REMAINING_20261105_v63_REVIEW')['actual_HF_Tnil_through2']==[5,5,3]
 L=load('AMBIENT_B18_T_EXACT_LOEWY11_20261105_v63');assert pins['reports/AMBIENT_B18_T_EXACT_LOEWY11_20261105_v63/report.json']=='26c005cd7d6958970791fc3267cf4f1cfa5bcb2a34cd778e7f010c83958d1fbd'and L['actual_T_length']==11 and L['actual_Tnil_length']==16
 first=load('AMBIENT_B18_HONEST_NIL_T_FIRST_LINEAR_20261104_v63','phase3-all849-minimal-first-relations-and-rank.json');B=load('AMBIENT_B18_HONEST_NIL_T_QUADRATIC_CLOSURE_20261104_v63','known-quadratic-closure-and-complete814-source-projection.json');A=load('AMBIENT_B18_HONEST_NIL_T_QUADRATIC_VALUES_20261105_v63','phase3-all661-actual-six-quadratic-values-and-rank.json')
 save('phase.json',dict(phase='all actual source/length gates pinned; independent homogeneous reconstruction'))
 mon=[list(itertools.combinations_with_replacement(range(8),d))for d in[1,2,3]]
 def exponent(m):return tuple(m.count(i)for i in range(8))
 exps=[[exponent(m)for m in level]for level in mon];idx=[{e:i for i,e in enumerate(level)}for level in exps]
 linear=[Q.pairs(first['all849_actual_minimal_first_columns40_sparse'][j])for j in first['rank_minor_columns']];assert len(linear)==35 and sum(map(len,linear))==197
 # Rebuild all280 quadratic products and the fixed174 source witness.
 q280=[]
 for c in linear:
  for var in range(8):
   z={}
   for row,a in c.items():
    target,j=divmod(row,8);e=list(exps[0][j]);e[var]+=1;Q.put(z,36*target+idx[1][tuple(e)],a)
   q280.append(z)
 assert q280==[Q.pairs(c)for c in B['all280_known_actual_variable_multiple_columns180_sparse']]
 cp=B['fixed_rank_pivot_columns'];rr=B['fixed_rank_pivot_rows'];inv=[Q.sparse(v)for v in B['fixed_minor_inverse']];small=[{j:q280[c][row]for j,c in enumerate(cp)if q280[c].get(row,F(0))}for row in rr]
 assert Q.mprod(small,inv)==Q.mprod(inv,small)==Q.eye(174)
 pi6=[Q.pairs(v)for v in B['exact_quotient_duals_rows_sparse']];assert all(not Q.mvec(pi6,c)for c in q280)
 nf6=[36*t+idx[1][exponent(m)]for t,m in B['quotient_coordinate_labels']];assert [[row.get(k,F(0))for k in nf6]for row in pi6]==[[F(i==j)for j in range(6)]for i in range(6)]
 qcols=[list(map(F,A['all661_actual_quadratic_six_values'][c]))for c in A['fixed_rank_columns']];piv=A['fixed_rank_rows'];qi=[Q.sparse(v)for v in A['fixed_minor_inverse']];free=[i for i in range(6)if i not in piv];assert free==[2,4,5]
 p36=[]
 for row in free:
  c=Q.mprod([{i:qcols[i][row]for i in range(3)if qcols[i][row]}],qi)[0]
  p36.append(Q.plus({row:F(1)},{piv[i]:-a for i,a in c.items()}))
 pi3=Q.mprod(p36,pi6);nf3=[nf6[i]for i in free];assert nf3==[35,143,179]
 assert [[row.get(k,F(0))for k in nf3]for row in pi3]==[[F(i==j)for j in range(3)]for i in range(3)]
 # ker(pi3) is exactly the known177 quadratic subspace:174 source lower
 # witness +three actual independent images, with quotient rank3.
 assert all(not Q.mvec(p36,Q.sparse(c))for c in qcols)
 quadratic=[{nf6[i]:a for i,a in enumerate(c)if a}for c in qcols]
 assert sum(map(len,quadratic))==9 and quadratic==[Q.pairs(c)for c in I['all3_actual_quadratic_NF_relations180_sparse']]
 assert pi3==[Q.pairs(c)for c in I['quadratic_projection3_by180_sparse']]
 raw=[];labels=[]
 for j,col in enumerate(linear):
  for e2,m2 in zip(exps[1],mon[1]):
   z={}
   for row,a in col.items():
    target,var=divmod(row,8);e=tuple(x+y for x,y in zip(e2,exps[0][var]));Q.put(z,120*target+idx[2][e],a)
   raw.append(z);labels.append(['first35',first['rank_minor_columns'][j],list(m2)])
 for j,col in enumerate(quadratic):
  for var in range(8):
   z={}
   for row,a in col.items():
    target,q=divmod(row,36);e=list(exps[1][q]);e[var]+=1;Q.put(z,120*target+idx[2][tuple(e)],a)
   raw.append(z);labels.append(['quadratic3',A['fixed_rank_columns'][j],var])
 assert len(raw)==1284 and sum(map(len,raw))==7164
 assert raw==[Q.pairs(v)for v in I['all1284_known_cubic_columns600_sparse']]and labels==I['all1284_source_labels']
 q=[];divisors=[]
 for target in range(5):
  for e in exps[2]:
   factors=[i for i,n in enumerate(e)for _ in range(n)];var=factors.pop();j=36*target+idx[1][exponent(factors)]
   q.append({8*k+var:row[j]for k,row in enumerate(pi3)if j in row});divisors.append([target,factors,var])
 assert q==[Q.pairs(v)for v in I['canonical_q24_by600_columns_sparse']]and divisors==I['all600_divisor_witnesses_for_Id_minus_Eq_in_known_cubic']
 E=[]
 for row in nf3:
  target,mi=divmod(row,36)
  for var in range(8):
   e=list(exps[1][mi]);e[var]+=1;E.append(120*target+idx[2][tuple(e)])
 assert E==I['E600_by24_unit_column_rows']
 def project(col):return Q.plus(*(Q.scaled(q[i],a)for i,a in col.items()))
 reduced=[project(c)for c in raw]+[Q.plus({i:F(1)},Q.scaled(q[row],-1))for i,row in enumerate(E)]
 assert reduced==[Q.pairs(v)for v in I['equivalent_known_cubic24_by1308_sparse_columns']]
 den=[1]*24
 for c in reduced:
  for row,a in c.items():den[row]=math.lcm(den[row],a.denominator)
 assert den==I['row_denominator_units24']
 integer=[{row:a*den[row]for row,a in c.items()}for c in reduced];assert integer==[Q.pairs(v)for v in I['equivalent_integer24_by1308_sparse_columns']]
 cp=R['fixed_pivot_columns'];rr=R['fixed_pivot_rows'];assert len(cp)==len(rr)==R['fixed_rank']==21
 minor=[{j:integer[c][row]for j,c in enumerate(cp)if integer[c].get(row,F(0))}for row in rr];inv=[Q.sparse(v)for v in R['fixed_minor_inverse']]
 assert minor==[Q.sparse(v)for v in R['fixed_minor']]and Q.mprod(minor,inv)==Q.mprod(inv,minor)==Q.eye(21)and Q.det_fixed(minor)==F(R['fixed_minor_determinant'])
 hs=hashlib.sha256()
 for c in integer:
  a=Q.mvec(inv,{i:c[row]for i,row in enumerate(rr)if c.get(row,F(0))});assert Q.plus(*(Q.scaled(integer[cp[j]],v)for j,v in a.items()))==c
  hs.update((json.dumps([str(a.get(i,F(0)))for i in range(21)],separators=(',',':'))+'\n').encode())
 assert hs.hexdigest()==R['all1308_factorization_stream_sha256']
 # Independent normalized quotient of the KNOWN cubic module, with no
 # assumption that there are no further actual cubic initial relations.
 free24=[i for i in range(24)if i not in rr];assert free24==[7,15,23]
 f3=[]
 for row in free24:
  a=Q.mprod([{j:integer[c][row]for j,c in enumerate(cp)if integer[c].get(row,F(0))}],inv)[0]
  f3.append(Q.plus({row:F(1)},{rr[j]:-v*F(den[rr[j]],den[row])for j,v in a.items()}))
 assert all(not Q.mvec(f3,c)for c in reduced)
 full3=[{j:c for j,col in enumerate(q)if(c:=sum((a*col.get(i,F(0))for i,a in row.items()),F(0)))}for row in f3];assert all(not Q.mvec(full3,c)for c in raw)
 rcols=[]
 for row in nf3:
  target,mi=divmod(row,36);e=list(exps[1][mi]);e[7]+=1;rcols.append(120*target+idx[2][tuple(e)])
 assert rcols==[119,479,599]
 assert [[row.get(k,F(0))for k in rcols]for row in full3]==[[F(i==j)for j in range(3)]for i in range(3)]
 # Thus multiplication by r maps the actual degree2 space onto actual
 # degree3, even when further actual cubic relations quotient it down.
 # m³=r*m²+m⁴, and Nakayama on m³/(r*m²) gives exact equality.
 partitions=[list(p)for p in itertools.combinations_with_replacement(range(1,7),3)if sum(p)==6];assert partitions==[[1,1,4],[1,2,3],[2,2,2]]
 proof=dict(independent_all280_quad_products_and174_fixed_source_inverse=True,ker_pi3_equals_known177_by_actual_rank3=True,
  all1284_original_homogeneous600_columns7164terms=True,all600_Id_minus_Eq_divisor_witnesses=True,all24_consistency_sources_retained=True,
  exact_known_quotient_proof='Id600-Eq has image in K3. Hence E and q induce mutually inverse maps between F600/K3 and F24/(qK3+(Id-qE)F24). This is not an injective projection assumption.',
  fixed_reduced_rank=21,both21_inverse_products_and1308_factor_hash=True,known_original_rank=597,known_cubic_dimension=3,actual_HF3_upper=3,
  actual_total_Tnil_length=16,actual_HF_through2=[5,5,3],remaining_actual_tail_length=3,remaining_tail_options=[[3],[2,1],[1,1,1]],
  known_cubic_quotient3_by24=f3,known_cubic_quotient3_by600=full3,degree2_NF_basis180=nf3,
  r_multiplication_degree3_representative_rows600=rcols,r_multiplication_to_known_quotient_is_I3=True,
  new_exact_module_corollary='m³Tnil=rA*m²Tnil by surjectivity modulo m⁴ and Nakayama',
  all_higher_filtrations='m^n Tnil=rA^(n-2)*m²Tnil for n>=2',
  m2_one_variable_length=6,m2_one_variable_minimal_generators=3,m2_one_variable_Jordan_options=[[4,1,1],[3,2,1],[2,2,2]],
  NOT_selecting_an_actual_tail_or_S_module_split=True,rA_NOT_DVR_ell=True,no_new_source_jet_GB_rank_search_or_actual_cubic_value=True)
 save('independent-homogeneous600-closure-and-r-reduction-proof.json',proof)
 for p in[Path(__file__),Path(Q.__file__),ROOT/f'research/notes/{ID}.md']:pins[str(p.relative_to(ROOT))]=sha(p)
 result=dict(id=ID,all_passed=True,frozen=True,independent_known_homogeneous_degree3_quotient_dimension=3,known_original_rank=597,actual_HF3_upper=3,actual_m3_equals_rA_m2=True,actual_tail_NOT_determined=True,no_new_source_jet_GB_or_rank_search=True,input_sha256=pins,evidence_sha256={p.name:sha(p)for p in OUT.iterdir()if p.is_file()and p.name!='report.json'},seconds=time.monotonic()-st,max_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
 assert sum(p.stat().st_size for p in OUT.iterdir())<32<<20
 save('report.json',result);signal.alarm(0);print(json.dumps({k:v for k,v in result.items()if k not in['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':main()
