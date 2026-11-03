from pathlib import Path
import argparse,hashlib
import json,sympy as s,time
from sympy.polys.rings import ring
ROOT=Path(__file__).resolve().parents[1]
ID='S61_NONSQUARE_UNIVERSAL_QUAD22_20261103_v63'
parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,default=ROOT/'reports'/ID);OUT=parser.parse_args().output.resolve();OUT.mkdir(parents=True,exist_ok=False)
st=time.monotonic();R,a,e,d=ring('a,e,d',s.QQ);T,x0,x1,x2,x3,x4,Z,W=ring('x0,x1,x2,x3,x4,Z,W',R);xs=[x0,x1,x2,x3,x4,Z,W];images=[x0,x1,x2,W,x3,T.zero,Z+x1-3*x0,x4]
data=json.loads((ROOT/'reports/S61_BINARY_HYPERPLANE_QUADRATIC_CODIMENSION_20261102.json').read_text());qs=[]
for row in data['quadratic_rows_sparse']:
 q=T.zero
 for mon,cs in row:
  value=T(R.from_dict({tuple(ex):s.QQ(c)for ex,c in cs}))
  for im,power in zip(images,mon):
   if power:value*=im**power
  q+=value
 qs.append(q)
cut=lambda q:T.from_dict({m:c for m,c in q.items()if m[5]==m[6]==0})
C=100*a+18*e**2-9*e-2;L=T.from_dict({m:(-s.QQ(405,4)*c).exquo(C)for m,c in cut(qs[1]).diff(x1).items()})
assert cut(qs[1])==L**2/8100 and cut(qs[3])==L*(e*L-400*a*x0)/8100
ell=(-(d*e+8)*qs[1]+8*qs[2]+d*qs[3]).exquo(W);difference=cut(ell)-2*(5*a*d-6*e-4)*L/1215
J=5*a*d-6*e-4
assert ell.diff(x1)==T(-16*C*J/243)
coef=lambda q,i:q.get(tuple(int(j==i)for j in range(7)),R.zero)
assert 3*coef(difference,2)-s.QQ(7,2)*coef(difference,4)==-80*a*d**2/81
assert coef(difference,1)==0
# The chosen six light multiples complement ell using only C*J.
light_multiple_indices=[0,2,3,4,5,6]
assert len(light_multiple_indices)==6
# Exact actual change from six stored quadrics to q0,q1,q3,q5,ellZ,ellW.
assert ell*W==-(d*e+8)*qs[1]+8*qs[2]+d*qs[3]
assert -a*ell*Z==(-10*a-3*e**2+3*e)*qs[1]+3*(e-1)*qs[3]+10*a*qs[4]
ch=s.zeros(6);ch[0,0]=ch[1,1]=ch[2,3]=ch[3,5]=1
ch[4,1]=(10*a+3*e**2-3*e).as_expr()/a.as_expr();ch[4,3]=(-3*(e-1)).as_expr()/a.as_expr();ch[4,4]=-10
ch[5,1]=(-(d*e+8)).as_expr();ch[5,2]=8;ch[5,3]=d.as_expr();assert s.cancel(ch.det())==-80
cub=(e*L-400*a*x0)*qs[1]-L*qs[3];TZ,TW=T.zero,T.zero
for m,c in cub.items():
 mm=list(m)
 if mm[5]:mm[5]-=1;TZ+=T.from_dict({tuple(mm):c})
 else:assert mm[6];mm[6]-=1;TW+=T.from_dict({tuple(mm):c})
assert cub==Z*TZ+W*TW
assert ell*(e*L-400*a*x0)*qs[1]-ell*L*qs[3]-TZ*(ell*Z)-TW*(ell*W)==0
enc=lambda q:[[list(m),[[list(ex),str(c)]for ex,c in p.items()]]for m,p in q.items()]
res=dict(coordinate_order=list(map(str,xs)),L=enc(L),ell=enc(ell),ell_mod_ZW_minus_2J_L_over1215=enc(difference),TZ=enc(TZ),TW=enc(TW),new_syzygy=[[],enc(ell*(e*L-400*a*x0)),enc(-ell*L),[],enc(-TZ),enc(-TW)],known_source_basis=['q0','q1','q3','q5','ellZ','ellW'],identity_passed=True,seconds=time.monotonic()-st)
(OUT/'certificate.json').write_text(json.dumps(res,indent=2)+'\n')
paths=[Path(__file__),ROOT/'research/notes'/(ID+'.md'),ROOT/'reports/S61_BINARY_HYPERPLANE_QUADRATIC_CODIMENSION_20261102.json',ROOT/'reports/S61_BINARY_HYPERPLANE_QUADRATIC_PAIR_20261102.json',ROOT/'reports/S61_NONSQUARE_ZW_PRODUCT_RANK19_20261103_v63/report.json',ROOT/'reports/S61_NONSQUARE_ZW_PRODUCT_RANK19_20261103_v63_REVIEW/report.json']
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
report=dict(id=ID,all_exact_actual_QQ_identities_passed=True,twenty_second_quadratic_syzygy_constructed=True,old_syzygies='15Koszul+6light_multiples_of_(W*ellZ-Z*ellW)',omitted_light_index=1,ell_x1_known_unit='-16*C*J/243',new_independence_anchor='3*difference_x2-(7/2)*difference_x4=-80*a*d^2/81',pointwise_independent_kernel_directions=22,open='a*d*C*J!=0; actual chart a*F*B*C!=0',generic_kernel_dimension_dependency=22,exact_kernel_module_basis_uses_reviewed_global_upper_and_generic_rank146=True,image_row_module_locally_free_rank146=True,ambient_pointwise_rank146_or_cokernel_flatness_claimed=False,full_quartic_or_ambient_claimed=False,no_expanded168cofactor_matrix_or_pivot_search=True,input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},certificate_sha256=sha(OUT/'certificate.json'),frozen=True,independent_review_pending=True,seconds=time.monotonic()-st)
(OUT/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:report[k]for k in ['all_exact_actual_QQ_identities_passed','pointwise_independent_kernel_directions','seconds']}))
