#!/usr/bin/env python3
"""SOURCE: actual seven K/ellK barcycles, after independent two-head values."""
import argparse,hashlib,json,re,signal,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_C_KBAR_ACTUAL_SEVEN_BASIS_SOURCE_20261105_v63'
PATHS={
 'values':'reports/AMBIENT_B18_C_CONNECTING_HEADS2930_CACHED_VALUES_20261105_v63/report.json',
 'oldvalues':'reports/AMBIENT_B18_C_CONNECTING_HEAD23_CACHED_VALUES_20261105_v63/phase2-full-coherent-head23-and-low-two-values.json',
 'oldpeer':'reports/AMBIENT_B18_C_CONNECTING_HEAD23_CACHED_VALUES_20261105_v63_REVIEW/report.json',
 'split':'reports/AMBIENT_B18_C_CONNECTING_KERNEL3_SPLIT_20261105_v63/actual-M7-stable-surjection-splitting-and-kernel-three-source-basis.json',
 'split_peer':'reports/AMBIENT_B18_C_CONNECTING_KERNEL3_SPLIT_20261105_v63_REVIEW/report.json',
 'uct':'reports/AMBIENT_B18_C_CONNECTING_SEVEN_SOURCE_20261105_v63/honest-seven-source-quotient-four-fixed-kernel-classes-and-edge-constraints.json',
 'uct_peer':'reports/AMBIENT_B18_C_CONNECTING_SEVEN_SOURCE_20261105_v63_REVIEW/report.json',
 'five':'reports/AMBIENT_B18_C_CONNECTING_FIVE_FORMAL_SOURCE_20261105_v63/actual-five-head-formal-RegDivell-operator-and-source-binding.json',
 'five_peer':'reports/AMBIENT_B18_C_CONNECTING_FIVE_FORMAL_SOURCE_20261105_v63_REVIEW/report.json',
 'namespace':'reports/AMBIENT_B18_C_CONNECTING_PLANE_ORDER_SOURCE_20261105_v63/all56-plane-orders-surviving-sigma2-adjoint-and-U9-namespace-correction.json',
}
MAXPRE=0
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def count(x):
 if isinstance(x,(int,float,Q)):return 1
 if isinstance(x,str):return int(bool(re.fullmatch(r'-?\d+(?:/\d+)?',x)))
 if isinstance(x,(list,tuple)):return sum(map(count,x))
 if isinstance(x,dict):return count(list(x.keys()))+count(list(x.values()))
 return 0
def h(q):return max(abs(q.numerator).bit_length(),q.denominator.bit_length())
def pre(n):
 global MAXPRE
 MAXPRE=max(MAXPRE,n);assert n<=1153
def add(a,b):pre(2*max(h(a),h(b))+1);q=a+b;assert h(q)<=384;return q
def mul(a,b):pre(h(a)+h(b));q=a*b;assert h(q)<=384;return q
def div(a,b):assert b;pre(h(a)+h(b));q=a/b;assert h(q)<=384;return q
def parse(s):
 a=str(s).split('/');n=int(a[0]);d=int(a[1])if len(a)>1 else 1
 assert d;pre(max(abs(n).bit_length(),abs(d).bit_length()));q=Q(n,d);assert h(q)<=384;return q
def mm(A,B):return[[sum_checked([mul(a,b)for a,b in zip(row,col)])for col in zip(*B)]for row in A]
def sum_checked(v):
 q=Q(0)
 for a in v:q=add(q,a)
 return q
def enc(A):return[[str(q)for q in row]for row in A]
def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--value-peer-report',required=True);args=ap.parse_args()
 start=time.monotonic();signal.alarm(60);paths=dict(PATHS,value_peer=args.value_peer_report)
 pins={p:sha(ROOT/p)for p in paths.values()};d={k:json.loads((ROOT/p).read_text())for k,p in paths.items()}
 input_count=count(d);prior_owned=input_count+30000;assert prior_owned<150000
 raw_read=sum((ROOT/p).stat().st_size for p in paths.values());assert raw_read+200000<33554432
 for k,v in d.items():
  if k.endswith('peer'):assert v['all_passed']
 assert d['values']['all_passed']and d['values']['source_peer_report_sha256']=='985e63285949c4503136288971c13494f898ae15e6e97325b0d17d5b87e6ec5c'
 # The already peered old split gives fixed canonical p1,b1,d,a0,a4.
 k=d['split']['known_values'];assert k
 p1,b1=parse(k['head1_rg_p1']),parse(k['head1_k2_b1']);cyc=parse(k['cyclic_top_d']);a0=parse(k['h0_k1_a0']);a4=parse(k['h4_k1_a4'])
 assert(p1,b1,cyc,a0,a4)==(Q(-55495,968256),Q(1331,9216),Q(-36,1681),Q(2,25),Q(-4248,210125))
 p2,b2=map(parse,d['values']['head29_A_B_rg_k2'][2:]);p3,b3=map(parse,d['values']['head30_A_B_rg_k2'][2:])
 a3=div(b3,b1);c3=div(add(p3,-mul(a3,p1)),cyc)
 as_=div(add(b2,-mul(Q(1,4),b3)),b1)
 cs=div(add(add(p2,-mul(Q(1,4),p3)),-mul(as_,p1)),cyc)
 # Old literal11 basis: five heads, rh3,z5h4,rh4,z5rh4,r²h4,r³h4.
 basis=[['h',j]for j in range(5)]+[['r*h3'],['z5*h4'],['r*h4'],['z5*r*h4'],['r²*h4'],['r³*h4']]
 K=[[Q(0)]*7 for _ in range(11)]
 for j,i in enumerate([6,8,9,10,5]):K[i][j]=Q(1)
 K[3][5]=1;K[1][5]=-a3;K[7][5]=-c3
 K[2][6]=1;K[3][6]=Q(-1,4);K[1][6]=-as_;K[7][6]=-cs
 # Known part of the actual edge. Its four unspecified entries occur only
 # in columns0,4, absent from K. They remain formal, never zero-filled.
 edge=[[Q(0)]*11 for _ in range(4)]
 edge[0][4]=cyc;edge[2][0]=a0;edge[2][4]=a4
 for j,p,b in[(1,p1,b1),(2,p2,b2),(3,p3,b3)]:edge[1][j]=p;edge[3][j]=b
 edge[1][7]=cyc
 assert all(K[i][j]==0 for i in[0,4]for j in range(7))
 assert mm(edge,K)==[[Q(0)]*7 for _ in range(4)]
 # Append the same accepted S=(h4,rh4,h0,h1). The inverse is explicit,
 # so no new rank or inverse search is performed.
 comp=[4,7,0,1];B=[row+[Q(i==j)for j in comp]for i,row in enumerate(K)]
 inv=[[Q(0)]*11 for _ in range(11)]
 for j,i in enumerate([6,8,9,10,5]):inv[j][i]=1
 inv[5][3]=1;inv[5][2]=Q(1,4);inv[6][2]=1
 inv[7][4]=1;inv[8][7]=1;inv[9][0]=1;inv[10][1]=1
 for row,coeffU,coeffS in[(8,c3,cs),(10,a3,as_)]:
  inv[row][3]=add(inv[row][3],coeffU);inv[row][2]=add(inv[row][2],add(mul(coeffU,Q(1,4)),coeffS))
 ident=[[Q(i==j)for j in range(11)]for i in range(11)]
 assert mm(B,inv)==ident==mm(inv,B)
 five=d['five'];assert five['matrix_roots']['five_formal_cycles_ELL0']==6
 monmap=[[[]for _ in range(11)]for _ in range(5)]
 for j in range(5):monmap[j][j]=[[[0]*8,'1']]
 for row,col,pairs in[(3,5,[(7,1)]),(4,6,[(5,1)]),(4,7,[(7,1)]),(4,8,[(5,1),(7,1)]),(4,9,[(7,2)]),(4,10,[(7,3)])]:
  exp=[0]*8
  for i,n in pairs:exp[i]=n
  monmap[row][col]=[[exp,'1']]
 nodes=[['external_actual_sigma56_by5','five',6,'ELL0_S_F9'],
  ['actual_Rbar_polynomial_matrix',monmap,'original8_z0_to_z6_r'],
  ['constant_QQ_matrix',enc(K)],['matrix_product',0,1],['matrix_product',3,2]]
 shapes=[[56,5],[5,11],[11,7],[56,11],[56,7]]
 data={'scope':'SOURCE actual seven barclasses of original zero-degree K/ellK and typed56x7 formal mod-ell cycles. Not full-R0 K lifts, not a full-K presentation or action matrix.',
  'actual_old11_basis':basis,'actual_Kbar11_by7':enc(K),'labels':['z5*h4','z5*r*h4','r²*h4','r³*h4','r*h3','U','S'],
  'actual_four_head_parameters':{'p2':str(p2),'b2':str(b2),'p3':str(p3),'b3':str(b3)},
  'actual_correction_coefficients':{'a3=b3/b1':str(a3),'c3=(p3-a3*p1)/d':str(c3),'aS=(b2-b3/4)/b1':str(as_),'cS=(p2-p3/4-aS*p1)/d':str(cs)},
  'edge_known_part4_by11':enc(edge),'unknown_edge_entries_preserved_formally':[[1,0,'p0'],[3,0,'b0'],[1,4,'p4'],[3,4,'b4']],
  'all28_edge_times_K_zero':True,'unknown_part_times_K_identically_zero_by_absent_rows0_4':True,
  'whole11_basis_K_then_S':enc(B),'explicit_whole11_inverse':enc(inv),'both121_inverse_products_verified':True,
  'source_operator_nodes':nodes,'source_operator_shapes':shapes,'actual_barcycles56_by7_root':4,
  'actual_old11_by_fivehead_polynomial_map5_by11':monmap,
  'cycle_and_descent_proof':'The accepted FIVE sigma6 is (I-nbar*mbar)E, so Vbar*sigma0 and full768 relation descent hold. Multiply by the exact old11 monomial map and the displayed K constants. The connecting edge vanishes; the accepted honest UCT identifies its full seven-dimensional kernel with K/ellK. The explicit11 inverse proves independence without discarding higherTor.',
  'literal_U9_namespace':'Same accepted source: U9 lastgA505938, oldphysicalf637790 only through acceptedtelescope. This packet changes no source clone/section.',
  'no_new_Physical_IFT_gauge_jet_or_edge_response':True,'no_original_H0_nil_cube_index_or_Kbar_action_claim':True,
  'input_sha256':pins,'code_sha256':sha(Path(__file__)),'maximum_exact_QQ_preop_bits':MAXPRE,'all_passed':True}
 data['before_QQ_SOURCE_preflight']={'owned_upper':prior_owned,'parsed_input_numeric':input_count,'raw_read':raw_read,'RAW_with_output_reserve':raw_read+200000,'committed_bits':384,'temporary_bits':1153,'seconds':60,'memory_bytes':2000000000}
 owned=count(d)+count(data)+count(K)+count(B)+count(inv)+count(edge)+count(nodes)+count(shapes)+8000;assert owned<150000
 data['owned_numeric_upper']=owned
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);p=out/'actual-seven-barbasis-typed56-cycle-operator-and-full11-unit-complement.json';p.write_text(json.dumps(data,indent=2)+'\n')
 rep={'id':STEM,'all_passed':True,'status':'SOURCE_ACTUAL_SEVEN_KBAR_BASIS_TYPED_BARCYCLES_NO_FULL_LIFTS',
  'input_sha256':pins,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/(STEM+'.md')),
  'evidence_sha256':{p.name:sha(p)},'maximum_exact_QQ_preop_bits':MAXPRE,'owned_numeric_upper':owned,'seconds':time.monotonic()-start}
 (out/'report.json').write_text(json.dumps(rep,indent=2)+'\n');print(json.dumps(rep));signal.alarm(0)
if __name__=='__main__':main()
