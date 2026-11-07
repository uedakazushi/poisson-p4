"""Complete leading-HL face metadata: exact extrema, positive cardinality bounds."""
from pathlib import Path
import hashlib,json,math,time
from s61_col53_projection_leaf_io_20261105_v64 import read_leaf,leaf_paths
from s61_sparse_generic5_source_view_profile_20261105_v64 import count_numeric
ROOT=Path(__file__).resolve().parents[1]
ID='S61_COL53_NEXT17_79_FULL_HL_FACE_SOURCE_20261105_v64'
FINAL='reports/S61_COL53_NEXT17_79_FINALIZED_SOURCE_20261105_v64'
VALUE='reports/S61_COL53_NEXT17_79_STREAMED_VALUE_20261105_v64'
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb')as f:
  for b in iter(lambda:f.read(65536),b''):h.update(b)
 return h.hexdigest()
def enc(d):return(json.dumps(d,separators=(',',':'))+'\n').encode()
def main():
 start=time.monotonic();out=ROOT/'reports'/ID;out.mkdir(exist_ok=False);pins={}
 rp=ROOT/FINAL/'report.json';assert sha(rp)=='548ebf5f4078707d6ee2bb6dd1867c7161c2d32bd553ebef4ce454600c7db877';r=json.loads(rp.read_text());pins[str(rp.relative_to(ROOT))]=sha(rp)
 gp=rp.parent/'compiled-full-SLP-and-fixed-projection-plan.json';assert sha(gp)==r['evidence_sha256'][gp.name];pins[str(gp.relative_to(ROOT))]=sha(gp)
 D=json.loads(gp.read_text());ns=D['nodes'];numeric=count_numeric(D);assert numeric==65701
 needed=set();children={};root=D['minor_B_root']
 def visit(k):
  if k in needed:return
  needed.add(k);n=ns[k];ch=[]
  if n['q']!='0':
   if n['op']in['copy','multiply']:ch=list(dict.fromkeys(n['args']))
   elif n['op']=='add':ch=list(dict.fromkeys(c for c in n['args']if ns[c]['q']!='0'and ns[c]['lead5'][:2]==n['lead5'][:2]))
  children[k]=ch
  for c in ch:visit(c)
 visit(root);order=sorted(needed);boxes={};cards={};peak=0;decodepeak=0;inputterms=0
 def own(transient=0,extra=0):
  nonlocal peak
  v=25000+numeric+len(needed)+len(order)+len(children)+sum(map(len,children.values()))+9*len(boxes)+transient+extra+4096
  peak=max(peak,v);assert v<150000,('metadata_only_owned',v)
 for k in order:
  own();n=ns[k];op=n['op'];ch=children[k];box=None;card=0
  if n['q']=='0':pass
  elif op=='input':
   own(16152)
   for f in leaf_paths(n):pins[f]=sha(ROOT/f)
   p,shift,dec=read_leaf(n);assert dec+2048<=16152;decodepeak=max(decodepeak,dec+2048)
   for t in p['terms']:
    e=tuple(p['shift'][j]+int(t[j+1])+shift[j]for j in range(5))
    if e[:2]!=tuple(n['lead5'][:2]):continue
    v=(e[2],e[3],e[4]);card+=1
    if box is None:box=[[a,a]for a in v]
    else:
     for j,a in enumerate(v):box[j][0]=min(box[j][0],a);box[j][1]=max(box[j][1],a)
   inputterms+=card;del p
  elif op=='copy':
   assert len(ch)==1;c=ch[0];assert boxes[c]is not None
   box=[[v+n['shift5'][j+2]for v in b]for j,b in enumerate(boxes[c])];card=cards[c]
  elif op=='multiply':
   a,b=n['args'];assert boxes[a]is not None and boxes[b]is not None
   box=[[boxes[a][j][v]+boxes[b][j][v]+n['shift5'][j+2]for v in range(2)]for j in range(3)];card=cards[a]*cards[b]
  else:
   assert op=='add'
   for c in ch:
    if boxes[c]is None:continue
    card+=cards[c]
    if box is None:box=[b[:]for b in boxes[c]]
    else:
     for j in range(3):box[j][0]=min(box[j][0],boxes[c][j][0]);box[j][1]=max(box[j][1],boxes[c][j][1])
  if box is not None:
   assert box[0][0]==n['lead5'][2]
   card=min(card,math.prod(hi-lo+1 for lo,hi in box));assert card>0
  else:assert card==0
  boxes[k]=box;cards[k]=card;own();assert time.monotonic()-start<50
 # Metadata of old committed cells only: no cell coefficients are recomputed.
 vp=ROOT/VALUE/'report.json';assert sha(vp)=='4b1ed0275771c6f13a8f250b14a71a125277d3b8e6f17fd5995142b8d4dc0c54';pins[str(vp.relative_to(ROOT))]=sha(vp)
 vr=json.loads(vp.read_text());mi=vp.parent/'node-values-index.jsonl';assert sha(mi)==vr['evidence_sha256'][mi.name];pins[str(mi.relative_to(ROOT))]=sha(mi)
 oldnodes=0;oldsegments=0
 with mi.open('rb')as f:
  for line in f:
   rec=json.loads(line);assert rec['node']in needed;oldnodes+=1;oldsegments+=rec['bytes']
 assert oldnodes==1884 and set(map(int,D['needed_node_degree']))<=needed
 details=[]
 for k in order:
  box=boxes[k];details.append({'node':k,'complete_HL_face_Sqt_box':box,'positive_support_cardinality_UPPER_not_exact':cards[k],'old_reusable_relative_s_cap':D['needed_node_degree'].get(str(k))})
 own(extra=4*len(details)+4096)
 rootbox=boxes[root];width=rootbox[0][1]-rootbox[0][0]
 d={'SOURCE_METADATA_only_no_support_convolution_or_signed_values':True,'full_leading_HL':ns[root]['lead5'][:2],'exact_positive_face_S_min_max':rootbox[0],'full_relative_s_degree':width,'complete_face_Sqt_coordinate_box':rootbox,'positive_root_support_UPPER_not_exact':cards[root],'full_HL_ancestor_nodes':len(needed),'old_nodes_exact_reusable':oldnodes,'new_ancestor_nodes':len(needed)-oldnodes,'old_committed_byteidentical_segments_reusable':oldsegments,'all_node_metadata':details,'old_prefix_values_NEVER_recomputed_or_mutated':True,'extension_rule':'retain old c_e under fixed qi for old s cap; compute only output exponents above that cap using complete required child bands; preserve SOURCE leads even if old prefix actualzero','no_complete_face_ZERO_or_NONZERO_claim':True,'no_SOURCE_profile_or_VALUE_GO_implied':True}
 (out/'complete-leading-HL-extrema-and-band-reuse.json').write_bytes(enc(d))
 for f in['code/s61_col53_next17_79_full_HL_face_source_20261105_v64.py','code/s61_col53_projection_leaf_io_20261105_v64.py','code/s61_sparse_generic5_primitive_binary_codec_20261105_v64.py','code/s61_sparse_generic5_source_view_profile_20261105_v64.py','code/s61_sparse_generic5_source_packed5_profile_v2_20261105_v64.py']:pins[f]=sha(ROOT/f)
 raw=sum((ROOT/f).stat().st_size for f in pins)+sum(p.stat().st_size for p in out.iterdir())+262144;assert raw<32*1024*1024
 report={'id':ID,'all_passed':True,'SOURCE_METADATA_only_no_new_chosen_VALUE':True,'complete_positive_HL_face_S_min_max':rootbox[0],'full_HL_ancestor_nodes':len(needed),'root_positive_cardinality_UPPER':cards[root],'metadata_owned_peak':peak,'decode_numeric_upper':decodepeak,'fresh_RAW_upper':raw,'full_F11_and_all_matrix_RHS_complements_preserved':True,'full_face_support_and_signed_result_STILL_UNRUN':True,'input_sha256':pins,'evidence_sha256':{p.name:sha(p)for p in out.iterdir()},'seconds':time.monotonic()-start}
 (out/'report.json').write_bytes(enc(report));print('report_sha256',sha(out/'report.json'));print(json.dumps({k:v for k,v in report.items()if k not in('input_sha256','evidence_sha256')}))
if __name__=='__main__':main()
