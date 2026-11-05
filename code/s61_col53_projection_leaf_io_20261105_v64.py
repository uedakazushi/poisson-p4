"""Accepted physical SOURCE leaf transforms for the compiled scalar graph."""
from pathlib import Path
from fractions import Fraction as Q
import hashlib,json
from s61_sparse_generic5_primitive_binary_codec_20261105_v64 import decode
from s61_sparse_generic5_source_view_profile_20261105_v64 import count_numeric
ROOT=Path(__file__).resolve().parents[1]
OLD='reports/S61_SPARSE_GENERIC5_CONSTRAINED83_NINTH_COL56_STAGED_SOURCE_V4_20261105_v64'
def leaf_paths(n):
 s=n['source'];s=s.get('legacy_source',s)
 if 'binding'in s:return[OLD+'/'+s['binding']['scope']]
 if'thin_payload_file'in s:return[s['thin_payload_file']]
 if'cell_file'in s:return[s['cell_file']]
 if'file'in s:return[s['file']]
 return[]
def read_leaf(n):
 src=n['source'];s=src.get('legacy_source',src);shift=[0]*5;scale=Q(1);decoded=0
 if 'binding'in s:
  b=s['binding'];doc=json.loads((ROOT/OLD/b['scope']).read_text());decoded=count_numeric(doc)
  v=(doc['common_operands']if'common_operands'in doc else doc['row_operands'])[b['label']]
  p=v['original_literal'];scale=Q(v['QQ_transform_scalar']);shift=v['transform_monomial_exponent5']
 elif'thin_payload_file'in s or('file'in s and'offset'in s):
  file=s.get('thin_payload_file',s.get('file'))
  with(ROOT/file).open('rb')as f:f.seek(s['offset']);body=f.read(s['bytes'])
  assert hashlib.sha256(body).hexdigest()==s['payload_sha256'];p=decode(body);decoded=count_numeric(p)
 elif'cell_file'in s or'file'in s:
  file=s.get('cell_file',s.get('file'));body=(ROOT/file).read_bytes();assert body[:9]==b'Q5CELL1\n\1';p=decode(body[9:]);decoded=count_numeric(p)
  if'legacy_source'in src:
   profile=src['literal_profile'];low=[p['shift'][j]+min(t[j+1]for t in p['terms'])for j in range(5)]
   shift=[profile['exponent_box5'][j][0]-low[j]for j in range(5)]
   scale=Q(profile['known_QQ_content'])/Q(p['content'])
   expected=(Q(353894400000),[1,0,0,0,1])if'col56-ninth-cleared'in file else(Q(1),[0]*5)
   assert(scale,shift)==expected
 elif s.get('QQ_monomial'):
  profile=src.get('literal_profile');shift0=s.get('exponent5',[b[0]for b in profile['exponent_box5']]if profile else[0]*5)
  p={'content':n['q'],'shift':shift0,'terms':[['1',0,0,0,0,0]]};decoded=count_numeric(p)
 else:
  assert s.get('known_zero')or s.get('known_source_zero')
  p={'content':'0','shift':[0]*5,'terms':[]};decoded=count_numeric(p)
 assert Q(p['content'])*scale==Q(n['q'])
 return p,shift,decoded
