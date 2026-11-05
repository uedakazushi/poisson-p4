#!/usr/bin/env python3
"""Restore exact K9270P02 active polynomials and operation cursor; never resume."""
import argparse,hashlib,json,struct
from fractions import Fraction
from pathlib import Path
VERSION='K9270P02'
def read_state(path):
 with Path(path).open('rb')as f:
  def get(fmt):
   n=struct.calcsize('<'+fmt);b=f.read(n);assert len(b)==n;return struct.unpack('<'+fmt,b)[0]
  def z(limit):
   n=get('H');sign=get('B');assert sign in[0,1]and n<=(limit+7)//8;b=f.read(n);assert len(b)==n
   if n:assert b[0]!=0
   value=int.from_bytes(b,'big');assert value.bit_length()<=limit
   return -value if sign else value
  assert f.read(8)==VERSION.encode();n=get('I');assert n<65536;header=json.loads(f.read(n));assert header['version']==VERSION;np=get('I');assert np<=512;polys={};total=0
  for _ in range(np):
   pid=get('I');assert pid not in polys;role=get('B');source=get('i');term=get('i');row=get('i');mask=get('i');complete=get('B');assert complete in[0,1];num=z(4096);den=z(4096);assert num>0 and den>0;factor=Fraction(num,den);assert factor.numerator==num and factor.denominator==den;nt=get('I');total+=nt;assert total<=135000;coefs=[];last=None
   for i in range(nt):
    mon=(get('i'),get('i'),get('i'));val=z(1024);assert val and(last is None or last<mon);coefs.append((mon,val));last=mon
   polys[pid]=dict(role=role,source=source,term=term,row=row,mask=mask,complete=bool(complete),factor=factor,coefficients=coefs)
  pending=None;flag=get('B');assert flag in[0,1]
  if flag:
   pending=dict(dst=get('i'),lhs=get('i'),rhs=get('i'),left=get('Q'),right=get('Q'));pending['key']=(get('q'),get('q'),get('q'));pending['scale']=z(16384);known=get('B');assert known in[0,1];pending['value_known']=bool(known);pending['uncommitted_value']=z(16384);assert pending['dst']in polys
  assert not f.read(1)
 for rr in header['matrix_ids']:
  assert all(i in polys for i in rr)
 assert all(i in polys for i in header['weight_ids'])
 if header['acc_id']>=0:assert header['acc_id']in polys
 return dict(header=header,polynomials=polys,pending=pending,total_terms=total)
def check_fixture(s):
 assert s['header']['stage']=='fixture'and s['header']['complete'];values=[]
 for p in s['polynomials'].values():values.append({m:Fraction(c)*p['factor']for m,c in p['coefficients']})
 expected=[{(1,0,0):Fraction(1,2),(0,1,0):Fraction(-1,3)},{(0,0,1):Fraction(1),(0,0,0):Fraction(1,5)},{(1,0,1):Fraction(1,2),(1,0,0):Fraction(1,10),(0,1,1):Fraction(-1,3),(0,1,0):Fraction(-1,15)},{(0,0,0):Fraction(7)}]
 assert values==expected and s['pending']is None
def check_pending(s):
 h=s['header'];ps=s['polynomials'];p=s['pending']
 if h['stage']!='fixture':
  assert 0<=h['P01_resume_pairs_completed']<=8094365
  if h['P01_term3_boundary_reached']:assert h['P01_resume_pairs_completed']==8094365
  if h['tail_id']>=0:assert h['tail_id']in ps and h['tail_start_term']in[0,4]and 0<=h['tail_terms_ready']<=10-h['tail_start_term']
 if p is None:return dict(all_passed=True,pending=False)
 if p['lhs']<0 or p['rhs']<0:return dict(all_passed=True,pending=True,kind='input record',full_input_binding_required=True)
 assert p['lhs']in ps and p['rhs']in ps and p['dst']in ps
 a=ps[p['lhs']];b=ps[p['rhs']];d=ps[p['dst']];i=p['left'];j=p['right'];assert i<len(a['coefficients'])and j<len(b['coefficients']);ma,ca=a['coefficients'][i];mb,cb=b['coefficients'][j]
 shift=-h['known_output_a_shift']if p['dst']==h['acc_id']else 0
 key=(ma[0]+mb[0]+shift,ma[1]+mb[1],ma[2]+mb[2]);assert key==p['key']
 scale=a['factor']*b['factor']/d['factor'];assert scale.denominator==1 and abs(p['scale'])==scale.numerator
 if p['value_known']:
  old=next((v for m,v in d['coefficients']if m==key),0);assert old+p['scale']*ca*cb==p['uncommitted_value']
 return dict(all_passed=True,pending=True,key_and_scale_verified=True,value_verified=p['value_known'],stale_scratch_ignored=not p['value_known'])
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('state',type=Path);p.add_argument('--fixture',action='store_true');p.add_argument('--check',action='store_true');a=p.parse_args();s=read_state(a.state)
 if a.fixture:check_fixture(s)
 print(json.dumps(dict(all_passed=True,version=VERSION,state_sha256=hashlib.file_digest(a.state.open('rb'),'sha256').hexdigest(),complete=s['header']['complete'],active_polynomials=len(s['polynomials']),exact_coefficient_terms=s['total_terms'],header=s['header'],pending=s['pending'],pending_cursor_check=check_pending(s),arithmetic_NOT_resumed=True),indent=2,default=str))
if __name__=='__main__':main()
