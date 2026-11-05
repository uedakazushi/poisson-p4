"""SOURCE support projection: leading H,L and relative s degree<=2; full5 keys."""
def inside(e,lead):
 return e[0]==lead[0]and e[1]==lead[1]and 0<=e[2]-lead[2]<=2
def input_project(exponents,lead,prospective):
 out=set()
 for e in exponents:
  assert len(e)==5 and tuple(e[:3])>=tuple(lead)
  if inside(e,lead)and e not in out:prospective(len(out)+1);out.add(tuple(e))
 return out
def copy_project(A,leadA,leadout,shift,prospective):
 assert tuple(leadout)==tuple(leadA[j]+shift[j]for j in range(3))
 out=set()
 for e in A:
  assert inside(e,leadA);v=tuple(e[j]+shift[j]for j in range(5))
  assert inside(v,leadout)
  if v not in out:prospective(len(out)+1);out.add(v)
 return out
def product_project(A,B,leadA,leadB,leadout,prospective,shift=(0,0,0,0,0)):
 assert tuple(leadout)==tuple(leadA[j]+leadB[j]+shift[j]for j in range(3))
 out=set();pairs=0
 for e in A:
  assert inside(e,leadA)
  for f in B:
   assert inside(f,leadB)
   if(e[2]-leadA[2])+(f[2]-leadB[2])>2:continue
   pairs+=1;v=tuple(e[j]+f[j]+shift[j]for j in range(5));assert inside(v,leadout)
   if v not in out:prospective(len(out)+1);out.add(v)
 return out,pairs
def sum_project(children,leads,leadout,prospective):
 assert tuple(leadout)==min(map(tuple,leads))
 out=set()
 for A,l in zip(children,leads):
  assert all(inside(e,l)for e in A)
  if tuple(l[:2])!=tuple(leadout[:2]):continue
  assert l[2]>=leadout[2]
  if l[2]-leadout[2]>2:continue
  for e in A:
   if inside(e,leadout)and e not in out:prospective(len(out)+1);out.add(e)
 return out
