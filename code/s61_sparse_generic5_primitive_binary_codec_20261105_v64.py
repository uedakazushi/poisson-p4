"""Uncompressed canonical primitive rational-polynomial integer encoding."""
from fractions import Fraction
MAGIC=b'Q5PRIM1\n'
def uint(n):
    assert isinstance(n,int)and n>=0
    out=bytearray()
    while n>=128:out.append((n&127)|128);n>>=7
    out.append(n);return bytes(out)
def sint(n):
    n=int(n);m=abs(n);raw=m.to_bytes((m.bit_length()+7)//8,'big')
    return uint(2*len(raw)+(n<0))+raw
def canonical(p):
    return {'content':str(Fraction(p['content'])),'shift':list(p['shift']),
            'terms':[[str(int(t[0])),*t[1:]]for t in p['terms']]}
def encode(p):
    q=Fraction(p['content']);assert q and len(p['shift'])==5 and p['terms']
    z=bytearray(sint(q.numerator)+sint(q.denominator))
    for n in p['shift']:z.extend(uint(2*n if n>=0 else-2*n-1))
    z.extend(uint(len(p['terms'])))
    for t in p['terms']:
        assert len(t)==6 and int(t[0])and min(t[1:])>=0
        z.extend(sint(int(t[0])))
        for n in t[1:]:z.extend(uint(n))
    return bytes(z)
def decode(data):
    at=0
    def u():
        nonlocal at
        n=0;s=0;start=at
        while True:
            assert at<len(data);v=data[at];at+=1;n|=(v&127)<<s
            if not(v&128):break
            s+=7
        assert data[start:at]==uint(n);return n
    def i():
        nonlocal at
        tag=u();n=tag//2;assert at+n<=len(data)
        raw=data[at:at+n];at+=n;assert not n or raw[0]
        value=int.from_bytes(raw,'big');assert value or not(tag&1)
        return-value if tag&1 else value
    num,den=i(),i();assert num and den>0
    q=Fraction(num,den);assert(q.numerator,q.denominator)==(num,den)
    shift=[]
    for _ in range(5):
        z=u();shift.append(z//2 if not z&1 else-(z//2)-1)
    n=u();assert n>0;terms=[]
    for _ in range(n):
        c=i();assert c;terms.append([str(c),*[u()for _ in range(5)]])
    assert at==len(data)
    return {'content':str(q),'shift':shift,'terms':terms}
