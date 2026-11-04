#!/usr/bin/env python3
"""One source-bound d-Euclid over the full Ka13 algebra; all LC units certified."""
import argparse,hashlib,json,math,resource,signal,subprocess,sys,time
from pathlib import Path
from s61_nonsquare_raw_H_B_K13_named_saturation_20261104_v63 import Ring,Stop,pa,pm
ROOT=Path(__file__).resolve().parents[1]
ID='S61_NONSQUARE_RAW_H_B_K13_D_EUCLID_20261104_v63'
sys.set_int_max_str_digits(0)
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def ip(p):return {int(i):int(v)for i,v in p}
def ser(p):return [[i,str(v)]for i,v in sorted(p.items(),reverse=True)]
def main(out,work):
    start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000))
    signal.signal(signal.SIGALRM,lambda *_:(_ for _ in()).throw(TimeoutError('60s total cap')));signal.alarm(60)
    out=out.resolve();work=work.resolve();out.mkdir(parents=True,exist_ok=False);work.mkdir(parents=True,exist_ok=False)
    source=ROOT/'reports/S61_NONSQUARE_RAW_H_B_K13_GK_INPUT_20261104_v63'
    assert sha(source/'report.json')=='632e472246977ecf915eb6eabb42c8949caee98f05f1ab2ba2f5c514920e3282'
    sr=read(source/'report.json');assert sr['all_passed']
    assert sha(source/'source-input.json')==sr['evidence_sha256']['source-input.json']
    c=read(source/'source-input.json');g=ip(c['g13']);ring=Ring(g);zero=ring.c(0);one=ring.c(1)
    cpp=ROOT/'code/s61_nonsquare_raw_H_B_K13_named_unit_20261104_v63.cpp';exe=work/'unit-check'
    with (out/'compile.log').open('w')as log:
        subprocess.run(['bash','-lc','source /workspace/cloud-setup/activate.sh\nc++ -O2 -std=c++17 "$1" -o "$2" $(pkg-config --cflags --libs flint)','compile',str(cpp),str(exe)],stdout=log,stderr=subprocess.STDOUT,check=True,timeout=max(1,55-(time.monotonic()-start)))
    def trim(p):
        p=list(p)
        while p and not p[-1][0]:p.pop()
        return p
    def encp(p):return [ring.enc(x)for x in p]
    def degree(p):return len(p)-1
    def padd(a,b,sign=1):
        p=[zero for _ in range(max(len(a),len(b)))]
        for i,x in enumerate(a):p[i]=ring.add(p[i],x)
        for i,x in enumerate(b):p[i]=ring.add(p[i],x if sign==1 else ring.neg(x))
        return trim(p)
    def pscale(p,s):return trim([ring.mul(x,s)for x in p])
    def pmul(a,b):
        if not a or not b:return []
        p=[zero for _ in range(len(a)+len(b)-1)]
        for i,x in enumerate(a):
            for j,y in enumerate(b):p[i+j]=ring.add(p[i+j],ring.mul(x,y))
        return trim(p)
    inputs=[trim([ring.dec(x)for x in c['H4']])]+[trim([ring.dec(x)for x in r['H4_remainder']])for r in c['gk']]
    assert list(map(degree,inputs))==[4,3,3] and inputs[0][-1]==one
    preflight=dict(coefficient_algebra='QQ[a]/g13, no field or etale assumption',input_degrees=[4,3,3],
        input_QQ_terms=sum(len(x[0])for p in inputs for x in p),input_max_bits=max(max([x[1].bit_length()]+[abs(v).bit_length()for v in x[0].values()])for p in inputs for x in p),
        unit_query_bound=9,hard_seconds=60,hard_AS_bytes=2000000000,hard_artifact_bytes=32*1024*1024,
        hard_explicit_polynomial_terms=150000,stop_on_any_proper_nonunit=True,no_other_B_evaluation_or_minor=True)
    (out/'preflight.json').write_text(json.dumps(preflight,indent=2)+'\n')
    units=[];events=[];nodes=[];polys=[];cofactor_slp=[];cofs=[];cache={};maxbits=preflight['input_max_bits'];stage='inputs';result=None
    def scalarbits(x):return max([x[1].bit_length()]+[abs(v).bit_length()for v in x[0].values()])
    def check_size(extra=0):
        live=sum(len(x[0])for p in polys for x in p)+sum(len(x[0])for cs in cofs for p in cs for x in p)
        if live>150000:raise Stop('explicit_support_guard',dict(terms=live))
        existing=sum(p.stat().st_size for p in out.iterdir()if p.is_file())
        if existing+extra+2*1024*1024>32*1024*1024:raise Stop('artifact_size_guard',dict(existing=existing,extra=extra))
    def save(status,extra=None):
        nonlocal maxbits
        maxbits=max([maxbits]+[scalarbits(x)for p in polys for x in p]+[scalarbits(x)for cs in cofs for p in cs for x in p])
        d=dict(status=status,stage=stage,g13=ser(g),units=units,events=events,polynomial_nodes=nodes,
            source_cofactor_slp=cofactor_slp,max_coefficient_bits_seen=maxbits,extra=extra,seconds=time.monotonic()-start)
        payload=json.dumps(d,indent=2)+'\n';existing=sum(p.stat().st_size for p in out.iterdir()if p.is_file()and p.name!='checkpoint.json')
        if len(payload.encode())+existing>31*1024*1024:raise Stop('checkpoint_size_guard',dict(payload_bytes=len(payload.encode()),existing=existing))
        (out/'checkpoint.json').write_text(payload)
    def node(p,cs,record,slp):
        j=len(nodes);polys.append(p);cofs.append(cs);nodes.append(dict(id=j,degree=degree(p),coefficients=encp(p),**record));cofactor_slp.append(dict(id=j,**slp));check_size();return j
    for j,(name,p)in enumerate(zip(['H4','gR','kR'],inputs)):
        cs=[[one]if i==j else []for i in range(3)]
        node(p,cs,dict(kind='source',name=name),dict(op='source_basis',index=j))
    def inv(tag,x):
        if not x[0]:raise Stop('zero_LC_not_inverted',dict(tag=tag))
        key=(tuple(sorted(x[0].items())),x[1])
        if key in cache:
            i,y=cache[key];events.append(dict(kind='unit_reuse',tag=tag,unit=i));return y,i
        if len(units)>=9:raise Stop('unit_query_guard',dict(tag=tag,completed=len(units)))
        content=math.gcd(*(abs(v)for v in x[0].values()));p={i:v//content for i,v in x[0].items()}
        m=max(p);B=max(abs(v).bit_length()for v in p.values());bg=max(abs(v).bit_length()for v in g.values())
        hadamard_bits=26*(max(B,bg)+5)+math.factorial(26).bit_length()+64
        predicted=(m+15)*(200+(hadamard_bits*30103+99999)//100000);check_size(predicted)
        i=len(units);save('before_LC_unit_certificate',dict(tag=tag,primitive_bits=B,predicted_certificate_bytes=predicted))
        if m==0:
            assert abs(p[0])==1;u=dict(gcd=[[0,'1']],S=[[0,str(p[0])]],T=[],R='1',unit=True)
        else:
            inp=work/f'unit-{i}.txt';answer=work/f'unit-{i}.json'
            with inp.open('w')as f:
                for z in[g,p]:
                    f.write(str(len(z))+'\n')
                    for k,v in sorted(z.items(),reverse=True):f.write(str(k)+' '+str(v)+'\n')
            subprocess.run(['bash','-lc','source /workspace/cloud-setup/activate.sh\nexec "$@"','Ka13-LC-unit',str(exe),str(inp),str(answer)],check=True,timeout=max(1,55-(time.monotonic()-start)))
            u=read(answer)
        rec=dict(tag=tag,input=ring.enc(x),primitive_numerator=ser(p),numerator_content=str(content),native_or_constant_certificate=u,predicted_certificate_bytes=predicted)
        if not u['unit']:
            units.append(rec);save('proper_base_component_retained',dict(tag=tag,gcd=u['gcd']));raise Stop('proper_base_component_retained',dict(tag=tag,gcd=u['gcd']))
        S,T=ip(u['S']),ip(u['T']);R=int(u['R']);assert pa(pm(S,p),pm(T,g))=={0:R}
        y=ring.make({j:v*x[1]for j,v in S.items()},content*R);assert ring.mul(x,y)==one
        rec['inverse']=ring.enc(y);units.append(rec);cache[key]=(i,y);save('LC_unit_certificate_saved',dict(tag=tag,unit=i));return y,i
    def monic(i):
        p=polys[i]
        if not p or p[-1]==one:return i
        v,u=inv('node_'+str(i)+'_LC',p[-1]);q=pscale(p,v);assert q[-1]==one
        cs=[pscale(z,v)for z in cofs[i]]
        assert pscale(q,p[-1])==p
        j=node(q,cs,dict(kind='unit_normalization',source=i,unit=u,old_LC=ring.enc(p[-1])),dict(op='scale',source=i,scalar=ring.enc(v)))
        save('unit_normalization_checked');return j
    def division(a,b):
        A=polys[a];B=polys[b];assert B and B[-1]==one
        R=list(A);Q=[zero for _ in range(max(0,len(A)-len(B)+1))]
        while R and len(R)>=len(B):
            power=len(R)-len(B);c=R[-1];Q[power]=ring.add(Q[power],c)
            for j,x in enumerate(B):R[power+j]=ring.sub(R[power+j],ring.mul(c,x))
            assert not R[-1][0];R=trim(R)
        Q=trim(Q);assert padd(pmul(Q,B),R)==A
        cs=[padd(x,pmul(Q,y),-1)for x,y in zip(cofs[a],cofs[b])]
        j=node(R,cs,dict(kind='monic_division_remainder',dividend=a,divisor=b,quotient=encp(Q),identity_checked=True),dict(op='subtract_product',left=a,polynomial=encp(Q),right=b))
        save('monic_division_checked');return j
    def gcd(a,b):
        while polys[b]:
            b=monic(b);r=division(a,b);a,b=b,r
        return monic(a)
    status='pending';extra={}
    try:
        save('inputs_checked');stage='gcd_H4_gR';first=gcd(0,1);stage='gcd_with_kR';result=gcd(first,2)
        D=polys[result];assert D and D[-1]==one
        stage='final_source_division_and_Bezout_checks';divisions=[]
        for i in range(3):
            rr=division(i,result);assert not polys[rr]
            divisions.append(dict(source=i,zero_remainder_node=rr))
        combined=[]
        for coeff,sourcepoly in zip(cofs[result],inputs):combined=padd(combined,pmul(coeff,sourcepoly))
        assert combined==D
        resultdata=dict(gcd_node=result,gcd_degree=degree(D),gcd_coefficients=encp(D),
            source_names=['H4','gR','kR'],source_Bezout_cofactors=[encp(z)for z in cofs[result]],
            source_Bezout_all_coefficients_checked=True,all_three_source_divisions=divisions,
            all_three_division_remainders_zero=True,all_inverted_LC_certified_units=True,
            no_field_or_irreducibility_or_etale_assumption=True)
        if degree(D)==1:
            s=ring.neg(D[0]);resultdata.update(root_s=ring.enc(s),equivalent_ideal='(H4,gR,kR)=(d-s) in Ka13[d]',commonroot_scheme='Ka13 via d=s',all_original_Bconditions_NOT_evaluated=True)
            status='complete_linear_commonroot_Ka13_retained'
        elif degree(D)==0:status='complete_unit_commonroot_cover_empty'
        else:status='complete_higher_degree_commonroot_cover_retained'
        payload=json.dumps(resultdata,indent=2)+'\n';check_size(len(payload.encode()));(out/'result.json').write_text(payload)
        extra=dict(gcd_degree=degree(D),result_node=result);save(status,extra)
    except Stop as e:status=e.status;extra=e.data;save(status,extra)
    complete=status.startswith('complete_')
    paths=[Path(__file__),ROOT/'code/s61_nonsquare_raw_H_B_K13_named_saturation_20261104_v63.py',cpp,ROOT/f'research/notes/{ID}.md',source/'report.json',source/'source-input.json',ROOT/'reports/S61_NONSQUARE_RAW_H_B_K13_G_BASE_GATE_20261104_v63/report.json']
    report=dict(id=ID,all_passed=complete,status=status,completed_gcd_degree=extra.get('gcd_degree'),
        unit_queries=len(units),all_required_unit_certificates_saved=complete,polynomial_nodes=len(nodes),field_integer_operations=ring.operations,max_coefficient_bits_seen=maxbits,
        no_automatic_Bmatrix_substitution_minor_Norm_newpoint_prime_or_GB=True,whole_Ka13_retained=True,
        actual_Bbad_existence_or_empty_NOT_inferred=status!='complete_unit_commonroot_cover_empty',independent_review_pending=True,
        input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},evidence_sha256={p.name:sha(p)for p in out.iterdir()if p.is_file()},seconds=time.monotonic()-start)
    assert sum(p.stat().st_size for p in out.iterdir()if p.is_file())<32*1024*1024
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');signal.alarm(0)
    print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--work',type=Path,required=True);a=p.parse_args()
    try:main(a.output,a.work)
    except Exception as e:
        if a.output.exists():(a.output/'failure.json').write_text(json.dumps(dict(all_passed=False,error=str(e),checkpoint_preserved=(a.output/'checkpoint.json').exists(),no_auto_retry=True),indent=2)+'\n')
        raise
