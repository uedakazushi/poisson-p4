#!/usr/bin/env python3
"""One bounded full Ka13 B21(root) and deterministic up-to35 minor gate."""
import argparse,hashlib,itertools,json,math,resource,signal,subprocess,sys,time
from pathlib import Path
from s61_nonsquare_raw_H_B_K13_named_saturation_20261104_v63 import Ring,Stop,pa,pm
ROOT=Path(__file__).resolve().parents[1]
ID='S61_NONSQUARE_RAW_H_B_K13_B21_ROOT_TEST_20261104_v63'
sys.set_int_max_str_digits(0)
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def ip(p):return {int(i):int(v)for i,v in p}
def ser(p):return [[i,str(v)]for i,v in sorted(p.items(),reverse=True)]
def main(out,work):
    start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000))
    signal.signal(signal.SIGALRM,lambda *_:(_ for _ in()).throw(TimeoutError('60s total cap')));signal.alarm(60)
    out=out.resolve();work=work.resolve();out.mkdir(parents=True,exist_ok=False);work.mkdir(parents=True,exist_ok=False)
    db=ROOT/'reports/S61_NONSQUARE_RAW_H_B_K13_D_EUCLID_20261104_v63'
    bb=ROOT/'reports/S61_NONSQUARE_RAW_H_B_K13_MIXED23_SCALARS_20261104_v63'
    mb=ROOT/'reports/S61_NONSQUARE_RAW_H_B_MINORS_COST_20261104_v63'
    bridge=ROOT/'reports/S61_NONSQUARE_RAW_H_B_K13_SOURCE_BRIDGE_20261104_v63_REVIEW'
    hb=ROOT/'reports/S61_NONSQUARE_RAW_H_B_K13_H4_INPUT_20261104_v63'
    ub=ROOT/'reports/S61_NONSQUARE_RAW_H_B_K13_NAMED_SATURATION_20261104_v63'
    assert sha(db/'report.json')=='004096adc77ddee31067f3f6bf58e804fc6af2057a7552d43462c88b2bb55f45'
    dr=read(db/'report.json');assert dr['all_passed']and dr['completed_gcd_degree']==1
    for name in ['result.json','checkpoint.json']:assert sha(db/name)==dr['evidence_sha256'][name]
    br=read(bb/'report.json');assert br['all_passed'];bp=bb/'native-certificate.json.B21.json'
    assert sha(bp)==br['evidence_sha256'][bp.name]
    assert read(bridge/'report.json')['all_passed']
    us=read(ub/'report.json');assert us['status']=='clopen_U_unit_cover_retained'and us['remaining_degree']==4
    result=read(db/'result.json');cp=read(db/'checkpoint.json');ring=Ring(ip(cp['g13']));zero=ring.c(0);one=ring.c(1);s=ring.dec(result['root_s'])
    Bc=[[ring.dec(v)for v in row]for row in read(bp)['entries']];assert len(Bc)==21 and all(len(p)==4 for p in Bc)
    cost=read(mb/'certificate.json');assert cost['literal_all21_unit_preconditioning_reconstructions_verified']
    pre=dict(input_B_field_coefficients=sum(bool(x[0])for p in Bc for x in p),input_B_QQ_terms=sum(len(x[0])for p in Bc for x in p),
        coefficient_algebra='QQ[a]/g13; no field or irreducibility assumption',B_shape=[7,3],input_d_degree_upper=3,
        Horner_ring_multiplication_count=63,minor_count_upper=35,minor_order='lexicographic3subsets ofrows0..6',
        determinant_ring_multiplications_per_minor=9,unit_query_bound=1,root_terms=len(s[0]),
        root_numerator_bits=max(abs(v).bit_length()for v in s[0].values()),root_denominator_bits=s[1].bit_length(),
        hard_seconds=60,hard_AS_bytes=2000000000,hard_artifact_bytes=32*1024*1024,hard_explicit_terms=150000,
        save_all_Bvalues_before_minors=True,no_automatic_next_minor_after_nonzero_nonunit=True)
    assert pre['input_B_field_coefficients']==69 and pre['input_B_QQ_terms']==789
    (out/'preflight.json').write_text(json.dumps(pre,indent=2)+'\n')
    sourcebind=dict(source_B_coefficient_artifact=str(bp.relative_to(ROOT)),source_B_coefficient_SHA=sha(bp),
        all84_source_images_accepted_by=str((bridge/'report.json').relative_to(ROOT)),
        original_prepared_B_source=str((mb/'certificate.json').relative_to(ROOT)),
        literal_preconditioning='rowunit_i*Bprepared_ij*colunit_j=Boriginal_ij*D_j; row/colunits only QQ,a,d,F,B,C,J. Thus maximalminor ideals coincide on FULLU.',
        original_rows14=cost['actual_original_block_rows14'],original_columns14=cost['actual_original_block_columns14'],
        integer_column_multipliers=cost['integer_column_multipliers_QQ_units'],
        whole_source_cover='G unit modulo actualI from ALLQ/L gates; H,h1,k1,m1 belongtoI bysource norm-cofactor identities; Ka13 e=t bridge and gcd d-s imply Ka13[d]/(d-s) localizedatU surjects onto R/I. Sourceg/k reduction uses H4, not removal of an unknown factor.',
        FULLU_unit_on_whole_H4_cover_from_named_saturation=True,
        all35zero_scope='Only candidate retained here; alloriginal22/sourceopen implications are left for separate explicit acceptance, not inferred as existence.')
    (out/'source-binding.json').write_text(json.dumps(sourcebind,indent=2)+'\n')
    records=[];values=[];minorrecords=[];unit=None;status='pending';maxbits=0
    def pbits(x):return max([x[1].bit_length()]+[abs(v).bit_length()for v in x[0].values()])
    def guard(extra=0):
        total=sum(p.stat().st_size for p in out.iterdir()if p.is_file())
        if total+extra+2*1024*1024>32*1024*1024:raise Stop('artifact_size_guard',dict(existing=total,extra=extra))
        if sum(len(x[0])for x in values)+sum(len(ip(x['value']['numerator']))for x in minorrecords)>150000:raise Stop('term_guard',{})
    def save(status,extra=None):
        d=dict(status=status,g13=ser(ring.g),root_s=ring.enc(s),minors=minorrecords,unit=unit,extra=extra,seconds=time.monotonic()-start)
        payload=json.dumps(d,indent=2)+'\n';existing=sum(p.stat().st_size for p in out.iterdir()if p.is_file()and p.name!='checkpoint.json')
        if existing+len(payload.encode())>31*1024*1024:raise Stop('checkpoint_size_guard',{})
        (out/'checkpoint.json').write_text(payload)
    try:
        for index,p in enumerate(Bc):
            value=p[3];steps=[];qs=[]
            for j in [2,1,0]:
                previous=value;product=ring.mul(previous,s);value=ring.add(product,p[j]);assert ring.sub(value,product)==p[j]
                steps.append(dict(d_coefficient=j,previous=ring.enc(previous),times_s=ring.enc(product),after_add=ring.enc(value)))
                if j>0:qs.append(value)
            quotient=[qs[1],qs[0],p[3]]
            assert quotient[0]==ring.dec(steps[1]['after_add'])and quotient[1]==ring.dec(steps[0]['after_add'])
            records.append(dict(entry=index,row=index//3,column=index%3,input_coefficients=[ring.enc(v)for v in p],
                value=ring.enc(value),Horner_steps=steps,division_by_d_minus_s_quotient=[ring.enc(v)for v in quotient],
                identity='B_entry(d)=(d-s)*quotient+B_entry(s)',all_coefficient_stages_checked=True))
            values.append(value);maxbits=max(maxbits,pbits(value),*(pbits(v)for v in quotient));guard()
        (out/'B21-values.json').write_text(json.dumps(dict(g13=ser(ring.g),root_s=ring.enc(s),entry_order='rowmajor7x3',entries=records),indent=2)+'\n')
        save('all21_Bvalues_saved_before_minor_search')
        A=[values[3*i:3*i+3]for i in range(7)]
        for rows in itertools.combinations(range(7),3):
            b=[A[i]for i in rows]
            co=[ring.sub(ring.mul(b[1][1],b[2][2]),ring.mul(b[1][2],b[2][1])),
                ring.sub(ring.mul(b[1][0],b[2][2]),ring.mul(b[1][2],b[2][0])),
                ring.sub(ring.mul(b[1][0],b[2][1]),ring.mul(b[1][1],b[2][0]))]
            det=ring.add(ring.sub(ring.mul(b[0][0],co[0]),ring.mul(b[0][1],co[1])),ring.mul(b[0][2],co[2]))
            minorrecords.append(dict(rows7=list(rows),value=ring.enc(det),Laplace2x2_cofactors=[ring.enc(v)for v in co],zero=not bool(det[0])))
            maxbits=max(maxbits,pbits(det),*(pbits(v)for v in co));save('minor_checked');guard()
            if not det[0]:continue
            content=math.gcd(*(abs(v)for v in det[0].values()));p={i:v//content for i,v in det[0].items()}
            m=max(p);B=max(abs(v).bit_length()for v in p.values());bg=max(abs(v).bit_length()for v in ring.g.values())
            hadamard=26*(max(B,bg)+5)+math.factorial(26).bit_length()+64
            estimate=(m+15)*(200+(hadamard*30103+99999)//100000);guard(estimate)
            save('before_first_nonzero_minor_unit_certificate',dict(rows=list(rows),predicted_unit_bytes=estimate))
            cpp=ROOT/'code/s61_nonsquare_raw_H_B_K13_named_unit_20261104_v63.cpp';exe=work/'unit-check';inp=work/'unit.txt';answer=out/'minor-unit-certificate.json'
            with inp.open('w')as f:
                for z in[ring.g,p]:
                    f.write(str(len(z))+'\n')
                    for i,v in sorted(z.items(),reverse=True):f.write(str(i)+' '+str(v)+'\n')
            with(out/'compile.log').open('w')as log:
                subprocess.run(['bash','-lc','source /workspace/cloud-setup/activate.sh\nc++ -O2 -std=c++17 "$1" -o "$2" $(pkg-config --cflags --libs flint)','compile',str(cpp),str(exe)],stdout=log,stderr=subprocess.STDOUT,check=True,timeout=max(1,55-(time.monotonic()-start)))
            subprocess.run(['bash','-lc','source /workspace/cloud-setup/activate.sh\nexec "$@"','Bminor-unit',str(exe),str(inp),str(answer)],check=True,timeout=max(1,57-(time.monotonic()-start)))
            u=read(answer);unit=dict(rows7=list(rows),input=ring.enc(det),primitive_numerator=ser(p),numerator_content=str(content),certificate=u)
            if u['unit']:
                S,T=ip(u['S']),ip(u['T']);R=int(u['R']);assert pa(pm(S,p),pm(T,ring.g))=={0:R}
                inverse=ring.make({i:v*det[1]for i,v in S.items()},content*R);assert ring.mul(det,inverse)==one
                unit['inverse']=ring.enc(inverse);status='unit_minor_actual_Bbad_EMPTY_given_source_gates'
            else:status='proper_nonunit_minor_base_components_retained'
            save(status);break
        else:status='all35_minors_zero_candidate_ONLY';save(status)
    except Stop as e:status=e.status;save(status,e.data)
    complete=status in ['unit_minor_actual_Bbad_EMPTY_given_source_gates','proper_nonunit_minor_base_components_retained','all35_minors_zero_candidate_ONLY']
    empty=status=='unit_minor_actual_Bbad_EMPTY_given_source_gates'
    proof=dict(status=status,whole_Ka13_algebra_used=True,source_map_explained='source-binding.json',
        if_unit='A unit prepared minor maps to a unit in R/I, but every actualBmaxminor is zero there; preconditioning factors are FULLU units. Thus R/I=0 as a scheme, after the accepted source-cover gates.',
        actual_Bbad_empty_given_accepted_gates=empty,all_FULLU_B_rank3_given_accepted_gates=empty,
        if_all35zero='No existence assertion here: this packet only retains the commonroot candidate. Alloriginal22/sourceopen gates require explicit final acceptance.',
        no_global_fullcubic57_quartic239_genericPoisson_exhaustion_inference=True,
        independent_review_of_this_new_gate_pending=True)
    (out/'scheme-proof.json').write_text(json.dumps(proof,indent=2)+'\n')
    paths=[Path(__file__),ROOT/'code/s61_nonsquare_raw_H_B_K13_named_saturation_20261104_v63.py',ROOT/'code/s61_nonsquare_raw_H_B_K13_named_unit_20261104_v63.cpp',ROOT/f'research/notes/{ID}.md',db/'report.json',db/'result.json',db/'checkpoint.json',bb/'report.json',bp,mb/'report.json',mb/'certificate.json',bridge/'report.json',ub/'report.json',ub/'scheme-proof.json',ROOT/'reports/S61_NONSQUARE_RAW_H_B_G0_REMOVED_SCHEME_BOUND_20261104_v63_REVIEW/report.json']
    report=dict(id=ID,all_passed=complete,status=status,all21_Bvalues_saved_before_minors=(out/'B21-values.json').exists(),
        minors_computed=len(minorrecords),unit_queries=int(unit is not None),first_nonzero_rows=unit['rows7']if unit else None,
        actual_Bbad_empty_given_accepted_gates=empty,all_FULLU_B_rank3_given_accepted_gates=empty,
        field_integer_operations=ring.operations,max_coefficient_bits_seen=maxbits,whole_Ka13_retained=not empty,
        no_other_solver_Norm_point_prime_GB_or_automatic_next_job=True,independent_review_pending=True,
        input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},evidence_sha256={p.name:sha(p)for p in out.iterdir()if p.is_file()},seconds=time.monotonic()-start)
    assert sum(p.stat().st_size for p in out.iterdir()if p.is_file())<32*1024*1024
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');signal.alarm(0)
    print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--work',type=Path,required=True);a=p.parse_args()
    try:main(a.output,a.work)
    except Exception as e:
        if a.output.exists():(a.output/'failure.json').write_text(json.dumps(dict(all_passed=False,error=str(e),B21_saved=(a.output/'B21-values.json').exists(),no_auto_retry=True),indent=2)+'\n')
        raise
