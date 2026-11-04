#!/usr/bin/env python3
"""READ-ONLY source/type/height preflight for the Ka13 full-cubic sufficient input."""
import argparse,hashlib,json,math,resource,signal,sys,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_RAW_H_B_K13_FULL_CUBIC_PREFLIGHT_20261104_v63'
sys.set_int_max_str_digits(0)
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main(out):
    start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out=out.resolve();out.mkdir(parents=True,exist_ok=False)
    cb=ROOT/'reports/S61_NONSQUARE_CUBIC_ZW25_UNITS_20261104_v63';cr=read(cb/'report.json');cp=read(cb/'checkpoint.json')
    assert sha(cb/'checkpoint.json')==cr['evidence_sha256']['checkpoint.json']and cp['completed_pivots']==8
    rr,cc=cp['remaining_rows'],cp['remaining_columns'];assert (len(rr),len(cc))==(17,24)
    B={(1,0,1):Q(50),(0,2,1):Q(18),(0,1,1):Q(-9),(0,1,0):Q(60),(0,0,1):Q(-2),(0,0,0):Q(-40)}
    monoms=set();scales=[];classes=[];types={};terms=nonzero=maxentry=coeffbits=0;flatlines=['17 24']
    for i in rr:
        for j in cc:
            p,q=cp['matrix25x32'][i][j];qd={tuple(m):Q(v)for m,v in q};assert all(v and min(m)>=0 and len(m)==3 for m,v in qd.items())
            if len(qd)==1 and (0,0,0)in qd:kind='1';scale=qd[(0,0,0)]
            elif len(qd)==1 and (1,0,0)in qd:kind='a';scale=qd[(1,0,0)]
            else:kind='B';scale=qd[(1,0,1)]/50;assert qd=={m:scale*v for m,v in B.items()}
            assert scale and scale.denominator==1
            flatlines.append(str({'1':0,'a':1,'B':2}[kind])+' '+str(scale.numerator)+' '+str(len(p)))
            for m,v in p:
                v=Q(v);assert v and v.denominator==1 and len(m)==3 and min(m)>=0;monoms.add(tuple(m));coeffbits=max(coeffbits,abs(v.numerator).bit_length())
                flatlines.append(str(v.numerator)+' '+' '.join(map(str,m)))
            types[kind]=types.get(kind,0)+1;scales.append(abs(scale.numerator));terms+=len(p);nonzero+=bool(p);maxentry=max(maxentry,len(p))
            classes.append(dict(source_row25=i,target_column32=j,denominator_unit=kind,denominator_integer_scale=str(scale),numerator_terms=len(p),numerator_coefficients_SHA256=hashlib.sha256(json.dumps(p,separators=(',',':')).encode()).hexdigest()))
    assert terms==84253 and nonzero==354 and len(monoms)==1113
    flat=ROOT/'reports/S61_NONSQUARE_THETA_P80_CUBIC_GENERIC_20261104_v63/source17x24-input.txt'
    assert flat.read_text()=='\n'.join(flatlines)+'\n'
    kb=ROOT/'reports/S61_NONSQUARE_RAW_H_B_K13_ACTUAL_SOURCE_SCHEME_20261104_v63';mp=read(kb/'monomial-map.json');nodes=[json.loads(s)for s in(kb/'monomial-nodes.jsonl').read_text().splitlines()]
    known={tuple(m):j for m,j in mp['monomial_nodes']}
    for axis,p in enumerate(mp['powers']):
        for k,j in enumerate(p):m=[0,0,0];m[axis]=k;known[tuple(m)]=j
    candidates=[m for m in known if m!=(0,0,0)];dp={(0,0,0):0};dec={};bounds=[max(m[j]for m in monoms)for j in range(3)]
    assert bounds==[13,23,6]
    for a in range(14):
        for e in range(24):
            for d in range(7):
                m=(a,e,d)
                if m==(0,0,0):continue
                if m in known:dp[m]=1;dec[m]=(m,(0,0,0));continue
                options=[(1+dp[(a-i,e-j,d-k)],(i,j,k),(a-i,e-j,d-k))for i,j,k in candidates if i<=a and j<=e and k<=d]
                value,left,right=min(options);dp[m]=value;dec[m]=(left,right)
    used=set()
    def factors(m):
        if m in known:used.add(known[m]);return [known[m]]
        left,right=dec[m];return factors(left)+factors(right)
    recipes=[dict(exponent=list(m),accepted_node_factors=factors(m))for m in sorted(monoms)]
    K=max(len(r['accepted_node_factors'])for r in recipes);assert K==3
    base_den=math.lcm(*(int(nodes[j]['value']['denominator'])for j in used));base_mag=max(max(abs(int(v)).bit_length()for a,v in nodes[j]['value']['numerator'])-int(nodes[j]['value']['denominator']).bit_length()+1 for j in used if nodes[j]['value']['numerator'])
    dg=ROOT/'reports/S61_NONSQUARE_RAW_H_B_K13_D_EUCLID_20261104_v63';g={int(i):int(v)for i,v in read(dg/'checkpoint.json')['g13']};lc=abs(g[13])
    companion=Q(sum(abs(v)for i,v in g.items()if i<13),lc);cbit=max(0,companion.numerator.bit_length()-companion.denominator.bit_length()+1)
    D=base_den**K*lc**(12*(K-1));mag=K*base_mag+(K-1)*(8+24*cbit);pmag=mag+coeffbits+math.ceil(math.log2(maxentry))
    rawnum=D.bit_length()+pmag;source_scale_lcm=math.lcm(*scales)
    sourceunits=read(kb/'source-open-units.json');Bvalue=sourceunits['FULLU_product_all6_factors'][3]
    Bden=int(Bvalue['denominator']);Bmag=max(abs(int(v)).bit_length()for i,v in Bvalue['numerator'])-Bden.bit_length()+1
    # Clearing by cLCM*a*B avoids all parameter inverses. Multiplication by B
    # costs <=12 leading-coefficient factors, and multiplication by a <=1.
    Dclear=D*Bden*lc**13;clearmag=pmag+Bmag+8+24*cbit+cbit+source_scale_lcm.bit_length()
    clear_num=Dclear.bit_length()+clearmag
    def bytesbound(nb,db):return 408*(13*(nb*30103//100000+40)+(db*30103//100000+40)+100)
    height=dict(accepted_base_nodes_used=len(used),max_known_factors=K,base_denominator_LCM=str(base_den),base_denominator_LCM_bits=base_den.bit_length(),
        base_coefficient_abs_log2_upper=base_mag,g_leading_coefficient_bits=lc.bit_length(),companion_l1_norm_power2_upper=cbit,
        uniform_monomial_denominator='baseD^3*lc(g13)^24',uniform_monomial_denominator_bits=D.bit_length(),
        source_QQ_integer_coefficient_bits=coeffbits,max_source_entry_terms=maxentry,
        numerator_matrix_common_denominator_bits_upper=D.bit_length(),numerator_matrix_numerator_bits_upper=rawnum,numerator_matrix_RAW_JSON_bytes_upper=bytesbound(rawnum,D.bit_length()),
        all_constant_denominators_LCM=str(source_scale_lcm),all_constant_denominators_LCM_bits=source_scale_lcm.bit_length(),
        B_value_denominator_bits=Bden.bit_length(),B_coefficient_abs_log2_upper=Bmag,
        cleared_matrix_common_denominator_bits_upper=Dclear.bit_length(),cleared_matrix_numerator_bits_upper=clear_num,
        cleared_matrix_RAW_JSON_bytes_upper=bytesbound(clear_num,Dclear.bit_length()),
        proof='Degree<13 products have degree<=24, requiring<=12 numeric lc(g13) pseudodivisions. l1 multiplication norm grows by<=13^2*M^24, M=max(1,sum|lowerg/lc|). Existing<=3factor recipes bound every monomial. Final a multiplication needs only1 lc factor. Canonical reduction cannot enlarge common denominators.',
        no_field_or_irreducibility_assumption=True)
    assert height['cleared_matrix_RAW_JSON_bytes_upper']<30*1024*1024
    binding=dict(source_checkpoint=str((cb/'checkpoint.json').relative_to(ROOT)),source_rows25=rr,target_columns32=cc,
        existing_flat_universal_input=str(flat.relative_to(ROOT)),flat_input_all84253_coefficients_byte_matched=True,flat_input_bytes=flat.stat().st_size,
        all408_scalar_types='numerator/denominator QQ[a,e,d] arrays; NO POWER or unspecialized source-variable atoms',denominator_class_counts=types,entries=classes,
        known_unit8_minor=cp['selected_minor_unit_determinant'],
        specialization_legal_on_entire_Ka13='All6 FULLU factors have units because source-open-units.json gives explicit productinverse. Thus existing8 unit operations and denominator1/a/B survive.',
        cleared_matrix_recipe='For denominator c*u, Cij=(cLCM/c)*(a*B/u)*numeratorij, u=1,a,B. Thus C=cLCM*a*B*M is a unit scalar multiple overKa13. Do NOT expand the parameter multipliers.',
        proof_scope='rank17 of thisactual17x24 impliesfullcubic57 byaccepted projectedcriterion/product41/Cdelta line. Smaller projectedrank is INCONCLUSIVE. Do NOTimportTheta-specificheavy6 or10x43criterion ontoKa13.',
        coefficient_evaluation_NOT_performed=True,rank_Gaussian_minors_NOT_performed=True,
        source_matrix_coefficients_are_external_pinned_inputs=True)
    plan=dict(only_if_separately_authorized='ONE60s/2GB/32MiB/150k actualcoefficientpreparation; no automaticrankjob.',
        materialize='408 canonicalKa13values ofclearedmatrix beforeanyrank. Reuse272acceptednodes; factor1113monomials<=3. Streamrows/checkpoint, no duplicationofallmonomialcofactors.',
        live_output_slots_upper=408*13,monomial_cache_slots_upper=1113*13,
        source_term_count_below150k=True,conservative_cleared_matrix_RAW_bytes=height['cleared_matrix_RAW_JSON_bytes_upper'],
        if_any_guard='Savecompletedrows. No automaticbudgetincrease or alternateinput.',
        future_distinct_rank_method='Aftercoeffs+denominatorpunitgate: restrictionofscalars Ka13^24->Ka13^17 gives221x312 QQmatrix. A positive fullrank coefficientminor mod an explicitlyauthorizedexistingprime provesQQsurjectivity onWHOLEKa13 withoutassumingitfield, hencefull57. No suchprimeor ranktest here.',
        if_modular_rank_deficient='Inconclusive; no full57counterclaim or automaticotherprime.',
        no_POWER207_newpoint_prime_or_coefficient_evaluation=True)
    for name,data in [('source-and-types.json',binding),('monomial-factor-recipes.json',recipes),('height-and-output-bounds.json',height),('next-coefficient-plan.json',plan)]:
        (out/name).write_text(json.dumps(data,indent=2)+'\n')
    paths=[Path(__file__),ROOT/f'research/notes/{ID}.md',flat,cb/'report.json',cb/'checkpoint.json',kb/'report.json',kb/'monomial-map.json',kb/'monomial-nodes.jsonl',kb/'source-open-units.json',dg/'report.json',dg/'checkpoint.json',ROOT/'reports/S61_NONSQUARE_CUBIC_ZW25_UNITS_20261104_v63_REVIEW/report.json',ROOT/'reports/S61_NONSQUARE_THETA_FULL_CUBIC_SOURCE_PREFLIGHT_20261104_v63/report.json',ROOT/'reports/S61_NONSQUARE_CUBIC_ZW25_AUGMENTED_INPUT_20261104_v63/input.json',ROOT/'reports/S61_NONSQUARE_PRODUCT41_FULL_OPEN_20261103_v63/report.json']
    report=dict(id=ID,all_passed=True,source_terms=terms,source_nonzero_entries=nonzero,source_monomials=len(monoms),known_base_nodes_used=len(used),factor_count_upper=K,
        denominator_class_counts=types,cleared_matrix_RAW_bytes_bound=height['cleared_matrix_RAW_JSON_bytes_upper'],
        no_Ka13_source_matrix_coefficient_evaluation_rank_point_prime_Gaussian_or_GB=True,
        input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},evidence_sha256={p.name:sha(p)for p in out.iterdir()if p.is_file()},seconds=time.monotonic()-start)
    assert sum(p.stat().st_size for p in out.iterdir()if p.is_file())<32*1024*1024
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
