"""ROOT direct three-variable factor pipeline SOURCE/code review; no factor call."""
from pathlib import Path
import hashlib,json,math
R=Path(__file__).resolve().parents[1];sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest();pins={}
def load(p,h):
 q=R/p;assert sha(q)==h;pins[p]=h;d=json.loads(q.read_text());assert d['all_passed'];return d
cp='reports/S61_TWO_TAIL_H_EXACT_QQ_CUSTODY_V2_20261105_v63/report.json'
c=load(cp,'239a327c413156523da1906895f4cd1958f253bd670684dd496c61f081564d01')
for p,h in c['runtime_input_sha256'].items():assert sha(R/p)==h;pins[p]=h
for p,h in c['input_sha256'].items():assert sha(R/p)==h;pins[p]=h
cpp='code/s61_two_tail_H_exact_QQ_factor_20261105_v63.cpp';driver='code/s61_two_tail_H_exact_QQ_entry_20261105_v63.py'
assert sha(R/cpp)=='9e0f779705fca66bc2fa6a39cdcdd9bbf29ce4a478a6212536265111881f8cd1'
assert sha(R/driver)==c['driver_sha256']=='31c2f20098468d96e41cbe4fe92cbf780691b7c0bef5623118981f52bf4b6ca1'
assert sha(R/c['native_executable'])==c['native_executable_sha256']=='c6249081922a8ed53749f1952a72f7e8dc067961140fd6296ea05d86edb55b11'
pins[cpp]=sha(R/cpp);pins[driver]=sha(R/driver);pins[c['native_executable']]=c['native_executable_sha256']
rows=[tuple(map(int,l.split()))for l in (R/c['fixed_H_input']).read_text().splitlines()[1:]];assert len(rows)==94
assert len({50*a+5*e+d for a in range(7)for e in range(10)for d in range(5)})==350
assert max(50*a+5*e+d for a,e,d,v in rows)==309
norm=sum(abs(v)for a,e,d,v in rows);assert norm<2**47
assert math.gcd(*(v for a,e,d,v in rows))==1
out=R/'reports/S61_TWO_TAIL_H_EXACT_QQ_ROOT_CODE_REVIEW_20261105_v63';out.mkdir(exist_ok=False)
r=dict(id=out.name,all_passed=True,SOURCE_CODE_only_no_factor_values=True,input_sha256=pins,
 manually_reviewed_complete_CPP_and_Python_driver=True,Kronecker_injective_factor_box350=True,Kronecker_input_degree=309,input_L1_bits=norm.bit_length(),
 factor_output357bit_bound_proof='Each primitive integer factor remains primitive under injective box substitution; degree of image factor<=309, Mahler measure<=input L2<=L1<2^47, binomial bound gives coefficient<2^356, conservatively357bits.',
 library_internal_factor_and_primitive_growth_OPAQUE_not_384bit_certified=True,native_AS_bytes=1900000000,whole_AS_bytes=2000000000,
 audit=['ONEdirect3variablefmpz_mpoly_factorSTATUSmust1; no univariatefactor decoding or newpointprime','nativeprimitivepositiveLEXfactors coeff357/term350/count11/aggregate3850/degree sums guarded; allinput94/40bit fixed','allnativephase/log/candidates retained; opaque internalstate unavailableonSTOP explicitlydisclosed','Pythoncanonicaljson/current/next/pending/cursors budget50000;stage actualpartialbits+nextfactorbits+18<=384 checked BEFOREmultiplication','mixedradix address50a+5e+d injective;componentdegree sums checked beforeaddition toavoidcarryalias','entireproduct comparedto94nonzerosandother350boxzeros;scalarunit restoration not assumed','1.5Mverificationpair/4MiB/shared45soft58driver60hard2GB unchanged; fallback/capraise/retry forbidden','multiple nonconstant exactfactors proveQQreducible; APIalone singlefactor doesnotprovegeometricirred; independentpostpeerneeded'],
 one_exact_factor_pipeline_authorized=True,output_stem='S61_TWO_TAIL_H_EXACT_QQ_ENTRY_20261105_v63',
 visible_identity_hard_bits=384,external_canonical_factor_bits=357,verification_pair_cap=1500000,accessible_owned_cap=50000,RAW_cap=4194304,
 independent_factor_identity_postpeer_required=True,no_factor_irreducibility_or_component_or_Fitting_conclusion_yet=True,code_sha256=sha(Path(__file__)))
p=out/'report.json';p.write_text(json.dumps(r,indent=2)+'\n');print(sha(p))
