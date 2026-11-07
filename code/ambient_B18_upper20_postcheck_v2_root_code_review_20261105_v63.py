"""ROOT full independent L90/denseOLD postchecker review; no post arithmetic."""
from pathlib import Path
import hashlib,json
R=Path(__file__).resolve().parents[1];sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest();pins={}
def load(p,h):
 q=R/p;assert sha(q)==h;pins[p]=h;d=json.loads(q.read_text());assert d['all_passed'];return d
stem='reports/AMBIENT_B18_UPPER20_';gate=stem+'POSTCHECK_INPUT_GATE_V2_20261105_v63/report.json'
g=load(gate,'be8b5dbca322484fb6f93f1f1da7d3f1db89a9403fdf30a750bf463d6e3e6ecc')
c=load(stem+'POSTCHECK_V2_PREPARED_CUSTODY_20261105_v63/report.json','f77060c0b00ec8cb4e9b61a8b8f4c5e125c6265c7d18ce270ba314631db7a029')
w='code/ambient_B18_upper20_independent_postcheck_v2_20261105_v63.py';wh='12e5adee5372368e6756d6a70d908912d8315909d29a5748aad1c855c4a2824c';assert sha(R/w)==c['worker_sha256']==wh;pins[w]=wh
for p,h in c['input_sha256'].items():assert sha(R/p)==h;pins[p]=h
assert g['derived_forecast57977']==57977 and g['private_numeric_cap']==60000 and g['RAW_cap']==4000000
out=R/'reports/AMBIENT_B18_UPPER20_POSTCHECK_V2_ROOT_CODE_REVIEW_20261105_v63';out.mkdir(exist_ok=False)
r=dict(id=out.name,all_passed=True,SOURCE_CODE_review_only_no_signed20_postcheck=True,input_sha256=pins,
 root_complete208line_math_and_allV2_control_diffs_reviewed=True,
 independent_algorithm=['L90=68pointfree+22formalpointL terms includes rawS0nonzero contributions beforeactualkernel restriction','pullGamma/D covectors throughknownN andmonicnu thencontractintrinsicSc; different from primary forwardNSc','full32source-mask vectors/basis18/nu1nu2zero andallGamma r2ell/Dquadraticdecompositions including rr once','FORMALraw16 values compared before15truefullpointkernel restrictions+10explicitSOURCEzero kernels','dense21component polynomialC0/typed multiplication and reverseFIRST contractions reconstruct1742factors/raw21partials','raw617+20raw563 carried once OLD/FIRST/outer; nativecoordinatesalreadycorrectedonce','all20fivepieces/upper2/Qsix exactmatched, fixedzero2x20/I20 certificate independentlychecked','remaining12opaque derivatives unqueried; no HF2/length/H0/nilconclusion'],
 prefix_and_resource_audit=['currentRz/nativeSc+kernelpoint/2piece/FIRSTrecipe/rawfour+FIRSTfv/corrected20vectors retainedinwork500','gc module no longer shadowed by OLDGamma gc alias; traceback/context/causeframes cleared+gc beforeexclusiveGUARD','R/N/G/nu/Sc retiredbeforedenseOLD; C/T/leaves retiredbeforeFIRST; no densephase overlap','source57977<60000/RAW4MB/3541063 unchanged; externalAS2GB andwall60 mandatory over55compute+cleanup+5guard'],
 committed_prior=354,temporary_prior=1063,whole_numeric_prior=57977,private_cap=60000,RAW_cap=4000000,
 one_fixed_read_only_signed20_postcheck_authorized=True,output_stem='AMBIENT_B18_UPPER20_INDEPENDENT_POSTCHECK_V2_20261105_v63',
 no_extra_cases_or_rank_search_or_new_physical_jets_authorized=True,code_sha256=sha(Path(__file__)))
p=out/'report.json';p.write_text(json.dumps(r,indent=2)+'\n');print(sha(p))
