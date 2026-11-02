from pathlib import Path
import json,hashlib,time,argparse
from fractions import Fraction as Q
parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True);args=parser.parse_args();st=time.time();root=Path('.');main=Path('reports/S61_FLAG_ZERO12_ZERO_NUISANCE_FORMAL_20261102.json');note=Path('research/notes/S61_FLAG_ZERO12_ZERO_NUISANCE_FORMAL_PROOF_REVIEW_20261102.md');j=json.loads(main.read_text());assert j['all_passed'];hashes=dict(j['input_sha256'])
for p,h in hashes.items():assert hashlib.sha256(Path(p).read_bytes()).hexdigest()==h,p
parse=lambda rows:{tuple(m):Q(c)for m,c in rows}
F=list(map(parse,j['polynomial_residuals']));C=list(map(parse,j['t_squared_cofactors']));t=parse(j['homogenizer_polynomial']);out={}
def addprod(f,g,sign=1):
 for m,c in f.items():
  for n,d in g.items():
   k=tuple(a+b for a,b in zip(m,n));z=out.get(k,Q(0))+sign*c*d
   if z:out[k]=z
   elif k in out:del out[k]
for f,c in zip(F,C):addprod(f,c)
addprod(t,t,-1);assert not out
assert (j['slice_dimension'],j['full_fixed_beta_dimension'])==(1,3);assert j['t_nonzero_normal_form']!='0';assert j['orbit_minor']=='2';assert not j['relative_nuisance_nilpotence_claimed'];assert not j['general_beta_removal_claimed']
for p in [main,note,Path(__file__).resolve().relative_to(Path.cwd())]:hashes[str(p)]=hashlib.sha256(p.read_bytes()).hexdigest()
report=dict(id='S61_FLAG_ZERO12_ZERO_NUISANCE_FORMAL_PROOF_REVIEW_20261102',all_passed=True,review_scope='Independent proof/source review; independent Fraction sparse multiplication of stored t-squared cofactors. No original135 reconstruction or native/primary rerun.',accepted_claims=['full fixed-beta completed local dimension3','t is nonzero and square-zero at the specified zero-nuisance point'],cofactor_identity_separate_multiplication=True,input_hashes=hashes,proof_note=str(note),seconds=time.time()-st)
outp=args.output;outp.open('x').write(json.dumps(report,indent=2)+'\n');print('PASS',len(hashes),time.time()-st,hashlib.sha256(outp.read_bytes()).hexdigest())
