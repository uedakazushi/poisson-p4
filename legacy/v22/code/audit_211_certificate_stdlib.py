"""Independent stdlib audit of 211 triple certificate, including all E8 jets."""
import json,argparse
from pathlib import Path
BASE = Path(__file__).resolve().parents[1]
from fractions import Fraction as Q
from itertools import combinations,permutations
import audit_certificates_stdlib as A
R=A.R;parse=A.parse

def deriv(p,ix):
 out={}
 for m,c in p.p.items():
  if m[ix]:
   mm=list(m);mm[ix]-=1;out[tuple(mm)]=c*m[ix]
 return R(out,p.d)
def monic_key(e):
 if e.iszero():return e.key()
 return (e/e.p[min(e.p)]).key()
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--model',type=Path,default=(BASE/'previous_v21/data/tree_jacobi_models/211_triple.json'));ap.add_argument('--jets',type=Path,default=(BASE/'data/e8_jets/211_triple.json'));ap.add_argument('--certificate',type=Path,default=(BASE/'data/211_triple_certificate.json'));ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
 model=json.loads(args.model.read_text());jets=json.loads(args.jets.read_text());cert=json.loads(args.certificate.read_text());checks=[]
 def ck(ok,label):
  if not ok:raise AssertionError(label)
  checks.append(label)
 allcoeff=[];t=A.SY['t'];n=[A.SY[f'n{i}'] for i in range(3)]
 for ci,c in enumerate(jets['components']):
  pivot=c['pivot'];tangent=c['tangent'];normal=c['normal'];rr=[parse(e) for e in c['parametrization']]
  xx=rr.copy();xx[pivot]=R.c(1);xx[tangent]=t
  for j,ix in enumerate(normal):xx[ix]+=n[j]
  env=dict(zip([f'x{i}' for i in range(5)],xx));P=[[R.c(0) for j in range(5)] for i in range(5)]
  for ij,e in model['linear_model']['P_entries'].items():
   i,j=map(int,ij);P[i][j]=parse(e,env);P[j][i]=-P[i][j]
  p=[[P[i][j]-xx[i]*P[pivot][j]+xx[j]*P[pivot][i] for j in range(5)] for i in range(5)]
  H=[[A.normal_linear(p[tangent][normal[j]],k) for j in range(3)] for k in range(3)]
  br=[]
  for i,j in [(1,2),(2,0),(0,1)]:
   a,b=normal[i],normal[j];br.append(p[a][b]-deriv(rr[b],A.TI)*p[a][tangent]-deriv(rr[a],A.TI)*p[tangent][b])
  B=[[A.normal_linear(br[j],k) for j in range(3)] for k in range(3)]
  for i in range(3):
   for j in range(3):ck((H[i][j]-parse(c['H'][i][j])).iszero(),f'H{ci}{i}{j}');ck((B[i][j]-parse(c['B'][i][j])).iszero(),f'B{ci}{i}{j}')
  res=parse(c['R']);ck((sum((H[i][i] for i in range(3)),R.c(0))+res).iszero(),f'trace{ci}')
  for ii in combinations(range(3),2):
   for jj in combinations(range(3),2):allcoeff+=A.coeffs(B[ii[0]][jj[0]]*B[ii[1]][jj[1]]-B[ii[0]][jj[1]]*B[ii[1]][jj[0]],A.TI)
  HH=A.matmul(H,H);allcoeff+=A.coeffs(sum((HH[i][i] for i in range(3)),R.c(0))-Q(7,18)*res*res,A.TI)
  det=R.c(0)
  for perm in permutations(range(3)):
   sign=(-1)**sum(perm[i]>perm[j] for i in range(3) for j in range(i+1,3));det+=sign*H[0][perm[0]]*H[1][perm[1]]*H[2][perm[2]]
  allcoeff+=A.coeffs(det+res**3/36,A.TI)
  C=A.matmul([[H[i][j]+(res/2 if i==j else R.c(0)) for j in range(3)] for i in range(3)],B)
  for row in C:
   for e in row:allcoeff+=A.coeffs(e,A.TI)
 ck({monic_key(e) for e in allcoeff if not e.iszero()}=={monic_key(parse(e)) for e in jets['all_e8_equations']},'all_e8_coefficients')
 jetchecks=len(checks);env={k:parse(v) for k,v in cert['initial_substitution'].items()};expected=[parse(e,env) for e in model['jacobi_equations']+jets['all_e8_equations']]
 for ri,rr in enumerate(cert['rounds']):
  inp=[parse(e) for e in rr['input_equations']];ck(A.keys(inp)==A.keys(expected),f'inputs{ri}')
  for k,c in enumerate(rr['certificates']):
   expr=R.c(0)
   for ix,v in c['combination'].items():
    factor=parse(v);ck(factor.d==(0,0) and set(factor.p)<={A.ZERO},f'rational_scalar{ri}:{k}:{ix}');expr+=factor*inp[int(ix)]
   ck((expr-parse(c['target'])).iszero(),f'combination{ri}:{k}')
  new={k:parse(v) for k,v in rr['new_substitution'].items()}
  for k,l in enumerate(rr['linear_equations']):
   if l!='1':ck(parse(l,new).iszero(),f'substitution{ri}:{k}')
  expected=[parse(e,new) for e in rr['input_equations']]
 ck(cert['unit'] and any(c['target']=='1' for c in cert['rounds'][-1]['certificates']),'unit_certificate_present')
 report={'all_passed':True,'checks':len(checks),'jet_reconstruction_checks':jetchecks,'certificate_checks':len(checks)-jetchecks,'rounds':len(cert['rounds']),'parameter_denominators':False,'scope':'algebraic certificate and exact E8 jets verified; geometric exhaustiveness of initial substitution is a separate proof','implementation':'Python stdlib ast/Fraction sparse exact arithmetic; no CAS or sampling'}
 args.output.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)
if __name__=='__main__':main()
