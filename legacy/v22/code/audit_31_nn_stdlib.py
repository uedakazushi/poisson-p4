"""Exact independent NN certificate audit over Q or Q(lambda), no CAS."""
import argparse,json,hashlib
from pathlib import Path
BASE = Path(__file__).resolve().parents[1]
from fractions import Fraction as Q
import univariate_fraction_stdlib as U
G=U.G

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--model',type=Path,default=(BASE/'previous_v21/data/tree_jacobi_models/31_external.json'));ap.add_argument('--jets',type=Path,default=(BASE/'data/e8_jets/31_external.json'));ap.add_argument('--certificate-directory',type=Path,default=(BASE/'data/31_nn'));ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
 model=json.loads(args.model.read_text());jets=json.loads(args.jets.read_text());source=model['jacobi_equations']+jets['all_e8_equations'];checks=[];reports=[]
 def ck(ok,label):
  if not ok:raise AssertionError(label)
  checks.append(label)
 cases=[('31_C_NN_k3_certificate',None,[]),('31_NN_k1_II_certificate','z21',[('z21',(0,1))]),('31_NN_k1_IA_certificate','z17',[('z17',(0,1)),('69*z17+28',(28,69))]),('31_NN_k1_IB_certificate','z17',[('z17',(0,1)),('3*z17+1',(1,3)),('621*z17**2+87*z17+52',(52,87,621))]),('31_NN_k1_IA_rational_certificate',None,[]),('31_NN_k1_IB_quadratic_field_certificate','theta',[]),('31_NN_k1_IB_rational_certificate',None,[])]
 for name,field,allowed in cases:
  p=args.certificate_directory/(name+'.json');d=json.loads(p.read_text());U.set_modulus((13511,0,1) if field=='theta' else None);start=len(checks);alias={field:U.SY['lam']} if field else {};used=set()
  def parse(s,env=None):
   v=U.parse(s,env,alias);used.update(U.supported_denominator(v.d,allowed));return v
  env={k:parse(v) for k,v in d['initial_substitution'].items()};expected=[parse(e,env) for e in source]
  for ri,rr in enumerate(d['rounds']):
   inp=[parse(e) for e in rr['input_equations']];ck(U.keys(expected)==U.keys(inp),name+f':inputs{ri}')
   for k,c in enumerate(rr['certificates']):
    expr=G.c(0)
    for ix,v in c['combination'].items():
     scalar=parse(v);U.uni(scalar.p);expr+=scalar*inp[int(ix)]
    ck((expr-parse(c['target'])).iszero(),name+f':combination{ri}:{k}')
   new={k:parse(v) for k,v in rr['new_substitution'].items()}
   for k,l in enumerate(rr['linear_equations']):
    if l!='1':ck(parse(l,new).iszero(),name+f':substitution{ri}:{k}')
   expected=[parse(e,new) for e in rr['input_equations']]
  if d['unit']:
   ck(any(c['target']=='1' for c in d['rounds'][-1]['certificates']),name+':unit')
  elif name=='31_C_NN_k3_certificate':
   ck(name=='31_C_NN_k3_certificate',name+':expected_nonunit_case')
   rem=[parse(e) for e in d['remaining_equations']];ck(U.keys(rem)==U.keys(expected),name+':remaining_set')
   terminal=json.loads((args.certificate_directory/'31_NN_k3_final_certificate.json').read_text());steps=terminal['steps'];z=U.SY
   ck((rem[steps[0]['equation_index']]-z['z16']**2).iszero(),'k3_square_equation')
   ck((rem[steps[1]['equation_index']]-(z['z10']+5*z['z16']*z['z9'])/5).iszero(),'k3_linear_equation')
   ck((rem[steps[2]['equation_index']]-z['z9']*(10*z['z12']-1)/10).iszero(),'k3_nonzero_parameter_equation')
   finalenv={'z16':G.c(0),'z10':G.c(0),'z12':G.c(Q(1,10))}
   for i,ix in enumerate(terminal['univariate_equation_indices']):ck((parse(d['remaining_equations'][ix],finalenv)-parse(terminal['univariate_equations'][i])).iszero(),f'k3_univariate{i}')
   ck((sum((parse(a)*parse(b) for a,b in zip(terminal['bezout_coefficients'],terminal['univariate_equations'])),G.c(0))-1).iszero(),'k3_Bezout_identity')
  else:
   ck(name=='31_NN_k1_IB_rational_certificate',name+':expected_family_case')
   rem=[parse(e) for e in d['remaining_equations']];ck(U.keys(rem)==U.keys(expected),name+':remaining_set')
   z=U.SY;ck((12*rem[19]-z['z11']+3*z['z18']*z['z20']).iszero(),name+':first_graph_relation')
   ck((4*rem[19]-4*rem[7]-z['z1']+3*z['z20']**2).iszero(),name+':second_graph_relation')
   graph={'z11':3*z['z18']*z['z20'],'z1':3*z['z20']**2}
   for i,e in enumerate(d['remaining_equations']):ck(parse(e,graph).iszero(),name+f':all_remaining_on_graph{i}')
   finalenv={k:parse(v,graph) for k,v in d['final_substitution'].items()}|graph
   for i,e in enumerate(source):ck(parse(e,finalenv).iszero(),name+f':full_source_on_family{i}')
   P=[[G.c(0) for j in range(5)] for i in range(5)]
   for ij,e in model['linear_model']['P_entries'].items():
    i,j=map(int,ij);P[i][j]=parse(e,finalenv);P[j][i]=-P[i][j]
   b=[]
   for omitted in range(5):
    i,j,k,l=[v for v in range(5) if v!=omitted];b.append(((-1)**omitted)*(P[i][j]*P[k][l]-P[i][k]*P[j][l]+P[i][l]*P[j][k]))
   F=sum((z[f'x{i}']*b[i]/5 for i in range(5)),G.c(0))
   ck(all(m[U.A.I['x2']]+m[U.A.I['x3']]>=2 for m in F.p),name+':F_in_plane_ideal_squared')
   coef={}
   for m,c in F.p.items():
    if tuple(m[U.A.I[f'x{i}']] for i in range(5))==(0,3,0,2,0):
     mm=list(m)
     for i in range(5):mm[U.A.I[f'x{i}']]=0
     coef[tuple(mm)]=c
   ck((G(coef,F.d)+Q(1,120)).iszero(),name+':F_nonzero_coefficient')
  if field=='theta':
   v=parse(d['initial_substitution']['z17']);ck((621*v*v+87*v+52).iszero(),name+':quadratic_exception_parameter')
   ck(not v.iszero(),name+':nonzero_exception_parameter')
  outcome='unit' if d['unit'] else ('closed-point contradiction via square and Bezout' if name=='31_C_NN_k3_certificate' else 'actual nonnormal quintic family with singular plane')
  reports.append({'name':name,'checks':len(checks)-start,'outcome':outcome,'denominator_factors':sorted(used),'passed':True});print(json.dumps(reports[-1]),flush=True)
 U.set_modulus(None)
 partition=json.loads((args.certificate_directory/'31_NN_partition_certificate.json').read_text());pstart=len(checks)
 for si,st in enumerate(partition['stages']):
  extra={k:U.parse(v) for k,v in partition['stage_substitutions'][si].items()};env={k:U.parse(v,extra) for k,v in partition['geometric_initial_substitution'].items()}|extra
  expect=[U.parse(e,env) for e in source];inp=[U.parse(e) for e in st['input_equations']];ck(U.keys(expect)==U.keys(inp),f'partition_inputs{si}')
  for ci,c in enumerate(st['certificates']):
   total=sum((U.parse(v)*inp[int(ix)] for ix,v in c['combination'].items()),G.c(0));ck((total-U.parse(c['target'])).iszero(),f'partition_target{si}:{ci}')
 reports.append({'name':'31_NN_partition_certificate','checks':len(checks)-pstart,'target_count':sum(len(st['certificates']) for st in partition['stages']),'passed':True})
 report={'all_passed':True,'checks':len(checks),'certificates':reports,'source_sha256':{label:{'filename':p.name,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for label,p in [('model',args.model),('jets',args.jets)]},'scope':'full NN algebraic certificate coverage conditional on geometric reductions; k3 uses closed-point square-zero and geometric z9!=0; every generic denominator exception covered; quadratic field checked modulo theta²+13511 for both conjugates; rational family has a singular plane','implementation':'stdlib ast/Fraction sparse arithmetic overQ(lambda), exact Euclidean polynomial normalization, no CAS/no sampling'}
 args.output.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)
if __name__=='__main__':main()
