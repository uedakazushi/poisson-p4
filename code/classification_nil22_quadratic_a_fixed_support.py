#!/usr/bin/env python3
"""Complete18variable torus-fixed truecut; no normaljet truncation."""
from pathlib import Path
import hashlib,json,os,subprocess,time
import sympy as s
from classification_rank3_mixed_component_review import COLUMNS,PAIRS,jacobi
ROOT=Path(__file__).resolve().parents[1]
def encode(poly,variables):
 terms=[]
 for ex,c in s.Poly(poly,*variables).terms():
  co=str(c.p)if c.q==1 else '('+str(c.p)+'/'+str(c.q)+')';terms.append('*'.join([co]+[str(v)+('^'+str(k)if k!=1 else'')for v,k in zip(variables,ex)if k]))
 return '('+'+'.join(terms)+')'if terms else'0'
def main():
 started=time.monotonic();path=ROOT/'reports/CLASSIFICATION_NIL22_QUADRATIC_A_WEIGHT_FRAME_20261107.json';frame=json.loads(path.read_text())['frame'];assert frame['true_normal_relative_weights']==[0] and frame['actual_full_source_F_rank']==48;block=next(b for b in frame['complete_blocks']if b['relative_weight']==0);cols=block['implicit_columns'];rows=block['implicit_rows'];assert len(cols)==17;zz=s.symbols('z0:17');n=s.Symbol('n');variables=zz+(n,);xx=s.symbols('h x1 x2 x3 x4');P=s.zeros(5)
 for ij,f in frame['canonical_point'].items():
  i,j=map(int,ij.split(','));P[i,j]=s.sympify(f)
 for col,z in zip(cols,zz):
  ij,ex=COLUMNS[col];P[ij]+=z*s.prod(x**k for x,k in zip(xx,ex))
 for col,a in frame['true_normal_vectors'][0].items():
  ij,ex=COLUMNS[int(col)];P[ij]+=n*s.Rational(a)*s.prod(x**k for x,k in zip(xx,ex))
 for i,j in PAIRS:P[j,i]=-P[i,j]
 sparse={ij:dict(s.Poly(s.expand(P[ij]),*xx).terms())for ij in PAIRS if P[ij]!=0};JJ=jacobi(sparse).applyfunc(s.expand);indices=[i for i,f in enumerate(JJ)if f!=0];equations=[JJ[i]for i in indices];assert len(equations)<=46;inverse=s.Matrix([[s.Rational(v)for v in r]for r in block['implicit_inverse']]);gauss=(inverse*s.Matrix([JJ[r]for r in rows])).applyfunc(s.expand)
 for i,f in enumerate(gauss):assert s.Poly(f,*variables).homogeneous_component(1).as_expr()==zz[i] if hasattr(s.Poly(f,*variables),'homogeneous_component')else s.expand(sum(c*s.prod(v**k for v,k in zip(variables,ex))for ex,c in s.Poly(f,*variables).terms()if sum(ex)==1))==zz[i]
 work=ROOT/'work/nil22-quadratic-a-fixed-support';work.mkdir(parents=True,exist_ok=True);script=work/'complete.sing';output=work/'certificate.txt';out=ROOT/'reports/CLASSIFICATION_NIL22_QUADRATIC_A_FIXED_SUPPORT_20261107.json'
 if output.exists()or out.exists():raise FileExistsError('Outputalreadyexists')
 lines=['ring r=0,('+','.join(map(str,variables))+'),ds;','ideal I='+','.join(encode(f,variables)for f in equations)+';','matrix T;','ideal G=liftstd(I,T);','int j; int k; int good; int chosen=0; intvec le;','for(j=1;j<=size(G);j++){ le=leadexp(G[j]); good=1; for(k=1;k<=17;k++){if(le[k]!=0){good=0;}} if(good==1 && le[18]>0){chosen=j;break;}}','write("'+str(output)+'","DIM "+string(dim(G)));','write("'+str(output)+'","CHOSEN "+string(chosen));','if(chosen>0){write("'+str(output)+'","POLY "+string(G[chosen]));for(j=1;j<=size(I);j++){write("'+str(output)+'","LIFT "+string(j)+" "+string(T[j,chosen]));}}','write("'+str(output)+'","COMPLETE");','exit;'];script.write_text('\n'.join(lines)+'\n');env=os.environ.copy();env['LD_LIBRARY_PATH']='/workspace/cloud-setup/prefix/usr/lib/x86_64-linux-gnu:/workspace/cloud-setup/prefix/usr/lib'+(':'+env['LD_LIBRARY_PATH']if env.get('LD_LIBRARY_PATH')else'');env['SINGULARPATH']='/workspace/cloud-setup/prefix/usr/share/singular/LIB';env['SINGULAR_PROCS_DIR']=str(ROOT/'work/nonsquare-v62/singular-runtime/prefix/usr/lib/x86_64-linux-gnu/singular/MOD');native_started=time.monotonic();run=subprocess.run([str(ROOT/'work/nonsquare-v62/singular-runtime/prefix/usr/bin/Singular'),'-q',str(script)],capture_output=True,text=True,env=env,timeout=85)
 assert run.returncode==0 and not any(t.lstrip().startswith('?')for t in run.stdout.splitlines()),run.stdout+'\n'+run.stderr;assert output.exists()and output.read_text().splitlines()[-1]=='COMPLETE';records=output.read_text().splitlines();dim=int(next(t[4:]for t in records if t.startswith('DIM ')));chosen=int(next(t[7:]for t in records if t.startswith('CHOSEN ')));result=dict(task_id='E_multiple_hyperplane',all_passed=True,input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in [Path(__file__),path,ROOT/'code/classification_rank3_mixed_component_review.py']},producer_bytes=Path(__file__).read_text(),complete_torus_fixed_truecut_variables=list(map(str,variables)),complete_all350_restricted_equations=list(map(str,JJ)),nonzero_equation_indices=indices,complete_local_ideal_generators=list(map(str,equations)),all17_Gauss_linear_polynomials=list(map(str,gauss)),Gauss_full350_rows=rows,Gauss_inverse=[[str(v)for v in r]for r in inverse.tolist()],local_standard_initial_dimension=dim,normal_purepower_member_found=chosen>0,
  mathematical_reduction='AlltrueNweights0, Tstartrivialonone normal. TheT-invariant102cut hasinvertible101implicitJacobian. UniquenessofequivariantformalIFTforcesall84nonzero-weightimplicitfunctions identically0. Hence thecompleteformaltruecut equals this EXACT18variable/<=46equation fixedcut; no normaldegree truncation.',scope='Ifpurepowerliftverified, all18 lowesthomogeneousforms arez_i,n^k sofixedcutlocaldim0. Krull48cuts+actualsourceF48 implieswhole7sourceclosurefull23component. Ifdim1/none, no componenthoodinferred.',native_input_sha256=hashlib.sha256(script.read_bytes()).hexdigest(),native_seconds=time.monotonic()-native_started,seconds=time.monotonic()-started)
 if chosen:
  envparse={str(v):v for v in variables};poly=s.Poly(s.sympify(next(t[5:]for t in records if t.startswith('POLY ')).replace('^','**'),locals=envparse),*variables);lifts=[s.Poly(0,*variables)for _ in equations]
  for line in records:
   if line.startswith('LIFT '):_,index,text=line.split(' ',2);lifts[int(index)-1]=s.Poly(s.sympify(text.replace('^','**'),locals=envparse),*variables)
  epol=[s.Poly(f,*variables)for f in equations];assert sum((l*f for l,f in zip(lifts,epol)),s.Poly(0,*variables))==poly;degree=min(sum(ex)for ex,c in poly.terms());initial=[(ex,c)for ex,c in poly.terms()if sum(ex)==degree];normalco=poly.coeff_monomial(n**degree);assert normalco!=0
  # Removeotherlowestdegree monomials using17explicit Gaussmembers. This
  # produces a pure n-power initialform, independentoflocaltieorder.
  corrected=poly;clifts=list(lifts);indexof={row:j for j,row in enumerate(indices)};gausspol=[s.Poly(f,*variables)for f in gauss]
  for ex,co in initial:
   if not any(ex[:17]):continue
   j=next(i for i,k in enumerate(ex[:17])if k);reducedex=list(ex);reducedex[j]-=1;mon=s.Poly.from_dict({tuple(reducedex):co},variables);corrected-=mon*gausspol[j]
   for k,row in enumerate(rows):clifts[indexof[row]]-=mon*s.Poly(inverse[j,k],*variables)
  assert sum((l*f for l,f in zip(clifts,epol)),s.Poly(0,*variables))==corrected;lowest=sum(c*s.prod(v**k for v,k in zip(variables,ex))for ex,c in corrected.terms()if sum(ex)==degree);assert s.expand(lowest-normalco*n**degree)==0
  result.update(normal_purepower_certificate=dict(polynomial=str(corrected.as_expr()),lowest_degree=degree,initial_coefficient=str(normalco),lowest_homogeneous_form=str(normalco*n**degree),original_ideal_lifts={str(i):str(v.as_expr())for i,v in enumerate(clifts)if v!=0}),all18_member_identities_verified=True,reduced_fixedcut_support_is_origin=True,whole7source_full23_component_supported=True)
 out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(dict(local_dimension=dim,purepower=chosen>0,degree=result.get('normal_purepower_certificate',{}).get('lowest_degree'),seconds=result['seconds'])))
if __name__=='__main__':main()
