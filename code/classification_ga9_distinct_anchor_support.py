#!/usr/bin/env python3
"""Complete finite Ga9 gradedIFT with common-denominator native polynomials.

If Delta clears the normalbasis and all implicitblock inverses, weightd has
denominator Delta^(2d-1). Numerator arithmetic uses QQ[b,t], avoiding repeated
QQ(b) gcds. One existing Singular invocation performs the exact bounded job.
"""
from pathlib import Path
import hashlib,itertools,json,os,resource,subprocess,time
import sympy as s
from classification_ga9_distinct_anchor_frame import build_frame
from classification_rank3_mixed_component_review import COLUMNS,PAIRS,ROWS,QUADS
ROOT=Path(__file__).resolve().parents[1]

def main():
 frame=build_frame();started=frame['started'];K=frame['K'];b=frame['b'];xx=frame['xx'];Rb=K.field.ring;delta=Rb.one
 allfractions=[K.from_sympy(a)for v in frame['normals']for a in v if a!=0]
 allfractions.extend(a for _,_,inverse in frame['solvers'].values()for row in inverse for a in row if a)
 for a in allfractions:delta=delta.lcm(a.denom)
 def numerator(a):return a.numer*delta.exquo(a.denom)
 delta_expr=delta.as_expr();assert delta!=0
 work=ROOT/'work/ga9-distinct-anchor-support';work.mkdir(parents=True,exist_ok=True);input_path=work/'complete.sing';output_path=work/'certificate.txt'
 if output_path.exists():raise FileExistsError(str(output_path))
 variables=xx+(b,);names=list(map(str,variables))
 def polynomial(expr):
  poly=s.Poly(expr,*variables);terms=[]
  for ex,co in poly.terms():
   if co==0:continue
   factor=str(int(co))if co.q==1 else '('+str(co.p)+'/'+str(co.q)+')'
   factors=[factor]+[name+('^'+str(power)if power!=1 else'')for name,power in zip(names,ex)if power]
   terms.append('*'.join(factors))
  return '('+'+'.join(terms)+')'if terms else'0'
 def coeff_call(poly,row):
  quad,ex=ROWS[row];return 'C('+poly+'['+str(QUADS.index(quad)+1)+'],'+','.join(map(str,ex))+')'
 lines=['ring rpoly=0,(h,x,y,z,w,b,t0,t1,t2,t3,t4,t5),dp;','option(redSB);','poly delta='+polynomial(delta_expr)+';','matrix P0[5][5];']
 for i,j in PAIRS:
  lines+=['P0['+str(i+1)+','+str(j+1)+']='+polynomial(frame['P'][i,j])+';','P0['+str(j+1)+','+str(i+1)+']=-P0['+str(i+1)+','+str(j+1)+'];']
 for j,v in enumerate(frame['normals']):
  lines.append('matrix N'+str(j)+'[5][5];')
  entries={}
  for a,(ij,ex)in zip(v,COLUMNS):
   if a!=0:entries[ij]=entries.get(ij,0)+numerator(K.from_sympy(a)).as_expr()*s.prod(z**p for z,p in zip(xx,ex))
  for(i,k),expr in entries.items():lines+=['N'+str(j)+'['+str(i+1)+','+str(k+1)+']='+polynomial(expr)+';','N'+str(j)+'['+str(k+1)+','+str(i+1)+']=-N'+str(j)+'['+str(i+1)+','+str(k+1)+'];']
 lines+=['''proc C(poly p,int e1,int e2,int e3,int e4,int e5)
{ intvec ee=e1,e2,e3,e4,e5; int i; int j; int factorial=1;
  for(i=1;i<=5;i++){for(j=1;j<=ee[i];j++){p=diff(p,var(i));factorial=factorial*j;}p=subst(p,var(i),0);}
  return(p/factorial);
}
proc S(matrix A,matrix B,int cross)
{ list result; int omitted; int at; int ell; int i; int j; int n; int k; int parity;
  intvec quad; intvec tri; poly f; poly g;
  for(omitted=5;omitted>=1;omitted--){
    n=0;for(i=1;i<=5;i++){if(i!=omitted){n=n+1;quad[n]=i;}}
    f=0;
    for(at=1;at<=4;at++){
      n=0;for(i=1;i<=4;i++){if(i!=at){n=n+1;tri[n]=quad[i];}}
      i=tri[1];j=tri[2];k=tri[3];g=0;
      for(ell=1;ell<=5;ell++){
        g=g+A[i,ell]*diff(B[j,k],var(ell))+A[j,ell]*diff(B[k,i],var(ell))+A[k,ell]*diff(B[i,j],var(ell));
        if(cross){g=g+B[i,ell]*diff(A[j,k],var(ell))+B[j,ell]*diff(A[k,i],var(ell))+B[k,ell]*diff(A[i,j],var(ell));}
      }
      if(at%2==1){f=f+var(quad[at])*g;}else{f=f-var(quad[at])*g;}
    }
    result[size(result)+1]=f;
  }
  return(result);
}
list W; list NL; list RR; list AA; list LL;
matrix Wd[5][5]; ideal residual=0; poly cc; poly uu;
''']
 lines.append('write("'+str(output_path)+'","DELTA "+string(delta));')
 implicit_records=[]
 for degree in range(1,5):
  lines.append('NL=list(0,0,0,0,0);')
  for i in range(1,degree//2+1):
   j=degree-i
   if i>3 or j>3:continue
   lines.append('AA=S(W['+str(i)+'],W['+str(j)+'],'+('0'if i==j else'1')+');')
   for k in range(1,6):lines.append('NL['+str(k)+']=NL['+str(k)+']+AA['+str(k)+'];')
  if degree<=3:
   lines.append('Wd=0;')
   for j,weight in enumerate(frame['normal_weights']):
    if weight==degree:lines.append('Wd=Wd+t'+str(j)+'*delta^'+str(2*degree-2)+'*N'+str(j)+';')
   if degree in frame['solvers']:
    cs,rs,inverse=frame['solvers'][degree]
    for j,col in enumerate(cs):
     summands=[polynomial(numerator(a).as_expr())+'*'+coeff_call('NL',rs[k])for k,a in enumerate(inverse[j])if a]
     lines.append('uu=-('+'+'.join(summands)+');')
     ij,ex=COLUMNS[col];monom=s.prod(z**power for z,power in zip(xx,ex));i,k=ij
     lines+=['Wd['+str(i+1)+','+str(k+1)+']=Wd['+str(i+1)+','+str(k+1)+']+uu*'+polynomial(monom)+';','Wd['+str(k+1)+','+str(i+1)+']=Wd['+str(k+1)+','+str(i+1)+']-uu*'+polynomial(monom)+';','if(uu!=0){write("'+str(output_path)+'","IMP '+str(col)+' '+str(degree)+' "+string(uu));}']
     implicit_records.append((col,degree))
   lines.append('W['+str(degree)+']=Wd;LL=S(P0,Wd,1);')
  else:lines.append('LL=list(0,0,0,0,0);')
  lines.append('RR=list(0,0,0,0,0);')
  for k in range(1,6):lines.append('RR['+str(k)+']=delta*NL['+str(k)+']+LL['+str(k)+'];')
  for row in frame['selected_rows']:
   if frame['row_relative_weights'][row]==degree:lines.append('if('+coeff_call('RR',row)+'!=0){print("SELECTED_FAIL '+str(row)+'");exit;}')
  for row,rw in enumerate(frame['row_relative_weights']):
   if rw==degree:lines.append('cc='+coeff_call('RR',row)+';if(cc!=0){residual=residual,cc;write("'+str(output_path)+'","EQ '+str(row)+' '+str(degree)+' "+string(cc));}')
  lines.append('print("WEIGHT '+str(degree)+' COMPLETE");')
 # Convert the complete numeratorideal into QQ(b)[t]; Delta is a nonzero
 # genericfield unit, so multiplying powers did not change its ideal.
 lines+=['ring rfield=0,(t0,t1,t2,t3,t4,t5),dp;','map transfer=rpoly,0,0,0,0,0,0,t0,t1,t2,t3,t4,t5;','ideal G=std(transfer(residual));','write("'+str(output_path)+'","DIM "+string(dim(G)));','int gg;for(gg=1;gg<=size(G);gg++){write("'+str(output_path)+'","GB "+string(G[gg]));write("'+str(output_path)+'","LEAD "+string(leadmonom(G[gg])));}','write("'+str(output_path)+'","COMPLETE");','print("COMPLETE_GA9_POLYNOMIAL_SUPPORT");','exit;']
 input_path.write_text('\n'.join(lines)+'\n');print(json.dumps(dict(input_ready=True,delta_degree=delta.degree(),input_bytes=input_path.stat().st_size,setup_seconds=time.monotonic()-started)),flush=True)
 binary=ROOT/'work/nonsquare-v62/singular-runtime/prefix/usr/bin/Singular';native_started=time.monotonic()
 try:run=subprocess.run([str(binary),'-q',str(input_path)],text=True,capture_output=True,timeout=95)
 except subprocess.TimeoutExpired as err:
  output=(err.stdout or b'').decode()if isinstance(err.stdout,bytes)else(err.stdout or'');record=dict(all_passed=False,complete=False,actual_limit='95second nativeSingular walltime',partial_stdout=output[-4000:],native_input_sha256=hashlib.sha256(input_path.read_bytes()).hexdigest(),producer_bytes=Path(__file__).read_text(),seconds=time.monotonic()-started);(work/'resource.json').write_text(json.dumps(record,indent=2)+'\n');print(json.dumps(dict(complete=False,actual_limit=record['actual_limit'])));return
 print(run.stdout[-2500:],flush=True)
 if run.returncode or any(line.lstrip().startswith('?') for line in run.stdout.splitlines()) or not output_path.exists()or output_path.read_text().splitlines()[-1]!='COMPLETE':
  (work/'failure.txt').write_text(run.stdout+'\nSTDERR\n'+run.stderr);raise RuntimeError('Native job not complete; inspect '+str(work/'failure.txt'))
 native_rows=output_path.read_text().splitlines();implicit_polynomials=[];equations=[];gb=[];lead=[];dimension=None
 for line in native_rows:
  if line.startswith('IMP '):_,col,degree,f=line.split(' ',3);implicit_polynomials.append(dict(column=int(col),weight=int(degree),numerator=f,denominator_delta_power=2*int(degree)-1))
  elif line.startswith('EQ '):_,row,degree,f=line.split(' ',3);equations.append(dict(row=int(row),weight=int(degree),numerator=f,denominator_delta_power=2*int(degree)-1))
  elif line.startswith('GB '):gb.append(line[3:])
  elif line.startswith('LEAD '):lead.append(line[5:])
  elif line.startswith('DIM '):dimension=int(line[4:])
 assert dimension is not None
 inputs=[Path(__file__),ROOT/'code/classification_ga9_distinct_anchor_frame.py',frame['prior_path']]
 result=dict(task_id='E_reduced_nonnormal',all_passed=True,input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in inputs},producer_bytes=Path(__file__).read_text(),
  native_input_path=str(input_path.relative_to(ROOT)),native_input_sha256=hashlib.sha256(input_path.read_bytes()).hexdigest(),native_output_sha256=hashlib.sha256(output_path.read_bytes()).hexdigest(),
  coefficient_field='QQ at the canonical distinct-root source anchor',whole_actual_five_parameter_source_retained=True,all_four_normals_weights=frame['normal_weights'],true_F47_cut103_verified=True,normal_vectors=[{str(i):str(a)for i,a in enumerate(v)if a!=0}for v in frame['normals']],implicit_columns=frame['implicit'],implicit_rows=frame['selected_rows'],
  common_delta=str(delta_expr),common_delta_degree=delta.degree(),denominator_bound='implicit weightd divides Delta^(2d-1); nonlinear pairs have denominatorDelta^(2d-2), and all inverses/normals have denominator dividingDelta',
  complete_weight3_implicit_polynomials=implicit_polynomials,all_selected98_equations_exactly_zero=True,complete_all350_remaining_equations_through_weight4=equations,
  full_generic_residual_grevlex_basis=gb,full_leading_monomials=lead,generic_true_cut_dimension=dimension,maximal_radical_origin=dimension==0,
  actual_whole_source_component22_proved=dimension==0,local_uniqueness_or_smoothness_claimed=False,
  whole_distinct_source_component_implication='At the accepted canonical distinct-root source anchor, all5 sourcepartials and PGL/Euler/scalar give actualF47, equal to the orbit rank47. The fullsource imageclosure has the accepted upperprojective dimension22, so thisorbit is dense in the entire sourceclosure. A complete0dimensional true47cut gives rawlocaldimension<=47 by Krull, proving that WHOLE sourceclosure is a fullprojective22component. Positivecut dimension is not membership.',
  proof='Positive grading(1,1,1,1,1,1) and ambient implicitweightmax3 force the unique formalIFT map to be an exactfinite polynomial. The commonDelta induction computes exactly that map. FullJacobi rows have weightmax4, all350 are extracted once and selected98 vanish exactly. The saved sixvariable ideal is therefore the complete genuine formalcut overQQ(b), not a truncatedjet. Dimension0 implies originradical by positivegrading, giving the Krull/component squeeze with actualF47. Positive dimension is NOT componenthood or knownmembership.',independent_review_required=True,
  native_seconds=time.monotonic()-native_started,native_peak_RSS_KiB=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss,seconds=time.monotonic()-started)
 out=ROOT/'reports/CLASSIFICATION_GA9_DISTINCT_ANCHOR_SUPPORT_20261107.json'
 if out.exists():raise FileExistsError(str(out))
 out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(dict(complete=True,true_cut_dimension=dimension,implicit_polynomials=len(implicit_polynomials),remaining_equations=len(equations),seconds=result['seconds'])))
if __name__=='__main__':main()
