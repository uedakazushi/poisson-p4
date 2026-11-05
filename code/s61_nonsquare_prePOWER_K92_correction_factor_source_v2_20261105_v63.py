#!/usr/bin/env python3
"""SOURCE ONLY: rank-support cancellation at26579 and finite F0/B0 profile gate."""
import argparse,hashlib,importlib.util,json,math,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_PREPOWER_K92_CORRECTION_FACTOR_SOURCE_V2_20261105_v63';W=(1,2244,68);INF=10**12
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def read(p):return json.loads(p.read_text())
def dump(p,x):p.write_text(json.dumps(x,indent=2,default=str)+'\n')
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(parents=True,exist_ok=False);pins={}
 def src(rel):p=ROOT/rel;pins[rel]=sha(p);return read(p)
 dbase='reports/S61_NONSQUARE_PREPOWER_K92_EXACT_D49_20261105_v63';dr=src(dbase+'/report.json');assert dr['all_passed']and pins[dbase+'/report.json']=='bd86df91dc0b68b5f44d25d9fd3530b5a0f5e18734fab14ec955abd4cbb3d527';ds=src(dbase+'/next-source-Schur-valuation.json');di=src(dbase+'/D49-index.json');assert ds['new_determinant_order_lower_bound']==26579
 support=src('reports/S61_NONSQUARE_PREPOWER_K92_SCHUR_GAP_SOURCE_PREFLIGHT_20261105_v63/minplus-proof.json');A=support['Aplus_min'];B=support['B_min'];C=support['C_min'];star=support['Ainverse_minplus_closure'];alpha=ds['row_dual_potentials'];beta=ds['column_dual_potentials'];assert sum(alpha)+sum(beta)==26579
 internal=[max(beta[j]-B[m][j]for j in range(7)if B[m][j]<INF)for m in range(4)];Fmin=[[min(C[i][k]+star[k][m]for k in range(4))for m in range(4)]for i in range(7)];assert internal==[-2243,-2175,-2243,-2243]
 assert all(B[m][j]>=beta[j]-internal[m]for m in range(4)for j in range(7));assert all(Fmin[i][m]>=alpha[i]+internal[m]for i in range(7)for m in range(4))
 Fs=[(i,m)for i in range(7)for m in range(4)if Fmin[i][m]==alpha[i]+internal[m]];Bs=[(m,j)for m in range(4)for j in range(7)if B[m][j]==beta[j]-internal[m]];assert len(Fs)==14 and len(Bs)==18
 # Every possible F0 coefficient is a direct C coefficient: no nonconstant
 # A^-1 path can reach this exact target degree.
 for i,m in Fs:
  target=alpha[i]+internal[m];assert C[i][m]==target and all(C[i][k]+star[k][m]>target for k in range(4)if k!=m)
 assert [j for m,j in Bs if m==1]==[0]
 Dlead=[[Q(0)for j in range(7)]for i in range(7)];Drefs=[]
 for e in di['entries']:
  i=e['residual_row'];j=e['residual_column'];threshold=alpha[i]+beta[j];assert e['valuation']>=threshold
  if e['valuation']==threshold:
   assert e['terms']and e['leading_QQ_coefficient']is not None;Dlead[i][j]=Q(e['leading_QQ_coefficient']);Drefs.append(dict(e,scaled_leading_coefficient=e['leading_QQ_coefficient']))
 pattern=[(i,j)for i in range(7)for j in range(7)if Dlead[i][j]];assert len(pattern)==9 and all(i==1 or j in[0,1]for i,j in pattern);assert Dlead[3][0]and all(not Dlead[3][j]for j in range(1,7))
 dump(out/'rank-support-zero-proof.json',dict(original_order_lower_bound=26579,outer_row_potentials=alpha,outer_column_potentials=beta,inner_potentials=internal,F_min=Fmin,B_min=B,
  factor_identity='Pscaled=t^-alpha * C*A^-1*B * t^-beta = (t^-alpha*C*A^-1*t^-internal)*(t^internal*B*t^-beta); BOTH factors formalpower series by all displayed inequalities',
  leading_factor_support_F7x4=Fs,leading_factor_support_B4x7=Bs,F0_uses_ONLY_direct_C_coefficients=True,D_leading_known_matrix=Dlead,D_leading_source_references=Drefs,
  D_rank_upper3_cover=dict(rows=[1],columns=[0,1]),D_row3_ONLYcol0_nonzero=str(Dlead[3][0]),B_innerrow1_ONLYpossiblecol0=True,
  proof=['rankD<=3, so dimkerD>=4','Drow3 forces x0=0 for xinkerD','B0internalrow1 kills xinkerD, so image(F0B0|kerD) hasdimension<=3','a nonzero x is killedbyD and F0B0; hence D-F0B0 singular'],
  resulting_coefficient_ZERO=dict(Schur_order=26579,normalized11_source_monomial=[16,91,21]),new_integer_order_lower_bound=26580,
  ZERO_is_determinant_layer_ONLY_not_all_entries_or_full_det=True,NO_new_correction_coefficient_products=True))
 base='reports/S61_NONSQUARE_PREPOWER_K92_ONE_ROW_PRECONDITION_20261105_v63';rt=src(base+'/retropicalized-source.json');row=src(base+'/actual-row11.json');R=[[Q(v)for v in rr]for rr in rt['row_transform']];Qmat=[[Q(v)for v in rr]for rr in rt['column_transform']];DR=math.lcm(*(v.denominator for rr in R for v in rr));DQ=math.lcm(*(v.denominator for rr in Qmat for v in rr));DD=DR*DQ;RZ=[[int(v*DR)for v in rr]for rr in R];QZ=[[int(v*DQ)for v in rr]for rr in Qmat];Ds=list(map(int,rt['constant_clearing_rowLCMs']));rp=rt['row_potentials'];cp=rt['column_potentials']
 hp=ROOT/'code/s61_nonsquare_prePOWER_K92_tropical_source_preflight_v2_20261105_v63.py';pins[str(hp.relative_to(ROOT))]=sha(hp);sp=importlib.util.spec_from_file_location('sourcegate',hp);h=importlib.util.module_from_spec(sp);sp.loader.exec_module(h);bind=src('reports/S61_NONSQUARE_PREPOWER_K92_DET11_GUARDED_ATTEMPT_20261104_v63/source-binding.json');state=ROOT/'reports/S61_NONSQUARE_PREPOWER_K92_DET11_GUARDED_ATTEMPT_20261104_v63/state.bin';pins[str(state.relative_to(ROOT))]=sha(state);assert pins[str(state.relative_to(ROOT))]=='cd5fb2431c2f65649153d82890d52b52ce026f9d0e478bc30592687f17b20f02';M,inventory,kept,skip=h.selected_inputs(state,bind['actual_pure_matrix_roots'])
 for p in row['columns']:M[1][p['column']]=dict(state_id=-1001-p['column'],c=Q(p['QQ_common_factor']),co=[(tuple(m),int(v))for m,v in p['primitive_integer_polynomial']])
 request_specs=[('F',i,m,4+i,m,alpha[i]+internal[m])for i,m in Fs]+[('B',m,j,m,4+j,beta[j]-internal[m])for m,j in Bs];physical={};requests=[];linear_products=0;Nmax=0
 for kind,oi,oj,u,v,g in request_specs:
  contributors=[];L=0
  for i in range(11):
   if not RZ[u][i]:continue
   for j in range(11):
    if not QZ[j][v]:continue
    p=M[i][j]
    for idx,(mon,z)in enumerate(p['co']):
     gap=rp[i]+cp[j]-sum(a*b for a,b in zip(mon,W))
     if gap!=g:continue
     zz=p['c']*z*Ds[i];assert zz.denominator==1;key=(i,j,g)
     if key not in physical:physical[key]=dict(id=len(physical),row=i,column=j,gap=g,state_polynomial_id=p['state_id'],source_index=idx,monomial=mon,source_integer=str(z),source_QQfactor=str(p['c']),row_LCM=str(Ds[i]),cleared_integer=str(zz.numerator))
     contributors.append(dict(physical_id=physical[key]['id'],R_integer=str(RZ[u][i]),Q_integer=str(QZ[j][v])));L+=abs(RZ[u][i]*zz.numerator*QZ[j][v]);linear_products+=2
  assert contributors;Nmax=max(Nmax,L);requests.append(dict(kind=kind,row=oi,column=oj,global_row=u,global_column=v,gap=g,contributors=contributors,integer_numerator_abs_bound=str(L)))
 Dnums=[[int(v*DD)for v in rr]for rr in Dlead];assert all(Q(v,DD)==q for rr,qq in zip(Dnums,Dlead)for v,q in zip(rr,qq));product_terms=[(i,m,j)for i,m in Fs for mm,j in Bs if mm==m];assert len(product_terms)==65
 Hbound=DD*max(abs(v)for rr in Dnums for v in rr)+4*Nmax*Nmax;bits=Hbound.bit_length();live=kept+row['actual_total_terms']+len(physical)+32+14+18+49+3*49+15000;raw_upper=(len(physical)+32+49+121)*((7*bits+20)//3+300)+500000
 dump(out/'next32-factor-profile-source.json',dict(fixed_weight=W,outer_potentials=[alpha,beta],inner_potentials=internal,F0_support=Fs,B0_support=Bs,F0_requires_NO_nonconstant_Ainverse=True,common_linear_denominator=str(DD),
  physical_atoms=list(physical.values()),requests32=requests,known_D_leading_integer_numerators=Dnums,F0B0_products65=product_terms,
  next_H_integer_recipe='Hnum=DD*Dnum - Fnum*Bnum; H=Hnum/DD^2',new_H_determinant_known_ZERO_by_support=True,next_goal='exactH rank/kernel/unitminor to enable next typedSchur; NOT nonzero7det attempt',
  next_job='ONE32finite integer coefficient sums +65products +7x7rank/kernel/unitminor. No series/weight/K/source572. RequiresROOTbudget.',
  linear_scalar_products=linear_products,quadratic_scalar_products=65,coefficient_numerator_bit_upper=Nmax.bit_length(),H_integer_bit_upper=bits,whole_live_upper=live,artifact_RAW_upper=raw_upper,
  NO_requested_F_or_B_coefficients_or_H_evaluated=True))
 pins[str(Path(__file__).relative_to(ROOT))]=sha(Path(__file__));note=ROOT/'research/notes'/f'{ID}.md';pins[str(note.relative_to(ROOT))]=sha(note);total=sum(p.stat().st_size for p in out.iterdir());external=sum((ROOT/p).stat().st_size for p in pins);assert live<150000 and raw_upper+external<32*1024*1024
 report=dict(id=ID,all_passed=True,frozen=True,source_only=True,correction_factor_support_and_D_kernel_prove_layer26579_ZERO=True,no_correction_coefficients_or_determinant_computed=True,determinant_order_now_atleast=26580,
  next_finite_linear_requests=32,physical_atoms=len(physical),linear_scalar_products=linear_products,quadratic_products=65,H_integer_bit_upper=bits,whole_live_upper=live,next_RAW_upper=raw_upper,
  no_full_det11_zero_or_global57_claim=True,total_RAW_before_report=total,external_read_RAW=external,input_sha256=pins,evidence_sha256={p.name:sha(p)for p in out.iterdir()},seconds=time.monotonic()-start)
 dump(out/'report.json',report);signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
