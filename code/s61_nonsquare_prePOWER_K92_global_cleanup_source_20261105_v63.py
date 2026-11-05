#!/usr/bin/env python3
"""Typedgeneric rowcleanup, all7 callbacks, fresh onecolumn specialization cost."""
import argparse,hashlib,json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_PREPOWER_K92_GLOBAL_CLEANUP_SOURCE_20261105_v63'
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def read(p):return json.loads(p.read_text())
def dump(p,x):p.write_text(json.dumps(x,indent=2,default=str)+'\n')
def bitsupport(pp):
 s=0
 for g,v in pp:
  if int(v):s|=1<<int(g)
 return s
def exponents(n):
 while n:b=n&-n;yield b.bit_length()-1;n-=b
def main(out):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));out.mkdir(parents=True,exist_ok=False);pins={}
 def src(rel):p=ROOT/rel;pins[rel]=sha(p);return read(p)
 actual='reports/S61_NONSQUARE_PREPOWER_K92_WHOLE_ROW_ENTRY_20261105_v63';ar=src(actual+'/report.json');assert ar['all_passed']and ar['first_nonzero_obstruction_column']==0;ap=src(actual+'/all-source-and-prefix-identities.json');N0=src(actual+'/whole-row-column0.json');assert not N0['is_ZERO']
 rowframe=src('reports/S61_NONSQUARE_PREPOWER_K92_ONE_ROW_PRECONDITION_20261105_v63/retropicalized-source.json');state=ROOT/'reports/S61_NONSQUARE_PREPOWER_K92_DET11_GUARDED_ATTEMPT_20261104_v63/state.bin';pins[str(state.relative_to(ROOT))]=sha(state);assert pins[str(state.relative_to(ROOT))]=='cd5fb2431c2f65649153d82890d52b52ce026f9d0e478bc30592687f17b20f02'
 src('reports/S61_NONSQUARE_PREPOWER_K92_DET11_GUARDED_ATTEMPT_20261104_v63/source-binding.json');hp=ROOT/'code/s61_nonsquare_prePOWER_K92_tropical_source_preflight_v2_20261105_v63.py';pins[str(hp.relative_to(ROOT))]=sha(hp)
 src('reports/S61_NONSQUARE_PREPOWER_K92_EXACT_ABC72_20261105_v63/report.json');src('reports/S61_NONSQUARE_PREPOWER_K92_EXACT_D49_20261105_v63/report.json');src('reports/S61_NONSQUARE_PREPOWER_FULL57_SOURCE_LIFT_20261104_v63/report.json')
 nodes=[]
 def node(name,shape,op,args,**more):nodes.append(dict(name=name,shape=shape,op=op,args=args,**more));return name
 node('M',[11,11],'external_exact_generic_source',[],state_path=str(state.relative_to(ROOT)),loader=str(hp.relative_to(ROOT)),literal_loader='selected_inputs(state, sourcebinding.actual_pure_matrix_roots)[0]',entry_semantics='p.c*sum p.co[mon]*a^mon0*e^mon1*d^mon2; exactstoredgenericpolynomials, NOT univariateliterals')
 node('Mprime',[11,11],'elementary_row_add',['M'],target=1,source=2,multiplier=2)
 node('Dr',[11,11],'diagonal_parameter_monomials',[],entries=[dict(QQ=str(c),a_power=-p)for c,p in zip(rowframe['constant_clearing_rowLCMs'],rowframe['row_potentials'])])
 node('Dc',[11,11],'diagonal_parameter_monomials',[],entries=[dict(QQ='1',a_power=-p)for p in rowframe['column_potentials']])
 node('R1',[11,11],'literal_QQ_matrix',[],value=rowframe['row_transform']);node('Q1',[11,11],'literal_QQ_matrix',[],value=rowframe['column_transform']);node('G',[11,11],'matrix_product',['R1','Dr','Mprime','Dc','Q1'])
 for name,rr,cc in [('P',[0,1],[0,1]),('Q',[2,3],[2,3]),('A21',[2,3],[0,1]),('Btop',[0,1],list(range(4,11))),('Bbottom',[2,3],list(range(4,11))),('C0',[4],list(range(4))),('C2',[6],list(range(4))),('D0',[4],list(range(4,11))),('D2',[6],list(range(4,11)))]:node(name,[len(rr),len(cc)],'submatrix',['G'],rows=rr,columns=cc)
 node('deltaP',[],'determinant2',['P']);node('deltaQ',[],'determinant2',['Q']);node('adjP',[2,2],'classical_adjugate2',['P']);node('adjQ',[2,2],'classical_adjugate2',['Q']);node('dc',[1,4],'matrix_linear_combination',['C2','C0'],coefficients=[1,-59375]);node('dd',[1,7],'matrix_linear_combination',['D2','D0'],coefficients=[1,-59375])
 node('dcL',[1,2],'submatrix',['dc'],rows=[0],columns=[0,1]);node('dcR',[1,2],'submatrix',['dc'],rows=[0],columns=[2,3]);node('v',[1,2],'matrix_product',['dcR','adjQ']);node('u',[1,2],'scaled_difference',['deltaQ','dcL','v','A21'],formula='deltaQ*dcL-v*A21');node('w',[1,2],'matrix_product',['u','adjP']);node('vp',[1,2],'scalar_matrix_product',['deltaP','v']);node('delta',[],'scalar_product',['deltaP','deltaQ']);node('N',[1,7],'scaled_difference',['delta','dd','w','Btop','vp','Bbottom'],formula='delta*dd-w*Btop-vp*Bbottom')
 for j in range(7):node('N'+str(j),[],'matrix_entry',['N'],row=0,column=j)
 node('E',[11,11],'identity_except_one_row',['delta','w','vp'],row=6,coefficients=dict(row6='delta',row4='-59375*delta',row0='-w[0]',row1='-w[1]',row2='-vp[0]',row3='-vp[1]'))
 node('E_inverse',[11,11],'identity_except_one_row',['delta','w','vp'],row=6,coefficients=dict(row6='1/delta',row4='59375',row0='w[0]/delta',row1='w[1]/delta',row2='vp[0]/delta',row3='vp[1]/delta'),chart_inverse_only='delta')
 node('cleaned',[11,11],'matrix_product',['E','G']);node('P5',[5,5],'submatrix',['cleaned'],rows=[0,1,2,3,6],columns=[0,1,2,3,4]);node('P5inverse',[5,5],'block_triangular_inverse',['P5','deltaP','deltaQ','N0'],blocks='[[A4,b],[0,N0]]; A4^-1=[[adjP/deltaP,0],[-adjQ*A21*adjP/(deltaPdeltaQ),adjQ/deltaQ]]',chart_inverse_only=['deltaP','deltaQ','N0'])
 node('residual6',[6,6],'Schur_complement',['cleaned','P5inverse'],pivot_rows=[0,1,2,3,6],pivot_columns=[0,1,2,3,4],formula='D6-C6*P5inverse*B6',all_seven_N_callbacks_retained=True)
 dump(out/'typed-generic-rowcleanup-and-pivot.json',dict(field='QQ(a,e,d), on originalFULLU plus displayed newchart',nodes=nodes,
  generic_lift='replace everypotentialt^r by a^-r, since aFULLUunit; otherparametersremainindependent',
  specialization='a=t^-1,e=t^-2244,d=t^-68; allgenericnodes then equal acceptedABC/D and savedN0 univariate specialization',
  do_NOT_use_univariate_N0_as_generic_RF=True,chart=['deltaP!=0','deltaQ!=0','N0!=0'],unknown_chart_complement_RETAINED=True,
  generic_nonzero_only='deltaP,Q specializewithconstant1; N0 specializesnonzero4131terms. TheseRFs are nonzero, NOTglobalFULLUunits.',
  row_identity='E*G row6 first4ZERO; last7=N. detE=delta. E_inverse*E=I universalrowalgebra.',pivot_identity='detP5=delta*N0; rankG=5+rankresidual6 onchart',
  source_scope='onefixedcandidate11submatrix, transportedthroughacceptedfull57source; NOTnewfullrank57/det11nonzero orfullUclaim',independent_N0_peer_pending=True))
 L=int(ap['uniform_integer_denominator']);prefix={k:pp for k,pp in ap['all_prefix_integer_polynomials'].items()if k in['delta','w0','w1','vp0','vp1']};PM={k:(bitsupport(pp),sum(abs(int(v))for g,v in pp))for k,pp in prefix.items()};prefixslots=sum(len(p)for p in prefix.values());literal={};input_meta={}
 for rel in['reports/S61_NONSQUARE_PREPOWER_K92_EXACT_ABC72_20261105_v63/ABC72.jsonl','reports/S61_NONSQUARE_PREPOWER_K92_EXACT_D49_20261105_v63/D49.jsonl']:
  p=ROOT/rel;pins[rel]=sha(p)
  with p.open()as f:
   for ln,line in enumerate(f):
    e=json.loads(line);b=e.get('block','D');i=e['block_row']if b!='D'else e['residual_row'];j=e['block_column']if b!='D'else e['residual_column']
    if j==0 or not((b=='B'and i<4)or(b=='D'and i in[0,2])):continue
    c=Q(e['QQ_content']);fac=c.numerator*(L//c.denominator);S=0;N=0
    for g,v in e['sparse_integer_polynomial']:S|=1<<int(g);N+=abs(int(v)*fac)
    literal[b,i,j]=(S,N);input_meta[b,i,j]=dict(path=rel,line=ln,entry_index=e['entry_index'],cleared_factor=str(fac),source_terms=len(e['sparse_integer_polynomial']))
 def prod(p,q):
  S,N=p;T,M=q
  if S.bit_count()>T.bit_count():S,T=T,S
  U=0
  for g in exponents(S):U|=T<<g
  return U,N*M,S.bit_count()*T.bit_count()
 profiles=[]
 for j in range(1,7):
  d2=literal['D',2,j];d0=literal['D',0,j];dd=d2[0]|d0[0],d2[1]+59375*d0[1];terms=[prod(PM['delta'],dd)]+[prod(PM['w'+str(k)],literal['B',k,j])for k in range(2)]+[prod(PM['vp'+str(k)],literal['B',2+k,j])for k in range(2)];S=0;N=0
  for s,n,c in terms:S|=s;N+=n
  row=S.bit_count();scratch=max(s.bit_count()for s,n,c in terms);litcount=sum(input_meta[b,i,j]['source_terms']for b,i in[('D',0),('D',2)]+[('B',k)for k in range(4)]);owned=2*prefixslots+2*litcount+3*row+2*scratch+15000;raw=row*(90+(N.bit_length()+2)//3)+1000000
  profiles.append(dict(column=j,purpose='fillONEactualgenericpivotrowcallback viaSAMEformal specialization, notresumeold7rowtest',prefix5_polynomials_source=actual+'/all-source-and-prefix-identities.json',literal_inputs=[dict(block=b,row=i,**input_meta[b,i,j])for b,i in[('D',0),('D',2)]+[('B',k)for k in range(4)]],prefix_coefficients=prefixslots,literal_input_coefficients=litcount,output_support_upper=row,maximum_product_support_upper=scratch,owned_upper=owned,new_RAW_upper=raw,coefficient_bits_upper=N.bit_length(),integer_abs_upper=str(N),integer_pairs_upper=sum(c for s,n,c in terms),common_denominator=str(L**5)))
 best=min(profiles,key=lambda x:(x['integer_pairs_upper'],x['owned_upper'],x['column']))
 dump(out/'fresh-single-column-profiles.json',dict(all6_callbacks_retained=profiles,selected_cheapest_single_column=best,NO_actual_N1toN6_values_evaluated=True,not_automatic_all6_budget=True))
 pins[str(Path(__file__).relative_to(ROOT))]=sha(Path(__file__));note=ROOT/'research/notes'/f'{ID}.md';pins[str(note.relative_to(ROOT))]=sha(note);new=sum(p.stat().st_size for p in out.iterdir());external=sum((ROOT/p).stat().st_size for p in pins)
 report=dict(id=ID,all_passed=True,frozen=True,source_only=True,no_new_signed_callback_values=True,genericRF_lift_and_formal_specialization_separated=True,elementary_row_cleanup_and_inverse_typed=True,
  generic_chart_only=['deltaP','deltaQ','N0'],all_complements_retained=True,all7_full_callbacks_retained=True,next_single_column=best['column'],next_owned_upper=best['owned_upper'],next_NEW_RAW_upper=best['new_RAW_upper'],next_integer_pairs_upper=best['integer_pairs_upper'],
  no_det11_nonzero_or_global57_claim=True,total_new_RAW_before_report=new,external_read_RAW=external,input_sha256=pins,evidence_sha256={p.name:sha(p)for p in out.iterdir()},seconds=time.monotonic()-start)
 assert new+external<33554432;dump(out/'report.json',report);signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
