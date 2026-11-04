#!/usr/bin/env python3
"""ONE e reconstruction; only input preparation/cost for H4 and22conditions."""
import argparse,hashlib,json,math,resource,signal,subprocess,sys,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_RAW_H_B_K13_E_RECONSTRUCTION_20261104_v63';sys.set_int_max_str_digits(0)
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def up(p):return {int(a):Q(v)for a,v in p}
def ip(p):return {int(a):int(v)for a,v in p}
def add(a,b,c=1):
 z=dict(a)
 for i,v in b.items():
  z[i]=z.get(i,0)+c*v
  if not z[i]:del z[i]
 return z
def mul(a,b):
 z={}
 for i,c in a.items():
  for j,d in b.items():z[i+j]=z.get(i+j,0)+c*d
 return{k:v for k,v in z.items()if v}
def main(out,work):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.signal(signal.SIGALRM,lambda *_:(_ for _ in()).throw(TimeoutError('60s total cap')));signal.alarm(60);out=out.resolve();work=work.resolve();out.mkdir(parents=True,exist_ok=False);work.mkdir(parents=True,exist_ok=False)
 b=ROOT/'reports/S61_NONSQUARE_RAW_H_B_MIXED_LENGTH_K13_PREFLIGHT_20261104_v63';wb=ROOT/'reports/S61_NONSQUARE_RAW_H_B_MIXED_MINOR_W_20261104_v63';assert read(wb/'report.json')['all_W_components_zero']
 p=read(b/'mixed-minor-W-source-preflight.json');s=read(b/'K13-source-input.json');g={a:int(v)for a,v in up(s['g13']).items()};V=list(map(up,p['normalized_V_modg13']));DV=math.lcm(*(v.denominator for z in V for v in z.values()));vi=[{a:int(v*DV)for a,v in z.items()}for z in V]
 H=list(map(up,s['coefficient_arrays_by_e_low_to_high'][0]));assert H[16]=={0:Q(1)};DH=math.lcm(*(v.denominator for z in H[:16]for v in z.values()));hi=[{a:int(v*DH)for a,v in z.items()}for z in H[:16]]
 inp=out/'native-input.txt'
 with inp.open('w')as f:
  for z in[g]+vi+hi:
   f.write(str(len(z))+'\n')
   for a,v in sorted(z.items(),reverse=True):f.write(str(a)+' '+str(v)+'\n')
  f.write(str(DH)+'\n'+str(DV)+'\n')
 cpp=ROOT/'code/s61_nonsquare_raw_H_B_K13_e_reconstruction_20261104_v63.cpp';exe=work/'K13-e'
 with(out/'compile.log').open('w')as log:subprocess.run(['bash','-lc','source /workspace/cloud-setup/activate.sh\nc++ -O2 -std=c++17 "$1" -o "$2" $(pkg-config --cflags --libs flint)','compile',str(cpp),str(exe)],stdout=log,stderr=subprocess.STDOUT,check=True,timeout=max(1,58-(time.monotonic()-start)))
 with(out/'native.log').open('w')as log:subprocess.run(['bash','-lc','source /workspace/cloud-setup/activate.sh\nexec "$@"','K13-e',str(exe),str(inp),str(out/'native-certificate.json')],stdout=log,stderr=subprocess.STDOUT,check=True,timeout=max(1,58-(time.monotonic()-start)))
 n=read(out/'native-certificate.json')
 if not n['pivot_unit']:raise RuntimeError('Proper pivot gcd retained; no field inverse promoted')
 j=n['selected'];pc=int(n['pivot_content']);prim=ip(n['pivot_primitive']);assert vi[j]=={a:pc*v for a,v in prim.items()};assert add(mul(ip(n['unit_S']),prim),mul(ip(n['unit_T']),g))=={0:int(n['unit_constant'])}
 tn=ip(n['e_numerator']);td=int(n['e_denominator']);L=g[13]
 for i,r in enumerate(n['EV_rows']):
  raw=add({a:DH*v for a,v in vi[i-1].items()}if i else{},mul(hi[i],vi[15]),-1);assert raw==ip(r['raw']);assert {a:v*L**r['L_power']for a,v in raw.items()}==add(mul(g,ip(r['quotient'])),ip(r['remainder']))
  diff=add({a:v*DH*L**r['L_power']for a,v in mul(tn,vi[i]).items()},{a:v*td for a,v in ip(r['remainder']).items()},-1);assert diff==mul(g,ip(r['EV_equals_tV_quotient']))
 assert ip(n['e_raw_numerator'])==mul(ip(n['EV_rows'][j]['remainder']),ip(n['unit_S']))
 assert {a:v*L**n['e_reduction_L_power']for a,v in ip(n['e_raw_numerator']).items()}==add(mul(g,ip(n['e_reduction_quotient'])),ip(n['e_reduction_remainder']))
 assert {a:v*int(n['e_raw_denominator'])for a,v in tn.items()}=={a:v*td for a,v in ip(n['e_reduction_remainder']).items()}
 proof=dict(exact_triple_projection='QQ[a,a^-1,e]/(h1,k1,m1) ISOMORPHIC QQ[a]/g13 viae=t(a)',no_irreducibility_or_kernel_rank_assumption_needed=True,reason='Companionh1(E)=0; acceptedNkV=0 andNm1V=0 implyk1(E)V=m1(E)V=0. AllEV=tV=>p(E)V=p(t)V. CertifiedunitVj=>allthreep(t)=0. ThusTsurjectsKa13; lengthT<=13=dimQQKa13=>iso.',e_numerator=n['e_numerator'],e_denominator=n['e_denominator'],exact_T_length=13,explicit_H4_cover_length=52,actual_bad_scope='OnlyD(G) withsourceunitmembership; NOTactualBbadnonempty orwhole52bound.',all_companion_eigen_identities_replayed_integer=True,DH=str(DH),DV=str(DV),no_e_Euclid_newnorm=True)
 (out/'projection-isomorphism.json').write_text(json.dumps(proof,indent=2)+'\n')
 # PREPARE ONLY: bind actual original22/U/H and cost the next finiteNF input.
 sb=ROOT/'reports/S61_NONSQUARE_RAW_H_B_35_MINORS_20261104_v63';hb=ROOT/'reports/S61_NONSQUARE_RAW_H_B_PAIR_MONIC_COMPLETION_20261104_v63';src=read(sb/'certificate.json');hsrc=read(hb/'certificate.json');Pi=src['localized_ideal_generators_monic_QQ_after_namedunit_and_QQ_only_scaling'];U=src['FULLU_product'];HH=hsrc['monic_H_laurent'];assert len(Pi)==22 and sum(map(len,Pi))==2963
 counts=dict(source22_terms=sum(map(len,Pi)),source22_degree_a_e_d=[max(m[i]for z in Pi for m,v in z)for i in range(3)],U_terms=len(U),U_degrees=[max(m[i]for m,v in U)for i in range(3)],H4_terms=len(HH),H4_min_a=min(m[0]for m,v in HH),H4_degrees=[max(m[i]for m,v in HH)for i in range(3)],t_terms=len(tn),t_numerator_bits=max(abs(v).bit_length()for v in tn.values()),t_denominator_bits=td.bit_length(),Ka13_basis_dimension=13,H4_extension_dimension_over_Ka13=4,model_matrix_shape_over_Ka13=[4,88],matrix_scalar_QQ_slots_upper=4*88*13,all22_source_images_scalar_slots_upper=22*4*13,source_coefficient_maxbits=max(max(abs(Q(v).numerator).bit_length(),Q(v).denominator.bit_length())for z in Pi+[U,HH]for m,v in z),no_finite_NF_images_matrix_or_Gaussian_computed=True)
 recipe=dict(coefficient_algebra='QQ[a]/g13; a_inversefromacceptedclosedformula; e=t(a)',quotient='Ka13[d]/H4(a,t,d), monic4',source_file=str((sb/'certificate.json').relative_to(ROOT)),source22_key='localized_ideal_generators_monic_QQ_after_namedunit_and_QQ_only_scaling',U_key='FULLU_product',H4_source=str((hb/'certificate.json').relative_to(ROOT)),H4_key='monic_H_laurent',columns='FororiginalPi_i, allnormalforms d^j*Pi_i, i0..21,j0..3; basis1,d,d²,d³ =>4x88Ka13.',target_NOT_computed='FutureU^52 membership unconditionalQQfinite52; U^4onlyifKa13reducedfieldconditionindependentlyaccepted.',next_budget_required=True,costs=counts)
 # A conservative, literal integer-height bound, NOT a computed normal form.
 def cleared_sum(polys):
  vals=[Q(v)for z in polys for m,v in z];den=math.lcm(*(v.denominator for v in vals));return den,sum(abs(int(v*den))for v in vals)
 DS,SS=cleared_sum(Pi);DHs,SH=cleared_sum([HH]);growth=abs(L)+sum(abs(v)for v in g.values());Tmax=max(sum(abs(v)for v in tn.values()),td);Amax=max(abs(g[0]),sum(abs(v)for a,v in g.items()if a>0))
 NS=SS*Tmax**11*growth**127;DSout=DS*td**11*abs(L)**127
 NH=SH*Amax**5*Tmax**8*growth**144;DHout=DHs*abs(g[0])**5*td**8*abs(L)**144;HM=max(NH,DHout)
 NP=(DHout+4*HM)**10*growth**120;DP=DHout**10*abs(L)**120;NM=NS*NP*growth**12;DM=DSout*DP*abs(L)**12
 bytebound=lambda slots,num,den:slots*(150+(num.bit_length()*30103+99999)//100000+(den.bit_length()*30103+99999)//100000)
 recipe['conservative_input_normalform_height_bounds']=dict(source_a_reduction_steps_upper=127,H4_a_reduction_steps_upper=144,H4_coefficient_slots52=52,H4_numerator_bits_upper=NH.bit_length(),H4_denominator_bits_upper=DHout.bit_length(),H4_RAW_bytes_upper=bytebound(52,NH,DHout),matrix_numerator_bits_upper=NM.bit_length(),matrix_denominator_bits_upper=DM.bit_length(),matrix_RAW_bytes_upper=bytebound(4576,NM,DM),matrix_worst_bound_exceeds32MiB=bytebound(4576,NM,DM)>32*1024*1024,integer_bound_reason='SubstitutionclearsQQden,tDen^e,g0^negative_a; primitiveg13pseudo-reductionmultipliesL1normbygrowth=absLC+L1(g)atmostonceperdegreeelimination. dcompanion10stepsboundsddegree7pluscolumn3. No NF coefficients computed.',next_small_guaranteed_input_target='H4coefficientsALONE,52QQslots anddisplayedRAWbound; actualheights may sharpen full4x88budget. No automaticjob.')
 (out/'finite-H4-source-preflight.json').write_text(json.dumps(recipe,indent=2)+'\n')
 paths=[Path(__file__),cpp,ROOT/f'research/notes/{ID}.md',b/'report.json',b/'mixed-minor-W-source-preflight.json',b/'K13-source-input.json',wb/'report.json',wb/'native-certificate.json.W-checkpoint.json',sb/'report.json',sb/'certificate.json',hb/'report.json',hb/'certificate.json']
 report=dict(id=ID,all_passed=True,unit_pivot_certified=True,all16_EV_tV_exact=True,exact_triple_projection_Ka13=True,exact_projection_length13_Hcover52=True,no_actual_bad_nonempty_empty_claim=True,finite_H4_all22_source_preflight_ONLY=True,no_e_Euclid_newnorm_Gaussian_point_prime=True,native_seconds=n['native_seconds'],input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},evidence_sha256={p.name:sha(p)for p in out.iterdir()if p.is_file()},seconds=time.monotonic()-start)
 assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2));print(json.dumps(counts,indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--work',type=Path,required=True);a=p.parse_args()
 try:main(a.output,a.work)
 except Exception as e:
  if a.output.exists():(a.output/'failure.json').write_text(json.dumps({'all_passed':False,'error':str(e),'completed_native_certificate_preserved':(a.output/'native-certificate.json').exists(),'no_auto_nextmethod':True},indent=2)+'\n')
  raise
