#!/usr/bin/env python3
"""Source-bound Fitting corollary and K13 INPUT preflight; NO e-Euclid."""
import argparse,gzip,hashlib,json,math,resource,signal,sys,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_RAW_H_B_MIXED_LENGTH_K13_PREFLIGHT_20261104_v63'
sys.set_int_max_str_digits(0)
def read(p):return json.loads(gzip.decompress(p.read_bytes())if p.suffix=='.gz'else p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def uni(p):return {int(k):Q(v)for k,v in p}
def add(a,b,c=Q(1)):
 z=dict(a)
 for m,v in b.items():
  z[m]=z.get(m,Q(0))+c*v
  if not z[m]:del z[m]
 return z
def mul(a,b):
 z={}
 for i,c in a.items():
  for j,d in b.items():z[i+j]=z.get(i+j,Q(0))+c*d
 return{k:v for k,v in z.items()if v}
def rem(p,g):
 z=dict(p);degree=max(g);lc=g[degree]
 while z and max(z)>=degree:
  k=max(z)-degree;c=z[max(z)]/lc;z=add(z,{i+k:v for i,v in g.items()},-c)
 return z
def encode(p):return [[k,str(v)]for k,v in sorted(p.items(),reverse=True)]
def decode2(p):
 z={}
 for m,c in p:
  assert m[2]==0 and min(m)>=0;z[m[0],m[1]]=Q(c)
 return z
def reduced2(p,g):
 out={}
 for e in sorted({j for a,j in p}):
  r=rem({a:v for(a,j),v in p.items()if j==e},g)
  for a,c in r.items():out[a,e]=c
 return out
def ser2(p):return [[[a,e,0],str(v)]for(a,e),v in sorted(p.items(),reverse=True)]
def stats(p):
 den=math.lcm(*(v.denominator for v in p.values()))
 return dict(terms=len(p),a_degree=max((a for a,e in p),default=-1),e_degree=max((e for a,e in p),default=-1),max_numerator_bits=max((abs(v.numerator).bit_length()for v in p.values()),default=0),max_denominator_bits=max((v.denominator.bit_length()for v in p.values()),default=0),common_denominator_bits=den.bit_length(),integer_clear_maxbits=max((abs(int(v*den)).bit_length()for v in p.values()),default=0))
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out=out.resolve();out.mkdir(parents=True,exist_ok=False)
 mb=ROOT/'reports/S61_NONSQUARE_RAW_H_B_MIXED_NORM16_20261104_v63';mi=ROOT/'reports/S61_NONSQUARE_RAW_H_B_MIXED_NORM16_INPUT_20261104_v63';oi=ROOT/'reports/S61_NONSQUARE_RAW_H_B_OLD_RESIDUAL_NORM16_INPUT_20261104_v63';fb=ROOT/'reports/S61_NONSQUARE_RAW_H_B_OLD_A117_FACTOR_20261104_v63';length=ROOT/'reports/S61_NONSQUARE_RAW_H_B_OLD_RESIDUAL_LENGTH_20261104_v63_REVIEW'
 assert sha(mb/'report.json')=='e3d05f0a5ecd42d735df280973d2e654c5ca6a161d869fc66f8918fd7a78b378';n=read(mb/'native-certificate.json');new=read(mi/'certificate.json');old=read(oi/'certificate.json');fac=read(fb/'certificate.json')
 g=uni(n['gcd']);assert max(g)==13 and g[0]!=0
 H=decode2(new['H_h1_monic16']);assert H==decode2(old['h1_e16_monic'])
 factors=[f for f in fac['factor_inventory']if f['degree_a']==13];assert len(factors)==1 and factors[0]['exponent']==1
 factor={m[0]:Q(v)for m,v in factors[0]['factor']};ratio=g[13]/factor[13];assert g=={a:ratio*v for a,v in factor.items()}
 ainv={a-1:-v/g[0]for a,v in g.items()if a>0};assert rem(add(mul({1:Q(1)},ainv),{0:Q(1)},Q(-1)),g)=={}
 def remcol(c):
  ans={}
  for e in range(16):
   for(a,j),v in decode2(c['matrix_rows_source_order'][e][0]).items():
    assert j==0;ans[a,e]=v
  return ans
 RK=remcol(old);RM=remcol(new)
 # Reuse exact source monic-division remainders; only reduce a modulo g13.
 inputs=[H,RK,RM];reduced=[reduced2(p,g)for p in inputs];assert all(a<13 for p in reduced for a,e in p)
 profiles=list(map(stats,reduced));assert profiles[0]['e_degree']==16 and all(s['e_degree']<16 for s in profiles[1:])
 coeffs=[]
 for p in reduced:
  coeffs.append([encode({a:v for(a,j),v in p.items()if j==e})for e in range(max(j for a,j in p)+1)])
 source=dict(known_field_not_new_beta_point=True,quotient_algebra='QQ[a]/(g13)',g13=encode(g),g13_source_factor_QQunit=str(ratio),matching_saved_degree13_factor=factors[0],irreducibility_source='PreviouslyexecutedOLD_A117_FACTOR exactFLINTQQfactorization. Reused receipt, no refactor. Independentirreducibilityreview isnotclaimedhere.',a_inverse=encode(ainv),a_inverse_identity='a*a_inverse-1 is -g13/g13_constant; checked byexact remainder',source_polynomial_names=['h1_monic','old_k1_mod_h1','mixed_m1_mod_h1'],coefficient_arrays_by_e_low_to_high=coeffs,source_profiles=profiles,all_reductions_only_numericQQleading_unit=True,no_e_Euclid_or_pivot_inverse_computed=True)
 (out/'K13-source-input.json').write_text(json.dumps(source,indent=2)+'\n')
 fitting=dict(base='R=QQ[a,a^-1]',common_free16_model='R[e]/(h1monic16), basis1,e,...,e15 sameexactsourceinbothinputs',presentation='T=coker([Mk1|Mm1])=R[e]/(h1,k1,m1)',two_square_subdeterminants='det(Dold*Mk1), det(Dnew*Mm1), Dold/Dnew nonzeroQQconstants',primitive_gcd13_in_Fitt0='SavedunivariateBezoutandknowncontent/aunits putg13 inidealof16maximalminors.',length_argument='T finiteQQ asquotientofold117algebra. OverPIDR itsorderpolynomiald_T generatesFitt0 and dividesg13. ThusQQlengthT<=13.',explicit_T_length_upper=13,monic_H4_cover_length_upper=52,NO_cokerB_module_identification=True,actual_bad_scope='D(G) only; G0Q/Lcasesretained.',actual_scheme_membership_additional_gate='OverFULLUcapD(G): mixedcoefficients are in(H,gsource,ksource) viaformal4x4normpenciladjugate. Sourceunitstrips and G inversion mustbeappliedliterallybefore a quotientmapfromT[d]/H toactualbad scheme isclaimed.',mixed_adjugate_operator='A(lambda)=s*(Mg+lambda*Mk); W=adj(A)e0. Lift4source monicdivisions: detA=s*(g+lambda*k)*W(e,lambda)-s*H*Q(e,lambda). Coefficientlambda1 givesmembership. No expandedmixedadjugate computed.',no_global_actual_Bbad_scheme_length_claim=True)
 (out/'fitting-length-lemma.json').write_text(json.dumps(fitting,indent=2)+'\n')
 plan=dict(actual_input_e_degrees=[s['e_degree']for s in profiles],actual_QQ_terms=sum(s['terms']for s in profiles),full_field_coefficient_slots_upper=3*17*13,max_live_field_coefficients_under1000=True,source_field_degree13=True,proposed_order='gcd_e(h1, mixed_m1_remainder), then gcdwithold_k1_remainder; no25degreeinputneeded',each_field_inverse='Forpivotp(a), computeexactQQ[a]xgcd(p,g13), storeS*p+T*g13=1. Allcoefficients mustbeverified before use.',zero_pivot='Exact13vectorzero check; neverinvert zero',nonunit_pivot='Ifgcd(p,g13)!=1, retain/split both quotientcomponents; no silentinverse/factoromission.',irreducibility_provenance_only_native_old_receipt=True,no_new_factor_or_field_Euclid_job=True,field_inverse_count_degree_bound=34,initial_maxQQ_coefficient_bits=max(max(s['max_numerator_bits'],s['max_denominator_bits'])for s in profiles),no_global_height_bound_for_Euclid_claimed=True,proposed_future_budget='ONE60s/2GB/32MiB checkpointednumberfieldEuclid onlyafterparentallocation; hard coefficient/output guards, noautomaticbudgetincrease',completion_degree0='ExplicitunitBezout givestripleT=0 (oldD(G)badabsent); G0casesretain',completion_degree1='e isrationalinK13; thenmonicH4andall22sourcefiniteNFchecksneedseparateinput/budget',completion_higher='Keepfinitequotient, nofalselinearorexistenceclaim')
 (out/'next-unit-inverse-plan.json').write_text(json.dumps(plan,indent=2)+'\n')
 # Better NEXT method: mixed16-minors from Nm1*adj(Nk)e0, input only.
 nb=ROOT/'reports/S61_NONSQUARE_RAW_H_B_OLD_RESIDUAL_NORM16_CANONICAL_20261104_v63';native=read(nb/'native-certificate.json.gz');V=[]
 for p in native['only16_solution_cofactors']:
  v={}
  for m,s in p:
   assert m[1:]==[0,0];q=Q(s);assert q.denominator==1;v[m[0]]=q.numerator
  V.append(v)
 common=math.gcd(*(abs(v)for p in V for v in p.values()));apow=min(a for p in V for a in p);assert common>0
 Vnorm=[{a-apow:Q(v//common)for a,v in p.items()}for p in V];Vr=[rem(p,g)for p in Vnorm]
 D=Q(read(mi/'cost-preflight.json')['matrix_uniform_QQ_denominator']);Nm=[]
 for row in new['matrix_rows_source_order']:
  rr=[]
  for p in row:
   u={a:v*D for(a,e),v in decode2(p).items()};assert all(v.denominator==1 for v in u.values());rr.append(rem(u,g))
  Nm.append(rr)
 DA=math.lcm(*(v.denominator for row in Nm for p in row for v in p.values()));DV=math.lcm(*(v.denominator for p in Vr for v in p.values()));BA=max(abs(int(v*DA)).bit_length()for row in Nm for p in row for v in p.values());BV=max((abs(int(v*DV)).bit_length()for p in Vr for v in p.values()),default=0)
 Lbits=abs(g[13].numerator).bit_length();gbits=max(abs(v.numerator).bit_length()for v in g.values());bound_num=BA+BV+8+12*(gbits+5);bound_den=DA.bit_length()+DV.bit_length()+12*Lbits
 Wplan=dict(source_Nm_constant_clear=str(D),V_common_QQ_content=str(common),V_common_known_a_power=apow,normalization='ONEcommonC*a^apowforALL16Ventries; noindependentcolumnscaling',normalized_V_before_reduction=[encode(p)for p in Vnorm],normalized_V_modg13=[encode(p)for p in Vr],Nm1_modg13=[[encode(p)for p in row]for row in Nm],V_nonzero_modg13=any(Vr),V_nonzero_components=sum(bool(p)for p in Vr),Nm1_reduced_terms=sum(len(p)for row in Nm for p in row),V_reduced_terms=sum(map(len,Vr)),common_numeric_denominator_bits=[DA.bit_length(),DV.bit_length()],integer_clear_coefficient_bits=[BA,BV],W_formula='W=Nm1*Vnormal modg13; NOTcomputedhere',W_scalar_convolution_products_upper=16*16*13*13,modular_reduction_coefficient_updates_upper=16*12*13,W_output_coefficient_slots_upper=16*13,W_integer_pseudoremainder_num_bits_upper=bound_num,W_numeric_denominator_bits_upper=bound_den,W_RAW_serialization_upper=16*13*(150+(bound_num+bound_den)*30103//100000),gcd_completion='Ifgcd(g13,allWi)=1, Fitt0(T)=R andT=0. SaveactualQQBezoutagainstg13 foranynonzeroWiifg13irreducibilityused, otherwisegcdallcomponents.',all_W_zero_completion='KeepTfinite, doNOTinferexactlength13orlinear_e withoutadditionalcertificate.',no_W_product_or_gcd_test_computed=True)
 (out/'mixed-minor-W-source-preflight.json').write_text(json.dumps(Wplan,indent=2)+'\n')
 paths=[Path(__file__),ROOT/f'research/notes/{ID}.md',mb/'report.json',mb/'native-certificate.json',mi/'report.json',mi/'certificate.json',oi/'report.json',oi/'certificate.json',fb/'report.json',fb/'certificate.json',length/'report.json',nb/'report.json',nb/'native-certificate.json.gz',ROOT/'reports/S61_NONSQUARE_RAW_H_B_MIXED_NORM16_20261104_v63_REVIEW/report.json']
 report=dict(id=ID,all_passed=True,source_Fitting_Tlength_atmost13_Hcover52=True,not_actualBmodule_or_global_scheme_length=True,exact_shared_h1_and_basis=True,saved_degree13_factor_binding=True,source_only_a_reductions=True,no_e_Euclid_or_new_norm_factorization_Gaussian_point_prime_GB=True,input_source_profiles=profiles,input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},evidence_sha256={p.name:sha(p)for p in out.iterdir()if p.is_file()},seconds=time.monotonic()-start)
 assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
