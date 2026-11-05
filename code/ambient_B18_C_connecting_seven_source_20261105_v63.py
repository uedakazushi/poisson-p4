#!/usr/bin/env python3
"""Honest7-dimensional source quotient for the original comparison edge, no values."""
import hashlib,itertools,json,signal,time
from fractions import Fraction as Q
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_C_CONNECTING_SEVEN_SOURCE_20261105_v63'
PATHS={
 'G':'reports/AMBIENT_B18_C_GENUINE_R0_TRANSPORT_SOURCE_20261105_v63/report.json',
 'G_peer':'reports/AMBIENT_B18_C_GENUINE_R0_TRANSPORT_SOURCE_20261105_v63_REVIEW/report.json',
 'old40':'reports/AMBIENT_B18_T_KNOWN_INITIAL_LENGTH_UPPER17_20261104_v63/known-initial40-and-length-upper17.json',
 'old40_peer':'reports/AMBIENT_B18_T_KNOWN_INITIAL_LENGTH_UPPER17_20261104_v63_REVIEW/report.json',
 'T':'reports/AMBIENT_B18_T_EXACT_LOEWY11_20261105_v63/exact-Loewy11-and-source-basis-corollary.json',
 'T_peer':'reports/AMBIENT_B18_T_EXACT_LOEWY11_20261105_v63_REVIEW/report.json',
 'heads':'reports/AMBIENT_B18_TOLD_TNIL_NATURAL_MAP_20261104_v63/actual-natural-Told-to-Tnil-source-map-and-origin.json',
 'heads_peer':'reports/AMBIENT_B18_TOLD_TNIL_NATURAL_MAP_20261104_v63_REVIEW/report.json',
}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def mm(A,B):return[[sum((A[i][k]*B[k][j]for k in range(len(B))),Q(0))for j in range(len(B[0]))]for i in range(len(A))]
def inverse(A):
 n=len(A);X=[list(map(Q,row))+[Q(i==j)for j in range(n)]for i,row in enumerate(A)]
 for j in range(n):
  k=next(k for k in range(j,n)if X[k][j]);X[j],X[k]=X[k],X[j];v=X[j][j];X[j]=[q/v for q in X[j]]
  for k in range(n):
   if k!=j and X[k][j]:v=X[k][j];X[k]=[a-v*b for a,b in zip(X[k],X[j])]
 return[row[n:]for row in X]
def det3(A):
 return sum(((-1)**sum(p[i]>p[j]for i in range(3)for j in range(i+1,3))*A[0][p[0]]*A[1][p[1]]*A[2][p[2]]for p in itertools.permutations(range(3))),Q(0))
def main():
 started=time.monotonic();signal.alarm(60);pins={p:sha(ROOT/p)for p in PATHS.values()}
 d={k:json.loads((ROOT/p).read_text())for k,p in PATHS.items()};assert all(v['all_passed']for k,v in d.items()if k.endswith('_peer'))
 assert d['G']['actual_genuine_C_ell_R0_mu']==3 and d['G']['actual_genuine_C_ell_m2_zero']
 assert d['T']['exact_length']==11 and d['T']['exact_Loewy_HF']==[5,3,2,1]
 assert d['heads']['old_minimal_generator_free56_rows']==[22,23,29,30,55]
 gens=[g for g in d['old40']['homogeneous_module_generators40']if all(sum(mon)==1 for h,mon,q in g)]
 assert len(gens)==37
 rmon=[0,0,0,0,0,0,0,1]
 A=[[sum((Q(q)for h,m,q in g if h==j and m==rmon),Q(0))for j in range(5)]for g in gens]
 B=[row[:]for row in A];labels=list(range(37));rank=0;pc=[];selected=[]
 for col in range(5):
  k=next((i for i in range(rank,37)if B[i][col]),None)
  if k is None:continue
  B[rank],B[k]=B[k],B[rank];labels[rank],labels[k]=labels[k],labels[rank];selected.append(labels[rank])
  v=B[rank][col];B[rank]=[x/v for x in B[rank]]
  for i in range(37):
   if i!=rank and B[i][col]:v=B[i][col];B[i]=[a-v*b for a,b in zip(B[i],B[rank])]
  pc.append(col);rank+=1
 assert rank==3 and pc==[0,1,2]
 full=[A[i]for i in selected];minor=[row[:3]for row in full];iv=inverse(minor);det=det3(minor);assert det
 I3=[[Q(i==j)for j in range(3)]for i in range(3)];assert mm(minor,iv)==I3==mm(iv,minor)
 factors=mm([row[:3]for row in A],iv);assert mm(factors,full)==A
 null=[[Q(0),Q(0)],[Q(0),Q(0)],[Q(1,4),Q(0)],[Q(1),Q(0)],[Q(0),Q(1)]]
 assert mm(A,null)==[[0,0]for _ in range(37)]
 R=[[Q(0)]*7 for _ in range(7)];R[5][2]=Q(1,4);R[5][3]=1;R[6][4]=1
 assert mm(R,R)==[[0]*7 for _ in range(7)]
 assert all(sum((row[j]*R[i][j]for j in range(5)),Q(0))==0 for row in A for i in range(7))
 # Sourcebasis [h3,rh3,h4,rh4,h0,h1,h2-h3/4].
 C=[[Q(0)]*7 for _ in range(7)]
 for j,i in enumerate([3,5,4,6,0,1,2]):C[i][j]=1
 C[3][6]=Q(-1,4);ci=inverse(C);I7=[[Q(i==j)for j in range(7)]for i in range(7)];assert mm(C,ci)==I7==mm(ci,C)
 J=[[Q(0)]*7 for _ in range(7)];J[1][0]=1;J[3][2]=1;assert mm(mm(ci,R),C)==J
 oldbasis=d['T']['actual_filtered_monomial_basis11'];deep=[r for r in oldbasis if r['degree']>=2];assert len(deep)==3
 z5h4=next(r for r in oldbasis if r['target']==4 and r['monomial_original8']==[0,0,0,0,0,1,0,0])
 data={
  'scope':'OriginalHONESTedge delta:Told→G=C[ell], sameRbar=R0/(ell). Noedgecoefficientvalues; nothomologyafteruncontrolledbasechange.',
  'source_module':'Told length11/HF5321, minimalheads h0..4=oldfree56rows22,23,29,30,55',
  'target_module':'GenuineG≅J2⊕k² withallz0/r²0, previoussourceandpeer; muRbar3/dim4',
  'I':'(zA0,...,zA6,rA²) inRbar, exactannihilatorofG',
  'all37_linear_relation_r_coefficients_rows37_columns5':[[str(x)for x in row]for row in A],
  'fixed_rank3_rows':selected,'fixed_rank3_columns':[0,1,2],'fixed_rank3_minor':str(det),
  'fixed3_inverse':[[str(x)for x in row]for row in iv],'both_inverse_products_verified':True,
  'all37_row_factorizations':[[str(x)for x in row]for row in factors],
  'cyclic_top_chi_kernel5_by2':[[str(x)for x in row]for row in null],
  'actual_edge_constraint':'Incanonicaltargetbasis(g,rg,k1,k2), letchi_jbethegcoordinateofdelta(h_j). Thenchi=(0,0,c/4,c,d); c,dremainUNKNOWNandnotboth0becauseactualedgesurjective. Thesearetargetmoduleconstraints, NOTaninvalidfree-polynomialM0L1identity.',
  'honest_source_quotient':'M=Told/I*Told hasdim7/minimal5, actionrR7below, z0/r²0. Completeness: allactualminimalrelationconstants0; modIallhighercoefficientsvanishandalllinearrelationsarethe37knownspan(exactoldHF1=3). Thusnoadditionalrelationismod-Ihidden.',
  'source7_basis':['h0','h1','h2','h3','h4','r*h3','r*h4'],
  'actual_r_action7':[[str(x)for x in row]for row in R],
  'source7_module_type':'J2² ⊕ k³','canonical_chain_basis':['h3','r*h3','h4','r*h4','h0','h1','h2-h3/4'],
  'source7_basis_change':[[str(x)for x in row]for row in C],'source7_basis_inverse':[[str(x)for x in row]for row in ci],
  'both_basis_products_identity':True,'canonical_r7':[[str(x)for x in row]for row in J],
  'five_head_edge_parameter_count':17,
  'parameter_count_proof':'Fiveheads×4targetcoords=20QQparameters. Thethreeindependentchi constraintsleave17. AllremainingheadimagesstillUNKNOWN; surjectivityimposesrank3onthe3targettoprows,notanewchosenvalue.',
  'fixed_I_Told_source_basis':{'first':z5h4,'other_three':deep},
  'fixed_four_kernel_classes_proof':'m²Told(length3)⊂I*Told andz5*h4isindependentdegree1. SourceMlength7givesI*Toldlength11-7=4,sothesefourknownclassesformbasis. AllareinactualK/ellKbecausetargetIannihilated.',
  'honest_exact_kernel_reduction':'0→I*Told(length4)→K/ellK(length7)→ker(M7→G4)(length3)→0. ThisisNOTKbar/I*KbarandNOTK/IK; noArtin-Rees/Torconnectingimageisdiscarded.',
  'allhigher_Told_classes_edge_zero':'G m²ann impliesoldm²Told maps0. Thefirstclassz5*h4alsomaps0; remainingfirstclassesr*h3,r*h4areRtimesheadimages.',
  'future_actual_source_gate':'ConstructthefiveallorderformalheadcycleliftsandprojectRegDivell(V_full*sigma)throughgenuineqT7/kernelcoordinates BEFOREanyfiniteNFdivision. Current17parametersareNOTcomputed.',
  'no_new_physical_source_gauge_IFT_jet_or_edge_values':True,
  'no_original_H0_or_nilcube_claim':True,
  'source_profile':{'fixed_QQ_matrix':[37,5],'numeric_envelope':50000,'hard_numeric':150000,'hard_seconds':60,'hard_memory_bytes':2147483648,'hard_RAW':33554432},
 }
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);ep=out/'honest-seven-source-quotient-four-fixed-kernel-classes-and-edge-constraints.json';ep.write_text(json.dumps(data,indent=2)+'\n')
 report={'id':STEM,'all_passed':True,'status':'SOURCE_EDGE_FACTORS_HONEST7_J2SQUARED_K3_FOUR_FIXED_KERNEL_CLASSES_NO_VALUES',
  'source7_dimension':7,'fixed_old_kernel_classes_dimension':4,'actual_five_head_edge_values_computed':False,
  'input_sha256':pins,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/(STEM+'.md')),
  'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-started}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report));signal.alarm(0)

if __name__=='__main__':main()
