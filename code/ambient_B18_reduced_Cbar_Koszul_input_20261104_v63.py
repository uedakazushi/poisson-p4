#!/usr/bin/env python3
"""Fixed polynomial resolution input, no source/gauge values or new kernel."""
import argparse,hashlib,itertools,json
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_REDUCED_CBAR_KOSZUL_INPUT_20261104_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def main(out):
 out.mkdir(parents=True,exist_ok=False)
 files={
 'reports/AMBIENT_B18_REDUCED_COMPARISON_CBAR_VALUES_20261104_v63_REVIEW/report.json':'fb4546b8667fc50c668ded9cc62c5597f3d19c08cbe1a6f116b6532209a09e6c',
 'reports/AMBIENT_B18_REDUCED_COMPARISON_GENERIC_RANK2_20261104_v63_REVIEW/report.json':'18300f976757ec9de2ad7d69f450db8cd794d235b0d8d65948c264156431c266',
 'reports/AMBIENT_B18_REDUCED_COMPARISON_CBAR_VALUES_20261104_v63/actual-eight-nil-images-and-reduced-Cbar-cokernel.json':'4e7436474f1fdda82e17313f54a8ff58b18940ee7b3343dba8a920352272e1e5',
 'reports/AMBIENT_B18_FULLV56_COMPARISON_UNIT_20261104_v63/actual-fullV56-m0-and-typed-comparison.json':'97393693b708dd834e0eeae3196306269ac2d6540f3ff303864f9aa5d5a66aff',
 'reports/AMBIENT_B18_FULLV56_COMPARISON_UNIT_20261104_v63_REVIEW/report.json':'01ffb1d9168132c5c3769b93f70d16c5060a333652fe2f1304af68d22779ff8b',
 'reports/AMBIENT_B18_RADICAL_COMPARISON_OPERATOR_20261104_v63/actual-radical-normalizer-comparison-operator.json':'376d5fd2b8f2ab8a45651b0374d5b983ebe2df8f8a507108675942e853adbba5',
 'reports/AMBIENT_B18_HONEST_T_EULER_OPERATOR_20261104_v63/actual-honest-T116-936-Euler-operator.json':'a91fe6a6a4a1bf105ec030d6b3b3b95353d02a60616fa94bc558e3cf0a6cae94',
 }
 for p,h in files.items():assert sha(ROOT/p)==h,p
 for p in files:
  if p.endswith('REVIEW/report.json'):assert read(ROOT/p)['all_passed']
 basis=[[(sec,list(I))for sec in range(3)for I in itertools.combinations(range(8),d)]for d in range(9)]
 ranks=[len(v)for v in basis];assert ranks==[3,24,84,168,210,168,84,24,3]
 maps=[];count=0
 for d in range(1,9):
  ix={(s,tuple(I)):j for j,(s,I)in enumerate(basis[d-1])};cols=[]
  for sec,I in basis[d]:
   col=[]
   for j,v in enumerate(I):
    m=[0]*8;m[v]=2 if sec==2 and v==7 else 1
    col.append([ix[(sec,tuple(I[:j]+I[j+1:]))],m,(-1)**j]);count+=1
   cols.append(col)
  maps.append(cols)
 assert count==3072
 for d in range(1,8):
  for col in maps[d]:
   sums={}
   for j,m,a in col:
    for i,n,b in maps[d-1][j]:
     key=(i,tuple(x+y for x,y in zip(m,n)));sums[key]=sums.get(key,0)+a*b
   assert not any(sums.values())
 P=[[Q(1),Q(0),Q(0)],[Q(11,12),Q(0),Q(1)],[Q(0),Q(1),Q(0)]]
 Pinv=[[Q(1),Q(0),Q(0)],[Q(0),Q(0),Q(1)],[Q(-11,12),Q(1),Q(0)]]
 assert [[sum(P[i][k]*Pinv[k][j]for k in range(3))for j in range(3)]for i in range(3)]==[[Q(int(i==j))for j in range(3)]for i in range(3)]
 result={
 'actual_S_variable_order':['zA0','zA1','zA2','zA3','zA4','zA5','zA6','rA'],
 'DVR_ell_distinct_from_rA':True,
 'new_gamma_order':['tauA+11/12*tauB','tauA*tauB','tauB'],
 'old_tau_target_P_gamma_columns':[[str(c)for c in row]for row in P],
 'P_inverse':[[str(c)for c in row]for row in Pinv],'P_determinant':'-1',
 'last_regular_sequence_powers':[1,1,2],
 'basis_sector_wedge':basis,'ranks':ranks,'all8_differentials_sparse_columns':maps,
 'all3072_polynomial_entries_and_all_chain_squares_checked':True,
 'exactness':'three direct Koszul resolutions of actual regular sequences',
 'minimality':'all entries in actual maximal ideal',
 'source_target_orientation':'P^-1*barVrad=d1_gamma*m; barVrad*n=P*d1_gamma',
 'old_to_new_source_map':'Z=diag(Res4(A^t),I8), barVrad*Z=barVold',
 'actual_new_nil_source_columns':[20,21,22,23,28,29,30,31],
 'proposed_m_n_reg_division_slots':1632,
 'm_n_h780_built':False,
 'conditional_without_complement_honest_Tnil_shape':[128,972],
 'conditional_after24unit_cancellations_shape':[104,948],
 'no_new_source_jet_GB_or_actual_nil_cycle_kernel':True,
 }
 ep=out/'actual-reduced-Cbar-three-sector-Koszul.json';ep.write_text(json.dumps(result,separators=(',',':'))+'\n')
 rels=list(files)+[f'research/notes/{STEM}.md',str(Path(__file__).relative_to(ROOT))]
 report={'id':STEM,'all_passed':True,'status':'PROVED_actual_reduced_Cbar_minimal_resolution_and_source_section_preflight',
 'actual_Cprime_bar_minimal_resolution_ranks':ranks,'all3072_entries_and_chain_squares_checked':True,
 'literal_d1_d2_d3_shapes':[[3,24],[24,84],[84,168]],
 'actual_m_n_h780_or_Tnil_presentation_constructed':False,
 'input_sha256':{p:sha(ROOT/p)for p in rels},'evidence_sha256':{ep.name:sha(ep)}}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:report[k]for k in ['id','all_passed','literal_d1_d2_d3_shapes']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);main(p.parse_args().output.resolve())
