#!/usr/bin/env python3
"""Actual projectedZW25x31 plus one full-Hpsi source-bound input; no pivots."""
import argparse,hashlib,json,resource,signal,time
from pathlib import Path
import sympy as s
from sympy.polys.rings import ring
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_CUBIC_ZW25_AUGMENTED_INPUT_20261104_v63'
def read(p):return json.loads(p.read_text())
def sha(p):
 h=hashlib.sha256()
 with p.open('rb')as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def matrix(p):
 lines=p.read_text().splitlines();r,c=map(int,lines[0].split());assert len(lines)==1+r*c;entries=[]
 for line in lines[1:]:
  t=line.split();n=int(t[0]);assert len(t)==1+4*n;entries.append([[list(map(int,t[2+4*i:5+4*i])),t[1+4*i]]for i in range(n)])
 return[[entries[i*c+j]for j in range(c)]for i in range(r)]
def main(out):
 began=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(parents=True,exist_ok=False)
 pre=ROOT/'reports/S61_NONSQUARE_ELL_LFRAME_PREFLIGHT_20261104_v63';pr=read(pre/'report.json');pc=read(pre/'certificate.json');n=read(pre/'native.json')
 for name,h in pr['artifact_sha256'].items():assert sha(pre/name)==h
 q10=matrix(pre/'quad10x26.txt');G31=matrix(pre/'genuine16x31.txt');G40=matrix(pre/'genuine16x40.txt');assert len(q10)==10 and len(G31)==16
 assert G31[6]==G31[1]and all(G40[6][j]==G40[1][j]for j in range(5,40))and any(G40[6][j]!=G40[1][j]for j in range(5))
 assert n['genuine_C_factors_cancelled'][1]==n['genuine_C_factors_cancelled'][6]==1
 sourcep=ROOT/'reports/S61_BINARY_NONSQUARE_CUBIC_REPLAY_20261103_v62/source-report.json';src=read(sourcep);mons=src['normal_monomials'];ids=[]
 for v in[6,3]:m=[0]*8;m[v]=m[5]=1;ids.append(mons.index(m))
 raw=ROOT/'reports/S61_NONSQUARE_ELL_SIX_BORDERS_20261104_v63/frozen57-source.txt';basis=read(ROOT/'reports/S61_BINARY_NONSQUARE_CUBIC_ROW_BASIS_20261103_v62/report.json');assert sha(raw)==basis['remaining_matrix_sha256']
 R,a,e,d=ring('a,e,d',s.QQ);dec=lambda ts:R.from_dict({tuple(m):s.QQ(v)for m,v in ts});enc=lambda p:[[list(m),str(v)]for m,v in p.items()]
 H=[[R.zero,R.zero]for _ in range(16)]
 with raw.open()as f:
  assert f.readline().split()==['57','91','57']
  for i in range(57):
   for j in range(91):
    line=f.readline();assert line
    if i<41 or j not in ids:continue
    t=line.split();count=int(t[0]);terms={tuple(map(int,t[2+4*k:5+4*k])):s.QQ(t[1+4*k])for k in range(count)};H[i-41][ids.index(j)]=R.from_dict(terms)
  assert not f.read().strip()
 oldh=ROOT/'reports/S61_NONSQUARE_CUBIC_PROJECTED15_HEAVY_PREFLIGHT_20261104_v63/certificate.json';old=read(oldh);C=100*a+18*e**2-9*e-2;delta=-64000*C**2
 assert H[6][0]-H[1][0]==dec(old['unscaled_G61_minus_G20_heavy_coefficients']['6'])and H[6][1]-H[1][1]==dec(old['unscaled_G61_minus_G20_heavy_coefficients']['3'])
 z=dec(old['delta_scaled_HZ_HW_coefficients']['Z']);w=dec(old['delta_scaled_HZ_HW_coefficients']['W']);scales=[(-40)**3*C**(3-k)for k in n['genuine_C_factors_cancelled']]
 assert scales[1]==scales[6]==delta
 psi=[scale*(-w*hz+z*hw)for scale,(hz,hw)in zip(scales,H)];assert psi[6]==psi[1]
 keepG=[i for i in range(16)if i!=6];projected=q10 and [row+[[]for _ in range(5)]for row in q10]+[G31[i]for i in keepG]
 assert len(projected)==25 and all(len(row)==31 for row in projected)
 aug=[row+[[]]for row in projected[:10]]+[row+[enc(psi[i])]for row,i in zip(projected[10:],keepG)]
 assert len(aug)==25 and all(len(row)==32 for row in aug)
 profile={'rows':25,'columns':32,'entries_nonzero':sum(bool(v)for row in aug for v in row),'core_terms':sum(len(v)for row in aug for v in row),'max_entry_terms':max(len(v)for row in aug for v in row),'max_core_total_degree':max((sum(m)for row in aug for v in row for m,c in v),default=-1),'psi_nonzero_rows':sum(bool(p)for p in psi),'psi_total_terms':sum(len(p)for p in psi),'projected31_terms':sum(len(v)for row in projected for v in row)}
 note=ROOT/'research/notes'/(ID+'.md');deps=[ROOT/'reports'/n/'report.json'for n in['S61_NONSQUARE_ELL_LFRAME_PREFLIGHT_20261104_v63_REVIEW','S61_NONSQUARE_ELL_COLON_RANK10_20261104_v63_REVIEW','S61_NONSQUARE_ELL_LFRAME_BORDERS_20261104_v63_REVIEW','S61_NONSQUARE_CUBIC_PROJECTED15_UNIT_PAIR_20261104_v63_REVIEW']];assert all(read(p)['all_passed']for p in deps)
 labels=[['ZW_quad10',i,pc['quad10_rows'][i]]for i in range(10)]+[['ZW_genuine',i,41+i,basis['remaining_original_row_indices'][41+i]]for i in keepG]
 evidence={'matrix25x32':aug,'source25_labels':labels,'target31_recipe':{'full_restricted40':'pure35inL,x0,y1,y2,y3 plusH5','deleted9_pure_columns':list(range(9)),'retained31_full40_columns':pc['remaining26_full35_column_indices']+list(range(35,40)),'NO_heavy5_Schur_has_been_applied':True},
 'added_column':'psi_i=lambda_i*(-w_delta*rawHZ_i+z_delta*rawHW_i); Q10psi=0','rawHZ_HW_indices_in_frozen91':ids,
 'raw16_HZ_HW_coefficients':[[enc(p)for p in row]for row in H],'genuine16_source_clearing_units':list(map(enc,scales)),'psi16_coefficients':list(map(enc,psi)),
 'removed_genuine_relative_row':6,'kept_genuine_relative_row':1,'duplicate_zero_in_both_pi_and_psi':True,'source_unit_line':'G6-G1 isCdeltaafterequalclearing; pi(Cdelta)=psi(Cdelta)=0, notliteralzeroinfulltarget',
 'generic25_independence_proof':'AcceptednonzeroDelta34 givesrestricted36x40rank34; nineexactmonomialpivots remove9rank; onequadzerooverlapandoneGduplicatearedeletionsofzerorows. Remaining25x31hasgenericrank25, henceaug25x32also25. No newpoint.',
 'full57_sufficient_condition':'rankaug25=25 makesallremaining15genuineindependentunder(pi,psi)modquad10; globalunitCdelta supplies16th fullimage, acceptedproduct41 supplies41. Fullcubicrank57. Notaconverse whenpureT27directionsare needed.',
 'rankdrop_coverage':'Ifprojectedgenuine rank14, one source line remainsafterCdelta; psi nonzeroonit suffices. Ifpsi0,pureT27maystillrescuefullrank. Ifprojectedrank<=13, onepsicolumnalonescannotcertifyallmissingdirections.',
 'Delta34_zero_not_same_as_projected_rankdrop':True,'no_rank_job_or_huge_Delta10_determinant_expanded':True,'profile':profile}
 (out/'input.json').write_text(json.dumps(evidence,separators=(',',':'))+'\n');bounded=profile['core_terms']<=100000 and(out/'input.json').stat().st_size<32*1024*1024
 report={'id':ID,'all_source_binding_and_duplicate_checks_passed':True,'input_profile':profile,'allowed_input100k_and32MiB_guard_passed':bounded,'no_Gaussian_rank_or_point_job_in_preflight':True,'knownunit_Gaussian_authorized_separately_once_if_guard_passes':True,
 'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in[Path(__file__),note,pre/'report.json',pre/'certificate.json',pre/'genuine16x31.txt',pre/'genuine16x40.txt',pre/'quad10x26.txt',sourcep,raw,oldh,*deps]},'input_json_sha256':sha(out/'input.json'),'seconds':time.monotonic()-began,'frozen':True}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k!='input_sha256'},indent=2),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);main(p.parse_args().output.resolve())
