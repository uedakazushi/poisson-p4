#!/usr/bin/env python3
"""Independent SymPy-ring replay of the ORIGINAL two-row cubic split block.

Does not import the producer or run a growing Schur/GB computation.
"""
import argparse,gzip,hashlib,json,time
from pathlib import Path
import sympy as sp
from sympy.polys.rings import ring
ROOT=Path(__file__).resolve().parents[1]
STEM='S61_BINARY_NONSQUARE_SCALAR_PRODUCT_ROWS_20261103_v62'
R,a,e,P=ring('a,e,P',sp.QQ)
def sha(p):
 h=hashlib.sha256()
 with p.open('rb')as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()
def decode(data):
 q={}
 for ex,c in data:
  t=tuple(ex);assert len(t)==3 and min(t)>=0 and t not in q;q[t]=sp.QQ(c);assert q[t]
 return R.from_dict(q)
def parse(s):
 z=s.split();n=int(z[0]);assert len(z)==1+4*n
 return decode([[list(map(int,z[2+4*k:5+4*k])),z[1+4*k]]for k in range(n)])
def dot(u,v):return sum((p*q for p,q in zip(u,v)),R.zero)
def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args();out=args.output.resolve();out.mkdir(parents=True,exist_ok=False);start=time.monotonic();folder=ROOT/'reports'/STEM;cert=folder/'certificate.json.gz'
 with gzip.open(cert,'rt')as f:d=json.load(f)
 parent=json.loads((folder/'report.json').read_text());assert parent['all_passed']
 status=ROOT/'research/handoff/20261103_v62/v9_short_unit_pivots/STATUS.json';v9=json.loads(status.read_text());main=next(c for c in v9['cases']if c['case']=='main')
 source=ROOT/d['source'];assert sha(source)==d['source_sha256']==main['residual_sha256']=='6d21ebfbeade3c5a0a7bd2dc091bde5ec274622ddf899add905242164d616ea3'
 rows={int(i):[decode(p)for p in row]for i,row in d['rows_sparse'].items()};assert set(rows)=={0,7}
 # Stream the actual full source, and compare every coefficient of the two
 # original physical rows. No entrywise primitive/scalar normalization occurs.
 with source.open()as f:
  assert list(map(int,f.readline().split()))==[26,60,1];f.readline();f.readline()
  for i in range(26):
   for j in range(60):
    line=f.readline();assert line
    if i in rows:assert parse(line)==rows[i][j]
  assert not f.read().strip()
 cofs={int(i):[decode(p)for p in row]for i,row in d['cofactors_sparse'].items()};assert set(cofs)=={0,7}
 scalar_stats={}
 for i in [0,7]:
  assert len(cofs[i])==60 and dot(rows[i],cofs[i])==a
  scalar_stats[str(i)]={'active':sum(bool(q)for q in cofs[i]),'terms':sum(map(len,cofs[i])),'degree':max(sum(ex)for q in cofs[i]for ex in q)}
 assert scalar_stats['0']['terms']==64 and scalar_stats['7']['terms']==34
 c=cofs[7][43];assert c and c.is_ground
 T=[[decode(p)for p in row]for row in d['two_row_right_inverse_sparse']];assert len(T)==60 and all(len(row)==2 for row in T)
 for i,r in enumerate([0,7]):
  for j in [0,1]:assert dot(rows[r],[row[j]for row in T])==(a if i==j else R.zero)
 assert sum(len(q)for row in T for q in row)==206
 # Independent symbolic check of transport through the stated column43
 # replacement, BEFORE any whole-matrix row scaling. Let E replace column43
 # by cofs[7]. E^{-1}T0 has the following polynomial entries (c is QQ*).
 t0=[row[0]for row in T];inverse_image=[q-cofs[7][j]*t0[43].exquo(c)if j!=43 else t0[43].exquo(c)for j,q in enumerate(t0)]
 for i,r in enumerate([0,7]):
  transformed=[rows[r][j]if j!=43 else dot(rows[r],cofs[7])for j in range(60)]
  assert dot(transformed,inverse_image)==(a if i==0 else R.zero)
 other_image=dot(rows[0],cofs[7]);schur=[a*rows[0][j]-other_image*rows[7][j]for j in range(60)if j!=43]
 projected=[q for j,q in enumerate(inverse_image)if j!=43];assert dot(schur,projected)==a*a
 # Audit the native helper's immutable provenance, without rerunning its large
 # output computation. The exact scalar identity + constant c prove the
 # stated column replacement invertible on the whole source module.
 helper=ROOT/'work/nonsquare-v62/product-field-v7/witness-image-stats.cpp';manifest=ROOT/'research/handoff/20261103_v62/v9_short_unit_pivots/manifest.json';md=json.loads(manifest.read_text());hs=sha(helper)
 assert hs in manifest.read_text()
 log=folder/'full-row7-pivot.log';text=log.read_text();assert'EXACT_WITNESS_REPLACE 43 TARGET 7'in text and'TERMS 1801170'in text and'PASS 1'in text
 generated=d['native_outputs']['row7_replacement'];residual=ROOT/generated['path'];assert sha(residual)==generated['sha256']==parent['generated_residual']['sha256']and residual.stat().st_size==generated['bytes']
 for path,wanted in parent['input_sha256'].items():assert sha(ROOT/path)==wanted
 inputs=[Path(__file__),cert,folder/'report.json',ROOT/'code/s61_binary_nonsquare_scalar_product_rows_20261103_v62.py',ROOT/'research/notes'/f'{STEM}.md',status,manifest,helper,log]
 report={'id':'S61_NONSQUARE_SCALAR_PRODUCT_ROWS_20261103_V63_REVIEW','all_passed':True,'independent_algebra':'SymPy QQ polynomial ring (producer uses sparse Fraction convolution)','source_sha256':sha(source),'actual_original_rows_compared_all_coefficients':[0,7],'scalar_identities_verified':2,'scalar_statistics':scalar_stats,'two_row_right_inverse_all4_identities_verified':True,'two_row_right_inverse_terms':206,'column43_cofactor_nonzero_QQ_constant':str(c),'abstract_second_row_transport_identity_verified':'unscaled59-column Schur row times transported vector = a^2','whole_source_native_replay_not_repeated':True,'native_helper_pinned_in_V9_manifest':True,'native_output_sha256_checked':generated['sha256'],'accepted_localization':'a*P*F*B0*C*J !=0; P=ad, J=5P-6e-4','split_projection_proved':'ORIGINAL2x60 block has right inverse T/a','constructive_native_pivots_executed':32,'no_allpoint_product41_or_fullcubic57_claim':True,'no_full_cokernel_freeness_or_quartic_or_ambient_claim':True,'input_sha256':{str(p.resolve().relative_to(ROOT)):sha(p)for p in inputs},'seconds':time.monotonic()-start}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:report[k]for k in ['all_passed','scalar_identities_verified','two_row_right_inverse_all4_identities_verified','abstract_second_row_transport_identity_verified','seconds']},indent=2))
if __name__=='__main__':main()
