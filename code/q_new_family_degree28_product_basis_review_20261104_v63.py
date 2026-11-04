#!/usr/bin/env python3
"""Independent sparse-dictionary universal Inv28 products and coefficient minor."""
import argparse,gzip,hashlib,json,resource,signal,time
from pathlib import Path
from q_new_family_degree16_invariant_probe_20261103_v63 import character
ROOT=Path(__file__).resolve().parents[1];ID='Q_NEW_FAMILY_DEGREE28_PRODUCT_BASIS_20261104_v63_REVIEW';P=32003
def read(p):
 with (gzip.open(p,'rt') if p.suffix=='.gz' else p.open())as f:return json.load(f)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def decode(a,d):
 mons=[(i,j)for i in range(d+1)for j in range(d+1-i)];assert len(mons)==len(a)
 return {e:c for e,c in zip(mons,a)if c}
def mul(a,b):
 r={}
 for (i,j),c in a.items():
  for (k,l),d in b.items():
   e=(i+k,j+l);r[e]=(r.get(e,0)+c*d)%P
 return {e:c for e,c in r.items()if c}
def det(a):
 a=[r[:]for r in a];v=1
 for i in range(len(a)):
  j=next(j for j in range(i,len(a))if a[j][i])
  if j!=i:a[i],a[j]=a[j],a[i];v=-v
  pivot=a[i][i];v=v*pivot%P;ip=pow(pivot,-1,P)
  for j in range(i+1,len(a)):
   c=a[j][i]*ip%P
   if c:a[j][i:]=[(x-c*y)%P for x,y in zip(a[j][i:],a[i][i:])]
 return v%P
def main(out):
 st=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));out.mkdir(parents=True,exist_ok=False)
 primary=ROOT/'reports/Q_NEW_FAMILY_DEGREE28_PRODUCT_BASIS_20261104_v63';pre=ROOT/'reports/Q_NEW_FAMILY_DEGREE28_PRODUCT_PREFLIGHT_20261104_v63'
 report=read(primary/'report.json');assert report['all_passed'] and report['complete_QQ_Inv28_basis']
 hp=ROOT/'reports/Q_NEW_FAMILY_FULL_DEGREE24_HESSIAN_BASIS_20261104_v63/native-certificate.json.gz';cp=ROOT/'reports/Q_NEW_FAMILY_FULL_DEGREE24_UNIVERSAL_BASIS_20261104_v63/native-certificate.json.gz';qp=ROOT/'reports/Q_NEW_FAMILY_DEGREE26_POISSON_BASIS_20261104_v63/native-certificate.json.gz'
 h,c,q=map(read,[hp,cp,qp]);f=read(pre/'saved-formal-low-columns.json.gz')
 assert f['F4']==c['F4'] and f['p1_p2_p3_p4']==c['E5'][1:5]
 assert f['primitive12']==[h['degree12_columns'][i]for i in h['degree12_selected_columns'][5:]]
 assert f['primitive16']==h['degree16_columns'][h['degree16_selected_columns'][-1]]
 assert f['Inv14']==[q['degree14_columns'][i]for i in q['degree14_selected_columns']]
 assert f['Inv24']==[h['degree24_columns'][i]for i in h['degree24_selected_columns']]
 F=decode(f['F4'],4);p=[decode(x,8)for x in f['p1_p2_p3_p4']];a=[decode(x,12)for x in f['primitive12']];b=decode(f['primitive16'],16);o=[decode(x,14)for x in f['Inv14']];v=[decode(x,24)for x in f['Inv24']]
 polys=[mul(F,x)for x in v]
 polys +=[mul(x,mul(p[j],p[k]))for x in a for j in range(4)for k in range(j,4)]
 polys +=[mul(b,x)for x in a]
 polys +=[mul(o[i],o[j])for i in range(3)for j in range(i,3)]
 mons=[(i,j)for i in range(29)for j in range(29-i)]
 cols=[[z.get(e,0)for e in mons]for z in polys];assert len(cols)==117
 native=read(primary/'native-certificate.json.gz');assert cols==native['columns117']
 minor=read(primary/'minor.json');assert len(minor['columns'])==len(minor['coefficient_rows'])==108
 matrix=[[cols[j][i]for j in minor['columns']]for i in minor['coefficient_rows']];assert matrix==minor['matrix'];dd=det(matrix);assert dd==minor['minor']==27379
 ip=ROOT/'reports/Q_NEW_FAMILY_HYPERSURFACE_RECOGNITION_20261103_v63/invariant-certificates.json'
 chars=character(list(map(tuple,read(ip)['Cartan_weights'])),28);assert chars[28]['invariant_dimension']==108
 proof=dict(all117_times435_coefficients_independently_equal=True,selected108_minor_mod32003=dd,character_dimension28=108,
  all80_entrance_arrays_bound_to_frozen_lower_degree_certificates=True,complete_QQ_Inv28_basis=True,
  no_source28_restriction_or_new_point_or_defining_equation=True)
 (out/'independent-proof.json').write_text(json.dumps(proof,indent=2)+'\n')
 paths=[Path(__file__),ROOT/'research/notes'/(ID+'.md'),primary/'report.json',primary/'native-certificate.json.gz',primary/'minor.json',primary/'definitions.json',pre/'saved-formal-low-columns.json.gz',hp,cp,qp,ip,ROOT/'code/q_new_family_degree16_invariant_probe_20261103_v63.py']
 rr=dict(id=ID,all_passed=True,complete_QQ_Inv28_basis=True,input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},evidence_sha256={'independent-proof.json':sha(out/'independent-proof.json')},seconds=time.monotonic()-st)
 (out/'report.json').write_text(json.dumps(rr,indent=2)+'\n');signal.alarm(0);print(json.dumps(dict(all_passed=True,seconds=rr['seconds'])))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
