#!/usr/bin/env python3
"""Two exact actual-source coefficient queries of the frozen lazy quartic DAG."""
from pathlib import Path
import hashlib,json,time
from s61_nonsquare_lazy_quartic_oracle_20261103_v63 import Oracle,ROOT

def main():
 st=time.monotonic();folder=ROOT/'reports/S61_NONSQUARE_LAZY_QUARTIC_ORACLE_20261103_v63';path=folder/'oracle.json';o,data=Oracle.load(path)
 oldpath=ROOT/'reports/S61_BINARY_HYPERPLANE_HEAVY_SQUARE_PIVOT_20261102.json';old=json.loads(oldpath.read_text());a,e,d=o.a,o.e,o.d
 P=125*a*a*d-450*a*e-200*a-54*e**3-27*e**2+33*e+6;R=15*a*d*e-10*a*d+200*a+18*e**2+6*e-12;B=(50*a+18*e**2-9*e-2)*d+60*e-40
 expected=[-3*o.K(P**2)/(500*o.K(B**2)),3*o.K(R)/(10*o.K(B))];roots=dict(zip(data['genuine_quartic_original_rows'],data['genuine_quartic_roots']));mon=(0,0,0,0,0,2,0,0);checks=[]
 for row,want in zip([77,131],expected):
  n=o.coefficient(roots[row],mon);got=o.core_value(n);assert got==want
  checks.append({'original_row':row,'normal_monomial':list(mon),'query_node':n,'exact_QQ_core_coefficient':str(got)})
 sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest();result={'all_passed':True,'two_actual_H2_coefficients_match_independent_frozen_source_formulas':True,'no_actual_core_point_or_finitefield':True,'no_large_parameter_matrix_expansion':True,'original_DAG_nodes':len(data['arithmetic']['nodes']),'nodes_after_two_queries':len(o.nodes),'checks':checks,'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in [Path(__file__),ROOT/'code/s61_nonsquare_lazy_quartic_oracle_20261103_v63.py',path,oldpath]},'seconds':time.monotonic()-st,'frozen':True}
 out=folder/'API_H2_CHECK.json';assert not out.exists();out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:result[k]for k in ['all_passed','original_DAG_nodes','nodes_after_two_queries','seconds']}))
if __name__=='__main__':main()
