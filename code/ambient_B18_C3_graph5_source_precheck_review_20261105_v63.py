#!/usr/bin/env python3
"""Reproduce two pre-emission SOURCE blockers, without reading back a root."""
import ast,hashlib,json,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[1]
ID='AMBIENT_B18_C3_GRAPH5_CACHED_SOURCE_PRECHECK_20261105_v63_REVIEW'
def sha(p):return hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(pathlib.Path(p).read_text())
def save(p,x):pathlib.Path(p).write_text(json.dumps(x,indent=2)+'\n')
def main():
    out=ROOT/'reports'/ID;out.mkdir(exist_ok=False)
    s=ROOT/'reports/AMBIENT_B18_C3_GRAPH5_CACHED_SOURCE_20261105_v63';r=read(s/'report.json')
    assert sha(s/'report.json')=='84de2029149eee7ddc187004d2f356ab3675f044e4a2198c3f75b3e56672a8df'
    for p,h in r['input_sha256'].items():assert sha(ROOT/p)==h
    D=read(s/'compiled-ten-even-roots-and-provided-joint105-plane3-cache.json')
    bad=[x for x in D['operator_records_originalID_op_args'] if x[1]=='c' and x[3]!=1]
    constants=[x for x in D['operator_records_originalID_op_args'] if x[1]=='c']
    assert len(bad)==32 and len(constants)==45 and [20084,'c',1,4] in bad and [1675,'c',-4,5] in bad
    source=ROOT/'code/ambient_B18_C3_graph5_cached_source_20261105_v63.py'
    worker=ROOT/'code/ambient_B18_C3_graph5_cached_values_20261105_v63.py'
    assert "v[0]=(Q(r[2]),1)" in source.read_text() and "v=[parse(r[2])]+[Q(0)]*5" in worker.read_text()
    C2=read(ROOT/'reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63/actual-physical12-and-C2-right-exact-quotient.json')
    raw=read(ROOT/'reports/AMBIENT_B18_ELL2_RAW33_QUOTIENT_20261105_v63/actual-ell2-product33-and-raw45-actions.json')
    basis=raw['actual33_basis'];mixed=[]
    for row in C2['physical12_coefficients'][7:]:
        vals=[[basis[i],c] for i,c in row['bA'] if basis[i][0] and basis[i][1]+basis[i][2]==1]
        assert vals;mixed.append({'physical_row':row['physical_row'],'actual_bA2_H_mixed_terms':vals})
    evidence={'source_literal_constant_arity_blocker':bad,'constant_count':45,'noninteger_constant_count':32,
      'source_majorant_and_prepared_worker_drop_denominator':True,
      'prepared_worker_not_source_faithful_before_fix':True,
      'five_full_columns_have_actual_nonzz_H_quadratics':mixed,
      'existing_source_carry_wording_requires_correction':True,
      'mathematical_repair':'For H mixed terms use r/ell*RH; its eB plane cubic carry is0 because RH2_eB=delta_i2_plane=0. The eB projection theorem may still hold after this correction.',
      'positive_quartic_witness':'For a nonzero r*ell^2 eB cubic, multiply by r to get r^2*ell^2 in CI(3,3); ell multiplication would vanish modulo ell^3.',
      'no_root_value_cubic_coefficient_or_rank_computed':True,'no_primary_bytes_changed':True,
      'status':'SOURCE_HOLD_before_ANY_value_emission; corrected version and profile need fresh audit'}
    save(out/'two-pre-emission-source-blockers.json',evidence)
    report={'id':ID,'all_passed':False,'audit_blockers_reproduced':True,'source_or_value_GO_fulfilled':False,
      'input_sha256':{str(p.relative_to(ROOT)):sha(p) for p in [s/'report.json',s/'compiled-ten-even-roots-and-provided-joint105-plane3-cache.json',source,worker]},
      'evidence_sha256':{'two-pre-emission-source-blockers.json':sha(out/'two-pre-emission-source-blockers.json')},
      'code_sha256':sha(__file__),'no_value_emission':True,'no_source_gauge_or_jet_solve':True}
    save(out/'report.json',report);print(json.dumps({'review_hold':True,'report':str(out/'report.json'),'sha256':sha(out/'report.json')}))
if __name__=='__main__':main()
