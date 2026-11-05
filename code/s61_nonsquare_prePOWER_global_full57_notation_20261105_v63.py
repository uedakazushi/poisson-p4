#!/usr/bin/env python3
"""Literal notation/provenance clarification; no polynomial calculation."""
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ID = 'S61_NONSQUARE_PREPOWER_GLOBAL_FULL57_NOTATION_20261105_v63'


def sha(p):
    with p.open('rb') as f:
        return hashlib.file_digest(f,'sha256').hexdigest()


def main():
    out=ROOT/'reports'/ID
    out.mkdir(exist_ok=False)
    refs=[
      'reports/S61_NONSQUARE_PREPOWER_GLOBAL_FULL57_SOURCE_20261105_v63/report.json',
      'reports/S61_NONSQUARE_PREPOWER_GLOBAL_FULL57_SOURCE_20261105_v63/global-operator-recipe.json',
      'reports/S61_NONSQUARE_PREPOWER_KAPPA_EXTRA2_ENTRY_20261104_v63/root239.json',
      'reports/S61_NONSQUARE_PREPOWER_KAPPA_EXTRA2_ENTRY_20261104_v63_REVIEW/report.json',
      'code/data/s61_binary_nonsquare_cubic_row_basis_20261103_v62.json']
    x=json.loads((ROOT/refs[4]).read_text())['relations'][5]
    assert x['deleted_original_row']==14
    rr=[v['original_row'] for v in x['coefficients']]
    assert rr==[7,8,10,12,14,15,19,21,22,24,26,31] and max(rr)<42
    own=Path(__file__)
    pins={p:sha(ROOT/p) for p in refs}
    pins[str(own.relative_to(ROOT))]=sha(own)
    d=dict(
      T='literal 7208-term QQ polynomial in KAPPA_EXTRA2_ENTRY/root239.json; T=B*(kappa dot raw pure column2) with the certified source normalization',
      T_NOT='the unrelated auxiliary 900*a-126*e-121 used in the earlier Q2 unit-completion proof',
      A='literal current15 row2 / original32 column22, not the g factor of P14',
      D12='the literal generic12 determinant in DIRECT_M12_FULL15_SOURCE with its QQ/a/L/gamma/B row/column normalizations; ONE a15 source coefficient has independent peer',
      old25_remaining_Q9_other8_pure_tail_columns=[1,2,3,4,6,7,8,9],
      direct15_consumes_tail2=True,
      seven_other_pure_tails_still_retained_in_full15x22=[1,3,4,6,7,8,9],
      literal_product_only_relation=x,
      product14_coefficient='12*a/5',
      all_nonzero_rows_in_this_relation_are_original_product_rows=True,
      no_old_frozen_bytes_changed=True,
      no_new_matrix_value_or_coefficient_calculation=True)
    p=out/'clarification.json';p.write_text(json.dumps(d,indent=2)+'\n')
    r=dict(id=ID,all_passed=True,frozen=True,source_only=True,input_sha256=pins,
           evidence_sha256={p.name:sha(p)})
    (out/'report.json').write_text(json.dumps(r,indent=2)+'\n')
    print(sha(out/'report.json'))


if __name__=='__main__':
    main()
