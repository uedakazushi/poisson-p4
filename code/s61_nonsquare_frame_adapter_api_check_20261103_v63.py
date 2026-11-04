#!/usr/bin/env python3
"""Small literal coefficient fixtures for the new frame adapter API."""
import argparse
import gzip
import json
from pathlib import Path
import sympy as s
from s61_nonsquare_oracle_frame_adapter_20261103_v63 import FrameAdapter,ROOT
from s61_nonsquare_structural_reductions_review_20261103_v63 import read,sha


def main(frame_path,output):
    frame_path=frame_path.resolve()
    frame=read(frame_path);f=FrameAdapter.load(frame_path);o=f.oracle;before=len(o.nodes)
    mono=(2,0,0,0,0,0,0,0)
    assert f.core_value(f.quad_coefficient(0,mono))==o.K(o.a)
    source=read(ROOT/'reports/S61_BINARY_NONSQUARE_CUBIC_REPLAY_20261103_v62/source-report.json')
    lines=gzip.decompress((ROOT/'reports/S61_BINARY_NONSQUARE_CUBIC_ROW_BASIS_20261103_v62/remaining.txt.gz').read_bytes()).decode().splitlines()
    checks=[]
    for pos in[0,41]:
        j=next(j for j in range(91)if int(lines[1+pos*91+j].split()[0]))
        ts=lines[1+pos*91+j].split();n=int(ts[0]);expected=o.K(o.R.from_dict({tuple(map(int,ts[2+4*k:5+4*k])):s.QQ(ts[1+4*k])for k in range(n)}))
        mon=tuple(source['normal_monomials'][j]);assert f.core_value(f.basis57_coefficient(pos,mon))==expected
        checks.append(dict(basis_position=pos,monomial=list(mon)))
    heavy=(0,0,0,0,0,2,0,0);old=read(ROOT/'reports/S61_BINARY_HYPERPLANE_HEAVY_SQUARE_PIVOT_20261102.json')
    assert f.core_value(f.quartic_coefficient(77,heavy))==o.K.from_expr(s.sympify(old['coefficient_formulas'][0]))
    result=dict(all_literal_API_fixtures_passed=True,quad0_v0square='a',basis57_product_and_genuine_fixtures=checks,
        encoded_scaled_quartic77_Hsquare_equals_old_T34_b77=True,initial_nodes=before,query_nodes_added=len(o.nodes)-before,
        input_sha256={str(p.relative_to(ROOT)):sha(p)for p in[Path(__file__),frame_path,ROOT/'code/s61_nonsquare_oracle_frame_adapter_20261103_v63.py']})
    with output.open('x')as stream:json.dump(result,stream,indent=2);stream.write('\n')
    print(json.dumps({k:v for k,v in result.items()if k!='input_sha256'},indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--frame',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();main(a.frame,a.output)
