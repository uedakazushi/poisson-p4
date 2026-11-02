"""Prepare the vertical certificate from the original B-free core matrix."""
from collections import defaultdict
from fractions import Fraction
from functools import lru_cache
from pathlib import Path
import json
import subprocess
import sys

from s61_binary_hyperplane_Bfree_horizontal_20261102 import build_source
from s61_binary_hyperplane_root_model_20261102 import ROOT,CODE,P,digest,hashes,residue,runtime_environment

HORIZONTAL=ROOT/'reports/S61_BINARY_HYPERPLANE_BFREE_HORIZONTAL_20261102.json'
HORIZONTAL_CODE=CODE/'s61_binary_hyperplane_Bfree_horizontal_20261102.py'
EXTRACT=CODE/'s61_binary_Bfree_horizontal_steps_20261102/binary_hyperplane_Bfree_targets_exact_core_input_20261102.py'


def obtain_core(wd,provided,expected):
    reference=json.loads(HORIZONTAL.read_text())
    assert reference['all_passed'] and reference['source']['original_equations']==135
    assert reference['modular_certificate']['generated_inputs_sha256']['core_matrix']==expected
    inherited=reference['source']['input_sha256'].copy()
    assert all(digest(ROOT/name)==value for name,value in inherited.items())
    static=hashes([HORIZONTAL,HORIZONTAL_CODE,EXTRACT])
    static.update(inherited)
    if provided is not None:
        core=Path(provided).resolve()
        assert digest(core)==expected
        return core,static,{'mode':'hash-verified output of BFREE_HORIZONTAL','core_matrix':expected}
    env=runtime_environment()
    def run(label,args):
        subprocess.run(list(map(str,args)),cwd=ROOT,env=env,check=True,
                       capture_output=True,text=True,timeout=150)
    source=build_source(wd/'poisson-hyperplane-Bfree-cubic-20261102',env,run)
    assert source['source_rows_sha256']==reference['source']['source_rows_sha256']
    local_script=wd/'rebuild_exact_core.py'
    local_script.write_text(EXTRACT.read_text().replace('/tmp/',str(wd)+'/'))
    run('exact_core',[sys.executable,local_script])
    core=wd/'hyperplane-Bfree-t-targets-exact-core-20261102.in'
    assert digest(core)==expected
    static.update(source['input_sha256'])
    return core,static,{'mode':'original135 regenerated','core_matrix':expected,
                        'normalized_source_rows':source['source_rows_sha256']}


@lru_cache(maxsize=65536)
def rational(text):
    return Fraction(text)


def emit_field_row(stream,row):
    entries=[]
    for j,poly in enumerate(row):
        while poly and poly[-1]==0:
            poly.pop()
        if poly:
            entries.append(f"{j} {len(poly)} "+' '.join(map(str,poly)))
    stream.write(str(len(entries))+'\n')
    for entry in entries:
        stream.write(entry+'\n')


def prepare(core,wd,cert):
    """Keep only delta0..2 after exact d=3 and construct two Fp[a] matrices."""
    row_bounds=cert['row_bounds']
    col_bounds=cert['column_bounds']
    full_path,lead_path=wd/'full0',wd/'leading0'
    modular_paths={d:wd/f'field_d{d}' for d in (3,4)}
    modular={d:path.open('w') for d,path in modular_paths.items()}
    for stream in modular.values():
        stream.write('32003 321 238 0\n0 1\n')
    count=0
    total_jet_terms=0
    max_entry_degree=-1
    with core.open() as source,full_path.open('w') as jets,lead_path.open('w') as leading:
        assert list(map(int,source.readline().split()))==[328,238,321]
        jets.write('321 238\n')
        leading.write('321 238\n')
        for i in range(328):
            field_rows={3:[],4:[]}
            tight=[]
            for j in range(238):
                words=source.readline().split()
                terms=int(words[0])
                assert len(words)==1+4*terms
                if i>=321:
                    continue
                exact=defaultdict(Fraction)
                field4=defaultdict(int)
                for k in range(terms):
                    c=rational(words[1+4*k])
                    aa,ee,dd=map(int,words[2+4*k:5+4*k])
                    assert min(aa,ee,dd)>=0
                    pc=residue(c)
                    exact[aa,ee]+=c*3**dd
                    field4[aa]=(field4[aa]+pc*pow(2,ee,P)*pow(4,dd,P))%P
                    count+=1
                exact={m:c for m,c in exact.items() if c}
                jet={}
                field3=defaultdict(int)
                if exact:
                    degree=max(a for a,e in exact)
                    assert degree<=row_bounds[i]+col_bounds[j]
                    tight.append(degree-col_bounds[j])
                    max_entry_degree=max(max_entry_degree,degree)
                for (aa,ee),c in exact.items():
                    delta=row_bounds[i]+col_bounds[j]-aa
                    assert delta>=0
                    if delta<=2:
                        jet[ee,delta]=c
                    field3[aa]=(field3[aa]+residue(c)*pow(2,ee,P))%P
                order=sorted(jet.items())
                jets.write(str(len(order))+''.join(f' {c} 0 {ee} {d}' for (ee,d),c in order)+'\n')
                lead=[(ee,c) for (ee,d),c in order if d==0]
                leading.write(str(len(lead))+''.join(f' {c} 0 {ee}' for ee,c in lead)+'\n')
                total_jet_terms+=len(jet)
                for d,values in [(3,field3),(4,field4)]:
                    field_rows[d].append([values.get(k,0) for k in range(max(values,default=-1)+1)])
            if i<321:
                assert max(tight,default=0)==row_bounds[i]
                for d,stream in modular.items():
                    emit_field_row(stream,field_rows[d])
        assert not source.read().strip()
    for stream in modular.values():
        stream.close()
    return full_path,lead_path,modular_paths,dict(source_terms=count,jet_terms=total_jet_terms,
                                                  maximum_entry_a_degree=max_entry_degree)


def coefficient_line(line):
    z=line.split()
    count=int(z[0])
    assert len(z)==1+3*count
    result={}
    for k in range(count):
        c=rational(z[1+3*k])
        x,e=map(int,z[2+3*k:4+3*k])
        assert x==0 and e not in result
        residue(c)
        result[e]=c
    return result


def at_two(poly):
    return sum(residue(c)*pow(2,e,P) for e,c in poly.items())%P
