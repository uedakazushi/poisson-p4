#!/usr/bin/env python3
"""Exact first-jet differential of a QQ rational circuit at the Sp4 identity.

Restored from the serialized circuit semantics. This is current verification
code, not a recovered byte-identical copy of an unavailable historical builder.
The ten Sp4 differential columns are retained, even though their values are
specialized to the identity; other source inputs remain indeterminates.
"""
from pathlib import Path
import json
import sympy as s
ROOT=Path(__file__).resolve().parents[1]

def exact_matrix(body):
    if body.get('field_relation'):
        raise ValueError('Use the native circuit verifier for algebraic constants')
    names=body['input_names'];groups=['g'+str(i) for i in range(10)]
    group_values=[1,0,0,1,0,0,0,0,0,0]
    variables=[n for n in names if n not in groups]
    active=list(names);K=s.QQ.frac_field(*variables)
    cert=json.loads((ROOT/'legacy/v02/data/symplectic_slice_certificate.json').read_text())
    outputs=[body['output_u1_through_u61'][i]for i,b in enumerate(cert['free_variable_basis'])if b['weight']==1]
    assert len(outputs)==16
    needed=set(outputs);stack=list(outputs)
    while stack:
        op,arg=body['nodes'][stack.pop()]
        deps=arg if op in ('+','*')else[arg]if op=='inv'else[]
        for j in deps:
            if j not in needed:needed.add(j);stack.append(j)
    vals={};diffs={};zero,one=K.zero,K.one
    for i in sorted(needed):
        op,arg=body['nodes'][i]
        if op=='q':v=K.from_sympy(s.Rational(arg));d={}
        elif op=='v':
            name=names[arg]
            v=K.convert(group_values[groups.index(name)])if name in groups else K.gens[variables.index(name)]
            d={arg:one}
        elif op=='+':
            v=zero;d={}
            for j in arg:
                v+=vals[j]
                for a,z in diffs[j].items():d[a]=d.get(a,zero)+z
        elif op=='*':
            v=one;d={}
            for j in arg:
                old=v;v*=vals[j]
                d={a:z*vals[j]for a,z in d.items()}
                for a,z in diffs[j].items():d[a]=d.get(a,zero)+old*z
        elif op=='inv':
            v=one/vals[arg];d={a:-v*v*z for a,z in diffs[arg].items()}
        else:raise ValueError(op)
        vals[i]=v;diffs[i]={a:z for a,z in d.items()if z}
    matrix=[[diffs[j].get(a,zero)for a in range(len(active))]for j in outputs]
    return K,active,variables,matrix,len(needed)
