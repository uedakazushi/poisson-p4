#!/usr/bin/env python3
"""Source-only graph redundancy audit; no higher coefficient evaluation."""
import argparse,hashlib,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
import sympy as s
from ambient_B18_actual_top_coefficient_preflight_20261104_v63 import DAG
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_ZERO_GRAPH_SOURCE_AUDIT_20261104_v63'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main(out):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False)
    pp=ROOT/'reports/AMBIENT_B18_ZERO_NORMALIZER_TENSOR_MODEL_20261104_v63/actual-normalizer-tensor-free4-DAG.json';P=read(pp)
    peer=ROOT/'reports/AMBIENT_B18_ZERO_NORMALIZER_TENSOR_MODEL_20261104_v63_REVIEW/report.json';assert read(peer)['all_passed']
    rows=P['coordinate_rows36'];graphs=[64,68,73,74];odd=[18,19,28,29,60,67,72,63]
    assert rows[-12:]==odd+graphs
    # Four canonical divergence-free ey-wedge directions invisible in the
    # selected17 coordinates, but independently visible in the four graphs.
    D4=s.Matrix(P['four_divergence_relations40']);rawrows=P['zero_raw_rows40'];recon=s.Matrix(P['eliminated_four_coordinate_reconstruction'])
    directions=[]
    for graph in graphs:
        coord=s.Matrix([int(r==graph) for r in rows]);dep=recon*coord
        raw=[F(0)]*40
        for r,c in zip(rows,coord):raw[rawrows.index(r)]=F(c)
        for r,c in zip(P['eliminated_zero_rows4'],dep):raw[rawrows.index(r)]=F(c)
        assert D4*s.Matrix(raw)==s.zeros(4,1)
        assert [raw[rawrows.index(r)] for r in odd+P['actual_R0_rows9']]==[0]*17
        assert [raw[rawrows.index(r)] for r in graphs]==[F(int(r==graph)) for r in graphs]
        directions.append(raw)
    # A universal two-row elimination identity, retaining its actual nonunit
    # determinant and its even defect instead of declaring either invertible/0.
    ai,bi,ci,aj,bj,cj,ak,bk,ck,t,u=s.symbols('ai bi ci aj bj cj ak bk ck t u')
    delta=bi*cj-bj*ci;alpha=bk*cj-bj*ck;beta=bi*ck-bk*ci
    residual=delta*ak-alpha*ai-beta*aj
    assert s.expand(delta*(ak+bk*t+ck*u)-alpha*(ai+bi*t+ci*u)-beta*(aj+bj*t+cj*u)-residual)==0
    D=DAG();D.nodes=[tuple(n) for n in P['DAG_scalar_nodes']];D.cache={n:i for i,n in enumerate(D.nodes)};D.values=[];D.divisions=[]
    for n in D.nodes:
        if n[0]=='c':v=F(n[1],n[2])
        elif n[0]=='v':v=F(0)
        elif n[0]=='a':v=D.values[n[1]]+D.values[n[2]]
        elif n[0]=='m':v=D.values[n[1]]*D.values[n[2]]
        elif n[0]=='n':v=-D.values[n[1]]
        elif n[0]=='i':assert D.values[n[1]];v=1/D.values[n[1]]
        else:raise ValueError(n[0])
        D.values.append(v)
    D.zero=D.cache[('c',0,1)];D.one=D.cache[('c',1,1)]
    outputs=P['all12_zero_physical_difference_columns'];assert [r['physical_row'] for r in outputs]==odd+graphs
    aa,bb,cc,_=outputs[0]['column4'];ad,be,cf,_=outputs[1]['column4']
    det=D.add(D.mul(bb,cf),D.neg(D.mul(be,cc)));assert D.values[det]==0
    defects=[]
    for record in outputs[8:]:
        ag,bg,cg,_=record['column4']
        al=D.add(D.mul(bg,cf),D.neg(D.mul(be,cg)))
        bt=D.add(D.mul(bb,cg),D.neg(D.mul(bg,cc)))
        rr=D.total([D.mul(det,ag),D.neg(D.mul(al,aa)),D.neg(D.mul(bt,ad))])
        defects.append({'physical_graph_row':record['physical_row'],'delta_node':det,
          'row18_cofactor_node':al,'row19_cofactor_node':bt,'uncomputed_even_defect_node':rr,
          'identity':'delta*d_graph=alpha*d18+beta*d19+even_defect*e0',
          'delta_origin':'0','delta_nonzero_or_unit_not_claimed':True})
    relation=P['normalizer_tensor_relation_columns4_by48'];M32=relation[:32];G4=[r['column4'] for r in outputs[8:]]
    assert len(M32)==32 and all(len(r)==4 for r in M32)
    data={'DAG_scalar_nodes':D.nodes,'source_model_file':str(pp.relative_to(ROOT)),
      'free4_basis':P['free4_basis'],'same_R0_nodes9':P['R0_variable_nodes'],
      'pure96_source_unit_unknowns':sum((r['unknown48'] for r in P['normalizer_sources48_each']),[]),
      'pure96_source_unit_equations':sum((r['equations48'] for r in P['normalizer_sources48_each']),[]) if all('equations48'in r for r in P['normalizer_sources48_each']) else 'See literal normalizer_sources48_each in source_model_file',
      'exact_zero_odd_ideal_module_columns4_by32':M32,'four_graph_target_columns4':G4,
      'actual_finite_completion_target':'Find four coefficient vectors c_j in R0^32 with M32*c_j=G4_j, retaining completed local R0; generic fraction-field solutions do not suffice.',
      'universal_nonunit_cofactor_identities':defects,
      'canonical_zero40_directions4':[[str(c) for c in row] for row in directions],
      'all4_universal_directions_have_zero_selected17_and_independent_graphs':True,
      'ambient_only_Euler_divergence_rank2_Jacobi_cannot_eliminate_graph4':True,
      'source_specific_CA_CB_graph_redundancy_proved_or_disproved':False,
      'individual_graph_membership_does_not_give_common_graph_extension':True,
      'no_source3_4_corecubic_GB_new_beta_or_generic_field_solve':True}
    p=out/'actual-four-graph-membership-target.json';p.write_text(json.dumps(data,separators=(',',':'))+'\n')
    paths=[Path(__file__),ROOT/'research/notes'/f'{STEM}.md',pp,peer]
    rep={'id':STEM,'all_passed':True,'status':'ambient_identity_obstruction_and_literal_actual4_by32_membership_frontier',
      'four_universal_independent_graph_directions_verified':True,'exact_actual_odd_ideal_module_columns':32,
      'actual_graph_membership_targets':4,'nonunit_elimination_even_defects_retained':4,
      'source_specific_graph_redundancy_decided':False,'full_core_socle_depth_or_CI3_claimed':False,
      'input_sha256':{str(p.relative_to(ROOT)):sha(p) for p in paths},'evidence_sha256':{p.name:sha(p)},
      'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    (out/'report.json').write_text(json.dumps(rep,indent=2)+'\n');print(json.dumps(rep))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
