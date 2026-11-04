#!/usr/bin/env python3
"""One degree-two scalar-unit contraction of the frozen R19 and N56 maps.

No new physical source jet is constructed.  The output records every scalar
unit-system stage and checks every equation through degree two.
"""
import argparse,hashlib,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
from ambient_B18_top_graph_telescoping_review_20261104_v63 import solve

ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_R_N_PROJECTION_NORM_20261104_v63'
ZERO=(F(0),)*3
def read(p):return json.loads(p.read_text())
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()
def enc(a):
    if isinstance(a,(list,tuple)):return [enc(x) for x in a]
    if isinstance(a,dict):return {str(k):enc(v) for k,v in a.items()}
    return str(a) if isinstance(a,F) else a
def closure(nodes,roots):
    found=set();todo=list(roots)
    while todo:
        i=todo.pop()
        if i in found:continue
        found.add(i)
        if nodes[i][0] not in ['v','c']:todo.extend(nodes[i][1:])
    return sorted(found)
def mul(a,b):return tuple(sum((a[j]*b[k-j] for j in range(k+1)),F(0)) for k in range(3))
def evaluate(nodes,needed,binding):
    vals={}
    for i in needed:
        n=nodes[i];op=n[0]
        if op=='c':z=(F(n[1],n[2]),F(0),F(0))
        elif op=='v':z=binding.get(i,ZERO)
        elif op=='a':z=tuple(a+b for a,b in zip(vals[n[1]],vals[n[2]]))
        elif op=='n':z=tuple(-a for a in vals[n[1]])
        elif op=='m':z=mul(vals[n[1]],vals[n[2]])
        elif op=='i':
            a=vals[n[1]];assert a[0],('nonunit',i)
            z=(1/a[0],-a[1]/a[0]**2,a[1]**2/a[0]**3-a[2]/a[0]**2)
        else:raise ValueError(op)
        vals[i]=z
    return vals
def unit_solve(nodes,unknown,equations,J,external,roots):
    inverse,det=solve(J,[[F(i==j) for j in range(len(J))] for i in range(len(J))])
    needed=closure(nodes,equations+roots)
    binding=dict(external);binding.update({v:ZERO for v in unknown});stages=[]
    for degree in [1,2]:
        val=evaluate(nodes,needed,binding)
        assert all(val[r][0]==0 for r in equations)
        residual=[val[r][degree] for r in equations]
        correction=[-sum((inverse[i][j]*residual[j] for j in range(len(J))),F(0)) for i in range(len(J))]
        for u,c in zip(unknown,correction):
            p=list(binding[u]);p[degree]+=c;binding[u]=tuple(p)
        val=evaluate(nodes,needed,binding)
        assert all(all(q==0 for q in val[r][:degree+1]) for r in equations)
        stages.append({'degree':degree,'raw_equation_residual':residual,'correction':correction,
                       'all_unknown_coefficients':[[u,binding[u]] for u in unknown],
                       'all_equation_coefficients_after_correction':[val[r] for r in equations]})
    return val,binding,stages,inverse,det
def main(out):
    start=time.monotonic();signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
    out.mkdir(parents=True,exist_ok=False)
    rp=ROOT/'reports/AMBIENT_B18_R_ZERO_SLICE_IMAGE_20261104_v63/actual-R-zero-slice-image-DAG.json'
    np=ROOT/'reports/AMBIENT_B18_N_SAME_R0_BINDING_20261104_v63/actual-N-R0-binding-DAG.json'
    gates=[rp.parent/'report.json',np.parent/'report.json',
      ROOT/'reports/AMBIENT_B18_N_SAME_R0_BINDING_20261104_v63_REVIEW/report.json',
      ROOT/'reports/AMBIENT_B18_T1_T2_QUADRATIC_PHI_20261104_v63_REVIEW/report.json',
      ROOT/'reports/AMBIENT_B18_N_LITERAL_MODULE_CONE_20261104_v63_REVIEW/report.json']
    assert all(read(p)['all_passed'] for p in gates)
    R,N=map(read,[rp,np]);rows=R['same_actual_R0_rows9']
    assert rows==N['actual_same_R0_rows9']
    ru=R['group19_unknown_nodes'];re=R['all19_cut_equations']
    rroots=R['actual_R0_readback_nodes9']
    rv,rb,rs,ri,rd=unit_solve(R['DAG_scalar_nodes'],ru,re,
      [[F(x) for x in row] for row in R['all19_cut_Jacobian']],
      {R['frozen_actual_R_curve_parameter_node']:(F(0),F(1),F(0))},rroots)
    assert rd==-32
    rcoef=[rv[i] for i in rroots]
    assert [p[1] for p in rcoef]==list(map(F,R['actual_R0_slice_curve_velocity9']))
    base_indices=N['chosen7_R0_indices']+[N['Cohen_even8_R0_index']]
    nu=N['nu_tracefree_square_node'];aa=N['ninth_actual_even_a_node'];bb=N['ninth_actual_odd_b_node']
    ne=N['all56_actual_coefficient_Hensel_equations'];ns=N['all56_coefficient_unknown_nodes']
    ext={v:rcoef[i] for v,i in zip(N['actual_Cohen_R8_parameter_nodes'],base_indices)}
    nv,nb,nst,ni,nd=unit_solve(N['DAG_scalar_nodes'],ns,ne,
      [[F(x) for x in row] for row in N['unit56_Jacobian']],ext,[nu,aa,bb])
    assert nd==F(N['unit56_minor'])
    ninth=N['ninth_R0_index'];v=tuple(x-y for x,y in zip(rcoef[ninth],nv[aa]))
    nc,bc,ac=nv[nu],nv[bb],nv[aa]
    assert nc[0]==bc[0]==v[0]==0
    assert nc[1]==v[1]==0,('earliest_order_unexpected',nc,v)
    norm4=v[2]**2-nc[2]*bc[1]**2
    assert norm4,('projection_norm_coefficient_zero',nc,bc,v)
    evidence={'R_parameter_convention':'same frozen raw T, before any reparameterization',
      'same_actual_R0_rows9':rows,'R_group19_unit_inverse':ri,'R_group19_unit_stages':rs,
      'all9_actual_R0_curve_coefficients_through2':rcoef,
      'N_Cohen8_R0_indices':base_indices,'N_Cohen8_parameter_series':list(ext.items()),
      'N56_unit_inverse':ni,'N56_unit_stages':nst,
      'N_tracefree_nu_coefficients_through2':nc,'N_ninth_even_a_coefficients_through2':ac,
      'N_z1_b61_coefficients_through2':bc,'actual_v_rho61_minus_a_coefficients_through2':v,
      'actual_projection_annihilator':'fN=v^2-nu*z1^2',
      'first_nonzero_projection_norm_coefficient_degree':4,
      'projection_norm_coefficient_T4':norm4,
      'all19_and56_equations_checked_through2':True,
      'no_physical_source_jet_or_degree_above2':True}
    ep=out/'actual-R-N-projection-norm.json';ep.write_text(json.dumps(enc(evidence),separators=(',',':'))+'\n')
    paths=[Path(__file__),ROOT/'research/notes'/f'{STEM}.md',rp,np,*gates,
      ROOT/'code/ambient_B18_top_graph_telescoping_review_20261104_v63.py']
    rep={'id':STEM,'all_passed':True,'status':'actual_R0_curve_misses_N_projection_and_D4_support',
      'same_R0_coordinate_binding_checked':True,'all19_and56_unit_equations_zero_through2':True,
      'nu_first_coefficient_zero':True,'v_first_coefficient_zero':True,
      'projection_norm_first_nonzero_degree':4,'projection_norm_coefficient_T4':str(norm4),
      'generic_D0_and_T0_R_length2_by_full_ordinary_double_and_positive_piece_vanishing':True,
      'T0_flatness_curve_support_socle_depth_or_annihilator_claimed':False,
      'no_new_physical_source_jets_GB_points_primes_or_higher_degree':True,
      'input_sha256':{str(p.relative_to(ROOT)):sha(p) for p in paths},
      'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-start,
      'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    (out/'report.json').write_text(json.dumps(rep,indent=2)+'\n')
    print(json.dumps({k:rep[k] for k in ['all_passed','projection_norm_coefficient_T4','seconds','max_rss_kib']}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
