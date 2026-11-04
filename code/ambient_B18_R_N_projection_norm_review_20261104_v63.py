#!/usr/bin/env python3
"""Independent saved R19/N56 unit-stage replay through degree2 only.

Strict sparse QQ[T]/T^3 interpreter: every variable is explicitly bound.
No author worker, unit solve, source jet, new point, or higher degree is run.
"""
import argparse,hashlib,json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
from ambient_B18_T0_cokernel_Z_review_20261104_v63 import add as dadd,mul as dmul,inv as dinv,scale as dscale

ROOT=Path(__file__).resolve().parents[1]
AUTHOR='AMBIENT_B18_R_N_PROJECTION_NORM_20261104_v63'
ID=AUTHOR+'_REVIEW'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def poly(a):return {i:Q(x) for i,x in enumerate(a) if Q(x)}
def dense(a):return [a.get(i,Q(0)) for i in range(3)]
def add(a,b):
    c=dict(a)
    for k,x in b.items():
        c[k]=c.get(k,Q(0))+x
        if not c[k]:del c[k]
    return c
def scale(a,c):return {k:x*c for k,x in a.items() if x*c}
def mul(a,b):
    c={}
    for i,x in a.items():
        for j,y in b.items():
            if i+j<=2:c[i+j]=c.get(i+j,Q(0))+x*y
    return {k:x for k,x in c.items() if x}
def inverse(a):
    c=a.get(0,Q(0));assert c
    h=scale({i:x for i,x in a.items() if i},-1/c)
    return scale(add(add({0:Q(1)},h),mul(h,h)),1/c)
def closure(nodes,roots):
    seen=set();todo=list(roots)
    while todo:
        i=todo.pop()
        if i in seen:continue
        seen.add(i);op,*a=nodes[i]
        if op in ('a','m'):todo.extend(a)
        elif op in ('n','i'):todo.append(a[0])
        else:assert op in ('c','v')
    return sorted(seen)
def evaluate(nodes,needed,bindings):
    values={}
    for i in needed:
        op,*a=nodes[i]
        if op=='v':assert i in bindings,('UNBOUND_VARIABLE',i,nodes[i]);v=bindings[i]
        elif op=='c':v=poly([Q(a[0],a[1])])
        elif op=='a':v=add(values[a[0]],values[a[1]])
        elif op=='m':v=mul(values[a[0]],values[a[1]])
        elif op=='n':v=scale(values[a[0]],Q(-1))
        else:assert op=='i';v=inverse(values[a[0]])
        values[i]=v
    return values
def jacobian(nodes,needed,unknown,equations,bindings):
    vals={};indices={v:j for j,v in enumerate(unknown)}
    for i in needed:
        op,*a=nodes[i]
        if op=='v':
            assert i in bindings
            v=(bindings[i].get(0,Q(0)),{indices[i]:Q(1)} if i in indices else {})
        elif op=='c':v=(Q(a[0],a[1]),{})
        elif op=='a':v=dadd(vals[a[0]],vals[a[1]])
        elif op=='m':v=dmul(vals[a[0]],vals[a[1]])
        elif op=='n':v=dscale(vals[a[0]],Q(-1))
        else:assert op=='i';v=dinv(vals[a[0]])
        vals[i]=v
    assert all(vals[i][0]==0 for i in equations)
    return [[vals[i][1].get(j,Q(0)) for j in range(len(unknown))] for i in equations]
def check_stages(nodes,u,eq,J,external,roots,stages,savedinverse):
    assert len(u)==len(eq)==len(J) and len(set(u))==len(u)
    needed=closure(nodes,eq+roots)
    bind=dict(external);bind.update({v:{} for v in u})
    assert {i for i in needed if nodes[i][0]=='v'}<=set(bind)
    assert jacobian(nodes,needed,u,eq,bind)==J
    n=len(u);inv=[list(map(Q,r)) for r in savedinverse]
    assert len(inv)==n and all(len(r)==n for r in inv)
    for i in range(n):
        for j in range(n):
            assert sum((J[i][k]*inv[k][j] for k in range(n)),Q(0))==int(i==j)
            assert sum((inv[i][k]*J[k][j] for k in range(n)),Q(0))==int(i==j)
    assert len(stages)==2
    for degree,stage in enumerate(stages,1):
        assert stage['degree']==degree
        before=evaluate(nodes,needed,bind)
        residual=[before[i].get(degree,Q(0)) for i in eq]
        assert residual==list(map(Q,stage['raw_equation_residual']))
        correction=list(map(Q,stage['correction']))
        for i in range(n):assert sum((J[i][j]*correction[j] for j in range(n)),Q(0))==-residual[i]
        new={i:poly(a) for i,a in stage['all_unknown_coefficients']}
        assert set(new)==set(u)
        for i,c in zip(u,correction):
            assert new[i]==add(bind[i],{degree:c} if c else {})
        bind.update(new)
        after=evaluate(nodes,needed,bind)
        for i,stored in zip(eq,stage['all_equation_coefficients_after_correction']):
            assert dense(after[i])==list(map(Q,stored))
            assert not any(c for k,c in after[i].items() if k<=degree)
    return after,needed

def main(out,Rreview):
    started=time.monotonic();signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
    out.mkdir(parents=True,exist_ok=False)
    folder=ROOT/'reports'/AUTHOR;primary=read(folder/'report.json');assert primary['all_passed']
    hashes={}
    for p,h in primary['input_sha256'].items():assert sha(ROOT/p)==h;hashes[p]=h
    for p,h in primary['evidence_sha256'].items():assert sha(folder/p)==h;hashes[str((folder/p).relative_to(ROOT))]=h
    hashes[str((folder/'report.json').relative_to(ROOT))]=sha(folder/'report.json')
    # Actual R/group/R0 binding is a separate REQUIRED gate, not inferred from
    # the norm numbers or from physical nonmembership in N.
    rg=read(Rreview);assert rg['all_passed']
    Rpath=ROOT/'reports/AMBIENT_B18_R_ZERO_SLICE_IMAGE_20261104_v63/actual-R-zero-slice-image-DAG.json'
    assert rg['input_sha256'][str(Rpath.relative_to(ROOT))]==sha(Rpath)
    gates=[Rreview,
      ROOT/'reports/AMBIENT_B18_N_SAME_R0_BINDING_20261104_v63_REVIEW/report.json',
      ROOT/'reports/AMBIENT_B18_N_LITERAL_MODULE_CONE_20261104_v63_REVIEW/report.json',
      ROOT/'reports/AMBIENT_B18_T1_T2_QUADRATIC_PHI_20261104_v63_REVIEW/report.json',
      ROOT/'reports/AMBIENT_B18_W1_TAIL_BINDING_20261104_v63_REVIEW/report.json',
      ROOT/'reports/AMBIENT_R_COMPLEMENT_20261103_v63_REVIEW.json']
    for p in gates:assert read(p)['all_passed'];hashes[str(p.relative_to(ROOT))]=sha(p)
    scope=ROOT/'reports/AMBIENT_B18_W1_POLYNOMIAL_LAURENT_SCOPE_20261104_v63.json'
    assert read(scope)['all_passed'];hashes[str(scope.relative_to(ROOT))]=sha(scope)
    R=read(Rpath)
    Npath=ROOT/'reports/AMBIENT_B18_N_SAME_R0_BINDING_20261104_v63/actual-N-R0-binding-DAG.json'
    N=read(Npath);data=read(folder/'actual-R-N-projection-norm.json')
    rows=R['same_actual_R0_rows9'];assert rows==N['actual_same_R0_rows9']==data['same_actual_R0_rows9']
    assert rows==[27,108,109,112,113,117,118,119,61]
    ru=R['group19_unknown_nodes'];re=R['all19_cut_equations'];rroots=R['actual_R0_readback_nodes9']
    assert len(ru)==19
    rv,rneeded=check_stages(R['DAG_scalar_nodes'],ru,re,
      [list(map(Q,a)) for a in R['all19_cut_Jacobian']],
      {R['frozen_actual_R_curve_parameter_node']:{1:Q(1)}},rroots,
      data['R_group19_unit_stages'],data['R_group19_unit_inverse'])
    rcoeff=[dense(rv[i]) for i in rroots]
    assert rcoeff==[list(map(Q,a)) for a in data['all9_actual_R0_curve_coefficients_through2']]
    assert [p[1] for p in rcoeff]==list(map(Q,R['actual_R0_slice_curve_velocity9']))
    indices=N['chosen7_R0_indices']+[N['Cohen_even8_R0_index']]
    assert indices==data['N_Cohen8_R0_indices']
    assert [rows[i] for i in indices]==[108,109,112,113,117,118,119,27]
    ext={v:rv[rroots[i]] for v,i in zip(N['actual_Cohen_R8_parameter_nodes'],indices)}
    assert ext=={i:poly(a) for i,a in data['N_Cohen8_parameter_series']}
    nu=N['nu_tracefree_square_node'];aa=N['ninth_actual_even_a_node'];z1=N['ninth_actual_odd_b_node']
    ns=N['all56_coefficient_unknown_nodes'];ne=N['all56_actual_coefficient_Hensel_equations'];assert len(ns)==56
    nv,nneeded=check_stages(N['DAG_scalar_nodes'],ns,ne,
      [list(map(Q,a)) for a in N['unit56_Jacobian']],ext,[nu,aa,z1],
      data['N56_unit_stages'],data['N56_unit_inverse'])
    v=add(rv[rroots[N['ninth_R0_index']]],scale(nv[aa],Q(-1)))
    assert dense(nv[nu])==list(map(Q,data['N_tracefree_nu_coefficients_through2']))==[0,0,256]
    assert dense(nv[aa])==list(map(Q,data['N_ninth_even_a_coefficients_through2']))==[0,0,-472]
    assert dense(nv[z1])==list(map(Q,data['N_z1_b61_coefficients_through2']))
    assert nv[z1].get(0,Q(0))==0 and nv[z1][1]==-23
    assert rcoeff[N['ninth_R0_index']]==[0,0,-240]
    assert dense(v)==list(map(Q,data['actual_v_rho61_minus_a_coefficients_through2']))==[0,0,232]
    assert nv[nu].get(0,Q(0))==nv[nu].get(1,Q(0))==v.get(0,Q(0))==v.get(1,Q(0))==0
    norm4=v[2]**2-nv[nu][2]*nv[z1][1]**2
    assert norm4==-81600==Q(data['projection_norm_coefficient_T4'])
    assert data['first_nonzero_projection_norm_coefficient_degree']==4
    proof=dict(required_actual_R_group_R0_binding_review_passed=True,
      same_physical_R0_nine_and_Cohen8_rows_verified=True,
      independent_strict_QQ_T_mod_T3_interpreter=True,
      no_unbound_variables_defaulted_to_zero=True,
      actual_19_and56_numeric_Jacobians_independently_rebuilt=True,
      both_saved_unit_inverses_checked_left_and_right=True,
      every_raw_residual_correction_and_all75_equations_through2_verified=True,
      R_needed_nodes=len(rneeded),N_needed_nodes=len(nneeded),
      actual_nu2=256,actual_z1_first=-23,actual_a61_second=-472,
      actual_rho61_second=-240,actual_v_second=232,
      exact_projection_norm_T4_coefficient=str(norm4),
      higher_unknown_coefficients_do_not_enter_T4_by_verified_orders=True,
      R0_generic_point_misses_D4_support=True,
      all_positive_pieces_vanish_at_that_R0_generic_point=True,
      actual_D0_and_T0_generic_R_length2=True,
      no_nil_character_or_different_base_assumption=True,
      no_T0_flatness_entire_support_finite_socle_or_depth_claim=True,
      no_author_worker_unit_solve_new_source_jet_or_degree_above2=True)
    (out/'independent-R-N-norm-proof.json').write_text(json.dumps(proof,indent=2)+'\n')
    for p in (Path(__file__),ROOT/f'research/notes/{ID}.md',ROOT/'code/ambient_B18_T0_cokernel_Z_review_20261104_v63.py'):
        hashes[str(p.relative_to(ROOT))]=sha(p)
    report=dict(id=ID,all_passed=True,input_sha256=hashes,
      evidence_sha256={'independent-R-N-norm-proof.json':sha(out/'independent-R-N-norm-proof.json')},
      actual_projection_norm_and_generic_R_length2_independently_accepted=True,
      finite_socle_or_full_depth_not_decided=True,
      seconds=time.monotonic()-started,max_rss_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
    print(json.dumps(report))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--R-binding-review',type=Path,required=True)
    a=p.parse_args();main(a.output,a.R_binding_review.resolve())
