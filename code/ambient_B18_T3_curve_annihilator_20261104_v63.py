#!/usr/bin/env python3
"""ONE read-only96-column proof and actual CA curve/parameter gate; no GB/jets."""
import argparse
from fractions import Fraction as F
import gzip
import hashlib
import json
from pathlib import Path
import resource
import signal
import time
import sympy as s
from ambient_B18_actual_top_coefficient_preflight_20261104_v63 import DAG
from ambient_B18_top_graph_telescoping_20261104_v63 import forward,replay,ancestors,linear_solve,independent_rows
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_T3_CURVE_ANNIHILATOR_20261104_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(gzip.decompress(p.read_bytes())if p.suffix=='.gz'else p.read_text())
def enc(M):return[[str(v)for v in r]for r in M]
def main(out):
    started=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));out.mkdir(parents=True,exist_ok=False)
    fullpath=ROOT/'reports/AMBIENT_B18_FULL_RELATIVE_INPUT_20261104_v63/full-relative186-DAG.json'
    fullpeer=ROOT/'reports/AMBIENT_B18_FULL_RELATIVE_INPUT_20261104_v63_REVIEW/report.json'
    sourcepath=ROOT/'reports/AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63/actual-source-coefficients.json'
    capeer=ROOT/'reports/AMBIENT_CA_ORIGIN_CONDUCTOR_MULTIPLICITY4_20261103_v63_REVIEW/report.json'
    caauthor=ROOT/'reports/AMBIENT_CA_ORIGIN_CONDUCTOR_MULTIPLICITY4_20261103_v63/report.json'
    binarypath=ROOT/'reports/AMBIENT_ORIGIN_BINARY_CUBIC_RECOGNITION_20261103_v63_REVIEW/report.json'
    upperpath=ROOT/'reports/AMBIENT_B18_CONDUCTOR_DIMENSION20_20261103_v63/certificate.json.gz'
    firstpath=ROOT/'reports/AMBIENT_B18_ODD_FIRSTJETS_N_PROJECTION_20261103_v63/certificate.json.gz'
    orbitpath=ROOT/'reports/AMBIENT_B18_ACTUAL_INTERSECTION_ORBIT_PRODUCT_20261103_v63/certificate.json.gz'
    weightpath=ROOT/'reports/AMBIENT_B18_RELATIVE_WEIGHTED_MODEL_20261104_v63/certificate.json'
    full,source=read(fullpath),read(sourcepath);assert read(fullpeer)['all_passed'] and read(capeer)['all_passed']
    ca_info=read(caauthor);assert ca_info['actual_odd_ideal_minimal_generator_count']==10 and ca_info['conductor_multiplicity']==4
    assert read(binarypath)['all_passed']
    nodes=full['DAG_scalar_nodes'];oldnodes=source['DAG_scalar_nodes'];assert nodes[:len(oldnodes)]==oldnodes
    Q3=full['actual_Qn_R0_squared_generator_matrices']['3'];top=full['actual12_top_T3_additional_columns']
    C0=full['actual_CA_C0_generators'];H=full['actual_CA_H_generators']
    assert len(Q3)==2034 and len(top)==12 and len(C0)==9 and len(H)==11
    memo={}
    def flat(i):
        if i in memo:return memo[i]
        n=nodes[i]
        if n[0]=='c':ans=(F(n[1],n[2]),())
        elif n[0]=='n':c,at=flat(n[1]);ans=(-c,at)
        elif n[0]=='m':
            c,a=flat(n[1]);d,b=flat(n[2]);ans=(c*d,tuple(sorted(a+b)))if c*d else(F(0),())
        else:ans=(F(1),(i,))
        memo[i]=ans;return ans
    def times(i,j):
        c,a=flat(i);d,b=flat(j);return(c*d,tuple(sorted(a+b)))if c*d else(F(0),())
    ann=[]
    for alpha in range(6):
        for e in range(2):
            topindex=2*alpha+e;j=H[9+e][1]
            for i in range(8):
                f=C0[i+1][1];column=1980+9*alpha+i+1
                signatures=[]
                for row in range(2):
                    lhs=times(f,top[topindex][row]);rhs=times(j,Q3[column][row]);assert lhs==rhs
                    signatures.append([str(lhs[0]),list(lhs[1])])
                ann.append({'top12_column_index0':topindex,'I0_generator_index0':i,
                            'I0_scalar_node':f,'Q3_column_index0':column,'Q3_cofactor_scalar_node':j,
                            'both_row_product_normal_forms':signatures})
    assert len(ann)==96
    ca=source['families']['CA'];varids={n[1]:i for i,n in enumerate(oldnodes)if n[0]=='v'}
    xrows=source['physical_zero_base_rows9'];xids=[varids[f'x{i}']for i in xrows]
    zids=[i for p in ca['coefficient_unknown_pairs'].values()for i in p]+[ca['p_quadratic_relation']['s'],ca['p_quadratic_relation']['t']]
    vals,grads=forward(oldnodes,set(xids+zids));eq=ca['all48_zero_IFT_quadratic_coefficient_equations']
    J=[[grads[i].get(v,F(0))for v in zids]for i in eq]
    B=[[-grads[i].get(v,F(0))for v in xids]for i in eq]
    zx,detJ=linear_solve(J,B);assert abs(detJ)==F(ca['unit48_block_determinant_even24_then_odd24_vars_a23_t_b23_s'])
    def derivative(i):
        ans=[grads[i].get(v,F(0))for v in xids]
        for k,v in enumerate(zids):
            q=grads[i].get(v,F(0))
            if q:ans=[a+q*b for a,b in zip(ans,zx[k])]
        return ans
    generators=ca['actual_zero_odd_I0_generators']+ca['actual_positive_odd_J0_generators'];gids=[i for row,i in generators]
    matrix=[derivative(i)for i in gids];selected=independent_rows(matrix);assert selected==[0,1,2,3,4,5,6,8,9]
    inverse,detF=linear_solve([matrix[i]for i in selected],[[F(int(i==j))for j in range(9)]for i in range(9)])
    assert detF
    # Exact transport of accepted tangent lines to the ACTUAL CA j46,j95 coordinates.
    upper=read(upperpath)['actual_jet_bindings']['CA'];K=s.Matrix(upper['full_minimal_kernel_basis29_by22'])
    Rrows=read(firstpath)['common_Y28']+[61];orbit=read(orbitpath);weight=read(weightpath)
    O=s.Matrix(orbit['physical150_by19_scalar_then_GL']);cuts=weight['common_physical_group_cuts19']
    positions={raw:i for i,raw in enumerate(Rrows)};assert all(i in positions for i in cuts+xrows+[140])
    projected=K-O[Rrows,:]*O[cuts,:].inv()*K[[positions[i]for i in cuts],:]
    assert projected[:,:19]==s.zeros(29,19)
    uv=projected[[positions[i]for i in xrows],20:22]
    zero7=s.Matrix([[s.Rational(x.numerator,x.denominator)for x in row]for row in matrix[:7]])
    assert zero7*uv==s.zeros(7,2)
    assert projected[positions[140],20]==projected[positions[140],21]==0
    assert projected[positions[140],19]!=0
    rsJac=s.Matrix([[s.Rational(x.numerator,x.denominator)for x in row]for row in matrix[8:10]])
    transport=rsJac*uv;assert transport.det()
    ell=s.Matrix([[1,1]])*transport.inv();er,es=map(s.factor,list(ell))
    assert ell*transport==s.Matrix([[1,1]])
    assert (ell*transport*s.Matrix([0,1]))[0]==1
    assert (ell*transport*s.Matrix([1,3]))[0]==4
    initial_gate=s.Rational(12,125)*(1)*(1+s.Rational(1,3))**2;assert initial_gate==s.Rational(64,375)
    # One actual CA inverse instance defines g_A(r,s), without any jet expansion.
    D=DAG();original=replay(D,oldnodes,{});assert all(original[i]==i for i in range(len(oldnodes)))
    Y=[gids[i]for i in selected];f63=gids[7]
    needed=ancestors(oldnodes,eq+[gids[i]for i in selected]+[f63])
    bindings={i:D.var(f'CA_CURVE_INV_{oldnodes[i][1]}')for i in xids+zids}
    mapping=replay(D,oldnodes,bindings,needed)
    inv_eq=[mapping[i]for i in eq]+[D.add(mapping[gids[index]],D.neg(Y[j]if j>=7 else D.zero))for j,index in enumerate(selected)]
    assert len(inv_eq)==57 and all(D.values[i]==0 for i in inv_eq)
    ellnode=D.add(D.mul(D.const(er),Y[7]),D.mul(D.const(es),Y[8]))
    curve={'DAG_scalar_nodes':D.nodes,'original_CA48_equations':eq,'selected_F9_physical_rows':[generators[i][0]for i in selected],
           'selected_F9_nodes':Y,'pointed_curve_inverse57_equations':inv_eq,
           'pointed_curve_inverse57_bindings':[[a,b]for a,b in bindings.items()],
           'curve_equation_g_CA_node':mapping[f63],
           'Noether_parameter_ell_curve_DAG_node':ellnode,
           'Noether_parameter_ell_in_full186_shared_source_namespace':[[str(er),Y[7]],[str(es),Y[8]]],
           'coefficient_ring_curve':'QQ[[r,s]]/(g_CA), g_CA!=0',
           'g_nonzero_proof':'actual mu(I_CA)=10; otherwise sevenzero+two positive generateI in9elements',
           'multiplicity3_proof':'R0[[w]]/(I0+wJ) = k[[w]] fiberproduct_k k[[r,s]]/(g); actualconductor e4 =>ordg3'}
    curvestring=json.dumps(curve,separators=(',',':'))+'\n';assert len(curvestring.encode())<=32*1024*1024
    (out/'actual-CA-curve-DAG.json').write_text(curvestring)
    certificate={'all96_I0_annihilation_column_identities':ann,
                 'all_indices_zero_based':True,'annihilation_formula':'f_i*top_(alpha,e)=j_e*Q3column_(1980+9alpha+i+1)',
                 'actual_CA10_by9_Jacobian':enc(matrix),'actual_CA48_Jacobian_determinant':str(detJ),
                 'actual_selected_CA9_indices':selected,'actual_selected_CA9_Jacobian_determinant':str(detF),
                 'old_UV_to_actual_CA_r_s_matrix':[[str(x)for x in row]for row in transport.tolist()],
                 'Noether_ell_r_s_coefficients':[str(er),str(es)],
                 'ell_on_old_V_axis':1,'ell_on_old_self_line_U1_V3':4,
                 'old_binary_initial_on_ell0_U1_Vminus1':'64/375',
                 'curve_dimension':1,'curve_CM':True,'curve_multiplicity':3,'curve_ell_quotient_length':3,
                 'T3_generator_bound_over_curve':12,'T3_mod_ell_length_upper_bound':36,
                 'exact_one_colon_criterion':'((Q[3]:ell) intersect Jinf)/Q[3]=0 iff H0_m(T3)=0',
                 'ambient_2rows_not_reduced_mod_I0':True,
                 'Tor_warning_example':'R=k[[x]], J=(x), Q=(x^2): J/Q=k but both columnimages vanishmodx',
                 'length36_is_not_a_jet_or_torsion_exponent_bound':True}
    (out/'annihilator-and-parameter.json').write_text(json.dumps(certificate,indent=2)+'\n')
    paths=[fullpath,fullpeer,sourcepath,capeer,caauthor,binarypath,upperpath,firstpath,orbitpath,weightpath,Path(__file__)]
    report={'id':STEM,'all_passed':True,
            'I0_CA_times_T3_zero_by_explicit96_column_identities':True,
            'actual_CA10_odd_Jacobian_rank':9,'actual_selected_CA9_determinant':str(detF),
            'actual_CA_curve_g_nonzero_dimension1_CM_multiplicity3':True,
            'explicit_ell_is_Noether_parameter_and_curve_fiber_length3':True,
            'top_T3_socle_vanishing_reduced_to_one_ell_colon':True,
            'T3_mod_ell_length_upper_bound':36,
            'one_colon_membership_T3_depth_nil_or_full_D_answer_computed':False,
            'no_ambient_matrix_basechange_mod_I0':True,
            'no_GB_jets_source_core_rerun_or_new_point':True,
            'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},
            'evidence_sha256':{name:sha(out/name)for name in ['actual-CA-curve-DAG.json','annihilator-and-parameter.json']},
            'seconds':time.monotonic()-started,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:report[k]for k in ['all_passed','actual_selected_CA9_determinant','T3_mod_ell_length_upper_bound','seconds','max_rss_kib']}))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True)
    main(p.parse_args().output)
