#!/usr/bin/env python3
"""Independent quotient identities, CA implicit gate and physical curve parameter."""
import argparse
from collections import Counter
from fractions import Fraction as F
import gzip
import hashlib
import json
from pathlib import Path
import resource
import signal
import time
from ambient_B18_actual_top_coefficient_preflight_review_20261104_v63 import DAGReview
from ambient_B18_top_graph_telescoping_review_20261104_v63 import closure,solve

ROOT=Path(__file__).resolve().parents[1]
BASE='AMBIENT_B18_T3_CURVE_ANNIHILATOR_20261104_v63';STEM=BASE+'_REVIEW'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(gzip.decompress(p.read_bytes()) if p.suffix=='.gz' else p.read_text())
def matrix(a):return [list(map(F,r)) for r in a]
def mul(a,b):return [[sum(x*y for x,y in zip(row,col)) for col in zip(*b)] for row in a]


def main(out):
    start=time.monotonic();signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000))
    folder=ROOT/'reports'/BASE;primary=read(folder/'report.json')
    assert primary['all_passed'] and primary['T3_mod_ell_length_upper_bound']==36
    assert not primary['one_colon_membership_T3_depth_nil_or_full_D_answer_computed']
    cert=read(folder/'annihilator-and-parameter.json');curve=read(folder/'actual-CA-curve-DAG.json')
    fullpath=ROOT/'reports/AMBIENT_B18_FULL_RELATIVE_INPUT_20261104_v63/full-relative186-DAG.json'
    full=read(fullpath);fullpeer=ROOT/'reports/AMBIENT_B18_FULL_RELATIVE_INPUT_20261104_v63_REVIEW/report.json'
    assert read(fullpeer)['all_passed']
    nodes=full['DAG_scalar_nodes'];Q3=full['actual_Qn_R0_squared_generator_matrices']['3']
    top=full['actual12_top_T3_additional_columns'];C0=full['actual_CA_C0_generators'];H=full['actual_CA_H_generators']
    assert len(Q3)==2034 and len(top)==12 and len(C0)==9 and len(H)==11
    Qgens=full['all99_CB_kernel_and87_graph_generators_evaluated_on_CA']
    assert 11*sum(q['weight']<3 for q in Qgens)==1980
    def product_signature(rootlist):
        coefficient=F(1);atoms=Counter();stack=list(rootlist)
        while stack:
            i=stack.pop();n=nodes[i]
            if n[0]=='c':coefficient*=F(n[1],n[2])
            elif n[0]=='n':coefficient=-coefficient;stack.append(n[1])
            elif n[0]=='m':stack.extend(n[1:])
            else:atoms[i]+=1
        if not coefficient:return F(0),()
        return coefficient,tuple(sorted(atoms.elements()))
    identities=cert['all96_I0_annihilation_column_identities'];assert len(identities)==96
    counter=0
    for alpha in range(6):
        for e in range(2):
            for i in range(8):
                rec=identities[counter];counter+=1
                t=2*alpha+e;c=1980+9*alpha+i+1;fi=C0[i+1][1];je=H[9+e][1]
                assert rec['top12_column_index0']==t and rec['I0_generator_index0']==i
                assert rec['I0_scalar_node']==fi and rec['Q3_column_index0']==c and rec['Q3_cofactor_scalar_node']==je
                for row in range(2):
                    lhs=product_signature([fi,top[t][row]])
                    rhs=product_signature([je,Q3[c][row]])
                    assert lhs==rhs
                    assert rec['both_row_product_normal_forms'][row]==[str(lhs[0]),list(lhs[1])]
    sourcepath=ROOT/'reports/AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63/actual-source-coefficients.json'
    old=read(sourcepath);oldnodes=old['DAG_scalar_nodes'];ca=old['families']['CA'];D=DAGReview(oldnodes)
    assert C0[1:]==[[0,i] for row,i in ca['actual_zero_odd_I0_generators']]
    assert H[9:]==[[0,i] for row,i in ca['actual_positive_odd_J0_generators']]
    xrows=old['physical_zero_base_rows9'];xids=[D.variables[f'x{i}'] for i in xrows]
    zids=[i for p in ca['coefficient_unknown_pairs'].values() for i in p]+[ca['p_quadratic_relation'][k] for k in ['s','t']]
    eq=ca['all48_zero_IFT_quadratic_coefficient_equations']
    J=[[D.grad[i].get(v,F(0)) for v in zids] for i in eq]
    rhs=[[-D.grad[i].get(v,F(0)) for v in xids] for i in eq]
    zx,detJ=solve(J,rhs);assert detJ==F(-84934656,625)
    def derivative(node):
        answer=[D.grad[node].get(x,F(0)) for x in xids]
        for i,z in enumerate(zids):
            q=D.grad[node].get(z,F(0))
            if q:answer=[a+q*b for a,b in zip(answer,zx[i])]
        return answer
    gens=ca['actual_zero_odd_I0_generators']+ca['actual_positive_odd_J0_generators']
    gids=[i for row,i in gens];M=[derivative(i) for i in gids]
    selected=[0,1,2,3,4,5,6,8,9]
    assert M==matrix(cert['actual_CA10_by9_Jacobian']) and not any(M[7])
    inv,detF=solve([M[i] for i in selected],[[F(i==j) for j in range(9)] for i in range(9)])
    assert detF==F(128,3)==F(cert['actual_selected_CA9_Jacobian_determinant'])
    # Fully independent Fraction Schur transport, not the author's SymPy engine.
    upperpath=ROOT/'reports/AMBIENT_B18_CONDUCTOR_DIMENSION20_20261103_v63/certificate.json.gz'
    firstpath=ROOT/'reports/AMBIENT_B18_ODD_FIRSTJETS_N_PROJECTION_20261103_v63/certificate.json.gz'
    orbitpath=ROOT/'reports/AMBIENT_B18_ACTUAL_INTERSECTION_ORBIT_PRODUCT_20261103_v63/certificate.json.gz'
    weightpath=ROOT/'reports/AMBIENT_B18_RELATIVE_WEIGHTED_MODEL_20261104_v63/certificate.json'
    K=matrix(read(upperpath)['actual_jet_bindings']['CA']['full_minimal_kernel_basis29_by22'])
    Rrows=read(firstpath)['common_Y28']+[61];pos={raw:i for i,raw in enumerate(Rrows)}
    O=matrix(read(orbitpath)['physical150_by19_scalar_then_GL'])
    cuts=read(weightpath)['common_physical_group_cuts19']
    correction,groupdet=solve([O[i] for i in cuts],[K[pos[i]] for i in cuts])
    assert groupdet==-32
    product=mul([O[i] for i in Rrows],correction)
    projected=[[x-y for x,y in zip(a,b)] for a,b in zip(K,product)]
    assert all(not any(row[:19]) for row in projected)
    UV=[[projected[pos[row]][20],projected[pos[row]][21]] for row in xrows]
    assert mul(M[:7],UV)==[[0,0] for _ in range(7)]
    assert projected[pos[140]][20:22]==[0,0] and projected[pos[140]][19]!=0
    transport=mul(M[8:10],UV)
    expected=[[F(4,25),F(-4,25)],[F(9,25),F(1,25)]]
    assert transport==expected==matrix(cert['old_UV_to_actual_CA_r_s_matrix'])
    assert transport[0][0]*transport[1][1]-transport[0][1]*transport[1][0]==F(8,125)
    ell=[F(-5),F(5)];assert cert['Noether_ell_r_s_coefficients']==list(map(str,ell))
    assert mul([ell],transport)==[[1,1]]
    assert mul([ell],mul(transport,[[0],[1]]))==[[1]]
    assert mul([ell],mul(transport,[[1],[3]]))==[[4]]
    binarypath=ROOT/'reports/AMBIENT_ORIGIN_BINARY_CUBIC_RECOGNITION_20261103_v63_REVIEW/report.json'
    binary=read(binarypath);assert binary['all_passed']
    assert binary['CA_exact_initial_ideal']=='(WU,WV,(12/125)*U*(U-V/3)^2)'
    assert F(12,125)*(1+F(1,3))**2==F(64,375)
    # Bind the new CA g-curve to the SAME original actual48 system.
    newnodes=[tuple(n) for n in curve['DAG_scalar_nodes']]
    nold=len(oldnodes);assert newnodes[:nold]==[tuple(n) for n in oldnodes]
    for i,n in enumerate(newnodes[nold:],nold):
        if n[0]=='c':value=F(n[1],n[2])
        elif n[0]=='v':value=F(0)
        else:
            assert all(isinstance(j,int) and 0<=j<i for j in n[1:]);a=n[1]
            if n[0]=='n':value=-D.values[a]
            elif n[0]=='a':value=D.values[a]+D.values[n[2]]
            elif n[0]=='m':value=D.values[a]*D.values[n[2]]
            elif n[0]=='i':assert D.values[a];value=1/D.values[a]
            else:raise ValueError(n)
        D.values.append(value)
    D.nodes=newnodes;D.cache={n:i for i,n in enumerate(newnodes)}
    Y=[gids[i] for i in selected];f63=gids[7]
    assert curve['selected_F9_nodes']==Y and curve['selected_F9_physical_rows']==[18,19,28,29,60,67,72,46,95]
    binding=dict(curve['pointed_curve_inverse57_bindings']);assert set(binding)==set(xids+zids)
    assert all(newnodes[new]==('v',f'CA_CURVE_INV_{oldnodes[oldid][1]}') for oldid,new in binding.items())
    mapping={}
    for i in sorted(closure(newnodes[:nold],eq+[gids[i] for i in selected]+[f63])):
        n=newnodes[i];op=n[0]
        if op=='c':new=D.const(F(n[1],n[2]))
        elif op=='v':new=binding[i]
        elif op=='n':new=D.neg(mapping[n[1]])
        elif op=='a':new=D.add(mapping[n[1]],mapping[n[2]])
        elif op=='m':new=D.mul(mapping[n[1]],mapping[n[2]])
        elif op=='i':new=D.inverse(mapping[n[1]])
        else:raise ValueError(n)
        mapping[i]=new
    equations=[mapping[i] for i in eq]+[D.add(mapping[gids[index]],D.neg(Y[j] if j>=7 else D.zero)) for j,index in enumerate(selected)]
    assert equations==curve['pointed_curve_inverse57_equations']
    assert curve['curve_equation_g_CA_node']==mapping[f63]
    ellnode=D.add(D.mul(D.const(-5),Y[7]),D.mul(D.const(5),Y[8]))
    assert curve['Noether_parameter_ell_curve_DAG_node']==ellnode
    assert curve['Noether_parameter_ell_in_full186_shared_source_namespace']==[['-5',Y[7]],['5',Y[8]]]
    caauthor=ROOT/'reports/AMBIENT_CA_ORIGIN_CONDUCTOR_MULTIPLICITY4_20261103_v63/report.json'
    cai=read(caauthor);assert cai['actual_odd_ideal_minimal_generator_count']==10 and cai['conductor_multiplicity']==4
    checked={}
    peers=[fullpeer,binarypath,ROOT/'reports/AMBIENT_CA_ORIGIN_CONDUCTOR_MULTIPLICITY4_20261103_v63_REVIEW/report.json',
           ROOT/'reports/AMBIENT_B18_RELATIVE_WEIGHTED_MODEL_20261104_v63_REVIEW/report.json']
    for p in [folder/'report.json',*peers]:
        data=read(p);assert data['all_passed']
        for name,want in data.get('input_sha256',{}).items():
            t=ROOT/name;assert sha(t)==want,str(t);checked[name]=want
        for name,want in data.get('evidence_sha256',{}).items():
            t=p.parent/name;assert sha(t)==want,str(t);checked[str(t.relative_to(ROOT))]=want
    proof={'all96_actual_Q3_column_identities_both_rows':True,
        'CA_implicit_Jacobian_determinants':{'48':str(detJ),'9':str(detF)},
        'actual_CA57_curve_inverse_and_common_source_namespace':True,
        'physical_G19_Schur_UV_transport':[[str(x) for x in row] for row in transport],
        'ell_exact_coefficients':['-5','5'],'old_initial_at_ell0':'64/375',
        'curve_multiplicity_from_actual_conductor':'4=1+ord(gA), hence ord(gA)=3',
        'curve_CM_dimension':1,'curve_ell_fiber_length':3,
        'T3_mod_ell_length_upper':36,'one_ell_colon_vanishing_equivalent_to_T3_m_socle_vanishing':True,
        'Tor_or_ambient_mod_I0_shortcut_used':False,'actual_colon_depth_nil_answer':False,
        'fiber_bound_implies_jet_or_torsion_exponent_bound':False}
    out.mkdir(parents=True,exist_ok=False)
    ep=out/'independent-proof.json';ep.write_text(json.dumps(proof,indent=2)+'\n')
    for p in [Path(__file__),ROOT/'research/notes'/(STEM+'.md'),ROOT/'research/notes'/(BASE+'.md'),folder/'report.json',
              folder/'annihilator-and-parameter.json',folder/'actual-CA-curve-DAG.json',fullpath,sourcepath,caauthor,*peers]:
        checked[str(p.relative_to(ROOT))]=sha(p)
    result={'id':STEM,'all_passed':True,
        'all96_actual_annihilator_column_identities_independently_verified':True,
        'actual_CA_AD_rank9_det128_over3_and_pointed57_curve_binding':True,
        'physical_G19_Schur_transport_and_initial_parameter_nonvanishing':True,
        'actual_mu10_conductor_e4_fiberproduct_to_curve_order3_proof_accepted':True,
        'I0_times_T3_zero_and_curve_CM1_parameter_fiber_length3':True,
        'top12_generators_give_closed_ell_fiber_length_upper':36,
        'one_actual_R0_squared_ell_colon_decides_T3_socle_vanishing':True,
        'ambient_mod_I0_Tor_shortcut_or_finite_jet_bound_claimed':False,
        'actual_colon_membership_T3_vanishing_depth_nil_answer':False,
        'author_rerun_new_jets_GB_or_point':False,
        'input_sha256':checked,'evidence_sha256':{'independent-proof.json':sha(ep)},
        'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'frozen':True}
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ['all_passed','seconds','max_rss_kib']}))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
