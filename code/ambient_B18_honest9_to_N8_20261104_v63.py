#!/usr/bin/env python3
"""ONE pure-unit CA96/N56/prefix160 actual honest9 -> same-R0 N8 map."""
import argparse
from fractions import Fraction as F
import hashlib,json
from pathlib import Path
import resource,signal,time
from ambient_B18_actual_top_coefficient_preflight_20261104_v63 import DAG,Quadratic
from ambient_B18_top_graph_telescoping_20261104_v63 import ancestors,forward,linear_solve,independent_rows
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_HONEST9_TO_N8_20261104_v63'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def enc(M):return [[str(v) for v in row] for row in M]
def main(out):
    started=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));out.mkdir(parents=True,exist_ok=False)
    cp=ROOT/'reports/AMBIENT_B18_N_LITERAL_MODULE_CONE_20261104_v63/actual-N-literal-cone-DAG.json'
    cr=cp.parent/'report.json';np=ROOT/'reports/AMBIENT_B18_N_SAME_R0_BINDING_20261104_v63/actual-N-R0-binding-DAG.json'
    sp=ROOT/'reports/AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63/actual-source-coefficients.json'
    up=ROOT/'reports/AMBIENT_B18_T3_HONEST10_UNIT_REDUCTION_20261104_v63/unit-reduced-presentation-DAG.json'
    cone,N,old,unit=map(read,[cp,np,sp,up]);assert read(cr)['all_passed']
    D=DAG();D.nodes=[tuple(n) for n in cone['DAG_scalar_nodes']];D.cache={n:i for i,n in enumerate(D.nodes)};D.values=[]
    for n in D.nodes:
        o=n[0]
        if o=='c':v=F(n[1],n[2])
        elif o=='v':v=F(0)
        elif o=='a':v=D.values[n[1]]+D.values[n[2]]
        elif o=='m':v=D.values[n[1]]*D.values[n[2]]
        elif o=='n':v=-D.values[n[1]]
        elif o=='i':assert D.values[n[1]];v=1/D.values[n[1]]
        else:raise ValueError(o)
        D.values.append(v)
    Q=object.__new__(Quadratic);Q.D=D;Q.S,Q.T=N['quadratic_trace_S_T_nodes'];Q.zero=(D.zero,D.zero);Q.one=(D.one,D.zero);Q.p=(D.zero,D.one)
    ca=old['families']['CA'];oldnodes=old['DAG_scalar_nodes'];var={n[1]:i for i,n in enumerate(oldnodes) if n[0]=='v'}
    xids=[var[f'x{r}'] for r in old['physical_zero_base_rows9']]
    az=[i for p in ca['coefficient_unknown_pairs'].values() for i in p]+[ca['p_quadratic_relation']['s'],ca['p_quadratic_relation']['t']]
    AA=[D.var(f'HN_CA_A_{i}') for i in az];BB=[D.var(f'HN_CA_B_{i}') for i in az];unknown96=AA+BB
    binding={i:tuple(p) for i,p in zip(xids,N['all9_R0_quadratic_p_readback_pairs'])}
    binding.update({i:(a,b) for i,a,b in zip(az,AA,BB)})
    sourcegens=ca['actual_zero_odd_I0_generators']+ca['actual_positive_odd_J0_generators'];selected=[0,1,2,3,4,5,6,8,9]
    Yrows=[sourcegens[i][0] for i in selected];retained=unit['remaining_original_free_generator_indices']
    honestrows=[None if i==0 else Yrows[i-1] for i in retained];assert honestrows==[None,18,19,28,29,67,72,46,95]
    evenroots={r:[ca['complete_normalized_physical_w_coefficients'][d][r][0] for d in range(4)] for r in honestrows[1:]}
    eq48=ca['all48_zero_IFT_quadratic_coefficient_equations'];needed=sorted(ancestors(oldnodes,eq48+[i for row in evenroots.values() for i in row]));mp={}
    for i in needed:
        n=oldnodes[i];o=n[0]
        if o=='c':v=Q.scalar(D.const(F(n[1],n[2])))
        elif o=='v':assert i in binding;v=binding[i]
        elif o=='a':v=Q.add(mp[n[1]],mp[n[2]])
        elif o=='m':v=Q.mul(mp[n[1]],mp[n[2]])
        elif o=='n':v=Q.neg(mp[n[1]])
        elif o=='i':v=Q.inv(mp[n[1]])
        else:raise ValueError(o)
        mp[i]=v
    E96=[i for e in eq48 for i in mp[e]];assert len(E96)==96 and all(D.values[i]==0 for i in E96)
    images=[]
    for row in honestrows[1:]:
        caeven=Q.total(mp[i] for i in evenroots[row]);p=Q.add(tuple(cone['all150_actual_quadratic_p_pairs'][row]),Q.neg(caeven))
        a=D.add(p[0],D.mul(D.const(F(1,2)),D.mul(Q.S,p[1])));b=p[1]
        assert D.values[a]==D.values[b]==0
        images.append({'honest_physical_row':row,'CA_even_scalar_basechanged_to_N_p_pair':list(caeven),'N_tracefree_even_odd':[a,b]})
    U56=cone['all56_coefficient_unknown_nodes'];base=cone['actual_R8_parameter_nodes8'];z=cone['z1_through_z8_actual_odd_coefficient_nodes'];EN=cone['all56_actual_coefficient_Hensel_equations']
    unknown160=U56+base+unknown96;E160=EN+z+E96;assert len(unknown160)==len(E160)==160
    _,grad=forward(D.nodes,set(unknown160))
    J=[[grad[e].get(v,F(0)) for v in unknown160] for e in E160]
    rhs=[[F(int(56<=i<64 and i-56==j)) for j in range(8)] for i in range(160)]
    dz,det160=linear_solve(J,rhs);assert det160
    J96=[[grad[e].get(v,F(0)) for v in unknown96] for e in E96]
    _,det96=linear_solve(J96,[[F(int(i==j)) for j in range(96)] for i in range(96)])
    expected48=F(ca['unit48_block_determinant_even24_then_odd24_vars_a23_t_b23_s']);assert det96==expected48**2
    assert det160==F(cone['unit56_minor'])*F(cone['unit8_odd_minor'])*det96
    origins=[[sum(grad[im['N_tracefree_even_odd'][1]].get(v,F(0))*dz[k][j] for k,v in enumerate(unknown160)) for j in range(8)] for im in images]
    # Prefix systems contain only specified pure-unit algebraic functions.
    roots=EN+z+E96+[im['N_tracefree_even_odd'][1] for im in images]
    need=sorted(ancestors(D.nodes,roots));assert {i for i in need if D.nodes[i][0]=='v'}<=set(unknown160)
    prefix=[];systems=[]
    def substitute(oldnodes,ids,binding):
        m={}
        for i in ids:
            n=oldnodes[i];o=n[0]
            if o=='v':assert i in binding;m[i]=binding[i]
            elif o=='c':m[i]=D.const(F(n[1],n[2]))
            elif o=='a':m[i]=D.add(m[n[1]],m[n[2]])
            elif o=='m':m[i]=D.mul(m[n[1]],m[n[2]])
            elif o=='n':m[i]=D.neg(m[n[1]])
            elif o=='i':m[i]=D.inverse(m[n[1]])
            else:raise ValueError(o)
        return m
    originalnodes=D.nodes[:]
    for level in range(1,8):
        binds={i:D.var(f'HN_PREFIX{level}_{D.nodes[i][1]}') for i in unknown160}
        m=substitute(originalnodes,need,binds);prefix.append(m)
        eq=[m[i] for i in EN]+[D.add(m[i],D.neg(z[j] if j<level else D.zero)) for j,i in enumerate(z)]+[m[i] for i in E96]
        assert len(eq)==160 and all(D.values[i]==0 for i in eq)
        systems.append({'level':level,'all160_variable_bindings':[[i,j] for i,j in binds.items()],'all160_unit_equations':eq,'unit160_determinant':str(det160)})
    telescope=[];columns=[[D.one]+[D.zero]*7];v=cone['actual_v_node']
    for j,im in enumerate(images):
        b=im['N_tracefree_even_odd'][1];previous=D.zero;c=[];div=[]
        for i in range(8):
            current=prefix[i][b] if i<7 else b;num=D.add(current,D.neg(previous))
            h=D.put(('r',num,z[i],f'HN_honest{j+1}:z{i+1}'),origins[j][i]);c.append(h)
            div.append({'node':h,'numerator':num,'coordinate':z[i],'current':current,'previous':previous,'origin':str(origins[j][i])});previous=current
        assert previous==b
        columns.append([D.add(im['N_tracefree_even_odd'][0],D.mul(c[0],v))]+c[1:])
        telescope.append({'honest_index':j+1,'odd_source_function':b,'all8_coordinate_divisions':div,'membership_coefficients8':c,
                          'exact_identity':'sum zi*ci telescopesfrom0 tob;unit160prefix8=original'})
    rank=len(independent_rows([[D.values[i] for i in c] for c in columns]));assert rank==8
    data={'DAG_scalar_nodes':D.nodes,'honest9_generator_physical_rows':honestrows,
          'N8_module_generator_order':cone['eight_N_module_generators'],'source_and_target_traces_not_identified':True,
          'CA_scalar48_unknown_nodes':az,'CA96_quadratic_N_basechange_unknown_nodes':unknown96,
          'CA96_coefficient_pairs_of_original48':[[i,[a,b]] for i,a,b in zip(az,AA,BB)],
          'all96_unit_basechange_equations':E96,'unit96_Jacobian':enc(J96),'unit96_determinant':str(det96),
          'all8_honest_nonunit_generators_on_N':images,
          'pure_unit160_variable_order':unknown160,'pure_unit160_equation_order':E160,
          'unit160_Jacobian':enc(J),'unit160_determinant':str(det160),
          'all7_prefix160_systems':systems,'all64_regular_coordinate_divisions':telescope,
          'actual_N_z8_coordinate_nodes':z,'actual_R0_v_node':v,
          'B_N_8_by9_columns':columns,'B_N_origin_matrix_columns':[[str(D.values[i]) for i in c] for c in columns],
          'B_N_origin_rank':rank,'first_honest_generator_maps_to_N_unit':True,
          'column_identity':'h=(a,b), b=sum zi ci; maps to(a+c1*v)n0+sum_i>=2 ci*n_(i-1),usingv=z1tau',
          'no_old_r_g_function_composition_or_opaque_branch_pullback':True,
          'nil_correction_row_transition9_and_T3_solution_not_constructed':True}
    encoded=json.dumps(data,separators=(',',':'))+'\n';assert len(encoded.encode())<32*1024*1024
    (out/'honest9-to-N8-DAG.json').write_text(encoded)
    paths=[Path(__file__),cp,cr,np,sp,up]
    rep={'id':STEM,'all_passed':True,'status':'actual_same_R0_Honest9_to_N8_map_no_nil_correction',
         'CA96_unit_determinant':str(det96),'prefix160_unit_determinant':str(det160),
         'pure_prefix160_system_count':7,'regular_coordinate_division_count':64,
         'B_N_shape':[8,9],'B_N_origin_rank':rank,'actual_generator1_maps_to1':True,
         'old_nonunit_r_g_composed_or_trace_bases_identified':False,
         'nil_correction_transition9_extension_or_T3_solved':False,
         'no_jets_GB_new_points_or_p_inverse':True,
         'input_sha256':{str(p.relative_to(ROOT)):sha(p) for p in paths},
         'evidence_sha256':{'honest9-to-N8-DAG.json':sha(out/'honest9-to-N8-DAG.json')},
         'DAG_nodes':len(D.nodes),'certificate_bytes':len(encoded.encode()),'seconds':time.monotonic()-started,
         'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    (out/'report.json').write_text(json.dumps(rep,indent=2)+'\n');print(json.dumps({k:rep[k] for k in ['all_passed','B_N_origin_rank','CA96_unit_determinant','prefix160_unit_determinant','DAG_nodes','certificate_bytes','seconds','max_rss_kib']}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
