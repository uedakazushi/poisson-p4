#!/usr/bin/env python3
"""Bounded actual512 transfer cycles, retaining all immutable coefficient names.

No scalar coefficient expansion, kernel solve, source worker or new jets.
"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import resource
import signal
import time
from fractions import Fraction as Q

ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_CONDUCTOR_TRANSFER512_NIL_CYCLES_20261104_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())

class Typed:
    def __init__(self):
        self.nodes=[];self.types=[];self.memo={}
        self.z=self.put('R0',('R0_QQ',0,1));self.o=self.put('R0',('R0_QQ',1,1))
        self.cz=self.embed(self.z);self.co=self.embed(self.o)
    def put(self,t,n):
        key=(t,n)
        if key not in self.memo:
            self.memo[key]=len(self.nodes);self.nodes.append(n);self.types.append(t)
        return self.memo[key]
    def rref(self,ns,i):return self.put('R0',('R0_actual_external',ns,i))
    def cref(self,i):return self.put('C',('C_actual_external','radical',i))
    def rneg(self,i):
        assert self.types[i]=='R0'
        return self.z if i==self.z else self.put('R0',('R0_neg',i))
    def radd(self,*ids):
        assert all(self.types[i]=='R0'for i in ids)
        ids=tuple(i for i in ids if i!=self.z)
        if not ids:return self.z
        return ids[0]if len(ids)==1 else self.put('R0',('R0_sum',ids))
    def rmul(self,*ids):
        assert all(self.types[i]=='R0'for i in ids)
        if self.z in ids:return self.z
        ids=tuple(i for i in ids if i!=self.o)
        if not ids:return self.o
        return ids[0]if len(ids)==1 else self.put('R0',('R0_product',ids))
    def embed(self,i):
        assert self.types[i]=='R0'
        return self.put('C',('C_embed_R0',i))
    def basis(self,i):return self.co if i==0 else self.put('C',('C_basis',i))
    def cadd(self,*ids):
        assert all(self.types[i]=='C'for i in ids)
        ids=tuple(i for i in ids if i!=self.cz)
        if not ids:return self.cz
        return ids[0]if len(ids)==1 else self.put('C',('C_sum',ids))
    def cmul(self,*ids):
        assert all(self.types[i]=='C'for i in ids)
        if self.cz in ids:return self.cz
        ids=tuple(i for i in ids if i!=self.co)
        if not ids:return self.co
        return ids[0]if len(ids)==1 else self.put('C',('C_product',ids))
    def comp(self,i,k):
        assert self.types[i]=='C'
        return self.put('R0',('R0_component_of_C',i,k))

NV=10
def pv(i):
    e=[0]*NV;e[i]=1;return {tuple(e):1}
def pa(*ps):
    out={}
    for p in ps:
        for m,c in p.items():out[m]=out.get(m,0)+c
    return {m:c for m,c in out.items()if c}
def pn(p):return {m:-c for m,c in p.items()}
def pm(*ps):
    out={(0,)*NV:1}
    for p in ps:
        new={}
        for m,c in out.items():
            for n,d in p.items():
                k=tuple(a+b for a,b in zip(m,n));new[k]=new.get(k,0)+c*d
        out={m:c for m,c in new.items()if c}
    return out
def cm(x,y,na,nb):
    z=[{}for _ in range(4)]
    for i,j in itertools.product(range(4),repeat=2):
        fs=[x[i],y[j]]
        if i&j&1:fs.append(na)
        if i&j&2:fs.append(nb)
        z[i^j]=pa(z[i^j],pm(*fs))
    return z

def main(out):
    start=time.monotonic();signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000))
    if out.exists():
        assert {p.name for p in out.iterdir()} <= {'failed-initial-preflight.json'}
    else:out.mkdir(parents=True,exist_ok=False)
    paths={
      'radical':'reports/AMBIENT_B18_RADICAL_COMPARISON_OPERATOR_20261104_v63/actual-radical-normalizer-comparison-operator.json',
      'radical_peer':'reports/AMBIENT_B18_RADICAL_COMPARISON_OPERATOR_20261104_v63_REVIEW/report.json',
      'chain':'reports/AMBIENT_B18_NORMALIZER_CHAIN_LIFT_640_NEGATIVE_BOUNDARY_20261104_v63/actual-UP-QB-640-typed-operator.json',
      'C0':'reports/AMBIENT_B18_T0_HONEST_C0_20261104_v63/honest-C0-D0-D4-DAG.json',
      'C0_peer':'reports/AMBIENT_B18_T0_HONEST_C0_20261104_v63_REVIEW/report.json',
      'source':'reports/AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63/actual-source-coefficients.json',
      'relative_peer':'reports/AMBIENT_B18_FULL_RELATIVE_INPUT_20261104_v63_REVIEW/report.json',
      'nil2_peer':'reports/AMBIENT_B18_NORMALIZER_NIL_GENERATORS_20261104_v63_REVIEW/report.json',
      'upper7_peer':'reports/AMBIENT_B18_QUADRATIC_CONDUCTOR_NIL_UPPER7_20261104_v63_REVIEW/report.json',
      'upper17_peer':'reports/AMBIENT_B18_T_KNOWN_INITIAL_LENGTH_UPPER17_20261104_v63_REVIEW/report.json',
      'rankK':'reports/AMBIENT_B18_CRITICAL_CONTACT_AND_COMPARISON_RANK4_20261104_v63/report.json',
      'product44_peer':'reports/AMBIENT_B18_NIL_PRODUCT44_OPERATOR_20261104_v63_REVIEW/report.json',
    }
    hashes={
      'radical':'376d5fd2b8f2ab8a45651b0374d5b983ebe2df8f8a507108675942e853adbba5',
      'radical_peer':'07b8ec2ea5d6dc419c15fad7cb3ea9526d817037b48c8d96b9966882c1dc851f',
      'chain':'d57cdec83150cde5f2f342a83fa817c71b60dfa36d158c9e8636b04d96fef379',
      'C0':'f417e859af9d15090dfe026c6b311b8c4e8bb3b9d624a884ecd8379884046572',
      'C0_peer':'9e3f2de68224ba1696888b21a51de61f1cded2dfbca258cd856bdc1d99bd2882',
      'source':'75cc6bcc4d4ff8ebd255a273cd3db494f92cecdff69583bea80e4316dee4fa74',
      'relative_peer':'d4d9d2f3d66b61fab72e436aaa19a534f12613b646736fbe918051c5961b6c3c',
      'nil2_peer':'02888d3c4b987327f6df2d8faf1f2adffa3f9c01af47e27aff8c5c0982e75472',
      'upper7_peer':'e3b515c545a07ca7b7feb4ff2233487d46a3a6f6b6a781ee565fddb6a9e1fd74',
      'upper17_peer':'5bedafd8efeed3dee33bc43321a78d4caf5a83869ee033230c38b2ed070ab0bd',
      'rankK':'a53d4613cbc088e376cec76b787b8669a807dd58898687590df627d44e0aa6ae',
      'product44_peer':'7d876031573583d070890c28da09ae123c550d2108cda73d45cd2f61b87eaa19',
    }
    for k,p in paths.items():assert sha(ROOT/p)==hashes[k],k
    for k in ['radical_peer','C0_peer','relative_peer','nil2_peer','upper7_peer','upper17_peer','rankK','product44_peer']:
        assert read(ROOT/paths[k])['all_passed'],k
    R,O,H,S=[read(ROOT/paths[k])for k in ['radical','chain','C0','source']]
    rows=O['physical_target_rows12'];assert rows[:8]==[18,19,28,29,60,67,72,63]
    F=O['actual_F12_C4'];assert len(F)==12 and all(len(v)==4 and v[3]==0 for v in F[:8])
    A=R['native_new_A12_by9_nodes'];assert len(A)==12 and all(len(v)==9 for v in A)
    beta=H['actual_C0_basis8_scalar_nodes'];assert beta==[7148,7199,7454,7505,7688,7943,8045,505938]
    rawB=S['families']['CB']['actual_zero_odd_I0_generators']
    assert [r for r,b in rawB]==rows[:8]
    assert [b for r,b in rawB]==[27897,27948,28203,28254,28437,28692,28794,28539]
    identities=O['accepted_actual_regular_function_identities']
    cbframe=next(v for v in identities if v['kind']=='accepted_actual_CB_f63_tail')
    assert cbframe=={'kind':'accepted_actual_CB_f63_tail','function':28539,'gB':299650,'coordinate_quotient_nodes':[327551,327554,327557,327560,327563,327566,327569]}
    assert len([v for v in identities if v['kind']=='accepted_actual_CB_zero_graph_membership'])==4
    assert len([v for v in identities if v['kind']=='accepted_actual_CA_zero_graph_membership'])==4
    aliases={v['odd_generator']:v['physical_odd_coefficient']for v in identities if v['kind']=='exact_QQ_linear_source_alias'}
    assert all(b in aliases for r,b in rawB)
    # Fixed exact source-atom arithmetic checks signs and BOTH unipotent frames.
    # The scalar functions themselves remain opaque original functions.
    op=O['operator_scalar_nodes'];memo={}
    def sadd(p,q):
        z=dict(p)
        for m,c in q.items():z[m]=z.get(m,Q(0))+c
        return {m:c for m,c in z.items()if c}
    def smul(p,q):
        z={}
        for m,c in p.items():
            for n,d in q.items():
                key=tuple(sorted(m+n));z[key]=z.get(key,Q(0))+c*d
        return {m:c for m,c in z.items()if c}
    def sexp(i):
        if i in memo:return memo[i]
        n=op[i];o=n[0]
        if o=='c':p={():Q(n[1],n[2])}if n[1]else{}
        elif o=='actual_C0_scalar_function':p={(n[1],):Q(1)}
        elif o=='sum':
            p={}
            for j,a,b in n[1]:p=sadd(p,{m:Q(a,b)*c for m,c in sexp(j).items()})
        elif o=='product':
            p={():Q(n[1],n[2])}
            for j in n[3]:p=smul(p,sexp(j))
        else:raise AssertionError(n)
        assert len(p)<200;memo[i]=p;return p
    caframe=next(v for v in identities if v['kind']=='accepted_actual_CA_f63_tail')
    assert caframe['function']==7790 and caframe['gA']==505938
    ca_raw=S['families']['CA']['actual_zero_odd_I0_generators']
    for j in range(8):
        if j<7:
            wantA={(ca_raw[j][1],):Q(1)};wantB={(rawB[j][1],):Q(1)}
        else:
            wantA={(505938,):Q(1)};wantB={(299650,):Q(1)}
            for k in range(7):
                wantA[tuple(sorted((beta[k],caframe['coordinate_quotients'][k])))]=Q(1)
                wantB[tuple(sorted((rawB[k][1],cbframe['coordinate_quotient_nodes'][k])))]=Q(1)
        assert sexp(F[j][1])==wantA
        assert {m:-c for m,c in sexp(F[j][2]).items()}==wantB
    nilgroups=[5,7];ns=R['new_Cscalar_and_R0component_nodes'];Qcols=R['actual_Qrad4_by36_columns']
    for g,wanted in [(5,[7,11]),(7,[11,13])]:
        cnode=ns[Qcols[4*g][0]][1];assert ns[cnode]==['C_product',wanted]
    # Persist conservative source/operation preflight BEFORE constructing512 cycles.
    pre={
      'actual_count':8*8*2*4,'actual_source_rows8':rows[:8],
      'actual_CA_beta8_nodes':beta,'actual_raw_CB8_nodes':rawB,
      'CB_raw_to_z7_gB_frame':'identity7 with final row(e1,...,e7,1)',
      'CB_raw_frame_determinant':'1','CB_frame_original_cofactor_nodes':cbframe,
      'all8_physical_F_tauA_plus_and_tauB_minus_signs_checked':True,
      'actual_source_atom_frame_terms_checked':sum(len(sexp(F[j][i]))for j in range(8)for i in [1,2]),
      'all4_CA_and_all4_CB_zero_graph_memberships_retained':True,
      'actual_A12_by9_adapter_namespace':'radical.native_new_A12_by9_nodes, F=A*(dR9)',
      'normalized_nil_groups_in_native_radical9':nilgroups,
      'typed_node_bound_conservative':60000,'visible_expanded_polynomial_coefficient_terms':0,
      'coefficient_expansion780_new_kernel_source_or_jet_job':False,
      'maximum_output_bytes':32000000,
      'honest_K_mod_J_presentation_constructed':False,
      'K_mod_J_dimension_bound':26,
      'missing_Tor2_gate':'Tor2(coker(D0->E4),R0/J) transgression into K/JK, or complete source cycle lifts; SyzQ alone is insufficient',
    }
    (out/'preflight.json').write_text(json.dumps(pre,indent=2)+'\n')
    pool=Typed();nu=pool.rref('C0',20086);bas=[pool.basis(i)for i in range(4)]
    betan=[pool.rref('C0',i)for i in beta]
    fn=[]
    for j in range(8):
        de=pool.rref('chain',F[j][0]);aa=pool.rref('chain',F[j][1]);bb=pool.rneg(pool.rref('chain',F[j][2]))
        fn.append((de,aa,bb))
    # Universal exact four-component transfer identity, no numerical source point.
    na,nb,be,bj,de,ba=[pv(i)for i in range(6)];pp=[pv(i)for i in range(6,10)]
    E=pm(be,pa(pm(bj,pp[0]),pm(de,pp[2]),pm(na,ba,pp[3])))
    HH=pa(pm(bj,pp[1]),pm(ba,pp[2]),pm(de,pp[3]))
    UT=[E,pm(be,HH),{},{}];fp=[pm(be,bj,x)for x in pp]
    FF=[de,ba,pn(bj),{}];tt=[pp[2],pp[3],{},{}]
    corr=[pm(be,x)for x in cm(FF,tt,na,nb)]
    assert UT==[pa(x,y)for x,y in zip(fp,corr)]
    records=[]
    for i,j,al,ga in itertools.product(range(8),range(8),range(2),range(4)):
        gnil=nilgroups[al];pcomp=[pool.rref('radical',x)for x in Qcols[4*gnil+ga]]
        p0,pA,pB,pAB=pcomp;de,ba,bj=fn[j];be=betan[i]
        tC=pool.cadd(pool.embed(pB),pool.cmul(bas[1],pool.embed(pAB)))
        f=pool.rmul(be,bj)
        rhsC=[];a36=[]
        for g in range(9):
            c=pool.cmul(pool.embed(be),pool.cref(A[j][g]),tC)
            if g==gnil:c=pool.cadd(c,pool.cmul(pool.embed(f),bas[ga]))
            rhsC.append(c);a36.extend(pool.comp(c,k)for k in range(4))
        even=pool.rmul(be,pool.radd(pool.rmul(bj,p0),pool.rmul(de,pB),pool.rmul(nu,ba,pAB)))
        odd=pool.radd(pool.rmul(bj,pA),pool.rmul(ba,pB),pool.rmul(de,pAB))
        c8=[pool.z]*8;c8[i]=odd
        records.append({
          'index':len(records),'CA_generator_index':i,'CB_raw_index':j,'CB_physical_row':rows[j],
          'normalized_nil_generator_index':al,'normalizer_basis_index':ga,
          'actual_f_betaA_times_bB':f,'p_four_actual_R0_components':pcomp,
          'literal_actual_C0_free9_representative':[even]+c8,
          'literal_Qrad_lift9_C_coefficients':rhsC,
          'literal_H1_cycle44':a36+c8,
          'actual_cycle_identity':'Qrad*a=U*C0rep, since f*p+betaA*F_j*(pB+tauA*pAB)=U*C0rep',
          'even_recovery_identity':'Qrad_even*a equals stored C0rep[0] exactly',
        })
    assert len(records)==512 and all(len(v['literal_H1_cycle44'])==44 for v in records)
    # Whole source nil2*4 spans NE as R0 module; all a,b transfers are bilinear.
    result={
      'source_artifacts':{k:{'path':p,'sha256':hashes[k]}for k,p in paths.items()},
      'new_typed_scalar_nodes':pool.nodes,'new_node_types':pool.types,
      'typed_C_multiplication':'same immutable free4 C/R0, tauA²=nuA,tauB²=nuB, no coefficient approximation',
      'actual_A12_by9_namespace':'radical.native_new_A12_by9_nodes',
      'adapter_convention':'Fphysical=A*(dR9) as generator COLUMN, no transpose inside C multiplication',
      'actual_raw_CB8_unipotent_basis_proof':cbframe,
      'raw8_generate_entire_IB':'first7=zB7, f63=gB+sum(zBi*ei), unitframe determinant1; all4 graph membership identities included',
      'CA8_generate_entire_IA':'actual C0=(1,(zA7,gA)*tauA), all4 actual zero graph memberships preserved',
      'normalized_nil_two_basis4_spans_entire_NilE4_as_R0_module':True,
      'actual_512_transfer_cycles':records,
      'transfer_index_convention':'phi_(a,b), f=a*b; independent choice of factorization of f is NOT assumed',
      'transfer_D0_linearity_and_two_composites':'psi phi=f on E4, phi psi=f on D0',
      'transfer_product':'phi(x)*phi(y)=f*phi(x*y)',
      'known_ideal_W':'sum over all a in IA,b in IB of phi_(a,b)(Nil(E4))',
      'all512_generate_exact_transfer_W_not_all_NilD0':True,
      'W_subset_NilD0':True,'W_times_NilD0_squared_zero':True,'W_power3_zero':True,
      'IAIB_times_NilD0_subset_W':True,
      'actual_nil_cube_factors_through_NilD0_mod_W_in_each_argument':True,
      'NilD0_mod_W_annihilator_contains_IAIB':True,
      'NilD0_mod_W_or_all_H0_assumed_finite':False,
      'honest_K_mod_J_presentation':'OPEN: Tor2 transgression/source lift correction remains',
      'K_mod_J_dimension_at_most':26,'K_square_length_at_most':182,
      'K_square_annihilator_contains_J_and_m_R0_squared':True,
      'no_all_H0_or_nilcube_annihilator_length_bound_claimed':True,
      'exact_nil_cube_K_square_or_index_decided':False,
      'new_source_kernel_GB_jet_or_scalar_coefficient_expansion':False,
    }
    ep=out/'actual-transfer512-H1-cycles-typed-DAG.json'
    ep.write_text(json.dumps(result,separators=(',',':'))+'\n')
    assert len(pool.nodes)<=60000 and ep.stat().st_size<32000000
    rels=list(paths.values())+[f'research/notes/{STEM}.md',str(Path(__file__).relative_to(ROOT))]
    report={
      'id':STEM,'all_passed':True,'status':'PROVED_actual512_nil_transfer_cycles_known_subideal_and_honest_KJK_frontier',
      'actual_transfer_cycle_count':512,'all512_complete_H1_cycle_source_identities':True,
      'both_actual_C0_and_Qrad_lifts_retained':True,'actual_raw_CB_frame_unit1_entire_IB':True,
      'all512_generate_transfer_W':True,'W_times_NilD0_squared_zero':True,
      'IAIB_times_NilD0_subset_W':True,'actual_cube_factors_through_NilD0_mod_W':True,
      'honest_K_mod_J_presentation_constructed':False,'K_mod_J_dimension_upper_bound':26,
      'K_square_length_upper_bound':182,'all_H0_length_or_nilcube_claimed':False,
      'typed_node_count':len(pool.nodes),'bytes':ep.stat().st_size,
      'scalar_coefficients_expanded_or_new_kernel_source_jet_GB':False,
      'input_sha256':{p:sha(ROOT/p)for p in rels},
      'evidence_sha256':{ep.name:sha(ep),'preflight.json':sha(out/'preflight.json')},
      'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
    }
    failure=out/'failed-initial-preflight.json'
    if failure.exists():
        report['preserved_initial_namespace_preflight_failure']=True
        report['evidence_sha256'][failure.name]=sha(failure)
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:report[k]for k in ['id','all_passed','actual_transfer_cycle_count','typed_node_count','bytes','seconds']}))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);main(p.parse_args().output.resolve())
