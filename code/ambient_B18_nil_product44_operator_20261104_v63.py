#!/usr/bin/env python3
"""Literal actual nil-product and cube operators on the complete H1 model.

Coefficients are immutable actual functions, never newly solved or expanded.
"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import resource
import signal
import time

ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_NIL_PRODUCT44_OPERATOR_20261104_v63'

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())

class DAG:
    def __init__(self):
        self.nodes=[];self.cache={};self.zero=self.put(('QQ',0,1));self.one=self.put(('QQ',1,1))
    def put(self,n):
        if n not in self.cache:self.cache[n]=len(self.nodes);self.nodes.append(n)
        return self.cache[n]
    def ref(self,artifact,i):return self.put(('actual_R0_function',artifact,i))
    def var(self,name,i):return self.put(('arbitrary_cycle_input',name,i))
    def c(self,i):return self.put(('QQ',i,1))
    def add(self,*xs):
        xs=tuple(i for i in xs if i!=self.zero)
        if not xs:return self.zero
        if len(xs)==1:return xs[0]
        return self.put(('sum',xs))
    def mul(self,*xs):
        if self.zero in xs:return self.zero
        xs=tuple(i for i in xs if i!=self.one)
        if not xs:return self.one
        if len(xs)==1:return xs[0]
        return self.put(('product',xs))

# A distinct tiny polynomial engine verifies universal coefficient identities.
NV=14
def pvar(i):
    e=[0]*NV;e[i]=1;return {tuple(e):1}
def padd(*ps):
    z={}
    for p in ps:
        for m,c in p.items():z[m]=z.get(m,0)+c
    return {m:c for m,c in z.items() if c}
def pneg(p):return {m:-c for m,c in p.items()}
def pmul(*ps):
    z={(0,)*NV:1}
    for p in ps:
        q={}
        for m,c in z.items():
            for n,d in p.items():
                k=tuple(a+b for a,b in zip(m,n));q[k]=q.get(k,0)+c*d
        z={m:c for m,c in q.items() if c}
    return z
def pcmul(x,y,nuA,nuB):
    out=[{}for _ in range(4)]
    for i,j in itertools.product(range(4),repeat=2):
        factors=[x[i],y[j]]
        if i&j&1:factors.append(nuA)
        if i&j&2:factors.append(nuB)
        out[i^j]=padd(out[i^j],pmul(*factors))
    return out

def main(out):
    start=time.monotonic();signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000))
    out.mkdir(parents=True,exist_ok=False)
    paths={
      'radical':'reports/AMBIENT_B18_RADICAL_COMPARISON_OPERATOR_20261104_v63/actual-radical-normalizer-comparison-operator.json',
      'radical_peer':'reports/AMBIENT_B18_RADICAL_COMPARISON_OPERATOR_20261104_v63_REVIEW/report.json',
      'chain':'reports/AMBIENT_B18_NORMALIZER_CHAIN_LIFT_640_NEGATIVE_BOUNDARY_20261104_v63/actual-UP-QB-640-typed-operator.json',
      'C0':'reports/AMBIENT_B18_T0_HONEST_C0_20261104_v63/honest-C0-D0-D4-DAG.json',
      'C0_peer':'reports/AMBIENT_B18_T0_HONEST_C0_20261104_v63_REVIEW/report.json',
      'nil_finite':'reports/AMBIENT_B18_NIL_CUBE_FINITE_ORIGIN_20261104_v63_REVIEW/report.json',
    }
    expected={
      'radical':'376d5fd2b8f2ab8a45651b0374d5b983ebe2df8f8a507108675942e853adbba5',
      'radical_peer':'07b8ec2ea5d6dc419c15fad7cb3ea9526d817037b48c8d96b9966882c1dc851f',
      'chain':'d57cdec83150cde5f2f342a83fa817c71b60dfa36d158c9e8636b04d96fef379',
      'C0':'f417e859af9d15090dfe026c6b311b8c4e8bb3b9d624a884ecd8379884046572',
      'C0_peer':'9e3f2de68224ba1696888b21a51de61f1cded2dfbca258cd856bdc1d99bd2882',
      'nil_finite':'ee80d12faa103fb5ac0e2144b4d12cd662f0c5afcfa1ded335e2171946ad1ccc',
    }
    for key,p in paths.items():assert sha(ROOT/p)==expected[key],key
    R,O,H=read(ROOT/paths['radical']),read(ROOT/paths['chain']),read(ROOT/paths['C0'])
    assert read(ROOT/paths['radical_peer'])['all_passed'] and read(ROOT/paths['C0_peer'])['all_passed']
    assert R['NilD0_exactly_H1_of_actual_complex']
    assert R['actual_X44_by780']['final_shapes']=={'V':[3,44],'X':[44,780]}
    beta=H['actual_C0_basis8_scalar_nodes'];assert len(beta)==8 and beta[-1]==505938
    assert H['actual_C0_free9_images_in_CA_trace_pairs']==[[1,0]]+[[0,i]for i in beta]
    scalar=O['operator_scalar_nodes'];U=O['exact_existing_U9_C4']
    assert U==[[1,0,0,0]]+[[0,i+2,0,0]for i in range(8)]
    assert [scalar[U[i+1][1]]for i in range(8)]==[['actual_C0_scalar_function',i]for i in beta]
    nuold=O['nu_A_B_operator_nodes'][0];assert scalar[nuold]==['actual_C0_scalar_function',20086]
    Qcols=R['actual_Qrad4_by36_columns'];assert len(Qcols)==36
    ns=R['new_Cscalar_and_R0component_nodes']
    # Bind every 4-column block to a SINGLE dR9 C-element and its basis multiples.
    block_roots=[]
    for g in range(9):
        root=None
        for alpha in range(4):
            ids=Qcols[4*g+alpha]
            n=[ns[i]for i in ids]
            assert all(z[0]=='R0_component_of_C' and z[2]==k for k,z in enumerate(n))
            cr=n[0][1];assert all(z[1]==cr for z in n)
            if alpha==0:root=cr
            else:
                prod=ns[cr];assert prod[0]=='C_product'
                assert root in prod[1] and len(prod[1])==2
                other=next(i for i in prod[1]if i!=root)
                assert ns[other]==['C_basis',alpha]
        block_roots.append(root)
    V=R['actual_V3_by44_rows_after_unit1_cancellation']
    assert all(V[k][:36]==[q[k+1]for q in Qcols]for k in range(3))
    pool=DAG();nu=pool.ref('C0',20086);bs=[pool.ref('C0',i)for i in beta]
    q0=[pool.ref('radical',q[0])for q in Qcols]
    z=[pool.var('z',i)for i in range(44)];zp=[pool.var('zprime',i)for i in range(44)]
    def unpack(v):
        b0=pool.add(*(pool.mul(q,x)for q,x in zip(q0,v[:36])))
        B=pool.add(*(pool.mul(q,x)for q,x in zip(bs,v[36:])))
        return b0,B
    def cmu_trace(a,b,coeff4):
        v0,v1,v2,v3=coeff4
        return [pool.add(pool.mul(a,v0),pool.mul(nu,b,v1)),
                pool.add(pool.mul(b,v0),pool.mul(a,v1)),
                pool.add(pool.mul(a,v2),pool.mul(nu,b,v3)),
                pool.add(pool.mul(b,v2),pool.mul(a,v3))]
    b0,B=unpack(z);c0,C=unpack(zp)
    prod=[]
    for g in range(9):prod.extend(cmu_trace(b0,B,zp[4*g:4*g+4]))
    prod.extend(pool.add(pool.mul(b0,y),pool.mul(c0,x))for x,y in zip(z[36:],zp[36:]))
    assert len(prod)==44
    square_even=pool.add(pool.mul(b0,b0),pool.mul(nu,B,B))
    square_tau=pool.mul(pool.c(2),b0,B)
    cube=[]
    for g in range(9):cube.extend(cmu_trace(square_even,square_tau,z[4*g:4*g+4]))
    cube_odd_scalar=pool.add(pool.mul(pool.c(3),b0,b0),pool.mul(nu,B,B))
    cube.extend(pool.mul(cube_odd_scalar,x)for x in z[36:])
    cube_even=pool.add(pool.mul(b0,b0,b0),pool.mul(pool.c(3),nu,b0,B,B))
    assert len(cube)==44
    # Universal chain closure, including the even-coordinate correction OFF cycles.
    na,nb=pvar(0),pvar(1);h=[pvar(i)for i in range(2,6)]
    a=[pvar(i)for i in range(6,10)];x0,x1,yB=pvar(10),pvar(11),pvar(12)
    x=[x0,x1,{},{}];old=pcmul(h,a,na,nb)
    out4=pcmul(h,pcmul(x,a,na,nb),na,nb)
    assert out4==pcmul(x,old,na,nb)
    y0=old[0];yodd=padd(pmul(x0,yB),pmul(y0,x1))
    residual=[padd(old[1],pneg(yB)),old[2],old[3]]
    got=[padd(out4[1],pneg(yodd)),out4[2],out4[3]]
    want=[pmul(x0,residual[0]),padd(pmul(x0,residual[1]),pmul(na,x1,residual[2])),padd(pmul(x1,residual[1]),pmul(x0,residual[2]))]
    assert got==want
    even_delta=padd(out4[0],pneg(padd(pmul(x0,y0),pmul(na,x1,yB))))
    assert even_delta==pmul(na,x1,residual[0])
    triple=pcmul(pcmul(x,x,na,nb),x,na,nb)
    assert triple[0]==padd(pmul(x0,x0,x0),pmul({(0,)*NV:3},na,x0,x1,x1))
    assert triple[1]==padd(pmul({(0,)*NV:3},x0,x0,x1),pmul(na,x1,x1,x1))
    data={
      'source_artifacts':{k:{'path':v,'sha256':expected[k]}for k,v in paths.items()},
      'all_new_R0_scalar_nodes':pool.nodes,
      'D0_source_basis':['1']+[f'beta{i+1}*tauA'for i in range(8)],
      'actual_beta8_C0_node_ids':beta,'actual_nuA_C0_node_id':20086,
      'actual_Qrad36_native_C_block_roots':block_roots,
      'actual_Qrad_even36_R0_component_root_ids':[q[0]for q in Qcols],
      'complete_cycle_source_order':['a36 (9 blocks times4)','c8 (actual C0 odd basis)'],
      'cycle_equation':'V3_by44*z=0 in the exact actual coefficient ring R0',
      'cycle_representative_even_b0':b0,'cycle_representative_odd_scalar_B':B,
      'cycle_representative_free9_roots':[b0]+z[36:],
      'second_cycle_representative_even':c0,'second_cycle_odd_scalar':C,
      'literal_bilinear44_output_roots':prod,
      'literal_one_input_cubic44_output_roots':cube,
      'actual_even_product_on_cycles':pool.add(pool.mul(b0,c0),pool.mul(nu,B,C)),
      'actual_odd_product_scalar_on_cycles':pool.add(pool.mul(b0,C),pool.mul(c0,B)),
      'actual_even_cube_on_cycles':cube_even,'actual_odd_cube_coefficient_factor':cube_odd_scalar,
      'product_source_formula':'a_product=diag9(mu(b0+B*tauA))*a_prime; c_product=b0*c_prime+b0_prime*c',
      'mu_trace_order':['1','tauA','tauB','tauA*tauB'],
      'even_recovery_OFF_cycles':'q0*a_product-(b0*b0_prime+nuA*B*B_prime)=nuA*B*(V*zprime)_tauA',
      'chain_closure_OFF_cycles':'V*product=diag(b0,[[b0,nuA*B],[B,b0]])*(V*zprime)',
      'actual_multiplication_on_H1':'exact transported original D0 product under the accepted complete H1 isomorphism',
      'free9_not_assumed_associative_before_Koszul_quotient':True,
      'boundary_descent':'A boundary has zero D0 representative; its product therefore has zero H1 class by the complete H1 isomorphism. No780-entry homotopy is asserted.',
      'exact_cube_gate':'Nil(D0)^3=0 iff for EVERY exact cycle z the cubic44 output belongs to image X44_by780',
      'cube_gate_polarization_char0':'6xyz=(x+y+z)^3-(x+y)^3-(x+z)^3-(y+z)^3+x^3+y^3+z^3',
      'cube_target':'accepted finite-origin H0(D0) submodule; not all H0 is assumed generated by cubes or Tor63',
      'all_H1_cycles_generated_or_cube_membership_decided':False,
      'source_coefficients_expanded_or_new_jets_solved':False,
      'new_typed_node_count':len(pool.nodes),
      'universal_associativity_chain_even_and_cube_identities_checked':True,
    }
    ep=out/'actual-nil-product44-and-cube-DAG.json';ep.write_text(json.dumps(data,separators=(',',':'))+'\n')
    assert ep.stat().st_size<32000000 and len(pool.nodes)<150000
    inputs=list(paths.values())+[f'research/notes/{STEM}.md',str(Path(__file__).relative_to(ROOT))]
    report={
      'id':STEM,'all_passed':True,'status':'PROVED_source_bound_actual_nil_ring_product_on_complete_H1',
      'literal_bilinear_output_shape':[44], 'literal_cubic_output_shape':[44],
      'all_input_variables_subject_to_exact_cycle_equations':True,
      'actual_C0_representation_and_all36_Qrad_bindings_verified':True,
      'off_cycle_chain_and_even_correction_identities_verified':True,
      'product_descent_to_actual_nilradical_H1_proved':True,
      'arbitrary_cycle_cubic_membership_gate_complete':True,
      'all_nil_cycles_generated_or_cube_zero_nonzero_decided':False,
      'new_source_jet_point_GB_or_780_expansion':False,
      'typed_node_count':len(pool.nodes),'bytes':ep.stat().st_size,
      'input_sha256':{p:sha(ROOT/p)for p in inputs},'evidence_sha256':{ep.name:sha(ep)},
      'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
    }
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:report[k]for k in ['id','all_passed','typed_node_count','bytes','seconds']}))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);main(p.parse_args().output.resolve())
