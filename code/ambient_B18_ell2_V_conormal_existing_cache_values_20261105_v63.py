#!/usr/bin/env python3
"""One cache-only composition of actual V conormal constants; no solve/jet."""
import argparse,collections,hashlib,json,math,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
from ambient_B18_ell2_raw33_quotient_20261105_v63 import add,scale,mv,srec
from ambient_B18_conductor_product_special_fiber20_20261104_v63 import compose
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_ELL2_V_CONORMAL_EXISTING_CACHE_VALUES_20261105_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def parse(v):return {i:F(c) for i,c in v}
def guard(v):
    assert all(max(abs(x.numerator).bit_length(),x.denominator.bit_length())<=384 for x in v.values()),'384bit cap'
    return v
class Jet56:
    def __init__(self,basis):
        self.basis=[tuple(m) for m in basis];self.idx={m:i for i,m in enumerate(self.basis)};self.one={0:F(1)}
        self.table={(i,j):self.idx[p] for i,m in enumerate(self.basis) for j,n in enumerate(self.basis) if (p:=tuple(x+y for x,y in zip(m,n))) in self.idx}
    def mul(self,a,b):
        o={}
        for i,c in a.items():
            for j,d in b.items():
                k=self.table.get((i,j))
                if k is not None:o=add(o,{k:c*d})
        return guard(o)
    def inv(self,a):
        c=a.get(0,F(0));assert c;n=dict(a);n.pop(0,None);n=scale(n,-1/c);v=self.one;o=dict(v)
        for _ in range(3):v=self.mul(v,n);o=add(o,v)
        o=guard(scale(o,1/c));assert self.mul(a,o)==self.one;return o
    def from_CA(self,p):
        # Original source CA variables(z7,r,s) -> (z7,r,ell), s=r+ell/5.
        o={}
        for m,c in p.items():
            for l in range(min(1,m[8])+1):
                n=tuple(m[:7])+(m[7]+m[8]-l,l)
                if n in self.idx:o=add(o,{self.idx[n]:c*math.comb(m[8],l)/F(5)**l})
        return guard(o)
def main(out,peer_path):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024))
    out=out.resolve();out.mkdir(parents=True,exist_ok=False);hs={}
    first=out.with_name(out.name+'_FIRST_LINEAR_READER_ATTEMPT')
    if first.exists():
        for p in first.iterdir():
            if p.is_file():hs[str(p.relative_to(ROOT))]=sha(p)
    def read(p):hs[p]=sha(ROOT/p);return json.loads((ROOT/p).read_text())
    def write(name,obj):
        p=out/name;p.write_text(json.dumps(obj,separators=(',',':'))+'\n');assert sum(x.stat().st_size for x in out.iterdir() if x.is_file())<32*1024*1024;return sha(p)
    write('phase.json',{'phase':'source peer/input gate','source_coefficients_evaluated':False})
    gate=read('reports/AMBIENT_B18_ELL2_V_CONORMAL_EXISTING_CACHE_GATE_20261105_v63/report.json');assert gate['all_passed'] and gate['whole_slot_bound']==92807
    peer=read(peer_path);assert peer['all_passed']
    ppins={}
    for k,v in peer.get('input_sha256',{}).items():
        p,h=(v['path'],v['sha256']) if isinstance(v,dict) else (k,v)
        assert sha(ROOT/p)==h;ppins[p]=h
    profilepath='reports/AMBIENT_B18_ELL2_V_CONORMAL_HADAMARD_SOURCE_PROFILE_20261105_v63/report.json'
    assert ppins.get(profilepath)==sha(ROOT/profilepath),'peer must bind actual conormal profile andprefixgA source'
    G=read('reports/AMBIENT_B18_ELL2_V_CONORMAL_EXISTING_CACHE_GATE_20261105_v63/whole-existing-cache-V-conormal-source-gate.json')
    adj=read('reports/AMBIENT_B18_ELL2_V_CONORMAL_HADAMARD_SOURCE_PROFILE_20261105_v63/actual26-conormal-and-symbolic-V-Hadamard-source-adjoints.json')
    full=read('reports/AMBIENT_B18_T3_FULL9_FIRSTJET_20261104_v63/full9-firstjet-certificate.json')
    joint=read('reports/AMBIENT_B18_T3_JOINT_CURVE_RESTRICTIONS_20261104_v63/joint-degree3-certificate.json')
    Y=read('reports/AMBIENT_B18_CONDUCTOR_PRODUCT_SPECIAL_FIBER20_20261104_v63/actual-special-fiber20-coefficients.json')
    oldP=read('reports/AMBIENT_B18_RAW_COKERNEL_ARTIN20_SOURCE_PREFLIGHT_20261104_v63/actual-Artin20-source-and-cost-preflight.json')
    S=read('reports/AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63/actual-source-coefficients.json')
    chain=read('reports/AMBIENT_B18_NORMALIZER_CHAIN_LIFT_640_NEGATIVE_BOUNDARY_20261104_v63/actual-UP-QB-640-typed-operator.json')
    for p,h in gate['input_sha256'].items():assert sha(ROOT/p)==h
    ring=Jet56(G['coefficient_ring_basis56']);zero=(0,)*9
    mons9=[]
    for inds in full['source2_monomials']:
        m=[0]*9
        for i in inds:m[i]+=1
        mons9.append(tuple(m))
    own={}
    for fam in ['CA','CB']:
        own[fam]={i:{m:F(c) for m,c in zip(mons9,cs) if F(c)} for i,cs in full['both57_source2'][fam]['all57_source2jets']}
    YB=[add({tuple(int(k==j) for k in range(9)):F(c) for j,c in enumerate(a) if F(c)},
            {tuple(m):F(c) for m,c in b}) for a,b in zip(Y['YB_in_YA_linear'],Y['YB_in_YA_quadratic'])]
    cb={i:compose(p,YB,2) for i,p in own['CB'].items()}
    shared=sorted(set(cb)&set(own['CA']));assert shared==list(range(2,11))
    assert all(cb[i]==own['CA'][i] for i in shared),'exact FULL9 nonlinear physical copy'
    second=dict(own['CA']);second.update({i:p for i,p in cb.items() if i not in second})
    ids=joint['all105_unknown_original_node_order'];unknown={i:ring.from_CA(second[i]) for i in ids}
    plane={i:{(0,)*7+tuple(m):F(c) for m,c in zip(joint['monomials'],cs) if F(c)} for i,cs in joint['all105_unknown_jets']}
    overlaps=[]
    for i in ids:
        p=ring.from_CA(plane[i]);old={j:c for j,c in unknown[i].items() if not any(ring.basis[j][:7])}
        low={j:c for j,c in p.items() if sum(ring.basis[j])<=2}
        assert low==old,('common plane cache overlap',i)
        unknown[i]=add(unknown[i],{j:c for j,c in p.items() if sum(ring.basis[j])==3})
        overlaps.append(i)
    write('phase1-actual105-B56-cache-composition.json',{'basis56':ring.basis,'unknown105_order':ids,
      'all105_cached_values_B56':[srec(unknown[i]) for i in ids],
      'all9_shared_physical_inverse_source2_exact':True,'all105_commonplane_lower_overlaps_exact':True,
      'no_unknown_coefficient_solved':True})
    nodes=S['DAG_scalar_nodes'];order=oldP['needed_closure_order'];roots=oldP['requested_root_nodes'];uses=collections.Counter(roots)
    for i in order:
        if nodes[i][0] in ['a','m','n','i']:uses.update(nodes[i][1:])
    val={};inv=[];peak=0;peakQQ=0;qq=0
    for i in order:
        n=nodes[i];op=n[0]
        if op=='v':v=unknown[i]
        elif op=='c':
            c=F(n[1],n[2]) if len(n)==3 else F(n[1]);v={0:c} if c else {}
        elif op=='a':v=add(val[n[1]],val[n[2]])
        elif op=='m':v=ring.mul(val[n[1]],val[n[2]])
        elif op=='n':v=scale(val[n[1]],-1)
        elif op=='i':v=ring.inv(val[n[1]]);inv.append([i,n[1],srec(val[n[1]]),srec(v)])
        else:raise ValueError(n)
        val[i]=guard(v);qq+=len(v);peak=max(peak,len(val));peakQQ=max(peakQQ,qq);assert peak<=414 and peakQQ<=23184
        if op in ['a','m','n','i']:
            for j in n[1:]:
                uses[j]-=1
                if not uses[j]:qq-=len(val[j]);del val[j]
    def var(i):m=[0]*9;m[i]=1;return {ring.idx[tuple(m)]:F(1)}
    targets=[var(i) for i in range(8)]+[add(var(7),var(8),F(1,5))]
    eqids=joint['all105_equation_original_node_order'];res=[dict(val[i]) for i in eqids]
    for k,t in enumerate(targets):res[48+k]=add(res[48+k],t,-1)
    assert all(not v for v in res),'all105 source equations in B56 must hold without a solve'
    nuA,nuB=[val[i] for i in oldP['nu_roots_CA_CB']];V=[];Fdata=[]
    assert chain['physical_target_rows12']==[x['physical_row'] for x in oldP['physical12_coefficient_roots']]
    assert chain['all48_existing_Q_frame_columns_verified']
    for q in oldP['physical12_coefficient_roots']:
        aA,aB,bA,bB=[val[q[k]] for k in ['aA_node','aB_node','bA_node','bB_node']];delta=add(aA,aB,-1)
        V.extend([[bA,scale(bB,-1),{}],[delta,{},scale(bB,-1)],[{},delta,bA],
          [scale(ring.mul(nuB,bB),-1),ring.mul(nuA,bA),delta]])
        Fdata.append({'row':q['physical_row'],'roots':q,'delta':srec(delta),'bA':srec(bA),'bB':srec(bB)})
    V.extend([[scale(var(i),-1),{},{}] for i in range(7)])
    ga={}
    for m,c in zip(joint['monomials'],joint['g_CA_all10_coefficients']):
        if F(c):ga[(0,)*7+tuple(m)]=F(c)
    gA=ring.from_CA(ga);assert all(sum(ring.basis[j])==3 for j in gA)
    V.append([scale(gA,-1),{},{}]);assert len(V)==56 and all(not v.get(0,0) for col in V for v in col)
    matrix=[]
    for col in V:
        outv={}
        for t in adj['scalar_monomial_adjoint_library']:
            assert tuple(t['monomial9']) in ring.idx
            c=col[t['physical_target']].get(ring.idx[tuple(t['monomial9'])],F(0))
            if c:outv=add(outv,parse(t['conormal26_value']),c)
        matrix.append(guard(outv))
    piv={};reprs={};receipts=[]
    for j,orig in enumerate(matrix):
        v=dict(orig);w={j:F(1)};steps=[]
        while v:
            p=min(v)
            if p not in piv:
                d=v[p];piv[p]=guard(scale(v,1/d));reprs[p]=guard(scale(w,1/d))
                receipts.append({'source_column':j,'pivot':p,'divide_by':str(d),'subtract':[[k,str(c)] for k,c in steps],'normalized':srec(piv[p]),'actual_source_column_combination':srec(reprs[p])});break
            d=v[p];steps.append((p,d));v=add(v,piv[p],-d);w=add(w,reprs[p],-d)
    assert len(piv)==26,'source complete V conormal matrix must have rank26'
    selected=[p['source_column'] for p in receipts];sidx={j:i for i,j in enumerate(selected)}
    factor=[]
    for orig in matrix:
        v=dict(orig);w={}
        while v:
            p=min(v);d=v[p];w=add(w,reprs[p],d);v=add(v,piv[p],-d)
        factor.append({sidx[j]:c for j,c in w.items()})
    actualbasis=[matrix[j] for j in selected];assert all(mv(actualbasis,v)==orig for v,orig in zip(factor,matrix))
    pivot_sequence=[p['pivot'] for p in receipts];inversions=sum(pivot_sequence[i]>pivot_sequence[j] for i in range(26) for j in range(i+1,26))
    det=F((-1)**inversions)
    for p in receipts:det*=F(p['divide_by'])
    guard({0:det});assert det
    write('phase2-all105-equations-and-primitive56-source-columns.json',{'basis56':ring.basis,
      'all105_equation_residuals':[srec(v) for v in res],'all105_equations_exact':True,'six_unit_inverse_products':inv,
      'all164_original_source_roots':[[i,srec(val[i])] for i in roots],'primitive_F12':Fdata,
      'nuA':srec(nuA),'nuB':srec(nuB),'actual_new_gA_plane_through3':srec(gA),
      'all56_primitive_V_columns':[[srec(v) for v in col] for col in V],
      'exact_allorder_primitive_frame_binding':'accepted chain48source equality, sourcephysical12rootsandprefix0gA505938=7790',
      'no_membership_quotient_derivatives_evaluated':True})
    evidence={'actual26_by56_columns':[srec(v) for v in matrix],'rank':26,'fixed_selected26_columns':selected,
      'fixed26_minor':[[str(matrix[j].get(i,F(0))) for j in selected] for i in range(26)],
      'fixed26_minor_determinant':str(det),'fixed_rank_receipts':receipts,
      'all56_columns_in_selected26_basis':[srec(v) for v in factor],'all56_column_factorizations_exact':True,
      'formal_unit_section_consequence':'the source-regular minimalcomparisonm hasthisconstantmatrix; selected26minor is a formalunit, so n=Eselected*(mselected)^-1, m*n=I26, h=0',
      'source_m_operator':'ordered matrix Hadamard+periodicaugmentation withactualregularV, thenknownunitC2resolutionprojection',
      'full_T2_coefficient_presentation_or_length_not_computed':True,'no_new_implicit_unknown_solve':True}
    write('actual-V-conormal26-by56-unit-minor-and-source-factorization.json',evidence)
    write('phase.json',{'phase':'complete','rank26':26})
    slots=sum(len(v) for v in unknown.values())+sum(len(v) for v in val.values())+sum(len(v) for col in V for v in col)+sum(map(len,matrix))+sum(map(len,factor))
    assert slots+18584+23184+7815+3136+1456+1224+672<92807
    report={'id':STEM,'all_passed':True,'actual_V_conormal_rank':26,'fixed26_minor_determinant':str(det),
      'all56_column_factorizations_exact':True,'all105_cached_source_equations_verified':True,
      'new_implicit_source_or_gauge_solve':False,'actual_T2_or_H0_length_computed':False,
      'source_live_nodes':peak,'source_live_QQ':peakQQ,'named_QQ_slots':slots,'whole_slot_bound':92807,
      'input_sha256':hs,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),
      'evidence_sha256':{p.name:sha(p) for p in out.iterdir() if p.is_file() and p.name!='phase.json'},'seconds':time.monotonic()-start}
    write('report.json',report)
    print(json.dumps({k:report[k] for k in ['all_passed','actual_V_conormal_rank','fixed26_minor_determinant','named_QQ_slots','seconds']}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);p.add_argument('--source-peer',required=True);a=p.parse_args();main(a.output,a.source_peer)
